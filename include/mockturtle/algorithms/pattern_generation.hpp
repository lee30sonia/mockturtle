/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file pattern_generation.hpp
  \brief Powerful Simulation Pattern Generation

  \author Siang-Yun Lee
*/

#pragma once

#include <random>
#include <mockturtle/networks/aig.hpp>
#include "../utils/progress_bar.hpp"
#include "../utils/stopwatch.hpp"
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/dont_cares.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/print.hpp>
#include <kitty/bit_operations.hpp>
#include "../utils/node_map.hpp"
#include "cnf.hpp"
#include "cleanup.hpp"
#include <percy/solvers/bsat2.hpp>

namespace mockturtle
{

struct patgen_params
{
  /*! \brief Number of initial random patterns to start with. */
  uint32_t num_random_pattern{1000};

  /*! \brief Whether to start with patterns in a file. */
  std::optional<std::string> patfile{};

  /*! \brief Whether to substitute constant nodes. */
  bool substitute_const{true};

  /*! \brief Number of patterns each node should have for both values. */
  uint32_t num_stuck_at{1};

  /*! \brief Fanout levels to consider for observability. -1 = no limit. */
  int observability_levels{-1};

  /*! \brief Whether to check and re-generate type 1 observable patterns. */
  bool observability_type1{false};

  /*! \brief Whether to check and re-generate type 2 observable patterns. */
  bool observability_type2{false};

  /*! \brief Whether to apply the distinguishing-node strategy. */
  bool distinguish_nodes{false};

  /*! \brief Whether to save generated patterns into file. */
  std::optional<std::string> write_pats{};

  /*! \brief Show progress. */
  bool progress{false};

  /*! \brief Be verbose. */
  bool verbose{false};

  /*! \brief Random seed. */
  std::default_random_engine::result_type random_seed{0};

  uint32_t conflict_limit{1000};
};

struct patgen_stats
{
  stopwatch<>::duration time_total{0};

  /*! \brief Time for simulations. */
  stopwatch<>::duration time_sim{0};

  /*! \brief Time for SAT solving. */
  stopwatch<>::duration time_sat{0};

  /* time for ODC computation */
  stopwatch<>::duration time_odc{0};

  /*! \brief Number of constant nodes found. */
  uint32_t num_constant{0};

  /*! \brief Number of total generated patterns (including random ones). */
  uint32_t num_total_patterns{0};

  /*! \brief Number of type1 unobservable nodes. */
  uint32_t unobservable_type1{0};

  /*! \brief Number of resolved type1 unobservable nodes. */
  uint32_t unobservable_type1_resolved{0};

  /*! \brief Number of type2 unobservable nodes. */
  uint32_t unobservable_type2{0};

  /*! \brief Number of resolved type2 unobservable nodes. */
  uint32_t unobservable_type2_resolved{0};

  /*! \brief Number of div0 distinguishing-node strategy generated patterns. */
  uint32_t num_div0_pats{0};

  void report() const
  {
    //std::cout << "[i] kernel: default_resub_functor\n";
    //std::cout << fmt::format( "[i]     constant-resub {:6d}                                   ({:>5.2f} secs)\n",
    //                          num_const_accepts, to_seconds( time_resubC ) );
    //std::cout << fmt::format( "[i]            0-resub {:6d}                                   ({:>5.2f} secs)\n",
    //                          num_div0_accepts, to_seconds( time_resub0 ) );
    //std::cout << fmt::format( "[i]            total   {:6d}\n",
    //                          (num_const_accepts + num_div0_accepts) );
  }
};

namespace detail
{

template<class Ntk>
class patgen_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using TT = unordered_node_map<kitty::partial_truth_table, Ntk>;

  /* without providing filename of random patterns */
  explicit patgen_impl( Ntk& ntk, patgen_params const& ps, patgen_stats& st )
    : ntk( ntk ), ps( ps ), st( st ), literals( node_literals( ntk ) ), 
      tts( ntk ), random( ps.random_seed ), sim( ntk.num_pis(), ps.num_random_pattern, ps.random_seed )
  {
    st.num_total_patterns = ps.num_random_pattern;
  }

  /* provide filename of (fixed) random patterns */
  explicit patgen_impl( Ntk& ntk, std::string const& patfile, patgen_params const& ps, patgen_stats& st )
    : ntk( ntk ), ps( ps ), st( st ), literals( node_literals( ntk ) ), 
      tts( ntk ), random( ps.random_seed ), sim( patfile, ps.num_random_pattern )
  {
    st.num_total_patterns = ps.num_random_pattern;
  }

  void run()
  {
    stopwatch t( st.time_total );

    /* start the managers */
    //progress_bar pbar{ntk.size(), "resub |{0}| node = {1:>4}   cand = {2:>4}   est. gain = {3:>5}", ps.progress};

    generate_cnf<Ntk>( ntk, [&]( auto const& clause ) {
      solver.add_clause( clause );
    }, literals );
    max_id = ntk.size() - 1;

    call_with_stopwatch( st.time_sim, [&]() {
      simulate_nodes<Ntk>( ntk, tts, sim );
    });

    if ( ps.num_stuck_at > 0 )
    {
      stuck_at_check();
      
      if ( ps.substitute_const )
      {
        for ( auto n : const_nodes )
        {
          if ( !ntk.is_dead( ntk.get_node( n ) ) )
            ntk.substitute_node( ntk.get_node( n ), ntk.get_constant( ntk.is_complemented( n ) ) );
        }
      }
    }

    if ( ps.observability_type2 )
    {
      observability_check();
    }

    if ( ps.distinguish_nodes )
    {
      distinguish_div0();
    }
  }

private:
  uint32_t new_lit()
  {
    solver.add_var();
    return make_lit( ++max_id );
  }

  void add_clauses_for_gate( node const& n, unordered_node_map<uint32_t, Ntk> const& lits, bool inverse = false )
  {
    /* currently only consider AIG */
    std::vector<uint32_t> l_fi;
    const auto c = lit_not_cond( lits[n], inverse );
    ntk.foreach_fanin( n, [&]( auto const& fi ){
      l_fi.emplace_back( lit_not_cond( lits.has( ntk.get_node(fi) )? lits[fi]: literals[fi], ntk.is_complemented( fi ) ) );
    });

    assert( ntk.is_and( n ) );
    /* AND gate */
    solver.add_clause( {l_fi[0], lit_not( c )} );
    solver.add_clause( {l_fi[1], lit_not( c )} );
    solver.add_clause( {lit_not( l_fi[0] ), lit_not( l_fi[1] ), c} );
  }

  void duplicate_fanout_cone_rec( node const& n, unordered_node_map<uint32_t, Ntk> const& lits, int level )
  {
    ntk.foreach_fanout( n, [&]( auto const& fo ){
      if ( ntk.visited( fo ) == ntk.trav_id() ) return true; /* skip */
      ntk.set_visited( fo, ntk.trav_id() );

      add_clauses_for_gate( fo, lits );

      if ( level == ps.observability_levels ) return true;

      duplicate_fanout_cone_rec( fo, lits, level+1 );
      return true; /* next */
    });
  }

  void make_lit_fanout_cone_rec( node const& n, unordered_node_map<uint32_t, Ntk>& lits, std::vector<uint32_t>& miter, int level )
  {
    ntk.foreach_fanout( n, [&]( auto const& fo ){
      if ( ntk.visited( fo ) == ntk.trav_id() ) return true; /* skip */
      ntk.set_visited( fo, ntk.trav_id() );

      lits[fo] = new_lit();

      if ( level == ps.observability_levels )
      {
        add_miter_clauses( ntk.make_signal( fo ), lits, miter );
        return true;
      }
      
      make_lit_fanout_cone_rec( fo, lits, miter, level+1 );
      return true; /* next */
    });
  }

  void add_miter_clauses( signal const& f, unordered_node_map<uint32_t, Ntk> const& lits, std::vector<uint32_t>& miter )
  {
    const auto a = lit_not_cond( literals[f], ntk.is_complemented( f ) );
    const auto b = lit_not_cond( lits[f], ntk.is_complemented( f ) );
    const auto c = new_lit();
    miter.emplace_back( c );

    solver.add_clause( {lit_not( a ), lit_not( b ), lit_not( c )} );
    solver.add_clause( {lit_not( a ), b, c} );
    solver.add_clause( {a, lit_not( b ), c} );
    solver.add_clause( {a, b, lit_not( c )} );
  }

  bool generate_observable_pattern( node const& n, bool value, std::vector<bool>& pattern )
  {
    solver.bookmark();

    /* for the original circuit, assert original n to be value */
    std::vector<pabc::lit> assumptions( 1 );
    assumptions[0] = lit_not_cond( literals[n], !value );

    /* literals for the duplicated part */
    unordered_node_map<uint32_t, Ntk> lits( ntk );
    std::vector<uint32_t> miter;

    ntk.foreach_fanin( n, [&]( auto const& fi ){
      lits[fi] = literals[fi];
    });
    lits[n] = new_lit();
    ntk.incr_trav_id();
    make_lit_fanout_cone_rec( n, lits, miter, 1 );

    /* the inverse version of n */
    add_clauses_for_gate( n, lits, true );

    /* the other copy of fanout cone */
    ntk.incr_trav_id();
    duplicate_fanout_cone_rec( n, lits, 1 );

    /* miter for POs */
    ntk.foreach_po( [&]( auto const& f ){
      if ( !lits.has( ntk.get_node( f ) ) ) return true; /* PO not in TFO, skip */

      add_miter_clauses( f, lits, miter );

      return true; /* next */
    });
    //for (auto i : miter) std::cout<<i<<" "; std::cout<<"\n";
    assert( miter.size() > 0 );
    const auto nlit = new_lit();
    miter.emplace_back( nlit );
    assumptions.emplace_back( lit_not( nlit ) );
    solver.add_clause( miter );

    /* solve and get answer */
    const auto res = call_with_stopwatch( st.time_sat, [&]() {
      return solver.solve( &assumptions[0], &assumptions[0] + 2, ps.conflict_limit );
    });
    if ( res == percy::synth_result::success )
    {
      pattern.clear();
      for ( auto j = 1u; j <= ntk.num_pis(); ++j )
        pattern.push_back(solver.var_value( j ));
    }
    else if ( res == percy::synth_result::failure )
    {
      //std::cout<<"UNSAT: node "<<unsigned(n)<<" is un-testable at value "<<value<<".\n";
      if ( false )//( ps.substitute_const )
      {
        ++st.num_constant;
        auto g = ntk.get_constant( !value );
        ntk.substitute_node( n, g );
        /* re-simulate */
        call_with_stopwatch( st.time_sim, [&]() {
          simulate_nodes<Ntk>( ntk, tts, sim );
        });
      }
    }
    else
    {
      //std::cout<<"solver timeout\n";
    }

    solver.rollback();
    return ( res == percy::synth_result::success );
  }

  void set_random_polarity()
  {
    std::vector<uint32_t> pols;
    for ( auto i = 1u; i <= ntk.num_pis(); ++i )
      if ( random() % 2 )
        pols.push_back( i );
    solver.set_polarity( pols );
  }

  void generate_more( node const& n, bool value, std::vector<std::vector<bool>> const& patterns )
  {
    solver.bookmark();
    std::vector<pabc::lit> assumptions( 1 );
    assumptions[0] = lit_not_cond( literals[n], !value );

    /* add blocking clauses */
    for ( auto i = 0u; i < patterns.size(); ++i )
    {
      auto const& pattern = patterns[i];
      std::vector<uint32_t> clause;
      for ( auto j = 0u; j < pattern.size(); ++j )
        clause.emplace_back( lit_not_cond( make_lit( j + 1 ), pattern[j] ) );
      solver.add_clause( clause );
    }

    auto num_generated = patterns.size();
    while ( num_generated < ps.num_stuck_at )
    {
      set_random_polarity();
      const auto res = call_with_stopwatch( st.time_sat, [&]() {
        return solver.solve( &assumptions[0], &assumptions[0] + 1, ps.conflict_limit );
      });
      if ( res == percy::synth_result::success )
      {
        std::vector<bool> pattern;
        for ( auto j = 1u; j <= ntk.num_pis(); ++j )
          pattern.push_back(solver.var_value( j ));

        if ( false ) //( ps.observability_type1 )
        {
          /* check if the found pattern is observable */ 
          bool observable = call_with_stopwatch( st.time_odc, [&]() { 
              return pattern_is_observable( ntk, n, pattern, ps.observability_levels );
            });
          if ( !observable )
          {
            //std::cout << "generated pattern is not observable! (" << unsigned(n) << " = " << value <<")" << std::endl;
            ++st.unobservable_type1;
            if ( generate_observable_pattern( n, value, pattern ) )
            {
              ++st.unobservable_type1_resolved;
              //std::cout << "after re-gen, now?? " << pattern_is_observable( ntk, n, pattern, ps.observability_levels ) <<"\n";
            }
          }
        }

        sim.add_pattern(pattern);
        ++st.num_total_patterns;
        ++num_generated;

        /* add blocking clauses */
        std::vector<uint32_t> clause;
        for ( auto j = 0u; j < pattern.size(); ++j )
          clause.emplace_back( lit_not_cond( make_lit( j + 1 ), pattern[j] ) );
        solver.add_clause( clause );
      }
      else break; /* can not generate more */
    }

    solver.rollback();
  }

private:
  void stuck_at_check()
  {
    std::vector<pabc::lit> assumptions( 1 );
    kitty::partial_truth_table zero = sim.compute_constant(false);

    ntk.foreach_gate( [&]( auto const& n ) 
    {
      //std::cout<<"processing node "<<unsigned(n)<<std::endl;
      if ( (tts[n] == zero) || (tts[n] == ~zero) )
      {
        bool value = !(tts[n] == ~zero); /* wanted value of n */

        assumptions[0] = lit_not_cond( literals[n], !value );
      
        set_random_polarity();
        const auto res = call_with_stopwatch( st.time_sat, [&]() {
          return solver.solve( &assumptions[0], &assumptions[0] + 1, ps.conflict_limit );
        });
        
        if ( res == percy::synth_result::success )
        {
          //std::cout << "SAT: add pattern. (" << n << " = " << value <<")" << std::endl;
          std::vector<bool> pattern;
          for ( auto j = 1u; j <= ntk.num_pis(); ++j )
            pattern.push_back(solver.var_value( j ));

          if ( ps.observability_type1 )
          {
            /* check if the found pattern is observable */ 
            bool observable = call_with_stopwatch( st.time_odc, [&]() { 
                return pattern_is_observable( ntk, n, pattern, ps.observability_levels );
              });
            if ( !observable )
            {
              //std::cout << "generated pattern is not observable! (" << unsigned(n) << " = " << value <<")" << std::endl;
              ++st.unobservable_type1;
              if ( generate_observable_pattern( n, value, pattern ) )
              {
                ++st.unobservable_type1_resolved;
                //std::cout << "after re-gen, now?? " << pattern_is_observable( ntk, n, pattern, ps.observability_levels ) <<"\n";
              }
            }
          }

          sim.add_pattern(pattern);
          ++st.num_total_patterns;

          if ( ps.num_stuck_at > 1 )
          {
            std::vector<std::vector<bool>> patterns;
            patterns.emplace_back( pattern );
            generate_more( n, value, patterns );
          }

          /* re-simulate */
          call_with_stopwatch( st.time_sim, [&]() {
            simulate_nodes<Ntk>( ntk, tts, sim );
            zero = sim.compute_constant(false);
          });
        }
        else if ( res == percy::synth_result::failure )
        {
          //std::cout << "UNSAT: this is a constant node. (" << n << ")" << std::endl;
          ++st.num_constant;
          const_nodes.emplace_back( value ? ntk.make_signal( n ) : !ntk.make_signal( n ) );
          return true; /* next gate */
        }
      }
      else if ( ps.num_stuck_at > 1 )
      {
        auto const& tt = tts[n];
        if ( kitty::count_ones( tt ) < ps.num_stuck_at )
        {
          /* collect the one patterns */
          std::vector<std::vector<bool>> patterns;
          for ( auto i = 0u; i < tt.num_bits(); ++i )
          {
            if ( kitty::get_bit( tt, i ) )
            {
              patterns.emplace_back();
              ntk.foreach_pi( [&]( auto const& pi ){ 
                patterns.back().emplace_back( kitty::get_bit( tts[pi], i ) );
              });
            }
          }

          generate_more( n, 1, patterns );

          /* re-simulate */
          call_with_stopwatch( st.time_sim, [&]() {
            simulate_nodes<Ntk>( ntk, tts, sim );
            zero = sim.compute_constant(false);
          });
        }
        else if ( kitty::count_zeros( tt ) < ps.num_stuck_at )
        {
          /* collect the zero patterns */
          std::vector<std::vector<bool>> patterns;
          for ( auto i = 0u; i < tt.num_bits(); ++i )
          {
            if ( !kitty::get_bit( tt, i ) )
            {
              patterns.emplace_back();
              ntk.foreach_pi( [&]( auto const& pi ){ 
                patterns.back().emplace_back( kitty::get_bit( tts[pi], i ) );
              });
            }
          }

          generate_more( n, 0, patterns );

          /* re-simulate */
          call_with_stopwatch( st.time_sim, [&]() {
            simulate_nodes<Ntk>( ntk, tts, sim );
            zero = sim.compute_constant(false);
          });
        }
      }

      return true; /* next gate */
    } );
  }

  void observability_check()
  {
    ntk.foreach_gate( [&]( auto const& n ) 
    {
      /* compute ODC */
      auto odc = call_with_stopwatch( st.time_odc, [&]() { 
          return observability_dont_cares<Ntk>( ntk, n, sim, tts, ps.observability_levels );
        });

      /* check if under non-ODCs n is always the same value */ 
      if ( ( tts[n] & ~odc ) == sim.compute_constant( false ) )
      {
        //std::cout << "under all observable patterns node "<< unsigned(n)<<" is always 0!" << std::endl;
        ++st.unobservable_type2;

        std::vector<bool> pattern(ntk.num_pis());
        if ( generate_observable_pattern( n, true, pattern ) )
        {
          sim.add_pattern(pattern);
          ++st.num_total_patterns;
          ++st.unobservable_type2_resolved;

          /* re-simulate */
          call_with_stopwatch( st.time_sim, [&]() {
            simulate_nodes<Ntk>( ntk, tts, sim );
          });

          //auto odc2 = call_with_stopwatch( st.time_odc, [&]() { return observability_dont_cares<Ntk>( ntk, n, sim, tts, ps.observability_levels ); });
          //std::cout<<"adding generated pattern, now? "<<( ( tts[n] & ~odc2 ) != sim.compute_constant( false ) )<<"\n";
        }
      }
      else if ( ( tts[n] | odc ) == sim.compute_constant( true ) )
      {
        //std::cout << "under all observable patterns node "<< unsigned(n)<<" is always 1!" << std::endl;
        ++st.unobservable_type2;

        std::vector<bool> pattern(ntk.num_pis());
        if ( generate_observable_pattern( n, false, pattern ) )
        {
          sim.add_pattern(pattern);
          ++st.num_total_patterns;
          ++st.unobservable_type2_resolved;

          /* re-simulate */
          call_with_stopwatch( st.time_sim, [&]() {
            simulate_nodes<Ntk>( ntk, tts, sim );
          });

          //auto odc2 = call_with_stopwatch( st.time_odc, [&]() { return observability_dont_cares<Ntk>( ntk, n, sim, tts, ps.observability_levels ); });
          //std::cout<<"adding generated pattern, now? "<<( ( tts[n] | odc2 ) != sim.compute_constant( true ) )<<"\n";
        }
      }

      return true; /* next gate */
    } );
  }

  void distinguish_div0()
  {
    ntk.foreach_gate( [&]( auto const& root ) 
    {
      ntk.foreach_gate( [&]( auto const& n ) 
      {
        if ( n <= root ) return true; /* only compare to 'later' nodes */

        if ( tts[root] == tts[n] || ~tts[root] == tts[n] )
        {
          auto nlit = new_lit();
          solver.add_clause( {literals[root], literals[n], nlit} );
          solver.add_clause( {literals[root], lit_not( literals[n] ), lit_not( nlit )} );
          solver.add_clause( {lit_not( literals[root] ), literals[n], lit_not( nlit )} );
          solver.add_clause( {lit_not( literals[root] ), lit_not( literals[n] ), nlit} );
          std::vector<pabc::lit> assumptions( 1, lit_not_cond( nlit, tts[root] == tts[n] ) );
        
          const auto res = call_with_stopwatch( st.time_sat, [&]() {
            return solver.solve( &assumptions[0], &assumptions[0] + 1, ps.conflict_limit );
          });
          
          if ( res == percy::synth_result::success )
          {
            std::vector<bool> pattern;
            for ( auto j = 1u; j <= ntk.num_pis(); ++j )
              pattern.push_back(solver.var_value( j ));

            sim.add_pattern(pattern);
            ++st.num_total_patterns;
            ++st.num_div0_pats;
            /* re-simulate */
            call_with_stopwatch( st.time_sim, [&]() {
              simulate_nodes<Ntk>( ntk, tts, sim );
            });
          }
        }

        return true; /* next */
      });
      return true; /* next */
    });
  }

private:
  Ntk& ntk;

  patgen_params const& ps;
  patgen_stats& st;

  node_map<uint32_t, Ntk> literals;
  percy::bsat_wrapper solver;
  uint32_t max_id;
  
  TT tts;
  std::vector<signal> const_nodes;

  std::default_random_engine random;

public:
  partial_simulator sim;
};

} /* namespace detail */

template<class Ntk>
partial_simulator pattern_generation( Ntk& ntk, patgen_params const& ps = {}, patgen_stats* pst = nullptr )
{
  /* TODO: check if basetype of ntk is aig */
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_clear_values_v<Ntk>, "Ntk does not implement the clear_values method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_make_signal_v<Ntk>, "Ntk does not implement the make_signal method" );
  static_assert( has_set_value_v<Ntk>, "Ntk does not implement the set_value method" );
  static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the has_size method" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement the has substitute_node method" );
  static_assert( has_value_v<Ntk>, "Ntk does not implement the has_value method" );
  static_assert( has_visited_v<Ntk>, "Ntk does not implement the has_visited method" );

  patgen_stats st;

  fanout_view<Ntk> fanout_view{ntk};

  if ( ps.patfile )
  {
    detail::patgen_impl p( fanout_view, *(ps.patfile), ps, st );
    p.run();

    if ( ps.write_pats )
    {
      p.sim.write_patterns( *(ps.write_pats) );
    }

    if ( ps.verbose )
      st.report();

    if ( pst )
      *pst = st;

    return p.sim;
  }
  else
  {
    detail::patgen_impl p( fanout_view, ps, st );
    p.run();

    if ( ps.write_pats )
    {
      p.sim.write_patterns( *(ps.write_pats) );
    }

    if ( ps.verbose )
      st.report();

    if ( pst )
      *pst = st;

    return p.sim;
  }
}

} /* namespace mockturtle */
