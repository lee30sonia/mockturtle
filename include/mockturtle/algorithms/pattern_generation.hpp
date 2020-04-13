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

#include <mockturtle/networks/aig.hpp>
#include "../utils/progress_bar.hpp"
#include "../utils/stopwatch.hpp"
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/dont_cares.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/print.hpp>
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

  /*! \brief Whether to substitute constant nodes. */
  bool substitute_const{true};

  /*! \brief Whether to check and re-generate type 1 observable patterns. */
  bool observability_type1{false};

  /*! \brief Whether to check and re-generate type 2 observable patterns. */
  bool observability_type2{false};

  /*! \brief Whether to save generated patterns into file. */
  std::optional<std::string> write_pats{};

  /*! \brief Show progress. */
  bool progress{false};

  /*! \brief Be verbose. */
  bool verbose{false};

  /*! \brief Random seed. */
  std::default_random_engine::result_type random_seed{0};

  uint32_t ODC_failure_limit{100};
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

  explicit patgen_impl( Ntk& ntk, patgen_params const& ps, patgen_stats& st )
    : ntk( ntk ), ps( ps ), st( st ), POs( ntk.num_pos() ), literals( node_literals( ntk ) ), 
      tts( ntk ), sim( ntk.num_pis(), ps.num_random_pattern, ps.random_seed )
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

    simulate_generate();
  }

private:
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

  void duplicate_fanout_cone_rec( node const& n, unordered_node_map<uint32_t, Ntk> const& lits )
  {
    ntk.foreach_fanout( n, [&]( auto const& fo ){
      if ( ntk.visited( fo ) == ntk.trav_id() ) return true; /* skip */
      ntk.set_visited( fo, ntk.trav_id() );

      add_clauses_for_gate( fo, lits );

      duplicate_fanout_cone_rec( fo, lits );
      return true; /* next */
    });
  }

  void make_lit_fanout_cone_rec( node const& n, unordered_node_map<uint32_t, Ntk>& lits )
  {
    ntk.foreach_fanout( n, [&]( auto const& fo ){
      if ( ntk.visited( fo ) == ntk.trav_id() ) return true; /* skip */
      ntk.set_visited( fo, ntk.trav_id() );

      solver.add_var();
      lits[fo] = make_lit( solver.nr_vars() - 1 );
      
      make_lit_fanout_cone_rec( fo, lits );
      return true; /* next */
    });
  }

  bool generate_observable_pattern( node const& n, bool value, std::vector<bool>& pattern )
  {
    solver.bookmark();

    /* for the original circuit, assert original n to be value */
    std::vector<pabc::lit> assumptions( 1 );
    assumptions[0] = lit_not_cond( literals[n], !value );

    /* literals for the duplicated part */
    unordered_node_map<uint32_t, Ntk> lits( ntk );
    ntk.foreach_fanin( n, [&]( auto const& fi ){
      lits[fi] = literals[fi];
    });
    solver.add_var();
    lits[n] = make_lit( solver.nr_vars() - 1 );
    ntk.incr_trav_id();
    make_lit_fanout_cone_rec( n, lits );

    /* the inverse version of n */
    add_clauses_for_gate( n, lits, true );

    /* the other copy of fanout cone */
    ntk.incr_trav_id();
    duplicate_fanout_cone_rec( n, lits );

    /* miter for POs */
    std::vector<uint32_t> miter;
    ntk.foreach_po( [&]( auto const& f ){
      if ( !lits.has( ntk.get_node( f ) ) ) return true; /* PO not in TFO, skip */

      const auto a = lit_not_cond( literals[f], ntk.is_complemented( f ) );
      const auto b = lit_not_cond( lits[f], ntk.is_complemented( f ) );

      solver.add_var();
      const auto c = make_lit( solver.nr_vars() - 1 );
      miter.emplace_back( c );

      solver.add_clause( {lit_not( a ), lit_not( b ), lit_not( c )} );
      solver.add_clause( {lit_not( a ), b, c} );
      solver.add_clause( {a, lit_not( b ), c} );
      solver.add_clause( {a, b, lit_not( c )} );

      return true; /* next */
    });
    solver.add_clause( miter );

    /* solve and get answer */
    const auto res = call_with_stopwatch( st.time_sat, [&]() {
      return solver.solve( &assumptions[0], &assumptions[0] + 1, 1000 );
    });
    if ( res == percy::synth_result::success )
    {
      pattern.clear();
      for ( auto j = 1u; j <= ntk.num_pis(); ++j )
        pattern.push_back(solver.var_value( j ));
    }
    else if ( res == percy::synth_result::failure ) 
      std::cout<<"UNSAT: node "<<unsigned(n)<<" is un-testable at value "<<value<<"\n";
    else
      std::cout<<"solver timeout\n";

    solver.rollback();
    return ( res == percy::synth_result::success );
  }

  void simulate_generate()
  {
    call_with_stopwatch( st.time_sim, [&]() {
      simulate_nodes<Ntk>( ntk, tts, sim );
    });

    std::vector<pabc::lit> assumptions( 1 );
    kitty::partial_truth_table zero = sim.compute_constant(false);

    ntk.foreach_gate( [&]( auto const& n ) 
    {
      //std::cout<<"processing node "<<unsigned(n)<<std::endl;
      if ( (tts[n] == zero) || (tts[n] == ~zero) )
      {
        bool value = !(tts[n] == ~zero); /* wanted value of n */

        assumptions[0] = lit_not_cond( literals[n], !value );
      
        const auto res = call_with_stopwatch( st.time_sat, [&]() {
          return solver.solve( &assumptions[0], &assumptions[0] + 1, 0 );
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
            ntk.foreach_po( [&]( auto const& f, auto i ){ POs.at(i) = ntk.get_node( f ); });
            //unordered_node_map<bool, Ntk> po_vals( ntk );
            bool observable = call_with_stopwatch( st.time_odc, [&]() { 
                return pattern_is_observable( ntk, n, pattern, POs );
              });
            if ( !observable )
            {
              std::cout << "generated pattern is not observable! (" << unsigned(n) << " = " << value <<")" << std::endl;
              if ( generate_observable_pattern( n, value, pattern ) )
                std::cout << "after re-gen, now?? " << pattern_is_observable( ntk, n, pattern, POs ) <<"\n";
            }
          }

          sim.add_pattern(pattern);
          ++st.num_total_patterns;

          /* re-simulate */
          call_with_stopwatch( st.time_sim, [&]() {
            simulate_nodes<Ntk>( ntk, tts, sim );
            zero = sim.compute_constant(false);
          });
        }
        else
        {
          //std::cout << "UNSAT: this is a constant node. (" << n << ")" << std::endl;
          ++st.num_constant;
          if ( ps.substitute_const )
          {
            auto g = ntk.get_constant( tts[n] == ~zero );
            ntk.substitute_node( n, g );
          }
          return true; /* next gate */
        }
      }

      else if ( ps.observability_type2 )
      {
        /* compute ODC */
        ntk.foreach_po( [&]( auto const& f, auto i ){ POs.at(i) = ntk.get_node( f ); });
        auto odc = call_with_stopwatch( st.time_odc, [&]() { 
            return observability_dont_cares_without_window<Ntk>( ntk, n, sim, tts, POs );
          });

        /* check if under non-ODCs n is always the same value */ 
        if ( ( tts[n] & ~odc ) == sim.compute_constant( false ) )
        {
          std::cout << "under all observable patterns node "<< unsigned(n)<<" is always 0!" << std::endl;

          std::vector<bool> pattern(ntk.num_pis());
          if ( generate_observable_pattern( n, true, pattern ) )
          {
            sim.add_pattern(pattern);
            ++st.num_total_patterns;

            /* re-simulate */
            call_with_stopwatch( st.time_sim, [&]() {
              simulate_nodes<Ntk>( ntk, tts, sim );
              zero = sim.compute_constant(false);
            });

            auto odc2 = call_with_stopwatch( st.time_odc, [&]() { 
              return observability_dont_cares_without_window<Ntk>( ntk, n, sim, tts, POs );
            });
            std::cout<<"adding generated pattern, now? "<<( ( tts[n] & ~odc2 ) != sim.compute_constant( false ) )<<"\n";
          }
        }
        else if ( ( tts[n] | odc ) == sim.compute_constant( true ) )
        {
          std::cout << "under all observable patterns node "<< unsigned(n)<<" is always 1!" << std::endl;

          std::vector<bool> pattern(ntk.num_pis());
          if ( generate_observable_pattern( n, false, pattern ) )
          {
            sim.add_pattern(pattern);
            ++st.num_total_patterns;

            /* re-simulate */
            call_with_stopwatch( st.time_sim, [&]() {
              simulate_nodes<Ntk>( ntk, tts, sim );
              zero = sim.compute_constant(false);
            });

            auto odc2 = call_with_stopwatch( st.time_odc, [&]() { 
              return observability_dont_cares_without_window<Ntk>( ntk, n, sim, tts, POs );
            });
            std::cout<<"adding generated pattern, now? "<<( ( tts[n] | odc2 ) != sim.compute_constant( true ) )<<"\n";
          }
        }
      }

      return true; /* next gate */
    } );
  }

private:
  Ntk& ntk;

  patgen_params const& ps;
  patgen_stats& st;

  std::vector<node> POs;

  node_map<uint32_t, Ntk> literals;
  percy::bsat_wrapper solver;
  
  TT tts;

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

} /* namespace mockturtle */
