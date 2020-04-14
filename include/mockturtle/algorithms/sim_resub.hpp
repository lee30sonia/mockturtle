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
  \file sim_resub.hpp
  \brief Simulation-Guided Resubstitution

  \author Siang-Yun Lee
*/

#pragma once

#include <mockturtle/networks/aig.hpp>
#include "../utils/progress_bar.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"
#include <mockturtle/algorithms/simulation.hpp>
#include <kitty/partial_truth_table.hpp>
//#include <kitty/print.hpp>
//#include <kitty/operators.hpp>
#include "../utils/node_map.hpp"
#include "cnf.hpp"
#include "cleanup.hpp"
#include <percy/solvers/bsat2.hpp>
#include "reconv_cut2.hpp"

namespace mockturtle
{

struct simresub_params
{
  /*! \brief Maximum number of divisors to consider. */
  uint32_t max_divisors{150};

  /*! \brief Maximum number of nodes added by resubstitution. */
  uint32_t max_inserts{2};

  /*! \brief Maximum fanout of a node to be considered as root. */
  uint32_t skip_fanout_limit_for_roots{1000};

  /*! \brief Maximum fanout of a node to be considered as divisor. */
  uint32_t skip_fanout_limit_for_divisors{100};

  /*! \brief Number of candidates found before every SAT solve validation. */
  uint32_t num_solve{1};

  /*! \brief Whether to save generated patterns into file. */
  std::optional<std::string> write_pats{};

  /*! \brief Show progress. */
  bool progress{false};

  /*! \brief Be verbose. */
  bool verbose{false};

  /*! \brief Maximum number of PIs of reconvergence-driven cuts. */
  uint32_t max_pis{8};
};

struct simresub_stats
{
  stopwatch<>::duration time_total{0};

  /* time for simulations */
  stopwatch<>::duration time_sim{0};

  /* time for SAT solving */
  stopwatch<>::duration time_sat{0};

  /* time for finding substitutions */
  stopwatch<>::duration time_eval{0};

  /* time for MFFC computation */
  stopwatch<>::duration time_mffc{0};

  /* time for divisor collection */
  stopwatch<>::duration time_divs{0};

  /* time for checking implications (containment) */
  stopwatch<>::duration time_collect_unate_divisors{0};

  /* time for executing the user-specified callback when candidates found */
  stopwatch<>::duration time_callback{0};

  /* time & number of r == d (equal) node substitutions */
  stopwatch<>::duration time_resub0{0};
  uint32_t num_div0_accepts{0};

  /* time & number of r == d1 &| d2 node substitutions */
  stopwatch<>::duration time_resub1{0};
  uint32_t num_div1_accepts{0};

  /*! \brief Initial network size (before resubstitution) */
  uint64_t initial_size{0};

  uint32_t num_cex{0};

  /*! \brief Total number of gain  */
  uint64_t estimated_gain{0};

  uint64_t num_total_divisors{0};

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

template<typename Ntk>
bool substitue_fn( Ntk& ntk, typename Ntk::node const& n, typename Ntk::signal const& g )
{
  ntk.substitute_node( n, g );
  //std::cout<<"substitute node "<<unsigned(n)<<" with node "<<unsigned(ntk.get_node(g))<<std::endl;
  return true;
};

template<typename Ntk>
bool undo_substitue_fn( Ntk& ntk, typename Ntk::node const& n, typename Ntk::signal const& g )
{
  /* create the deleted node n */
  std::vector<typename Ntk::signal> fanins;
  ntk.foreach_fanin( n, [&]( const auto& f ){
    fanins.emplace_back( f );
  });
  auto const new_n = ntk.create_and( fanins[0], fanins[1] );

  ntk.substitute_node( ntk.get_node( g ), ntk.is_complemented( g ) ? !new_n : new_n );
  //std::cout<<"undo: substitute node "<<unsigned(n)<<" with node "<<unsigned(ntk.get_node(g))<<std::endl;
  return true;
};

template<typename Ntk>
bool report_fn( Ntk& ntk, typename Ntk::node const& n, typename Ntk::signal const& g )
{
  //(void)ntk; (void)n; (void)g;
  std::cout<<"substitute node "<<unsigned(n)<<" with node "<<unsigned(ntk.get_node(g))<<std::endl;
  return false;
};

/* based on abcRefs.c */
template<typename Ntk>
class node_mffc_inside
{
public:
  using node = typename Ntk::node;

public:
  explicit node_mffc_inside( Ntk const& ntk )
    : ntk( ntk )
  {
  }

  int32_t run( node const& n, std::vector<node> const& leaves, std::vector<node>& inside )
  {
    /* increment the fanout counters for the leaves */
    for ( const auto& l : leaves )
      ntk.incr_fanout_size( l );

    /* dereference the node */
    auto count1 = node_deref_rec( n );

    /* collect the nodes inside the MFFC */
    node_mffc_cone( n, inside );

    /* reference it back */
    auto count2 = node_ref_rec( n );
    (void)count2;
    assert( count1 == count2 );

    for ( const auto& l : leaves )
      ntk.decr_fanout_size( l );

    return count1;
  }

private:
  /* ! \brief Dereference the node's MFFC */
  int32_t node_deref_rec( node const& n )
  {
    if ( ntk.is_pi( n ) )
      return 0;

    int32_t counter = 1;
    ntk.foreach_fanin( n, [&]( const auto& f ){
        auto const& p = ntk.get_node( f );

        ntk.decr_fanout_size( p );
        if ( ntk.fanout_size( p ) == 0 )
          counter += node_deref_rec( p );
      });

    return counter;
  }

  /* ! \brief Reference the node's MFFC */
  int32_t node_ref_rec( node const& n )
  {
    if ( ntk.is_pi( n ) )
      return 0;

    int32_t counter = 1;
    ntk.foreach_fanin( n, [&]( const auto& f ){
        auto const& p = ntk.get_node( f );

        auto v = ntk.fanout_size( p );
        ntk.incr_fanout_size( p );
        if ( v == 0 )
          counter += node_ref_rec( p );
      });

    return counter;
  }

  void node_mffc_cone_rec( node const& n, std::vector<node>& cone, bool top_most )
  {
    /* skip visited nodes */
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;
    ntk.set_visited( n, ntk.trav_id() );

    if ( !top_most && ( ntk.is_pi( n ) || ntk.fanout_size( n ) > 0 ) )
      return;

#if 0
    if ( !top_most && ntk.is_pi( n ) )
      return;
    if ( !top_most ) /* if all the fanouts are dangling, it is still in MFFC */
    {
      bool all_dangling = true;
      ntk.foreach_fanout( n, [&]( const auto& fo ){
        if ( ntk.fanout_size( fo ) > 0 )
        {
          all_dangling = false;
          return false; /* break loop */
        }
        return true;
      });
      if ( !all_dangling )
        return;
    }
#endif

    /* recurse on children */
    ntk.foreach_fanin( n, [&]( const auto& f ){
        node_mffc_cone_rec( ntk.get_node( f ), cone, false );
      });

    /* collect the internal nodes */
    cone.emplace_back( n );
  }

  void node_mffc_cone( node const& n, std::vector<node>& cone )
  {
    cone.clear();
    ntk.incr_trav_id();
    node_mffc_cone_rec( n, cone, true );
  }

private:
  Ntk const& ntk;
};

template<class NtkBase, class Ntk>
class simresub_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using TT = unordered_node_map<kitty::partial_truth_table, Ntk>;
  using resub_callback_t = std::function<bool( NtkBase&, node const&, signal const& )>;

  struct unate_divisors
  {
    using signal = typename Ntk::signal;

    std::vector<signal> positive_divisors;
    std::vector<signal> negative_divisors;

    void clear()
    {
      positive_divisors.clear();
      negative_divisors.clear();
    }
  };

  explicit simresub_impl( NtkBase& ntkbase, Ntk& ntk, simresub_params const& ps, simresub_stats& st, partial_simulator& sim, resub_callback_t const& callback = substitue_fn<NtkBase>, resub_callback_t const& undo_callback = undo_substitue_fn<NtkBase> )
    : ntkbase( ntkbase ), ntk( ntk ), ps( ps ), st( st ), callback( callback ), undo_callback( undo_callback ),
      tts( ntk ), sim( sim ), literals( node_literals( ntkbase ) )
  {
    st.initial_size = ntk.num_gates(); 

    auto const update_level_of_new_node = [&]( const auto& n ){
      ntk.resize_levels();
      update_node_level( n );
    };
    
    auto const update_level_of_existing_node = [&]( node const& n, const auto& old_children ){
      //std::cout<<unsigned(n)<<" is modified"<<std::endl;
      (void)old_children;
      update_node_level( n );
    };
    
    auto const update_level_of_deleted_node = [&]( const auto& n ){
      /* update fanout */
      //std::cout<<unsigned(n)<<" is deleted"<<std::endl;
      ntk.set_level( n, -1 );
    };
    
    ntk._events->on_add.emplace_back( update_level_of_new_node );
    
    ntk._events->on_modified.emplace_back( update_level_of_existing_node );
    
    ntk._events->on_delete.emplace_back( update_level_of_deleted_node );
  }

  void run()
  {
    stopwatch t( st.time_total );

    /* start the managers */
    progress_bar pbar{ntk.size(), "resub |{0}| node = {1:>4}   cand = {2:>4}   est. gain = {3:>5}", ps.progress};

    generate_cnf<NtkBase>( ntkbase, [&]( auto const& clause ) {
      solver.add_clause( clause );
    }, literals );

    //std::vector<node> PIs( ntk.num_pis() );
    //ntk.foreach_pi( [&]( auto const& n, auto i ){ PIs.at(i) = n; });

    call_with_stopwatch( st.time_sim, [&]() {
      simulate_nodes<Ntk>( ntk, tts, sim );
    });

    /* iterate through all nodes and try to replace it */
    auto const size = ntk.num_gates();
    ntk.foreach_gate( [&]( auto const& n, auto i ){
        if ( i >= size )
          return false; /* terminate */

        pbar( i, i, num_candidates, st.estimated_gain );

        /* skip nodes with many fanouts */
        if ( ntk.fanout_size( n ) > ps.skip_fanout_limit_for_roots )
          return true; /* next */

        /* use all the PIs as the cut */
        //auto const leaves = PIs;
        cut_manager<Ntk> mgr( ps.max_pis );
        auto const leaves = reconv_driven_cut( mgr, ntk, n );
        
        /* evaluate this cut */
        auto const g = call_with_stopwatch( st.time_eval, [&]() {
            return evaluate( n, leaves );
          });
        if ( !g ) return true; /* next */
        
        /* update progress bar */
        num_candidates++;
        st.estimated_gain += last_gain;
        
        return true; /* next */
      });

    if ( ps.num_solve > 1 && candidates.size() > 0 )
      validate();
  }

private:
  void update_node_level( node const& n, bool top_most = true )
  {
    uint32_t curr_level = ntk.level( n );

    uint32_t max_level = 0;
    ntk.foreach_fanin( n, [&]( const auto& f ){
        auto const p = ntk.get_node( f );
        auto const fanin_level = ntk.level( p );
        if ( fanin_level > max_level )
        {
          max_level = fanin_level;
        }
      });
    ++max_level;

    if ( curr_level != max_level )
    {
      ntk.set_level( n, max_level );

      /* update only one more level */
      if ( top_most )
      {
        ntk.foreach_fanout( n, [&]( const auto& p ){
            update_node_level( p, false );
          });
      }
    }
  }

  void collect_divisors_rec( node const& n, std::vector<node>& internal )
  {
    /* skip visited nodes */
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;
    ntk.set_visited( n, ntk.trav_id() );

    ntk.foreach_fanin( n, [&]( const auto& f ){
        collect_divisors_rec( ntk.get_node( f ), internal );
      });

    /* collect the internal nodes */
    if ( ntk.value( n ) == 0 && n != 0 ) /* ntk.fanout_size( n ) */
      internal.emplace_back( n );
  }

  bool collect_divisors( node const& root, std::vector<node> const& leaves, uint32_t required )
  {
    /* add the leaves of the cuts to the divisors */
    divs.clear();

    ntk.incr_trav_id();
    for ( const auto& l : leaves )
    {
      divs.emplace_back( l );
      ntk.set_visited( l, ntk.trav_id() );
    }

    /* mark nodes in the MFFC */
    for ( const auto& t : temp )
      ntk.set_value( t, 1 );

    /* collect the cone (without MFFC) */
    collect_divisors_rec( root, divs );

    /* unmark the current MFFC */
    for ( const auto& t : temp )
      ntk.set_value( t, 0 );

    /* check if the number of divisors is not exceeded */
    if ( divs.size() - leaves.size() + temp.size() >= ps.max_divisors - ps.max_pis )
      return false;

    /* get the number of divisors to collect */
    int32_t limit = ps.max_divisors - ps.max_pis - ( uint32_t( divs.size() ) + 1 - uint32_t( leaves.size() ) + uint32_t( temp.size() ) );

    /* explore the fanouts, which are not in the MFFC */
    int32_t counter = 0;
    bool quit = false;

    /* NOTE: this is tricky and cannot be converted to a range-based loop */
    auto size = divs.size();
    for ( auto i = 0u; i < size; ++i )
    {
      auto const d = divs.at( i );

      if ( ntk.fanout_size( d ) > ps.skip_fanout_limit_for_divisors )
        continue;

      /* if the fanout has all fanins in the set, add it */
      ntk.foreach_fanout( d, [&]( node const& p ){
          if ( ntk.fanout_size( p ) == 0 )
            return true; /* don't add dangling nodes */
          if ( ntk.visited( p ) == ntk.trav_id() || ntk.level( p ) > required )
            return true; /* next fanout */

          bool all_fanins_visited = true;
          ntk.foreach_fanin( p, [&]( const auto& g ){
              if ( ntk.visited( ntk.get_node( g ) ) != ntk.trav_id() )
              {
                all_fanins_visited = false;
                return false; /* terminate fanin-loop */
              }
              return true; /* next fanin */
            });

          if ( !all_fanins_visited )
            return true; /* next fanout */

          bool has_root_as_child = false;
          ntk.foreach_fanin( p, [&]( const auto& g ){
              if ( ntk.get_node( g ) == root )
              {
                has_root_as_child = true;
                return false; /* terminate fanin-loop */
              }
              return true; /* next fanin */
            });

          if ( has_root_as_child )
            return true; /* next fanout */

          divs.emplace_back( p );
          ++size;
          ntk.set_visited( p, ntk.trav_id() );

          /* quit computing divisors if there are too many of them */
          if ( ++counter == limit )
          {
            quit = true;
            return false; /* terminate fanout-loop */
          }

          return true; /* next fanout */
        });

      if ( quit )
        break;
    }

    /* get the number of divisors */
    num_divs = uint32_t( divs.size() );

    /* add the nodes in the MFFC */
    for ( const auto& t : temp )
    {
      divs.emplace_back( t );
    }

    assert( root == divs.at( divs.size()-1u ) );
    assert( divs.size() - leaves.size() <= ps.max_divisors - ps.max_pis );

    return true;
  }

  std::optional<signal> evaluate( node const& root, std::vector<node> const &leaves )
  {
    uint32_t const required = std::numeric_limits<uint32_t>::max();

    last_gain = 0;

    /* collect the MFFC */
    int32_t num_mffc = call_with_stopwatch( st.time_mffc, [&]() {
        node_mffc_inside collector( ntk );
        auto num_mffc = collector.run( root, leaves, temp );
        assert( num_mffc > 0 );
        return num_mffc;
      });

    /* collect the divisor nodes in the cut */
    bool div_comp_success = call_with_stopwatch( st.time_divs, [&]() {
        return collect_divisors( root, leaves, required );
      });

    if ( !div_comp_success )
      return std::nullopt;
    
    /* update statistics */
    st.num_total_divisors += num_divs;

    /* consider equal nodes */
    auto g = call_with_stopwatch( st.time_resub0, [&]() {
        return resub_div0( root, required );
      } );
    if ( g )
    {
      ++st.num_div0_accepts;
      last_gain = num_mffc;
      return g; /* accepted resub */
    }

    if ( ps.max_inserts == 0 || num_mffc <= 1 )
    {
      return std::nullopt;
    }

    /* collect level one divisors */
    call_with_stopwatch( st.time_collect_unate_divisors, [&]() {
        collect_unate_divisors( root, required );
      });

    /* consider equal nodes */
    g = call_with_stopwatch( st.time_resub1, [&]() {
        return resub_div1( root, required );
      } );
    if ( g )
    {
      ++st.num_div1_accepts;
      last_gain = num_mffc - 1;
      return g; /* accepted resub */
    }

    if ( ps.max_inserts == 1 || num_mffc == 2 )
    {
      return std::nullopt;
    }

    return std::nullopt;
  }

  void found_cex()
  {
    std::vector<bool> pattern;
    for ( auto j = 1u; j <= ntk.num_pis(); ++j )
      pattern.push_back(solver.var_value( j ));
    sim.add_pattern(pattern);
    ++st.num_cex;

    /* re-simulate */
    call_with_stopwatch( st.time_sim, [&]() {
      simulate_nodes<Ntk>( ntk, tts, sim );
    });
  }

  std::optional<signal> resub_div0( node const& root, uint32_t required ) 
  {
    (void)required;
    auto const& tt = tts[root];

    //for ( auto i = 0u; i < num_divs; ++i )
    for ( int i = num_divs-1; i >= 0; --i )
    {
      auto const d = divs.at( i );
      if ( tt == tts[d] || ~tt == tts[d] )
      {
        bool phase = ( ~tt == tts[d] );

        candidates.emplace_back( std::make_pair( root, phase ? !ntk.make_signal( d ) : ntk.make_signal( d ) ) );

        solver.add_var();
        auto nlit = make_lit( solver.nr_vars()-1 );
        solver.add_clause( {literals[root], literals[d], nlit} );
        solver.add_clause( {literals[root], lit_not( literals[d] ), lit_not( nlit )} );
        solver.add_clause( {lit_not( literals[root] ), literals[d], lit_not( nlit )} );
        solver.add_clause( {lit_not( literals[root] ), lit_not( literals[d] ), nlit} );
        lits_neq.emplace_back( lit_not_cond( nlit, !phase ) );
        
        if ( ps.num_solve == 1 )
        {
          if ( validate_one() )
            return phase ? !ntk.make_signal( d ) : ntk.make_signal( d );
        }
        else 
        {
          /* update network */
          call_with_stopwatch( st.time_callback, [&]() {
            callback( ntkbase, root, phase ? !ntk.make_signal( d ) : ntk.make_signal( d ) );
          });

          if ( lits_neq.size() == ps.num_solve )
            validate();
          return phase ? !ntk.make_signal( d ) : ntk.make_signal( d );
        }
      }
    }

    return std::nullopt;
  }

  bool validate_one()
  {
    std::vector<pabc::lit> assumptions( 1 );
    assumptions[0] = lits_neq[0];
    const auto res = call_with_stopwatch( st.time_sat, [&]() {
      return solver.solve( &assumptions[0], &assumptions[0] + 1, 0 );
    });

    lits_neq.clear();

    if ( res == percy::synth_result::success ) /* CEX found */
    {
      //std::cout<<"cex found for substituting node "<<unsigned(candidates[0].first)<<" with node "<<unsigned(ntk.get_node(candidates[0].second))<<std::endl;
      found_cex();
      candidates.clear();
      return false;
    }
    /* update network */
    call_with_stopwatch( st.time_callback, [&]() {
      callback( ntkbase, candidates[0].first, candidates[0].second );
    });
    candidates.clear();
    return true;
  }

  void validate()
  {
    if ( candidates.size() == 0u ) return;

    std::vector<uint32_t> pols;
    for ( auto i = 0u; i < lits_neq.size(); ++i )
    {
      if ( !lit_is_complemented( lits_neq.at( i ) ) )
        pols.push_back( lit2var( lits_neq.at( i ) ) );
    }
    solver.set_polarity( pols );

    solver.add_var();
    auto nlit = make_lit( solver.nr_vars()-1 );
    lits_neq.push_back( nlit );
    solver.add_clause( lits_neq );
    lits_neq.pop_back();
    std::vector<pabc::lit> assumptions( 1 );
    assumptions[0] = lit_not( nlit );

    const auto res = call_with_stopwatch( st.time_sat, [&]() {
      return solver.solve( &assumptions[0], &assumptions[0] + 1, 0 );
    });

    if ( res == percy::synth_result::success ) /* CEX found */
    {
      found_cex();
      refine();
    }
    else
    {
      candidates.clear();
      lits_neq.clear();
    }
  }

  void refine()
  {
    assert( candidates.size() == lits_neq.size() );
    auto const size_before = lits_neq.size();
    int j = 0;
    for ( int i = 0; i < int( candidates.size() ); ++i )
    {
      if ( solver.var_value( lit2var( lits_neq[j] ) ) ^ lit_is_complemented( lits_neq[j] ) )
      {
        /* undo false substitution */
        call_with_stopwatch( st.time_callback, [&]() {
          undo_callback( ntkbase, candidates[i].first, candidates[i].second );
        });
        lits_neq.erase( lits_neq.begin() + j );
        candidates.erase( candidates.begin() + i );
        --j; --i;
      }
      ++j;
    }
    assert( lits_neq.size() < size_before );
    validate();
  }

  void collect_unate_divisors( node const& root, uint32_t required )
  {
    udivs.clear();

    auto const& tt = tts[root];
    for ( auto i = 0u; i < num_divs; ++i )
    {
      auto const d = divs.at( i );

      if ( ntk.level( d ) > required - 1 )
        continue;

      auto const& tt_d = tts[d];

      /* check positive containment */
      if ( kitty::implies( tt_d, tt ) )
      {
        udivs.positive_divisors.emplace_back( ntk.make_signal( d ) );
        continue;
      }
      if ( kitty::implies( ~tt_d, tt ) )
      {
        udivs.positive_divisors.emplace_back( !ntk.make_signal( d ) );
        continue;
      }

      /* check negative containment */
      if ( kitty::implies( tt, tt_d ) )
      {
        udivs.negative_divisors.emplace_back( ntk.make_signal( d ) );
        continue;
      }
      if ( kitty::implies( tt, ~tt_d ) )
      {
        udivs.negative_divisors.emplace_back( !ntk.make_signal( d ) );
        continue;
      }
    }
  }

  std::optional<signal> resub_div1( node const& root, uint32_t required )
  {
    (void)required;
    auto const& tt = tts[root];

    /* check for positive unate divisors */
    for ( auto i = 0u; i < udivs.positive_divisors.size(); ++i )
    {
      auto const& s0 = udivs.positive_divisors.at( i );

      for ( auto j = i + 1; j < udivs.positive_divisors.size(); ++j )
      {
        auto const& s1 = udivs.positive_divisors.at( j );

        auto const& tt_s0 = ntk.is_complemented(s0)? ~(tts[s0]): tts[s0];
        auto const& tt_s1 = ntk.is_complemented(s1)? ~(tts[s1]): tts[s1];

        if ( ( tt_s0 | tt_s1 ) == tt )
        {
          auto l_r = literals[root];
          auto l_s0 = lit_not_cond( literals[ntk.get_node(s0)], ntk.is_complemented(s0));
          auto l_s1 = lit_not_cond( literals[ntk.get_node(s1)], ntk.is_complemented(s1));

          //std::cout<<"found substitution "<< unsigned(root)<<" = "<<(ntk.is_complemented(s0)?"~":"")<<unsigned(ntk.get_node(s0))<<" OR "<<(ntk.is_complemented(s1)?"~":"")<<unsigned(ntk.get_node(s1))<<"\n";
          auto g = ntk.create_or( s0, s1 );
          /* update CNF */
          literals.resize();
          solver.add_var();
          literals[ntk.get_node(g)] = make_lit( solver.nr_vars()-1 );
          auto l_g = lit_not_cond( literals[ntk.get_node(g)], ntk.is_complemented(g) );
          solver.add_clause( {lit_not( l_s0 ), l_g} );
          solver.add_clause( {lit_not( l_s1 ), l_g} );
          solver.add_clause( {l_s0, l_s1, lit_not( l_g )} );
          /* re-simulate */
          call_with_stopwatch( st.time_sim, [&]() {
              simulate_nodes<Ntk>( ntk, tts, sim );
            });
          candidates.emplace_back( std::make_pair( root, g ) );

          solver.add_var();
          auto nlit = make_lit( solver.nr_vars()-1 );
          solver.add_clause( {l_r, l_g, nlit} );
          solver.add_clause( {lit_not( l_r ), lit_not( l_g ), nlit} );
          lits_neq.emplace_back( lit_not( nlit ) );

          if ( ps.num_solve == 1 )
          {
            if ( validate_one() )
              return g;
          }
          else 
          {
            /* update network */
            call_with_stopwatch( st.time_callback, [&]() {
              callback( ntkbase, root, g );
            });

            if ( lits_neq.size() == ps.num_solve )
              validate();
            return g;
          }
        }
      }
    }

    /* check for negative unate divisors */
    for ( auto i = 0u; i < udivs.negative_divisors.size(); ++i )
    {
      auto const& s0 = udivs.negative_divisors.at( i );

      for ( auto j = i + 1; j < udivs.negative_divisors.size(); ++j )
      {
        auto const& s1 = udivs.negative_divisors.at( j );

        auto const& tt_s0 = ntk.is_complemented(s0)? ~(tts[s0]): tts[s0];
        auto const& tt_s1 = ntk.is_complemented(s1)? ~(tts[s1]): tts[s1];

        if ( ( tt_s0 & tt_s1 ) == tt )
        {
          auto l_r = literals[root];
          auto l_s0 = lit_not_cond( literals[ntk.get_node(s0)], ntk.is_complemented(s0));
          auto l_s1 = lit_not_cond( literals[ntk.get_node(s1)], ntk.is_complemented(s1));

          //std::cout<<"found substitution "<< unsigned(root)<<" = "<<(ntk.is_complemented(s0)?"~":"")<<unsigned(ntk.get_node(s0))<<" AND "<<(ntk.is_complemented(s1)?"~":"")<<unsigned(ntk.get_node(s1))<<"\n";
          auto g = ntk.create_and( s0, s1 );
          /* update CNF */
          literals.resize();
          solver.add_var();
          literals[ntk.get_node(g)] = make_lit( solver.nr_vars()-1 );
          auto l_g = lit_not_cond( literals[ntk.get_node(g)], ntk.is_complemented(g) );
          solver.add_clause( {lit_not( l_g ), l_s0} );
          solver.add_clause( {lit_not( l_g ), l_s1} );
          solver.add_clause( {lit_not( l_s0 ), lit_not( l_s1 ), l_g} );
          /* re-simulate */
          call_with_stopwatch( st.time_sim, [&]() {
              simulate_nodes<Ntk>( ntk, tts, sim );
            });
          candidates.emplace_back( std::make_pair( root, g ) );

          solver.add_var();
          auto nlit = make_lit( solver.nr_vars()-1 );
          solver.add_clause( {l_r, l_g, nlit} );
          solver.add_clause( {lit_not( l_r ), lit_not( l_g ), nlit} );
          lits_neq.emplace_back( lit_not( nlit ) );

          if ( ps.num_solve == 1 )
          {
            if ( validate_one() )
              return g;
          }
          else 
          {
            /* update network */
            call_with_stopwatch( st.time_callback, [&]() {
              callback( ntkbase, root, g );
            });

            if ( lits_neq.size() == ps.num_solve )
              validate();
            return g;
          }
        }
      }
    }

    return std::nullopt;
  }

private:
  NtkBase& ntkbase;
  Ntk& ntk;

  simresub_params const& ps;
  simresub_stats& st;

  /* callback functions to substitute and to undo wrong substitution */
  resub_callback_t const& callback;
  resub_callback_t const& undo_callback;

  /* temporary statistics for progress bar */
  uint32_t num_candidates{0};
  uint32_t last_gain{0};

  TT tts;
  partial_simulator& sim;

  node_map<uint32_t, NtkBase> literals;
  std::vector<uint32_t> lits_neq;
  std::vector<std::pair<node, signal>> candidates;
  percy::bsat_wrapper solver;

  unate_divisors udivs;

  std::vector<node> temp;
  std::vector<node> divs;
  uint32_t num_divs{0};
};

} /* namespace detail */

template<class Ntk>
void sim_resubstitution( Ntk& ntk, partial_simulator& sim, simresub_params const& ps = {}, simresub_stats* pst = nullptr )
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

  using resub_view_t = fanout_view<depth_view<Ntk>>;
  depth_view<Ntk> depth_view{ntk};
  resub_view_t resub_view{depth_view};

  simresub_stats st;

  detail::simresub_impl<Ntk, resub_view_t> p( ntk, resub_view, ps, st, sim, detail::substitue_fn<Ntk> );
  p.run();

  if ( ps.write_pats )
  {
    sim.write_patterns( *(ps.write_pats) );
  }

  if ( ps.verbose )
    st.report();

  if ( pst )
    *pst = st;
}

} /* namespace mockturtle */
