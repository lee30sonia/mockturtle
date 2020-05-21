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
#include <mockturtle/networks/xag.hpp>
#include "../utils/progress_bar.hpp"
#include "../utils/stopwatch.hpp"
#include "../utils/abc_resub.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/circuit_validator.hpp>
#include <kitty/partial_truth_table.hpp>
#include <bill/sat/interface/abc_bsat2.hpp>
#include "cleanup.hpp"
#include "reconv_cut2.hpp"
#include <bill/sat/interface/z3.hpp>
#include <bill/sat/interface/abc_bsat2.hpp>

namespace mockturtle
{

struct simresub_params
{
  /*! \brief Maximum number of divisors to consider. */
  uint32_t max_divisors{150};

  /*! \brief Maximum number of divisors to consider in computing k-resub. */
  uint32_t max_divisors_k{50};

  /*! \brief Maximum number of trials to call the engine for computing k-resub. */
  uint32_t num_trials_k{10};

  /*! \brief Maximum number of nodes added by resubstitution. */
  uint32_t max_inserts{2};

  /*! \brief Maximum fanout of a node to be considered as root. */
  uint32_t skip_fanout_limit_for_roots{1000};

  /*! \brief Maximum fanout of a node to be considered as divisor. */
  uint32_t skip_fanout_limit_for_divisors{100};

  /*! \brief Whether to save generated patterns into file. */
  std::optional<std::string> write_pats{};

  /*! \brief Whether to scan and substitute constant nodes first.
             Only safe if the provided patterns are stuck-at checked. */
  bool check_const{false};

  /*! \brief Whether to utilize ODC, and how many levels. 0 = no. -1 = Consider TFO until PO. */
  int odc_levels{0};

  /*! \brief Show progress. */
  bool progress{false};

  /*! \brief Be verbose. */
  bool verbose{false};

  /*! \brief Maximum number of PIs of reconvergence-driven cuts. */
  uint32_t max_pis{8};

  uint32_t conflict_limit{1000};

  uint32_t random_seed{0};
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
  stopwatch<>::duration time_div1_compare{0};
  uint32_t num_div1_accepts{0};

  /* time & number of k-resub */
  stopwatch<>::duration time_resubk{0};
  stopwatch<>::duration time_compute_function{0};
  uint32_t num_divk_accepts{0};

  /* time & number of (~)r == d1 ^ d2 node substitutions */
  stopwatch<>::duration time_xor{0};
  uint32_t num_xor_accepts{0};

  /*! \brief Initial network size (before resubstitution) */
  uint64_t initial_size{0};

  uint32_t num_cex{0};

  uint32_t num_cex_div0{0};
  uint32_t num_cex_div1{0};
  uint32_t num_cex_divk{0};
  uint32_t num_cex_xor{0};

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
bool substitute_fn( Ntk& ntk, typename Ntk::node const& n, typename Ntk::signal const& g )
{
  ntk.substitute_node( n, g );
  //std::cout<<"substitute node "<<unsigned(n)<<" with node "<<unsigned(ntk.get_node(g))<<std::endl;
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
  using validator_t = circuit_validator<Ntk, bill::solvers::z3, true, true>;
  using vgate = typename validator_t::gate;
  using fanin = typename vgate::fanin;
  using gtype = typename validator_t::gate_type;

  struct unate_divisors
  {
    using signal = typename Ntk::signal;

    std::vector<std::pair<signal, uint32_t>> positive_divisors;
    std::vector<std::pair<signal, uint32_t>> negative_divisors;

    void clear()
    {
      positive_divisors.clear();
      negative_divisors.clear();
    }

    void sort()
    {
      std::sort(positive_divisors.begin(), positive_divisors.end(), [](std::pair<signal, uint32_t> a, std::pair<signal, uint32_t> b) {
          return a.second > b.second;   
      });
      std::sort(negative_divisors.begin(), negative_divisors.end(), [](std::pair<signal, uint32_t> a, std::pair<signal, uint32_t> b) {
          return a.second > b.second;   
      });
    }
  };

  explicit simresub_impl( NtkBase& ntkbase, Ntk& ntk, simresub_params const& ps, simresub_stats& st, partial_simulator& sim, validator_params const& vps, resub_callback_t const& callback = substitute_fn<NtkBase> )
    : ntkbase( ntkbase ), ntk( ntk ), ps( ps ), st( st ), callback( callback ),
      tts( ntk ), sim( sim ), validator( ntk, vps )
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
    ntk._events->on_add.emplace_back( [&]( const auto& n ){
      validator.add_node( n );
      call_with_stopwatch( st.time_sim, [&]() {
        simulate_nodes<Ntk>( ntk, tts, sim );
      });
    });
    
    ntk._events->on_modified.emplace_back( update_level_of_existing_node );
    
    ntk._events->on_delete.emplace_back( update_level_of_deleted_node );
  }

  void run()
  {
    stopwatch t( st.time_total );

    /* start the managers */
    progress_bar pbar{ntk.size(), "resub |{0}| node = {1:>4}   cand = {2:>4}   est. gain = {3:>5}", ps.progress};

    //std::vector<node> PIs( ntk.num_pis() );
    //ntk.foreach_pi( [&]( auto const& n, auto i ){ PIs.at(i) = n; });

    call_with_stopwatch( st.time_sim, [&]() {
      simulate_nodes<Ntk>( ntk, tts, sim );
    });

    if ( ps.max_inserts > 1u )
    {
      abcresub::Abc_ResubPrepareManager( sim.compute_constant( false ).num_blocks() );
    }

    if ( ps.check_const )
    {
      auto const zero = sim.compute_constant( false );
      auto const one = sim.compute_constant( true );

      ntk.foreach_gate( [&]( auto const& n ){
        if ( tts[n] == zero )
          ntk.substitute_node( n, ntk.get_constant( false ) );
        else if ( tts[n] == one )
          ntk.substitute_node( n, ntk.get_constant( true ) );
      });

      call_with_stopwatch( st.time_sim, [&]() {
        simulate_nodes<Ntk>( ntk, tts, sim );
      });
    }

    /* iterate through all nodes and try to replace it */
    auto const size = ntk.num_gates();
    ntk.foreach_gate( [&]( auto const& n, auto i ){
        if ( i >= size )
          return false; /* terminate */

        pbar( i, i, candidates, st.estimated_gain );

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
        candidates++;
        st.estimated_gain += last_gain;

        /* update network */
        call_with_stopwatch( st.time_callback, [&]() {
            if ( ps.odc_levels != 0 ) validator.update();
            return callback( ntkbase, n, *g );
          });

        return true; /* next */
      });
    if ( ps.max_inserts > 1u )
    {
      abcresub::Abc_ResubPrepareManager( 0 );
    }
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

    /* collect the internal nodes */
    if ( ntk.value( n ) == 0 && n != 0 ) /* ntk.fanout_size( n ) */
      internal.emplace_back( n );

    ntk.foreach_fanin( n, [&]( const auto& f ){
        collect_divisors_rec( ntk.get_node( f ), internal );
      }); 
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
          if ( ntk.is_dead( p ) || ntk.visited( p ) == ntk.trav_id() || ntk.level( p ) > required )
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

  void collect_unate_divisors( node const& root, uint32_t required )
  {
    udivs.clear();

    auto const& tt = get_tt( root );;
    for ( auto i = 0u; i < num_divs; ++i )
    {
      auto const d = divs.at( i );

      if ( ntk.level( d ) > required - 1 )
        continue;

      auto const& tt_d = get_tt( d );;

      /* check positive containment */
      if ( kitty::implies( tt_d, tt ) )
      {
        udivs.positive_divisors.emplace_back( std::make_pair( ntk.make_signal( d ), kitty::count_ones( tt_d & tt ) ) );
        continue;
      }
      if ( kitty::implies( ~tt_d, tt ) )
      {
        udivs.positive_divisors.emplace_back( std::make_pair( !ntk.make_signal( d ), kitty::count_ones( ~tt_d & tt ) ) );
        continue;
      }

      /* check negative containment */
      if ( kitty::implies( tt, tt_d ) )
      {
        udivs.negative_divisors.emplace_back( std::make_pair( ntk.make_signal( d ), kitty::count_zeros( tt_d & tt ) ) );
        continue;
      }
      if ( kitty::implies( tt, ~tt_d ) )
      {
        udivs.negative_divisors.emplace_back( std::make_pair( !ntk.make_signal( d ), kitty::count_zeros( ~tt_d & tt ) ) );
        continue;
      }
    }

    udivs.sort();
  }

  void found_cex()
  {
    ++st.num_cex;
    sim.add_pattern( validator.cex );

    /* re-simulate */
    call_with_stopwatch( st.time_sim, [&]() {
      simulate_nodes<Ntk>( ntk, tts, sim );
    });

    if ( ps.max_inserts > 1u )
    {
      abcresub::Abc_ResubPrepareManager( sim.compute_constant( false ).num_blocks() );
    }
  }

  kitty::partial_truth_table get_tt( node const& n, bool inverse = false )
  {
    if ( ps.odc_levels == 0 )
      return inverse? ~tts[n]: tts[n];

    return ( inverse? ~tts[n]: tts[n] ) | observability_dont_cares( ntk, n, sim, tts, ps.odc_levels );
  }

  bool is_and( kitty::partial_truth_table const& tt1, kitty::partial_truth_table const& tt2, kitty::partial_truth_table const& tt )
  {
    for ( auto i = 0u; i < tt.num_blocks(); ++i )
    {
      if ( ( tt1._bits[i] & tt2._bits[i] ) != tt._bits[i] )
        return false;
    }
    return true;
  }
  bool is_or( kitty::partial_truth_table const& tt1, kitty::partial_truth_table const& tt2, kitty::partial_truth_table const& tt )
  {
    for ( auto i = 0u; i < tt.num_blocks(); ++i )
    {
      if ( ( tt1._bits[i] | tt2._bits[i] ) != tt._bits[i] )
        return false;
    }
    return true;
  }

  bool is_xor( kitty::partial_truth_table const& tt1, kitty::partial_truth_table const& tt2, kitty::partial_truth_table const& tt )
  {
    for ( auto i = 0u; i < tt.num_blocks(); ++i )
    {
      if ( ( tt1._bits[i] ^ tt2._bits[i] ) != tt._bits[i] )
        return false;
    }
    return true;
  }

private:
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
#if 1
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

    if ( ps.max_inserts < 1 || num_mffc <= 1 )
    {
      return std::nullopt;
    }

    /* collect level one divisors */
    call_with_stopwatch( st.time_collect_unate_divisors, [&]() {
        collect_unate_divisors( root, required );
      });

    /* consider single-gate resub */
    g = call_with_stopwatch( st.time_resub1, [&]() {
        return resub_div1( root, required );
      } );
    if ( g )
    {
      ++st.num_div1_accepts;
      last_gain = num_mffc - 1;
      return g; /* accepted resub */
    }

    if constexpr ( std::is_same<NtkBase, xag_network>::value )
    {
      /* consider X(N)OR resub */
      g = call_with_stopwatch( st.time_xor, [&]() {
          return resub_xor( root, required );
        } );
      if ( g )
      {
        ++st.num_xor_accepts;
        last_gain = num_mffc - 1;
        return g; /* accepted resub */
      }
    }

    if ( ps.max_inserts < 2 || num_mffc <= 2 )
    {
      return std::nullopt;
    }
#endif
    /* try k-resub */
    uint32_t size = 0;
    g = call_with_stopwatch( st.time_resubk, [&]() {
        return resub_divk( root, std::min( num_mffc - 1, int( ps.max_inserts ) ), size );
      } );
    if ( g )
    {
      ++st.num_divk_accepts;
      last_gain = num_mffc - size;
      return g; /* accepted resub */
    }

#if 0
    if ( ps.max_inserts < 3 || num_mffc <= 3 )
    {
      return std::nullopt;
    }

    /* consider X(N)OR resub */
    g = call_with_stopwatch( st.time_xor, [&]() {
        return resub_xor( root, required );
      } );
    if ( g )
    {
      ++st.num_xor_accepts;
      last_gain = num_mffc - 3;
      return g; /* accepted resub */
    }
#endif

    return std::nullopt;
  }

  std::optional<signal> resub_div0( node const& root, uint32_t required ) 
  {
    (void)required;
    auto tt = get_tt( root );
    auto ntt = get_tt( root, true );

    for ( auto i = 0u; i < num_divs; ++i )
    //for ( int i = num_divs-1; i >= 0; --i )
    {
      auto const d = divs.at( i );
      auto ttd = get_tt( d );

      if ( tt == ttd || ntt == ttd )
      {
        auto const g = ( tt == ttd ) ? ntk.make_signal( d ) : !ntk.make_signal( d );

        const auto valid = call_with_stopwatch( st.time_sat, [&]() {
          return validator.validate( root, g );
        });

        if ( !valid ) /* timeout */
        {
          continue;
        }
        else if ( *valid )
        {
          return g;
        }
        else
        {
          ++st.num_cex_div0;
          found_cex();
        }
      }
    }

    return std::nullopt;
  }

  std::optional<signal> resub_div1( node const& root, uint32_t required )
  {
    (void)required;
    auto const& tt = get_tt( root );;
    auto const& w = kitty::count_ones_fast( tt );
    auto const& nw = tt.num_bits() - w;

    /* check for positive unate divisors */
    for ( auto i = 0u; i < udivs.positive_divisors.size(); ++i )
    {
      auto const& s0 = udivs.positive_divisors.at( i ).first;
      auto const& tt_s0 = get_tt( ntk.get_node( s0 ), ntk.is_complemented(s0) );
      auto const& w_s0 = udivs.positive_divisors.at( i ).second;
      if ( w_s0 < uint32_t( w / 2 ) )
        break;

      for ( auto j = i + 1; j < udivs.positive_divisors.size(); ++j )
      {
        if ( w_s0 + udivs.positive_divisors.at( j ).second < w )
          break;

        auto const& s1 = udivs.positive_divisors.at( j ).first;
        auto const& tt_s1 = get_tt( ntk.get_node( s1 ), ntk.is_complemented(s1) );

        const auto isor = call_with_stopwatch( st.time_div1_compare, [&]() {
            return is_or( tt_s0, tt_s1, tt);
          });

        if ( isor )
        {
          fanin fi1; fi1.idx = 0; fi1.inv = !ntk.is_complemented( s0 );
          fanin fi2; fi2.idx = 1; fi2.inv = !ntk.is_complemented( s1 );
          vgate gate; gate.fanins = {fi1, fi2}; gate.type = gtype::AND;

          const auto valid = call_with_stopwatch( st.time_sat, [&]() {
            return validator.validate( root, {ntk.get_node( s0 ), ntk.get_node( s1 )}, {gate}, true );
          });

          if ( !valid ) /* timeout */
          {
            continue;
          }
          else if ( *valid )
          {
            return ntk.create_or( s0, s1 );
          }
          else
          {
            ++st.num_cex_div1;
            found_cex();
          }
        }
      }
    }

    /* check for negative unate divisors */
    for ( auto i = 0u; i < udivs.negative_divisors.size(); ++i )
    {
      auto const& s0 = udivs.negative_divisors.at( i ).first;
      auto const& tt_s0 = get_tt( ntk.get_node( s0 ), ntk.is_complemented(s0) );
      auto const& w_s0 = udivs.negative_divisors.at( i ).second;
      if ( w_s0 < uint32_t( nw / 2 ) )
        break;

      for ( auto j = i + 1; j < udivs.negative_divisors.size(); ++j )
      {
        if ( w_s0 + udivs.negative_divisors.at( j ).second < nw )
          break;

        auto const& s1 = udivs.negative_divisors.at( j ).first;
        auto const& tt_s1 = get_tt( ntk.get_node( s1 ), ntk.is_complemented(s1) );

        const auto isand = call_with_stopwatch( st.time_div1_compare, [&]() {
            return is_and( tt_s0, tt_s1, tt);
          });

        if ( isand )
        {
          fanin fi1; fi1.idx = 0; fi1.inv = ntk.is_complemented( s0 );
          fanin fi2; fi2.idx = 1; fi2.inv = ntk.is_complemented( s1 );
          vgate gate; gate.fanins = {fi1, fi2}; gate.type = gtype::AND;

          const auto valid = call_with_stopwatch( st.time_sat, [&]() {
            return validator.validate( root, {ntk.get_node( s0 ), ntk.get_node( s1 )}, {gate}, false );
          });

          if ( !valid ) /* timeout */
          {
            continue;
          }
          else if ( *valid )
          {
            return ntk.create_and( s0, s1 );
          }
          else
          {
            ++st.num_cex_div1;
            found_cex();
          }
        }
      }
    }

    return std::nullopt;
  }

  std::optional<signal> resub_divk( node const& root, uint32_t num_inserts, uint32_t& size ) 
  {    
    for ( auto j = 0u; j < ps.num_trials_k; ++j )
    {
      abc_resub rs( 2ul + num_divs, tts[root].num_blocks(), ps.max_divisors_k );
      rs.add_root( root, tts );
      rs.add_divisors( std::begin( divs ), std::begin( divs ) + num_divs, tts );

      auto const res = call_with_stopwatch( st.time_compute_function, [&]() {
        if constexpr ( std::is_same<NtkBase, xag_network>::value )
          return rs.compute_function( num_inserts, true );
        else
          return rs.compute_function( num_inserts, false );
      });
      if ( res )
      {
        auto const& index_list = *res;
        if ( index_list.size() == 1u ) /* div0 or constant */
        {
          const auto valid = call_with_stopwatch( st.time_sat, [&]() {
            if ( index_list[0] < 2 ) return validator.validate( root, ntk.get_constant( bool( index_list[0] ) ) );
            assert( index_list[0] >= 4 );
            return validator.validate( root, bool( index_list[0] % 2 ) ? !ntk.make_signal( divs[(index_list[0] >> 1u) - 2u] ) : ntk.make_signal( divs[(index_list[0] >> 1u) - 2u] ) );
          });

          if ( !valid ) /* timeout */
          {
            break;
          }
          else
          {
            if ( *valid )
            {
              size = 0u;
              if ( index_list[0] < 2 ) return ntk.get_constant( bool( index_list[0] ) );
              else return bool( index_list[0] % 2 ) ? !ntk.make_signal( divs[(index_list[0] >> 1u) - 2u] ) : ntk.make_signal( divs[(index_list[0] >> 1u) - 2u] );
            }
            else
            {
              ++st.num_cex_divk;
              found_cex();
              continue;
            }
          }
        }

        uint64_t const num_gates = ( index_list.size() - 1u ) / 2u;
        std::vector<vgate> gates( num_gates );
        size = 0u;
        for ( auto i = 0u; i < num_gates; ++i )
        {
          fanin f0; f0.idx = uint32_t( ( index_list[2*i] >> 1u ) - 2u ); f0.inv = bool( index_list[2*i] % 2 );
          fanin f1; f1.idx = uint32_t( ( index_list[2*i + 1u] >> 1u ) - 2u ); f1.inv = bool( index_list[2*i + 1u] % 2 );
          gates[i].fanins = { f0, f1 };
          gates[i].type = f0.idx < f1.idx ? gtype::AND : gtype::XOR;
          
          if constexpr ( std::is_same<NtkBase, xag_network>::value )
            ++size;
          else
            size += ( gates[i].type == gtype::AND )? 1u: 3u;
        }
        bool const out_neg = bool( index_list.back() % 2 );

        if ( size > num_inserts )
        {
          std::cout<<"circuit size exceed limit\n";
          return std::nullopt;
        }

        const auto valid = call_with_stopwatch( st.time_sat, [&]() {
          return validator.validate( root, std::begin( divs ), std::begin( divs ) + num_divs, gates, out_neg );
        });

        if ( !valid ) /* timeout */
        {
          break;
        }
        else
        {
          if ( *valid )
          {
            std::vector<signal> ckt;
            for ( auto n : divs )
            {
              ckt.emplace_back( ntk.make_signal( n ) );
            }

            for ( auto g : gates )
            {
              auto const f0 = g.fanins[0].inv ? !ckt[g.fanins[0].idx] : ckt[g.fanins[0].idx];
              auto const f1 = g.fanins[1].inv ? !ckt[g.fanins[1].idx] : ckt[g.fanins[1].idx];
              if ( g.type == gtype::AND )
              {
                ckt.emplace_back( ntk.create_and( f0, f1 ) );
              }
              else if ( g.type == gtype::XOR )
              {
                ckt.emplace_back( ntk.create_xor( f0, f1 ) );
              }
            }
            
            return out_neg ? !ckt.back() : ckt.back();
          }
          else
          {
            ++st.num_cex_divk;
            found_cex();
          }
        }
      }
      else /* loop until no result can be found by the engine */
      {
        return std::nullopt;
      }
    }

    return std::nullopt;
  }

  std::optional<signal> resub_xor( node const& root, uint32_t required ) 
  {
    (void)required;
    auto tt = get_tt( root );
    auto ntt = get_tt( root, true );

    for ( auto i = 0u; i < num_divs - 1; ++i )
    {
      auto const& s0 = divs.at( i );
      auto const& tt_s0 = get_tt( s0 );

      for ( auto j = i + 1; j < num_divs; ++j )
      {
        auto const& s1 = divs.at( j );
        auto const& tt_s1 = get_tt( s1 );

        const auto isxor = call_with_stopwatch( st.time_div1_compare, [&]() {
            return is_xor( tt_s0, tt_s1, tt);
          });
        const auto isxnor = call_with_stopwatch( st.time_div1_compare, [&]() {
            return is_xor( tt_s0, tt_s1, ntt);
          });

        if ( isxor || isxnor )
        {
          fanin fi1; fi1.idx = 0; fi1.inv = false;
          fanin fi2; fi2.idx = 1; fi2.inv = false;
          vgate gate; gate.fanins = {fi1, fi2}; gate.type = gtype::XOR;

          const auto valid = call_with_stopwatch( st.time_sat, [&]() {
            return validator.validate( root, {s0, s1}, {gate}, isxnor );
          });

          if ( !valid ) /* timeout */
          {
            continue;
          }
          else if ( *valid )
          {
            return isxor ? ntk.create_xor( ntk.make_signal( s0 ), ntk.make_signal( s1 ) ) : !ntk.create_xor( ntk.make_signal( s0 ), ntk.make_signal( s1 ) );
          }
          else
          {
            ++st.num_cex_xor;
            found_cex();
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

  resub_callback_t const& callback;

  /* temporary statistics for progress bar */
  uint32_t candidates{0};
  uint32_t last_gain{0};

  TT tts;
  partial_simulator& sim;

  validator_t validator;

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

  validator_params vps;
  vps.odc_levels = ps.odc_levels;
  vps.conflict_limit = ps.conflict_limit;
  vps.randomize = true;
  vps.random_seed = ps.random_seed;

  simresub_stats st;

  detail::simresub_impl<Ntk, resub_view_t> p( ntk, resub_view, ps, st, sim, vps, detail::substitute_fn<Ntk> );
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