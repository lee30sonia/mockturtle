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
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"
#include "../views/cnf_view.hpp"
#include <mockturtle/algorithms/simulation.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/operators.hpp>
#include "../utils/node_map.hpp"
#include "cnf.hpp"
#include "cleanup.hpp"
#include <percy/solvers/bsat2.hpp>
#include <bill/sat/interface/abc_bsat2.hpp>

namespace mockturtle
{

struct patgen_params
{
  /*! \brief Number of initial simulation patterns = 2^num_pattern_base. */
  uint32_t num_pattern_base{7};

  /*! \brief Number of reserved blocks(64 bits) for generated simulation patterns. */
  uint32_t num_reserved_blocks{1};

  /*! \brief Maximum number of divisors to consider. */
  uint32_t max_divisors{150};

  /*! \brief Maximum number of nodes added by resubstitution. */
  uint32_t max_inserts{2};

  /*! \brief Maximum fanout of a node to be considered as root. */
  uint32_t skip_fanout_limit_for_roots{1000};

  /*! \brief Maximum fanout of a node to be considered as divisor. */
  uint32_t skip_fanout_limit_for_divisors{100};

  /*! \brief Show progress. */
  bool progress{false};

  /*! \brief Be verbose. */
  bool verbose{false};

  /*! \brief Maximum number of PIs of reconvergence-driven cuts. */
  uint32_t max_pis{8};

  std::default_random_engine::result_type random_seed{0};
};

struct patgen_stats
{
  stopwatch<>::duration time_total{0};

  /* total time for initial simulation and complete pattern generation */
  stopwatch<>::duration time_simgen{0};

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

  /* time for doing substitution */
  stopwatch<>::duration time_substitute{0};

  /* time & number of r == d (equal) node substitutions */
  stopwatch<>::duration time_resub0{0};
  uint32_t num_div0_accepts{0};

  /* time & number of r == d1 &| d2 node substitutions */
  stopwatch<>::duration time_resub1{0};
  uint32_t num_div1_accepts{0};

  /*! \brief Initial network size (before resubstitution) */
  uint64_t initial_size{0};

  uint32_t num_constant{0};
  uint32_t num_generated_patterns{0};
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

template<class NtkBase, class Ntk>
class patgen_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using TT = unordered_node_map<kitty::partial_truth_table, Ntk>;

  explicit patgen_impl( NtkBase& ntkbase, Ntk& ntk, simresub_params const& ps, simresub_stats& st )
    : ntkbase( ntkbase ), ntk( ntk ), ps( ps ), st( st ), 
      tts( ntk ), phase( ntkbase ), sim( ntk.num_pis(), ps.num_pattern_base, ps.num_reserved_blocks, ps.random_seed )
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

    /* simulate all nodes and generate complete test patterns, finding out and replace constant nodes at the same time */
    call_with_stopwatch( st.time_simgen, [&]() {
      simulate_generate();
    });
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
  
  void simulate_generate()
  {
    call_with_stopwatch( st.time_sim, [&]() {
      simulate_nodes<Ntk>( ntk, tts, sim );
    });

    std::vector<bill::lit_type> assumptions( 1 ); /* bill::result::clause_type */
    kitty::partial_truth_table zero = sim.compute_constant(false);
  
    ntk.foreach_gate( [&]( auto const& n ) 
    {
      //std::cout<<"processing node "<<unsigned(n)<<std::endl;
      if ( (tts[n] == zero) || (tts[n] == ~zero) )
      {
        assumptions[0] = (tts[n] == ~zero)? ~( ntk.lit( n ) ): ntk.lit( n ); //lit_not_cond( literals[n], (tts[n] == ~zero) );
      
        const auto res = call_with_stopwatch( st.time_sat, [&]() {
          return ntk.solve( assumptions );
        });
        
        if ( res )
        {
          if ( *res )
          {
            //std::cout << "SAT: add pattern. (" << n << ")" << std::endl;
            std::vector<bool> pattern = ntk.pi_model_values();
            sim.add_pattern(pattern);
            ++st.num_generated_patterns;
  
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
            auto g = ntk.get_constant( tts[n] == ~zero );
            /* update network */
            call_with_stopwatch( st.time_substitute, [&]() {
              ntk.substitute_node( n, g );
            });
  
          }
        }
      }

      return true; /* next gate */
    } );

    normalizeTT();
  }

  void normalizeTT()
  {
    ntk.foreach_gate( [&]( auto const& n ){
      if ( kitty::get_bit( tts[n], 0 ) )
      {
        phase[n] = true;
        tts[n] = ~tts[n];
      }
      else
        phase[n] = false;
    });
  }

  void un_normalizeTT()
  {
    ntk.foreach_gate( [&]( auto const& n ){
      if ( phase.has(n) && phase[n] )
      {
        tts[n] = ~tts[n];
      }
    });
  }
  

private:
  NtkBase& ntkbase;
  Ntk& ntk;

  simresub_params const& ps;
  simresub_stats& st;

  /* temporary statistics for progress bar */
  uint32_t candidates{0};
  uint32_t last_gain{0};

  TT tts;
  unordered_node_map<bool, NtkBase> phase;
  partial_simulator<kitty::partial_truth_table> sim;

  unate_divisors udivs;
  binate_divisors bdivs;

  std::vector<node> temp;
  std::vector<node> divs;
  uint32_t num_divs{0};
};

} /* namespace detail */

template<class Ntk>
void pattern_generation( Ntk& ntk, patgen_params const& ps = {}, patgen_stats* pst = nullptr )
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

  using resub_view_t = cnf_view<fanout_view<depth_view<Ntk>>, true>;//, bill::solvers::bsat2>;
  using depth_view_t = depth_view<Ntk>;
  depth_view_t depth_view{ntk};
  fanout_view<depth_view_t> fanout_view{depth_view};
  resub_view_t resub_view(fanout_view, {.auto_update = false});
  simresub_stats st;

  detail::patgen_impl<Ntk, resub_view_t> p( ntk, resub_view, ps, st );
  p.run();

  if ( ps.verbose )
    st.report();

  if ( pst )
    *pst = st;
}

} /* namespace mockturtle */
