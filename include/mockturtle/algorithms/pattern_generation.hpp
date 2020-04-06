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
#include "../views/cnf_view.hpp"
#include <mockturtle/algorithms/simulation.hpp>
#include <kitty/partial_truth_table.hpp>
//#include <kitty/operators.hpp>
#include "../utils/node_map.hpp"
#include "cnf.hpp"
#include "cleanup.hpp"
#include <percy/solvers/bsat2.hpp>
#include <bill/sat/interface/abc_bsat2.hpp>

namespace mockturtle
{

struct patgen_params
{
  /*! \brief Number of initial random patterns to start with. */
  uint32_t num_random_pattern{1000};

  /*! \brief Whether to substitute constant nodes. */
  bool substitute_const{true};

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
    : ntk( ntk ), ps( ps ), st( st ), 
      tts( ntk ), sim( ntk.num_pis(), ps.num_random_pattern, ps.random_seed )
  {
    st.num_total_patterns = ps.num_random_pattern;
  }

  void run()
  {
    stopwatch t( st.time_total );

    /* start the managers */
    //progress_bar pbar{ntk.size(), "resub |{0}| node = {1:>4}   cand = {2:>4}   est. gain = {3:>5}", ps.progress};

    simulate_generate();
  }

private:
  /* `me` is one of `fo`'s fanin, get the other one. assuming 2-input gates */
  signal get_other_fanin( node const& fo, node const& me )
  {
    signal ret;
    ntk.foreach_fanin( fo, [&]( auto const& foi ){
      if ( ntk.get_node( foi ) != me )
      {
        ret = foi;
        return false;
      }
      return true;
    });
    return ret;
  }
#if 0
  void DFS( node const& n_ori, node const& n, std::vector<bool>& pattern, std::vector<pabc::lit>& assumptions, bool& fReturn, uint32_t& nFail )
  {
    //std::cout<<"DFS: node "<<unsigned(n)<<std::endl;
    ntk.foreach_fanout( n, [&]( auto const& fo ){
      if ( fReturn ) /* whether we have gotten an observable pattern */
        return false; /* terminate foreach_fanout, also terminating recursive DFS */

      /* add assumption if the other fanin of this node is controlling value */
      auto foi = get_other_fanin( fo, n );
      //std::cout<<"add assumption for fanout "<<unsigned(fo)<<", the other fanin is "<<(ntk.is_complemented(foi)?"~":"")<<unsigned(ntk.get_node(foi))<<std::endl;
      assumptions.emplace_back( lit_not_cond( literals[ntk.get_node(foi)], ntk.is_complemented(foi) ) );

      bool fEnd = false; /* whether has reached the end of a path */
      ntk.foreach_po( [&]( auto const& f ){ 
        if ( ntk.get_node(f) == fo ) /* DFS reached one of POs */
        {
          //std::cout<<"reach PO "<<unsigned(ntk.get_node(f))<<" try solving with "<<assumptions.size()<<" assumptions"<<std::endl;
          //for (auto j=1u; j<assumptions.size(); ++j) std::cout<<int(assumptions[j]/2)<<" "; std::cout<<"\n";
          const auto res = call_with_stopwatch( st.time_sat, [&]() {
            return solver.solve( &assumptions[0], &assumptions[0] + assumptions.size(), 0 );
          });
          if ( res == percy::synth_result::success )
          {
            for ( auto i = 1u; i <= ntk.num_pis(); ++i )
              pattern[i-1] = (solver.var_value( i ));
            fReturn = call_with_stopwatch( st.time_odc, [&]() { 
              return pattern_is_observable( ntk, n_ori, pattern, POs );
            });
            std::cout<<"SAT, now observable? "<<fReturn<<"\n";
          }
          if ( !fReturn ) /* UNSAT or SAT but still not observable */
          {
            if ( ++nFail >= ps.ODC_failure_limit )
            {
              std::cout<<"conflict limit reached!"<<std::endl;
              fReturn = true;
            }
          }

          fEnd = true;
          return false; /* terminate foreach_po */
        }
        return true; /* check the next PO */
      });

      if ( !fEnd )
        DFS( n_ori, fo, pattern, assumptions, fReturn, nFail );
      
      //std::cout<<"pop assumption for fanout "<<unsigned(fo)<<std::endl;
      assumptions.pop_back();

      return true; /* try next fanout (another path) */
    });   
  }

  bool generate_observable_pattern( node const& n, bool value, std::vector<bool>& pattern )
  {
    //std::cout<<"generating pattern for n to be "<<value<<" and observable.\n";
    std::vector<pabc::lit> assumptions( 1 );
    assumptions[0] = lit_not_cond( literals[n], !value );

    bool fReturn = false; /* whether we have gotten an observable pattern */
    uint32_t nFail = 0; /* number of failed paths */
    DFS( n, n, pattern, assumptions, fReturn, nFail );
    return fReturn;
  }
#endif
  bool generate_observable_pattern( node const& n, bool value, std::vector<bool>& pattern, unordered_node_map<bool, Ntk> const& po_vals )
  {
    //std::cout<<"generating pattern for n to be "<<value<<" and observable.\n";
    ntk.deactivate( n );

    std::vector<bill::lit_type> assumptions;
    if ( value )
    {
      ntk.foreach_fanin( n, [&]( auto const& f ){
        assumptions.emplace_back( ntk.lit( f ) );
        return true;
      });
    }
    else
    {
      std::vector<bill::lit_type> clause;
      auto nlit = bill::lit_type( ntk.add_var(), bill::lit_type::polarities::positive );
      ntk.foreach_fanin( n, [&]( auto const& f ){
        clause.emplace_back( ~( ntk.lit( f ) ) );
        return true;
      });
      clause.emplace_back( nlit );
      ntk.add_clause( clause );
      assumptions.emplace_back( ~nlit );
    }
    assumptions.emplace_back( value? ~( ntk.lit( n ) ): ntk.lit( n ) );
    ntk.foreach_po( [&]( auto const& f ){ 
      assumptions.emplace_back( po_vals[f] ? ~( ntk.lit( ntk.get_node(f) ) ): ntk.lit( ntk.get_node(f) ) );
      return false; /* only try for the first PO now */
    });

    const auto res = call_with_stopwatch( st.time_sat, [&]() {
      return ntk.solve( assumptions );
    });

    if ( res )
    {
      std::cout<<"generation "<< ((*res)? "SAT": "UNSAT") ;//<<"\n";
      if ( *res )
        pattern = ntk.pi_model_values();
      ntk.activate( n );
      return *res;
    }

    ntk.activate( n );
    return false;
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
        bool value = !(tts[n] == ~zero); /* wanted value of n */

        assumptions[0] = value? ntk.lit( n ): ~( ntk.lit( n ) ); //lit_not_cond( literals[n], (tts[n] == ~zero) );
      
        const auto res = call_with_stopwatch( st.time_sat, [&]() {
          return ntk.solve( assumptions );
        });
        
        if ( res )
        {
          if ( *res )
          {
            //std::cout << "SAT: add pattern. (" << n << ")" << std::endl;
            std::vector<bool> pattern = ntk.pi_model_values();
#if 0
            /* check if the found pattern is observable */ 
            ntk.foreach_po( [&]( auto const& f, auto i ){ POs.at(i) = ntk.get_node( f ); });
            unordered_node_map<bool, Ntk> po_vals( ntk );
            bool observable = call_with_stopwatch( st.time_odc, [&]() { 
                return pattern_is_observable( ntk, n, pattern, POs, po_vals );
              });
            if ( !observable )
            {
              std::cout << "generated pattern is not observable! (" << unsigned(n) << " = " << value <<")" << std::endl;
              generate_observable_pattern( n, value, pattern, po_vals );
              unordered_node_map<bool, Ntk> po_vals2( ntk );
              std::cout << " after re-gen, now? " << pattern_is_observable( ntk, n, pattern, POs, po_vals2 ) <<"\n";
            }
#endif
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
      }
#if 0
      else
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
          #if 0
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
          }
          #endif
        }
        else if ( ( tts[n] | odc ) == sim.compute_constant( true ) )
        {
          std::cout << "under all observable patterns node "<< unsigned(n)<<" is always 1!" << std::endl;
          #if 0
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
          }
          #endif
        }
      }
#endif
      return true; /* next gate */
    } );
  }

private:
  Ntk& ntk;

  patgen_params const& ps;
  patgen_stats& st;

  TT tts;
public:
  partial_simulator<kitty::partial_truth_table> sim;
};

} /* namespace detail */

template<class Ntk>
partial_simulator<kitty::partial_truth_table> pattern_generation( Ntk& ntk, patgen_params const& ps = {}, patgen_stats* pst = nullptr )
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

  using view_t = cnf_view<Ntk, true>; //, bill::solvers::bsat2>
  view_t view( ntk, {.auto_update = false} );

  patgen_stats st;

  detail::patgen_impl p( view, ps, st );
  p.run();

  if ( ps.verbose )
    st.report();

  if ( pst )
    *pst = st;

  return p.sim;
}

} /* namespace mockturtle */
