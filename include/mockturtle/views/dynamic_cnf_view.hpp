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
  \file cnf_view.hpp
  \brief Creates a CNF while creating a network

  \author Siang-Yun Lee (based on cnf_view.hpp by Mathias Soeken)
*/

#pragma once

#include <cstdint>
#include <vector>

#include "../algorithms/cnf.hpp"
#include "../traits.hpp"

#include <fmt/format.h>
#include <percy/solvers/bsat2.hpp>

namespace mockturtle
{

/* Differences in comparison to cnf_view:
 *** Should be created with existing network.
 *** Can add custom variables to the solver (which is not associated with any node).
 *** Can modify or delete nodes. Clauses are modified accordingly automatically.
 */

/*
template<typename Ntk, bool has_cnf_interface = has_solve_v<Ntk>&& has_value_v<Ntk>>
class dynamic_cnf_view
{
};

template<typename Ntk>
class dynamic_cnf_view<Ntk, true> : public Ntk
{
public:
  dynamic_cnf_view( Ntk const& ntk ) : Ntk( ntk )
  { }
};
*/

struct dynamic_cnf_view_params
{
  bool allow_modify{false};
};

template<typename Ntk>
class dynamic_cnf_view : public Ntk
{
public:
  using storage = typename Ntk::storage;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  /*dynamic_cnf_view( dynamic_cnf_view_params const& ps = {} ) 
  : ps_( ps ), literals_( node_literals<Ntk>() ), switches_( std::nullopt )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
    static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
    static_assert( has_make_signal_v<Ntk>, "Ntk does not implement the make_signal method" );
    static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
    static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
    static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
    static_assert( has_node_function_v<Ntk>, "Ntk does not implement the node_function method" );

    register_events();
  }*/

  dynamic_cnf_view( Ntk const& ntk, dynamic_cnf_view_params const& ps = {} )
  : ntk_( ntk ), ps_( ps ), literals_( node_literals( ntk ) ), switches_( ntk.size() )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
    static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
    static_assert( has_make_signal_v<Ntk>, "Ntk does not implement the make_signal method" );
    static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
    static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
    static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
    static_assert( has_node_function_v<Ntk>, "Ntk does not implement the node_function method" );

    register_events();

    /* unit clause for constant-0 */
    solver_.add_clause( {1} );

    solver_.set_nr_vars( ntk_.size() );
    ntk_.foreach_gate( [&]( auto const& n ) {
      on_add( n, false );
    } );
  }

  /* \brief Returns the variable associated to a node. */
  inline uint32_t var( node const& n ) const
  {
    return lit2var( literals_[n] );
  }

  /*! \brief Returns the literal associated to a node. */
  inline uint32_t lit( node const& n ) const
  {
    return literals_[n];
  }

  /*! \brief Returns the literal associated to a signal. */
  inline uint32_t lit( signal const& f ) const
  {
    return lit_not_cond( literals_[f], ntk_.is_complemented( f ) );
  }

  /*! \brief Returns the switching literal associated to a node. */
  inline uint32_t switch_lit( node const& n ) const
  {
    assert( ps_.allow_modify && "ps_.allow_modify is not turned on" );
    assert( !ntk_.is_pi( n ) && !ntk_.is_constant( n ) && "PI and constant node are not switch-able" );
    return switches_[ntk_.node_to_index( n )];
  }

  /*! \brief Solves the network with a set of custom assumptions.
   *
   * This function does not assert any primary output, unless specified
   * explicitly through the assumptions.
   *
   * The function returns `nullopt`, if no solution can be found (due to a
   * conflict limit), or `true` in case of SAT, and `false` in case of UNSAT.
   *
   * \param assumptions Vector of literals to be assumped when solving
   * \param limit Conflict limit (unlimited if 0)
   */
  inline std::optional<bool> solve( std::vector<int> const& assumptions_, int limit = 0 )
  {
    std::vector<int> assumptions = assumptions_;
    for ( auto i = 1u; i < switches_.size(); ++i )
    {
      if ( ntk_.is_pi( ntk_.index_to_node( i ) ) ) continue;
      assumptions.push_back( switches_[i] );
    }
    const auto res = solver_.solve( &assumptions[0], &assumptions[0] + assumptions.size(), limit );

    switch ( res )
    {
    case percy::success:
      return true;
    case percy::failure:
      return false;
    default:
    case percy::timeout:
      return std::nullopt;
    }
  }

  /*! \brief Solves the network by asserting all primary outputs to be true
   *
   * The function returns `nullopt`, if no solution can be found (due to a
   * conflict limit), or `true` in case of SAT, and `false` in case of UNSAT.
   *
   * \param limit Conflict limit (unlimited if 0)
   */
  inline std::optional<bool> solve( int limit = 0 )
  {
    std::vector<int> assumptions;
    ntk_.foreach_po( [&]( auto const& f ) {
      assumptions.push_back( lit( f ) );
    } );
    return solve( assumptions, limit );
  }

  /*! \brief Return model value for a node. */
  inline bool value( node const& n )
  {
    return solver_.var_value( var( n ) );
  }

  /*! \brief Return model value for a node (takes complementation into account). */
  inline bool value( signal const& f )
  {
    return value( ntk_.get_node( f ) ) != ntk_.is_complemented( f );
  }

  /*! \brief Whether a node is currently activated (included in CNF). */
  inline bool is_activated( node const& n ) const
  {
    return switch_lit( n ) & 0x1; 
    /* clauses are activated if switch literal is complemented */
  }

  /* \brief Returns all model values for all primary inputs. */
  std::vector<bool> pi_values()
  {
    std::vector<bool> values( ntk_.num_pis() );
    ntk_.foreach_pi( [&]( auto const& n, auto i ) {
      values[i] = value( n );
    } );
    return values;
  }

  /*! \brief Blocks last model for primary input values. */
  void block()
  {
    std::vector<uint32_t> blocking_clause( ntk_.num_pis() );
    ntk_.foreach_pi( [&]( auto const& n, auto i ) {
      blocking_clause[i] = lit_not_cond( literals_[n], value( n ) );
    });
    add_clause( blocking_clause );
  }

  /*! \brief Deactivates the clauses for a node. */
  void deactivate( node const& n )
  {
    if ( is_activated( n ) )
      switches_[ntk_.node_to_index( n )] ^= 0x1;
  }

  /*! \brief (Re-)activates the clauses for a node. */
  void activate( node const& n )
  {
    if ( !is_activated( n ) )
      switches_[ntk_.node_to_index( n )] ^= 0x1;
  }

  /*! \brief Number of variables. */
  inline uint32_t num_vars()
  {
    return solver_.nr_vars();
  }

  /*! \brief Number of clauses. */
  inline uint32_t num_clauses()
  {
    return solver_.nr_clauses();
  }

  /*! \brief Adds a clause to the solver. */
  void add_clause( std::vector<uint32_t> const& clause )
  {
    solver_.add_clause( clause );
  }

  /*! \brief Adds a clause from signals to the solver. */
  void add_clause( std::vector<signal> const& clause )
  {
    std::vector<uint32_t> lits( clause.size() );
    std::transform( clause.begin(), clause.end(), lits.begin(), [&]( auto const& s ) { return lit( s ); } );
    solver_.add_clause( lits );
  }

  /*! \brief Adds a clause to the solver.
   *
   * Entries are either all literals or network signals.
   */
  template<typename... Lit, typename = std::enable_if_t<
    std::disjunction_v<
      std::conjunction<std::is_same<Lit, uint32_t>...>,
      std::conjunction<std::is_same<Lit, signal>...>>>>
  void add_clause( Lit... lits )
  {
    if constexpr ( std::conjunction_v<std::is_same<Lit, uint32_t>...> )
    {
      solver_.add_clause( std::vector<uint32_t>{{lits...}} );
    }
    else
    {
      solver_.add_clause( std::vector<uint32_t>{{lit( lits )...}} );
    }
  }

private:
  void register_events()
  {
    ntk_.events().on_add.push_back( [this]( auto const& n ) { on_add( n ); } );
    ntk_.events().on_modified.push_back( [this]( auto const& n, auto const& previous ) {
      (void)previous; 
      on_modified( n ); 
    } );
    ntk_.events().on_delete.push_back( [this]( auto const& n ) { on_delete( n ); } );
  }

  void on_add( node const& n, bool add_var = true )
  {
    uint32_t node_lit;
    if ( add_var )
    {
      solver_.add_var();
      node_lit = make_lit( solver_.nr_vars()-1 );
      literals_.resize();
      literals_[n] = node_lit;
    }
    else
      node_lit = literals_[n];

    uint32_t switch_lit;
    if ( ps_.allow_modify )
    {
      solver_.add_var();
      switch_lit = make_lit( solver_.nr_vars()-1 );
      switches_.resize( ntk_.size() );
      switches_[ntk_.node_to_index( n )] = lit_not( switch_lit );
    }

    const auto _add_clause = [&]( std::vector<uint32_t> const& clause ) {
      if ( ps_.allow_modify )
      {
        std::vector<uint32_t> clause_ = clause;
        clause_.push_back( switch_lit );
        add_clause( clause_ );
      }
      else
        add_clause( clause );
    };

    std::vector<uint32_t> child_lits;
    ntk_.foreach_fanin( n, [&]( auto const& f ) {
      child_lits.push_back( lit( f ) );
    } );

    if constexpr ( has_is_and_v<Ntk> )
    {
      if ( ntk_.is_and( n ) )
      {
        detail::on_and( node_lit, child_lits[0], child_lits[1], _add_clause );
        return;
      }
    }

    if constexpr ( has_is_or_v<Ntk> )
    {
      if ( ntk_.is_or( n ) )
      {
        detail::on_or( node_lit, child_lits[0], child_lits[1], _add_clause );
        return;
      }
    }

    if constexpr ( has_is_xor_v<Ntk> )
    {
      if ( ntk_.is_xor( n ) )
      {
        detail::on_xor( node_lit, child_lits[0], child_lits[1], _add_clause );
        return;
      }
    }

    if constexpr ( has_is_maj_v<Ntk> )
    {
      if ( ntk_.is_maj( n ) )
      {
        detail::on_maj( node_lit, child_lits[0], child_lits[1], child_lits[2], _add_clause );
        return;
      }
    }

    if constexpr ( has_is_ite_v<Ntk> )
    {
      if ( ntk_.is_ite( n ) )
      {
        detail::on_ite( node_lit, child_lits[0], child_lits[1], child_lits[2], _add_clause );
        return;
      }
    }

    if constexpr ( has_is_xor3_v<Ntk> )
    {
      if ( ntk_.is_xor3( n ) )
      {
        detail::on_xor3( node_lit, child_lits[0], child_lits[1], child_lits[2], _add_clause );
        return;
      }
    }

    detail::on_function( node_lit, child_lits, ntk_.node_function( n ), _add_clause );
  }

  void on_modified( node const& n )
  {
    on_delete( n );
    on_add( n, false ); 
    /* reuse literals_[n] (so that the fanout clauses are still valid),
    but create a new switches_[n] to control a new set of gate clauses */
  }

  void on_delete( node const& n )
  {
    if ( !ps_.allow_modify )
    {
      assert( false && "dynamic_cnf_view_params.allow_modify is not turned on" );
      std::abort();
    }

    deactivate( n );
  }

private:
  Ntk const& ntk_;
  dynamic_cnf_view_params ps_;
  percy::bsat_wrapper solver_;
  node_map<uint32_t, Ntk> literals_;
  std::vector<uint32_t> switches_;
};

} /* namespace mockturtle */
