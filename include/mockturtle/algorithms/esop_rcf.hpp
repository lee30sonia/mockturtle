/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2020  EPFL
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
  \file esop_rcf.hpp
  \brief ESOP exact synthesis with Reduced Covering Functions

  \author Siang-Yun Lee
*/

#pragma once

#include <vector>
#include <cmath>
#include <fstream>

#include <kitty/cube.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/bit_operations.hpp>

#include <bill/dd/zdd.hpp>
#include <fmt/format.h>

#include "../traits.hpp"
#include "../utils/stopwatch.hpp"

namespace mockturtle
{

struct esop_synthesis_params
{
  uint32_t best{INT_MAX};
  uint32_t r{2};
};

namespace detail
{
/* a `g_var` object is 2^r covering variables in the RCF, all corresponding to a cube `x1` on X1 PIs
   assignment to the 2^r covering variables (corresponding to different value combinations of X2 PIs) decides the cost of a `g_var`
   hence data member `x2` stores the current assignment to the 2^r covering variables of this `g_var` */
struct g_var
{
  g_var( uint32_t const& n, uint32_t const& r, kitty::cube const& x1 )
    : n( n ), r( r ), x1( x1 ), x2( 0u )
  { }
  
  g_var& operator=( g_var other )
  {
    (void)other;
    return *this;
  }

  void reset_x2() { x2 = 0u; }

  void set_x2( uint32_t const& offset ) { x2 |= 1 << offset; }

  uint32_t cost( std::vector<uint32_t> const& costs ) const
  {
    return costs[x2];
  }

  std::vector<kitty::cube> esop_cubes()
  {
    /* variable order is incorrect */
    if ( x2 == 0u )
    {
      return std::vector<kitty::cube>();
    }
    
    auto c = x1;
    auto c2 = x1;
    switch( r )
    {
      case 0u:
        return std::vector<kitty::cube>( {x1} );
      case 1u:
      {
        switch( x2 )
        {
          case 3u: return std::vector<kitty::cube>( {x1} );
          
          case 1u: 
            c.add_literal( n-r, false );
            return std::vector<kitty::cube>( {c} );
          case 2u: 
            c.add_literal( n-r, true );
            return std::vector<kitty::cube>( {c} );
          default:
            assert( false );
        }
      }
      case 2u:
      {
        switch( x2 )
        {
          case 15u: return std::vector<kitty::cube>( {x1} );
          
          case 1u: // !x!y
            c.add_literal( n-r, false );
            c.add_literal( n-r+1, false );
            return std::vector<kitty::cube>( {c} );
          case 2u: // !xy
            c.add_literal( n-r, false );
            c.add_literal( n-r+1, true );
            return std::vector<kitty::cube>( {c} );
          case 3u: // !x
            c.add_literal( n-r, false );
            return std::vector<kitty::cube>( {c} );
          case 4u: // x!y
            c.add_literal( n-r, true );
            c.add_literal( n-r+1, false );
            return std::vector<kitty::cube>( {c} );
          case 5u: // !y
            c.add_literal( n-r+1, false );
            return std::vector<kitty::cube>( {c} );
          case 8u: // xy
            c.add_literal( n-r, true );
            c.add_literal( n-r+1, true );
            return std::vector<kitty::cube>( {c} );
          case 10u: // y
            c.add_literal( n-r+1, true );
            return std::vector<kitty::cube>( {c} );
          case 12u: // x
            c.add_literal( n-r, true );
            return std::vector<kitty::cube>( {c} );

          case 6u: // x^y = x ^ y
            c.add_literal( n-r, true );
            c2.add_literal( n-r+1, true );
            return std::vector<kitty::cube>( {c, c2} );
          case 7u: // !x+!y = 1 ^ xy
            c2.add_literal( n-r, true );
            c2.add_literal( n-r+1, true );
            return std::vector<kitty::cube>( {c, c2} );
          case 9u: // x^!y = x ^ !y
            c.add_literal( n-r, true );
            c2.add_literal( n-r+1, false );
            return std::vector<kitty::cube>( {c, c2} );
          case 11u: // !x+y = 1 ^ x!y
            c2.add_literal( n-r, true );
            c2.add_literal( n-r+1, false );
            return std::vector<kitty::cube>( {c, c2} );
          case 13u: // x+!y = 1 ^ !xy
            c2.add_literal( n-r, false );
            c2.add_literal( n-r+1, true );
            return std::vector<kitty::cube>( {c, c2} );
          case 14u: // x+y = 1 ^ !x!y
            c2.add_literal( n-r, false );
            c2.add_literal( n-r+1, false );
            return std::vector<kitty::cube>( {c, c2} );
          default:
            assert( false );
        }
      }
      default:
        assert( false );
    }

    assert( false );
    return std::vector<kitty::cube>();
  }

  uint32_t const& n;
  uint32_t const& r;

  kitty::cube const x1; /* (n-r) variables {0, 1, 2} */
  uint32_t x2; /* assignment of {0, 1} to 2^r variables (use 2^r bits, max value 2^2^r-1) */
};

struct xor_constraint
{
  xor_constraint( std::vector<uint32_t>& vids, uint32_t offset, bool even )
    : vids( vids ), offset( offset ), even( even )
  {}

  uint32_t zdd_node( bill::zdd_base& zdd, uint32_t const& block_size )
  {
    //std::cout<<"building ZDD for: "; print( block_size );

    uint32_t every = zdd.bottom();
    for ( auto& vid : vids )
    {
      every = zdd.union_( every, zdd.elementary( vid * block_size + offset ) );
    }
    uint32_t rest = zdd.nonsupersets( zdd.tautology(), every );

    uint32_t k = even ? 0 : 1;
    uint32_t zid = zdd.bottom();
    while ( k <= vids.size() )
    {
      zid = zdd.union_( zid, zdd.choose( every, k ) );
      k += 2;
    }
    zid = zdd.join( zid, rest );

    //zdd.print_sets( zid );
    return zid;
  }

  void print( uint32_t const& block_size ) const
  {
    for ( auto& vid : vids )
    {
      std::cout << vid * block_size + offset << " ^ ";
    }
    std::cout << "\b\b= " << !even << "\n";
  }

  bool valid( std::vector<g_var> const& vars ) const
  {
    bool res = even;
    for ( auto& vid : vids )
    {
      auto const& var = vars[vid];
      res ^= bool( var.x2 & ( 1 << offset ) );
    }
    return res;
  }

  std::vector<uint32_t> vids; /* indices of `g_var`s in `vars` */
  uint32_t offset; /* offset to all the variables (value of the assignment to X2 variables) */
  bool even; /* 1: make even number of vars 1, 0: make odd number of vars 1 */
};

/* one object that encodes 2^r constraints together */
struct xor_constraints_block
{
  xor_constraints_block( std::vector<uint32_t>& vids, uint32_t x2 )
    : vids( vids ), x2( x2 )
  {}

  void print( uint32_t block_size ) const
  {
    for ( auto i = 0u; i < block_size; ++i )
    {
      for ( auto& vid : vids )
      {
        std::cout << vid * block_size + i << " ^ ";
      }
      std::cout << "\b\b= " << ((x2 >> i) & 0x1) << "\n";
    }
  }

  void print_block() const
  {
    for ( auto& vid : vids )
    {
      std::cout << vid << " ^ ";
    }
    std::cout << "\b\b= " << x2 << "\n";
  }

  bool valid( std::vector<g_var> const& vars ) const
  {
    uint32_t tmp = 0u;
    for ( auto& vid : vids )
    {
      auto const& var = vars[vid];
      tmp ^= var.x2;
    }
    return tmp == x2;
  }

  /* given the assignments to all vars except for the last one, conclude on the x2 of the last g_var */
  uint32_t last_var_x2( std::vector<g_var> const& vars ) const
  {
    uint32_t tmp = 0u;
    for ( auto i = 0u; i < vids.size() - 1u; ++i )
    {
      auto const& var = vars[vids[i]];
      tmp ^= var.x2;
    }
    return tmp ^ x2;
  }

  /* given the assignments to `g_var`s from `vars[id]` to the end, calculate the right-hand-side value for the remaining equation */
  uint32_t rhs_current( std::vector<g_var> const& vars, uint32_t const id ) const
  {
    uint32_t tmp = x2;
    for ( auto i = 0u; i < vids.size(); ++i )
    {
      if ( vids[i] < id )
      {
        //std::cout << vids[i] << " ^ ";
        continue;
      }
      auto const& var = vars[vids[i]];
      tmp ^= var.x2;
      //std::cout << vids[i] << "(" << var.x2 << ") ^ ";
    }
    //std::cout << "\b\b= " << tmp << "\n";
    return tmp;
  }

  std::vector<uint32_t> vids; /* indices of `g_var`s in `vars` */
  uint32_t x2; /* the `x2`s of all `g_var`s XORed together should equal this value */
};

class esop_synthesis_impl
{
public:
  esop_synthesis_impl( kitty::dynamic_truth_table const& func, esop_synthesis_params const& ps )
    : func( func ), n( func.num_vars() ), r( ps.r ), block_size( 1u << r ),
      zdd( pow( 2, r ) * pow( 3, n - r ) ), best( ps.best ), ps( ps )
  {
    assert( r <= 2u && r <= n && "r parameter too large" );
    switch( r )
    {
      case 0u:
        costs = std::vector<uint32_t>( {0,1} );
        break;
      case 1u:
        costs = std::vector<uint32_t>( {0,1,1,1} );
        break;
      case 2u:
        costs = std::vector<uint32_t>( {0,1,1,1, 1,1,2,2, 1,2,1,2, 1,2,2,1} );
        break;
      default:
        assert( false );
    }
  }

  std::vector<kitty::cube> run()
  {
    create_covering_variables();
    construct_RCF();
    //write_maxixor();
    solve_RCF();

    std::cout << esops.size() << " min ESOPs of cost " << best << " found.\n";
    //for ( auto& esop : esops )
    //{
    //  for ( auto& c : esop )
    //  {
    //    c.print( n );
    //    std::cout << " ";
    //  }
    //  std::cout << "\n";
    //}

    return std::vector<kitty::cube>();//esops[0];
  }

private:
  void create_covering_variables()
  {
    vars.reserve( pow( 3, n - r ) );

    //for ( int mask1 = ( 1u << (n-r) ) -1; mask1 >=0 ; --mask1 ) 
    for ( auto mask1 = 0u; mask1 < ( 1u << (n-r) ); ++mask1 )
    {
      for ( auto bits1 = 0u; bits1 < ( 1u << (n-r) ); ++bits1 )
      {
        if ( ~mask1 & bits1 )
        {
          continue;
        }
        vars.emplace_back( n, r, kitty::cube( bits1, mask1 ) );
      }
    }
    std::reverse( vars.begin(), vars.end() );
    assert( vars.size() == pow( 3, n - r ) );
  }

  void construct_RCF()
  {
    for ( auto i = 0u; i < ( 1u << (n-r) ); ++i )
    {
      std::vector<uint32_t> v;
      for ( auto vid = 0u; vid < vars.size(); ++vid )
      {
        if ( contain( vars[vid].x1, i ) )
        {
          v.emplace_back( vid );
        }
      }

      uint32_t x2 = 0u;
      for ( auto bits2 = 0u; bits2 < block_size; ++bits2 )
      {
        RCF.emplace_back( xor_constraint( v, bits2, !kitty::get_bit( func, i * block_size + bits2 ) ) );
        if ( kitty::get_bit( func, i * block_size + bits2 ) )
        {
          x2 |= 0x1 << bits2;
        }
      }
      RCF_.emplace_back( xor_constraints_block( v, x2 ) );
    }
  }

  void solve_RCF()
  {
    /*for ( auto i = 0u; i < vars.size(); ++i )
    {
      std::cout<<i<<" ";
      vars[i].x1.print(n-r);
      std::cout<<"\n";
    }*/
    //for ( auto& c : RCF_ ) { c.print_block(); }
    //for ( auto& c : RCF ) { c.print( block_size ); }

    //write_maxixor();

    /* a complete process: try every 2^(2^r) assignment of `x2` to every `g_var` */
    //naive_solve_rec( 0u, 0u );
    //n_r = n - r; solve_phase1( vars.size() - 1u, 0u );

    build_ZDD();
  }

private: /* ZDD */
  void build_ZDD()
  {
    stopwatch<>::duration time_U(0);
    uint32_t zid = call_with_stopwatch( time_U, [&]() {
      return upper_bound();
    });
    std::cout<<"upper bound function built in " << to_seconds(time_U) << " sec.\n";
    for ( auto& c : RCF ) 
    { 
      zid = zdd.intersection( zid, c.zdd_node( zdd, block_size ) ); 
    }
    std::cout<<"constraints built\n";
    //zdd.print_sets( zid );
    zdd.foreach_set( zid, [&]( auto const& set ){
      //std::cout << fmt::format("{{ {} }}", fmt::join(set, ", "));
      for ( auto& v : vars )
      {
        v.reset_x2();
      }
      for ( auto& z : set )
      {
        vars.at( int( z / block_size ) ).set_x2( z % block_size );
      }
      auto cost = get_cost();
      if ( cost < best )
      {
        best = cost;
        esops.clear();
      }
      if ( cost <= best )
      {
        esops.emplace_back( make_esop() );
      }
      //std::cout << " -- " << get_cost() << ( valid() ? " (valid)" : " (invalid!!)" ) << "\n";
      return true;
    });
  }

  uint32_t upper_bound()
  {
    if ( ps.best == INT_MAX )
    {
      return zdd.tautology();
    }

    if ( r > 1 ) /* build tautology nodes for blocks of ZDD variables in each g_var */
    {
      for ( auto v = 0u; v < vars.size() - 1; ++v )
      {
        uint32_t zid = zdd.top();
        for ( int i = block_size - 1; i >= 0; --i )
        {
          zdd.ref( zid, 2 );
          zid = zdd.unique( v * block_size + i, zid, zid );
        }
        tautologies.push_back( zid );
        //zdd.print_sets( zid ); std::cout<<"-----\n";
      }
      tautologies.push_back( zdd.tautology( ( vars.size() - 1 ) * block_size ) );
    }

    uint32_t m = r == 2 ? 2 : 1; /* max. cost of a g_var */
    return upper_bound_rec( vars.size() - 1, m, ps.best );
  }

  uint32_t upper_bound_rec( uint32_t v, uint32_t const& m, uint32_t k )
  {
    if ( v == 0 )
    {
      return upper_bound_single( 0, k );
    }
    uint32_t zid = zdd.bottom();
    for ( auto i = 0u; i <= m; ++i )
    {
      if ( k < i ) break;
      zid = zdd.union_( zid, zdd.join( upper_bound_single( v, i ), upper_bound_rec( v - 1, m, k - i ) ) );
    }
    return zid;
  }

  uint32_t upper_bound_single( uint32_t v, uint32_t k )
  {
    uint32_t zid;
    uint32_t b = v * block_size;
    switch( r )
    {
      case 0u: // ~g0 = {{}}, 1 = {{}, {g0}}
        return ( k == 0u ) ? zdd.top() : zdd.union_( zdd.top(), zdd.elementary( v ) );
      case 1u: // ~g0 & ~g1 = {{}}, 1 = {{}, {g0}, {g1}, {g0, g1}}
        return ( k == 0u ) ? zdd.top() : tautologies.at( v );
      case 2u:
        if ( k == 0u ) return zdd.top();
        else if ( k >= 2u ) return tautologies.at( v );
        else // k == 1u
        {
          zid = zdd.join( zdd.elementary( b + 1 ), zdd.elementary( b + 2 ) ); // {g1, g2}
          zid = zdd.union_( zid, zdd.join( zdd.join( zdd.elementary( b ), zdd.elementary( b + 1 ) ), zdd.elementary( b + 2 ) ) ); // {g0, g1, g2}
          zid = zdd.union_( zid, zdd.join( zdd.elementary( b ), zdd.elementary( b + 3 ) ) ); // {g0, g3}
          zid = zdd.union_( zid, zdd.join( zdd.join( zdd.elementary( b ), zdd.elementary( b + 1 ) ), zdd.elementary( b + 3 ) ) ); // {g0, g1, g3}
          zid = zdd.union_( zid, zdd.join( zdd.join( zdd.elementary( b ), zdd.elementary( b + 2 ) ), zdd.elementary( b + 3 ) ) ); // {g0, g2, g3}
          zid = zdd.union_( zid, zdd.join( zdd.join( zdd.elementary( b + 1 ), zdd.elementary( b + 2 ) ), zdd.elementary( b + 3 ) ) ); // {g1, g2, g3}
          return zdd.difference( tautologies.at( v ), zid );
        }
      default: assert( false ); return 0u;
    }
  }

private: /* direct enumeration algorithms & maxixor write-out */
  /* try all possible assignments for g_vars (3^n_r - 1) to (3^n_r - 2^n_r) */
  /* first phase of a subproblem with n - r = n_r */
  bool solve_phase1( uint32_t vid, uint32_t cost )
  {
    if ( cost >= best ) /* use > if want all optimal solutions */
    {
      return true;
    }

    /* except for the first phase1, we need to check the other half of identical constraints */
    /* if the paired constraints have different values, conclude UNSAT */
    bool check_pairs = (n - r != n_r);

    if ( check_pairs )
    {
      /* the indexes of the two associated constraint blocks */
      auto bid1 = ( 1u << n_r ) - ( pow( 3, n_r ) - vid );
      // auto& bid2 = bid1 + ( 1u << n_r );
      auto x2 = RCF_[bid1].rhs_current( vars, vid + 1 );
      if ( x2 != RCF_[bid1 + ( 1u << n_r )].rhs_current( vars, vid + 1 ) )
      {
        /* UNSAT */
        //std::cout<<"UNSAT at vid "<<vid<<"\n";
        return false;
      }
    }

    bool is_last = ( vid == pow( 3, n_r ) - (1 << n_r) );

    for ( auto x2 = 0u; x2 < ( 1u << block_size ); ++x2 )
    {
      vars[vid].x2 = x2;
      //std::cout<<"assign "<<x2<<" to var "<<vid<<"\n";
      if ( is_last )
      {
        solve_phase2( vid - 1, cost + costs[x2] );
      }
      else if ( !solve_phase1( vid - 1, cost + costs[x2] ) )
      {
        return false;
      }
    }
    return true;
  }

  /* try all possible assignments for g_vars (3^n_r - 2^n_r - 1) to 3^(n_r - 1) */
  /* in the end, enter again phase1 for a smaller subproblem */
  /* in phase2, we never need to check for paired constraints */
  void solve_phase2( uint32_t vid, uint32_t cost )
  {
    if ( cost >= best ) /* use > if want all optimal solutions */
    {
      return;
    }

    /* termial case */
    if ( vid == 0u )
    {
      assert( n_r == 1u );
      auto x2 = RCF_[0].rhs_current( vars, 1 );
      if ( x2 != RCF_[1].rhs_current( vars, 1 ) )
      {
        /* UNSAT */
        return;
      }

      vars[vid].x2 = x2;
      //std::cout<<"assign "<<x2<<" to var "<<vid<<"\n";

      assert( cost + costs[x2] == get_cost() );
      if ( cost + costs[x2] < best ) /* use <= if want all optimal solutions */
      {
        esops.clear();
        best = cost + costs[x2];
        std::cout<<"new solution found with cost "<<cost + costs[x2]<<"\n";
        assert( valid() );
        esops.emplace_back( make_esop() );
      }
      return;
    }

    bool is_last = ( vid == pow( 3, n_r - 1 ) );

    for ( auto x2 = 0u; x2 < ( 1u << block_size ); ++x2 )
    {
      vars[vid].x2 = x2;
      //std::cout<<"assign "<<x2<<" to var "<<vid<<"\n";
      if ( is_last )
      {
        --n_r;
        solve_phase1( vid - 1, cost + costs[x2] );
        ++n_r;
      }
      else
      {
        solve_phase2( vid - 1, cost + costs[x2] );
      }
    }
  }
  
  /* vid: index of the currently deciding g_var. */
  /* cost: current accumulated cost. */
  void naive_solve_rec( uint32_t vid, uint32_t cost )
  {
    if ( cost >= best )
    {
      return;
    }

    /* all g_vars are assigned */
    if ( vid == vars.size() )
    {
      if ( cost < best )
      {
        esops.clear();
        best = cost;
      }
      assert( cost == get_cost() );
      std::cout<<"new solution found with cost "<<cost<<"\n";
      esops.emplace_back( make_esop() );
      return;
    }

    /* the last 2^(n-r) g_vars only appear in one block of constraints
       so at this point this is the last un-assigned g_var in a block
       we can decide its x2 directly */
    if ( vid >= vars.size() - ( 1u << (n-r) ) )
    {
      /* the associated constraint block */
      auto& block = RCF_[vid - ( vars.size() - ( 1u << (n-r) ) )];
      auto x2 = block.last_var_x2( vars );
      vars[vid].x2 = x2;
      assert( block.valid( vars ) );
      naive_solve_rec( vid + 1, cost + costs[x2] );
      return;
    }

    for ( auto x2 = 0u; x2 < ( 1u << block_size ); ++x2 )
    {
      vars[vid].x2 = x2;
      naive_solve_rec( vid + 1, cost + costs[x2] );
    }
  }

  void write_maxixor( const std::string file_name = "dump.txt" )
  {
    std::ofstream os( file_name.c_str(), std::ofstream::out );
    for ( auto& c : RCF_ ) 
    {
      for ( auto i = 0u; i < block_size; ++i )
      {
        os << ((c.x2 >> i) & 0x1) << ": ";
        for ( auto& vid : c.vids )
        {
          os << vid * block_size + i << " ";
        }
        os << "\n";
      }
    }
  }

private: /* helper functions */
  /* check if minterm `m` is contained in cube `c` */
  bool contain( kitty::cube const& c, uint32_t const& m )
  {
    return ( ( c._bits ^ m ) & (c._mask) ) == 0u;
  }

  uint32_t get_cost()
  {
    uint32_t cost = 0u;
    for ( auto& g : vars )
    {
      cost += g.cost( costs );
    }
    return cost;
  }

  /* translate into ESOP expression based on current g_var assignments */
  std::vector<kitty::cube> make_esop()
  {
    assert( valid() );
    std::vector<kitty::cube> esop;

    for ( auto& g : vars )
    {
      auto cubes = g.esop_cubes();
      esop.insert( esop.end(), cubes.begin(), cubes.end() );
    }

    return esop;
  }

  bool valid()
  {
    for ( auto& c : RCF )
    {
      if ( !c.valid( vars ) )
      {
        return false;
      }
    }
    return true;
  }

private:
  kitty::dynamic_truth_table const& func;
  std::vector<std::vector<kitty::cube>> esops;
  std::vector<g_var> vars;
  std::vector<xor_constraint> RCF; // not used in naive algorithm
  std::vector<xor_constraints_block> RCF_;

  uint32_t const n;
  uint32_t const r;
  uint32_t const block_size; // == 2^r
  uint32_t n_r; // (n - r) of current subproblem

  bill::zdd_base zdd;
  std::vector<uint32_t> tautologies; /* tautology ZDD nodes for each g_var */

  std::vector<uint32_t> costs;
  uint32_t best; /* cost of the best solution so far */

  esop_synthesis_params ps;
};

} // namespace detail

/*! \brief Performs ESOP exact synthesis
 */
std::vector<kitty::cube> esop_synthesis( kitty::dynamic_truth_table const& func, esop_synthesis_params const& ps = {} )
{
  detail::esop_synthesis_impl p( func, ps );
  return p.run();
}
} // namespace mockturtle
