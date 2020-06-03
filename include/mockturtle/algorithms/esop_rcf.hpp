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

  uint32_t cost( std::vector<uint32_t> const& costs ) const
  {
    return costs[x2];
  }

  std::vector<kitty::cube> esop_cubes()
  {
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

  void print( uint32_t block_size ) const
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
      res ^= ( var.x2 == offset );
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

  std::vector<uint32_t> vids; /* indices of `g_var`s in `vars` */
  uint32_t x2; /* the `x2`s of all `g_var`s XORed together should equal this value */
};

class esop_synthesis_impl
{
public:
  esop_synthesis_impl( kitty::dynamic_truth_table const& func, esop_synthesis_params const& ps )
    : func( func ), n( func.num_vars() ), r( ps.r ), block_size( 1u << r ), best( ps.best ), ps( ps )
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

public:
  std::vector<kitty::cube> run()
  {
    create_covering_variables();
    construct_RCF();
    write_maxixor();
    //solve_RCF();

    /*std::cout << esops.size() << " min ESOPs of cost " << best << " found.\n";
    for ( auto& esop : esops )
    {
      for ( auto& c : esop )
      {
        c.print( n );
        std::cout << " ";
      }
      std::cout << "\n";
    }*/

    return std::vector<kitty::cube>();//esops[0];
  }

private:
  void create_covering_variables()
  {
    vars.reserve( pow( 3, n - r ) );

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

    assert( vars.size() == pow( 3, n - r ) );
    //for ( auto g : vars ) { g.x1.print( n - r ); std::cout<<" "; } std::cout<<"\n";
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
    //for ( auto& c : RCF_ ) { c.print( block_size ); }

    /* a complete process: try every 2^(2^r) assignment of `x2` to every `g_var` */
    naive_solve_rec( 0u, 0u );
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
      //std::cout<<"new solution found with cost "<<cost<<"\n";
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

    for ( auto x2 = 0u; x2 < ( 1u <<block_size ); ++x2 )
    {
      vars[vid].x2 = x2;
      naive_solve_rec( vid + 1, cost + costs[x2] );
    }
  }

  /* translate into ESOP expression based on current g_var assignments */
  std::vector<kitty::cube> make_esop()
  {
    assert( valid() );
    std::vector<kitty::cube> esop;

    /*for ( auto& g : vars )
    {
      g.x1.print(n-r);
      std::cout<<" "<<g.x2<<"\n";
    }*/

    for ( auto& g : vars )
    {
      auto cubes = g.esop_cubes();
      esop.insert( esop.end(), cubes.begin(), cubes.end() );
    }

    return esop;
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

private:
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

  bool valid()
  {
    for ( auto& c : RCF_ )
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
  uint32_t const block_size; // 2^r

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
