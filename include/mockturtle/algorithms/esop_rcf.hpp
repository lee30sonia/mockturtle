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

#include <kitty/cube.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/bit_operations.hpp>

#include "../traits.hpp"

namespace mockturtle
{

struct esop_synthesis_params
{
  uint32_t r{2};
};

namespace detail
{
struct g_var
{
  g_var( uint32_t const& n, uint32_t const& r, kitty::cube const& x1 )
    : n( n ), r( r ), x1( x1 )
  { }

  uint32_t const& n;
  uint32_t const& r;

  kitty::cube x1; /* (n-r) variables {0, 1, 2} */
  //kitty::cube x2; /* r variables {0, 1} */
};

struct xor_constraint
{
  xor_constraint( std::vector<uint32_t>& vars, bool even )
    : vars( vars ), even( even )
  {}

  std::vector<uint32_t> vars;
  bool even; /* 1: make even number of vars 1, 0: make odd number of vars 1 */
};

class esop_synthesis_impl
{
public:
  esop_synthesis_impl( kitty::dynamic_truth_table const& func, esop_synthesis_params const& ps )
    : func( func ), n( func.num_vars() ), r( ps.r ), ps( ps )
  {
    assert( r <= 4u && r <= n && "r parameter too large" );
  }

public:
  std::vector<kitty::cube> run()
  {
    create_covering_variables();
    construct_RCF();

    return esop;
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
    for ( auto g : vars )
    {
      g.x1.print( n - r );
      std::cout<<" ";
    }
    std::cout<<"\n";
  }

  void construct_RCF()
  {
    auto block_size = 1u << r;
    for ( auto i = 0u; i < ( 1u << (n-r) ); ++i )
    {
      std::vector<uint32_t> v;
      for ( auto j = 0u; j < vars.size(); ++j )
      {
        if ( contain( vars[j].x1, i ) )
        {
          v.emplace_back( j );
        }
      }
      std::cout<<"bits "<< i*block_size <<" ~ "<<(i+1)*block_size<<": ";
      for (auto varId : v){ vars[varId].x1.print(n-r); std::cout<<" ";} std::cout<<"\n";
      for ( auto bits2 = 0u; bits2 < block_size; ++bits2 )
      {
        std::vector<uint32_t> ids;
        ids.resize( v.size() );
        for ( auto j = 0u; j < v.size(); ++j )
        {
          ids[j] = v[j] * block_size + bits2;
        }
        RCF.emplace_back( xor_constraint( ids, !kitty::get_bit( func, i * block_size + bits2 ) ) );
      }
    }
  }

private:
  /* check if minterm `m` is contained in cube `c` */
  bool contain( kitty::cube const& c, uint32_t const& m )
  {
    return ( ( c._bits ^ m ) & (c._mask) ) == 0u;
  }

private:
  kitty::dynamic_truth_table const& func;
  std::vector<kitty::cube> esop;
  std::vector<g_var> vars;
  std::vector<xor_constraint> RCF;

  uint32_t const n;
  uint32_t const r;

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
