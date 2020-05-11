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
  \file abc_resub.hpp
  \brief Interface to `abc_resub`.

  \author Heinz Riener
*/

#pragma once

#include <kitty/kitty.hpp>
#include <abcresub/abcresub.hpp>

namespace mockturtle
{

class abc_resub
{
public:
  explicit abc_resub( uint64_t num_divisors, uint64_t num_blocks_per_truth_table )
    : num_divisors( num_divisors )
    , num_blocks_per_truth_table( num_blocks_per_truth_table )
    , counter(0)
  {
    alloc();
  }

  virtual ~abc_resub()
  {
    release();
  }

  template<class node_type, class truth_table_storage_type>
  void add_divisor( node_type const& node, truth_table_storage_type const& tts, bool complement = false )
  {
    assert( abc_tts != nullptr && "assume that memory for truth tables has been allocated" );
    assert( abc_divs != nullptr && "assume that memory for divisors has been allocated" );

    assert( tts[node].num_blocks() == num_blocks_per_truth_table );
    for ( uint64_t i = 0ul; i < num_blocks_per_truth_table; ++i )
    {
      auto const tt = complement ? ~tts[node] : tts[node];
      Vec_WrdPush( abc_tts, tt._bits[i] );
      Vec_PtrPush( abc_divs, Vec_WrdEntryP( abc_tts, counter * num_blocks_per_truth_table + i ) );
    }
    ++counter;
  }

  template<class iterator_type, class truth_table_storage_type>
  void add_divisors( iterator_type begin, iterator_type end, truth_table_storage_type const& tts )
  {
    assert( abc_tts != nullptr && "assume that memory for truth tables has been allocated" );
    assert( abc_divs != nullptr && "assume that memory for divisors has been allocated" );

    while ( begin != end )
    {
      add_divisor( *begin, tts );
      ++begin;
    }
  }

  void compute_function( abcresub::Gia_ResbMan_t * p )
  {
    abcresub::Gia_ManResubPerform( p, abc_divs, num_blocks_per_truth_table, 100, /* fDebug = */1, /* fVerbose = */1 );
  }

protected:
  void alloc()
  {
    assert( abc_tts == nullptr );
    assert( abc_divs == nullptr );
    abc_tts = abcresub::Vec_WrdAlloc( num_divisors * num_blocks_per_truth_table );
    abc_divs = abcresub::Vec_PtrAlloc( num_divisors );
  }

  void release()
  {
    assert( abc_divs != nullptr );
    assert( abc_tts != nullptr );
    Vec_PtrFree( abc_divs );
    Vec_WrdFree( abc_tts );
    abc_divs = nullptr;
    abc_tts = nullptr;
  }

protected:
  uint64_t num_divisors;
  uint64_t num_blocks_per_truth_table;
  uint64_t counter;

  abcresub::Vec_Wrd_t * abc_tts{nullptr};
  abcresub::Vec_Ptr_t * abc_divs{nullptr};
}; /* abc_resub */

} /* namespace mockturtle */
