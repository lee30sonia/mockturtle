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

#include <string>
#include <vector>
#include <iostream>

#include <fmt/format.h>
#include <kitty/kitty.hpp>
#include <kitty/partial_truth_table.hpp>
#include "mockturtle/utils/stopwatch.hpp"

int main()
{
  using namespace mockturtle;
  using namespace kitty;

  std::vector<partial_truth_table> tts;
  for ( auto i = 0u; i < (1u<<20); ++i )
  {
    tts.emplace_back( 4096 );
  }
  uint64_t max64 = 1;
  max64 <<= 63;

  float t = 0;
  float t_fast = 0;

  for ( auto l = 0u; l < 20; ++l)
  {

    for ( auto i = 0u; i < tts.size(); ++i )
      create_random( tts.back() );
    
    stopwatch<>::duration t_count_ones( 0 );
    call_with_stopwatch( t_count_ones, [&]() {
      for ( uint64_t k = 0; k < max64; ++k )
      for ( auto j = 0u; j < (1u<<31); ++j )
      {
        for ( auto i = 0u; i < tts.size(); ++i )
          count_ones( tts[i] );
      }
    });
    std::cout << to_seconds( t_count_ones ) << " ";
    t += to_seconds( t_count_ones );


    for ( auto i = 0u; i < tts.size(); ++i )
      create_random( tts.back() );

    stopwatch<>::duration t_count_ones_fast( 0 );
    call_with_stopwatch( t_count_ones_fast, [&]() {
      for ( uint64_t k = 0; k < max64; ++k )
      for ( auto j = 0u; j < (1u<<31); ++j )
      {
        for ( auto i = 0u; i < tts.size(); ++i )
          count_ones_fast( tts[i] );
      }
    });
    std::cout << to_seconds( t_count_ones_fast ) << "\n";
    t_fast += to_seconds( t_count_ones_fast );

  }

  std::cout << t_fast << " " << t << "\n";
  


  return 0;
}
