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

#include <vector>

#include "kitty/dynamic_truth_table.hpp"
#include "kitty/constructors.hpp"
#include "kitty/print.hpp"
#include "kitty/esop.hpp"
#include "mockturtle/algorithms/esop_rcf.hpp"
#include "mockturtle/utils/stopwatch.hpp"

int main()
{
  using namespace mockturtle;

  kitty::dynamic_truth_table tt( 6u );
  kitty::create_from_hex_string( tt, "688c802028222222" );

  esop_synthesis_params ps;
  ps.r = 2u;
  //ps.best = 7u;
  stopwatch<>::duration runtime(0);

  esop_synthesis( tt, ps );
  call_with_stopwatch( runtime, [&]() {
    system("./maxixor dump.txt");
  });
  std::cout << "runtime: " << to_seconds( runtime ) << " sec.\n";

  return 0;


  for ( auto i = 0u; i < (1<<16); ++i )
  {
    tt._bits[0] = i;
    kitty::print_hex(tt); std::cout<<": ";
    auto esop = call_with_stopwatch( runtime, [&]() {
      //return esop_synthesis( tt, ps );
      return esop_from_pprm( tt );
    });
    std::cout << "ESOP of cost " << esop.size() << " found.\n";
  }
  std::cout << "runtime: " << to_seconds( runtime ) << " sec.\n";
  
  return 0;
}
