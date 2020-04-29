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
#include <fstream>

#include <fmt/format.h>
#include <mockturtle/algorithms/cnf.hpp>
#include <percy/solvers/bsat2.hpp>

using namespace mockturtle;

inline std::vector<std::string> split( const std::string& str, const std::string& sep )
{
  std::vector<std::string> result;

  size_t last = 0;
  size_t next = 0;
  std::string substring;
  while ( ( next = str.find( sep, last ) ) != std::string::npos )
  {
    substring = str.substr( last, next - last );
    if ( substring.length() > 0 )
    {
      std::string sub = str.substr( last, next - last );
      sub.erase( std::remove( sub.begin(), sub.end(), ' ' ), sub.end() );
      result.push_back( sub );
    }
    last = next + 1;
  }

  substring = str.substr( last );
  substring.erase( std::remove( substring.begin(), substring.end(), ' ' ), substring.end() );
  result.push_back( substring );

  return result;
}

void read_dimacs(percy::bsat_wrapper& solver, std::string file, int skip = 0)
{
  std::ifstream in( file, std::ifstream::in );
  std::string line;

  std::getline( in, line ); // header
  for (int i=0; i<skip; ++i) std::getline( in, line );
  while ( std::getline( in, line ) )
  {
    const auto tokens = split( line, " " );
    std::vector<uint32_t> clause;
    for ( auto t : tokens )
    {
      auto v = std::atoi(t.c_str());
      if ( v<0 ) clause.push_back( (-v-1) * 2 + 1 );
      else if ( v>0 ) clause.push_back( (v-1) * 2 );
    }
    if (clause.size() > 0)
      solver.add_clause( clause );
  }
}

void solve_assump(percy::bsat_wrapper& solver, int assump)
{
  std::vector<pabc::lit> assumptions( 1, (assump<0)? (-assump*2+1): (assump*2) );
  auto res = solver.solve( &assumptions[0], &assumptions[0] + 1, 0 );
  std::cout<<(( res == percy::synth_result::success )? "SAT" : "UNSAT")<<"\n";
}

int main()
{
  percy::bsat_wrapper solver;
  read_dimacs( solver, "initial.dimacs" );

  solve_assump(solver, -276);
  solver.bookmark();
  read_dimacs( solver, "first.dimacs", 3499 );
  solve_assump(solver, -276);
  solver.rollback();

  solve_assump(solver, 277);
  solver.bookmark();
  read_dimacs( solver, "second.dimacs", 3499 );
  solve_assump(solver, 277);
  solver.rollback();

  solve_assump(solver, 510);

  return 0;
}
