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

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/pattern_generation.hpp>
//#include <mockturtle/algorithms/pattern_generation_cnfview.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, float, float, float, bool> exp( "pattern_generation", "benchmark", "#PI", "size", "#pat", "#pat gen", "#const", "t_total", "t_sim", "t_SAT", "cec" );

  //for ( auto const& benchmark : epfl_benchmarks( ~hyp & ~mem_ctrl & ~experiments::log2 & ~experiments::div & ~experiments::sqrt & ~multiplier ) )
  for ( auto const& benchmark : iwls_benchmarks() )
  {
    //if ( benchmark != "iwls2005/mem_ctrl" ) continue;

    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) );
    auto size_before = aig.num_gates();

    patgen_params ps;
    patgen_stats st;

    ps.num_random_pattern = 256;
    ps.observability_type1 = true;
    ps.observability_type2 = true;
    ps.observability_levels = 5;
    ps.write_pats = "256sa1obs/" + benchmark + ".pat";
    //ps.patfile = "test.pat";
    ps.random_seed = 1689;
    ps.progress = false;

    pattern_generation( aig, ps, &st );
    aig = cleanup_dangling( aig );

    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    exp( benchmark, aig.num_pis(), size_before, st.num_total_patterns, st.num_total_patterns - ps.num_random_pattern, st.num_constant, to_seconds( st.time_total ), to_seconds( st.time_sim ), to_seconds( st.time_sat ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}
