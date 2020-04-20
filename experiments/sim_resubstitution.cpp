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
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/pattern_generation.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, float, float, float, float, bool> exp( "sim_resubstitution", "benchmark", "#PI", "size", "gain", "#pat", "#cex", "#const", "#div0", "#div1", "t_patgen", "t_resub", "t_sim", "t_SAT", "cec" );

  for ( auto const& benchmark : epfl_benchmarks( ~hyp & ~mem_ctrl & ~experiments::log2 & ~experiments::div & ~experiments::sqrt ) )
  {
    //if ( benchmark != "multiplier" ) continue;

    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig, orig;
    lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) );
    lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( orig ) );

    simresub_params ps;
    simresub_stats st;

    ps.max_pis = 10u; //100u; //8u;
    ps.max_divisors = 500u;
    ps.max_inserts = 1u;
    ps.progress = false;
    ps.use_odc = true;
    ps.odc_solve_limit = 10u;

    bool useExternal = true;
    auto pat_path = "patCEX/"; // "patABC/" "patgen/" "patCEX/" "stuck_at_10/" "stuck_at_10_obs/" 
    //ps.write_pats = "patCEX/" + benchmark + ".pat";

    patgen_stats st_pat;
    partial_simulator sim(1,1);
    if ( useExternal )
    {
      sim = partial_simulator( pat_path + benchmark + ".pat" );
      st_pat.num_total_patterns = sim.compute_constant( false ).num_bits();
    }
    else
    {
      sim = pattern_generation( aig, {.random_seed = 1689, .num_random_pattern = 1000, .num_stuck_at = 10, .observability_type1 = true, .write_pats = "stuck_at_10_obs/" + benchmark + ".pat"}, &st_pat );
      aig = cleanup_dangling( aig );
    }

    sim_resubstitution( aig, sim, ps, &st );
    aig = cleanup_dangling( aig );

    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    //std::cout << "num_total_divisors = " << st.num_total_divisors << std::endl;
    exp( benchmark, aig.num_pis(), orig.num_gates(), orig.num_gates() - aig.num_gates(), st_pat.num_total_patterns, st.num_cex, st_pat.num_constant, st.num_div0_accepts, st.num_div1_accepts, to_seconds( st_pat.time_total ), to_seconds( st.time_total ), to_seconds( st.time_sim ), to_seconds( st.time_sat ), cec );
  }

  exp.save();
  exp.table();
  //exp.compare();

  return 0;
}
