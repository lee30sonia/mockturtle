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
#include <lorina/verilog.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/pattern_generation.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/networks/aig.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, float, float, float, float, bool> exp( "sim_resubstitution", "benchmark", "#PI", "size", "gain", "#pat", "#cex", "#cex0", "#cex1", "#div0", "#div1", "t_patgen", "t_resub", "t_sim", "t_SAT", "cec" );

  //for ( auto const& benchmark : epfl_benchmarks( ~hyp & ~mem_ctrl & ~experiments::log2 & ~experiments::div & ~experiments::sqrt & ~multiplier ) )
  for ( auto const& benchmark : iwls_benchmarks() )
  {
    //if ( benchmark != "iwls2005/ac97_ctrl" ) continue;

    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) );

    simresub_params ps;
    simresub_stats st;

    ps.max_pis = 10u; //100u; //8u;
    ps.max_divisors = 200u;
    ps.max_inserts = 1u;
    ps.progress = false;
    //ps.odc_levels = 5;
    ps.check_const = true;

    bool useExternal = true;
    auto pat_path = "1024sa1/"; // "patABC/" "patgen/" "patCEX/" "stuck_at_10/" "stuck_at_10_obs/" 
    //ps.write_pats = "patCEX/" + benchmark + ".pat";

    patgen_stats st_pat;
    partial_simulator sim(1,1);
    if ( useExternal )
    {
      //sim = partial_simulator( pat_path + benchmark + ".pat", "rand/" + benchmark + ".pat", 4096 );
      sim = partial_simulator( pat_path + benchmark + ".pat" );
      st_pat.num_total_patterns = sim.compute_constant( false ).num_bits();
    }
    else
    {
      patgen_params ps_pat;
      ps_pat.random_seed = 1689;
      ps_pat.num_random_pattern = 256;
      ps_pat.num_stuck_at = 1;
      //ps_pat.distinguish_nodes = true;
      //ps_pat.observability_type1 = true;
      //ps_pat.observability_type2 = true;
      //ps_pat.write_pats = "sa5/" + benchmark + ".pat";
      //ps_pat.patfile = "rand/" + benchmark + ".pat";
      sim = pattern_generation( aig, ps_pat, &st_pat );
      aig = cleanup_dangling( aig );
    }

    const uint32_t size0 = aig.num_gates();
    sim_resubstitution( aig, sim, ps, &st );
    aig = cleanup_dangling( aig );

    const auto cec = abc_cec( aig, benchmark );
    exp( benchmark, aig.num_pis(), size0, size0 - aig.num_gates(), st_pat.num_total_patterns, st.num_cex, st.num_cex_div0, st.num_cex_div1, st.num_div0_accepts, st.num_div1_accepts, to_seconds( st_pat.time_total ), to_seconds( st.time_total ), to_seconds( st.time_sim ), to_seconds( st.time_sat ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}
