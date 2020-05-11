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

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, float, float, float, float, bool> exp( "sim_resubstitution", "benchmark", "#PI", "size", "gain", "#pat", "#cex", "#cex0", "#cex1", "#div0", "#div1", "t_patgen", "t_resub", "t_sim", "t_SAT", "cec" );

  //for ( auto const& benchmark : epfl_benchmarks( ~hyp & ~mem_ctrl & ~experiments::log2 & ~experiments::div & ~experiments::sqrt & ~multiplier ) )
  for ( auto const& benchmark : iwls_benchmarks() )
  {
    //if ( benchmark != "cavlc" ) continue;

    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig, orig;
    lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) );
    lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( orig ) );

    simresub_params ps;
    simresub_stats st;

    ps.max_pis = 10u;
    ps.max_divisors = 200u;
    ps.max_inserts = 1u;
    ps.progress = false;
    ps.check_const = false;

    bool useExternal = false;
    auto pat_path = "pats/"; 
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
      patgen_params ps_pat;
      ps_pat.random_seed = 1689;
      ps_pat.num_random_pattern = 1024;
      ps_pat.num_stuck_at = 1;
      ps_pat.write_pats = "pats/" + benchmark + ".pat"; 
      /* NOTE: you have to manually create directories build/pats/iwls2005/, if you want to use them later with `useExternal` */
      sim = pattern_generation( aig, ps_pat, &st_pat );
      aig = cleanup_dangling( aig );
    }

    sim_resubstitution( aig, sim, ps, &st );
    aig = cleanup_dangling( aig );

    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    exp( benchmark, aig.num_pis(), orig.num_gates(), orig.num_gates() - aig.num_gates(), st_pat.num_total_patterns, st.num_cex, st.num_cex_div0, st.num_cex_div1, st.num_div0_accepts, st.num_div1_accepts, to_seconds( st_pat.time_total ), to_seconds( st.time_total ), to_seconds( st.time_sim ), to_seconds( st.time_sat ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}
