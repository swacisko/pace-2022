//
// Created by sylwe on 08/04/2026.
//

#include "ExpConfig.h"

void ExpConfig::setParametersForAlgorithm() {
    if ( alg == HS ) {
        ihs_single_iteration_sec = 10;
        max_cycles_for_hs = 1e8;
        unhit_cycle_enumeration_type = 1;
    }
    if (alg == IHS) {
        ihs_single_iteration_sec = 3;
        unhit_cycle_enumeration_type = 1;
    }
    if (alg == MTZ) {
        ihs_iterations_in_mtz = 15;
        ihs_single_iteration_sec = 2;
        max_cycles_for_hs = 1e8;
        unhit_cycle_enumeration_type = 1;
    }
}

string parseAlgorithm(int alg)  {
    if(alg == HS) return "HS";
    if(alg == IHS) return "IHS";
    if(alg == MTZ) return "MTZ";
    if(alg == DIVERSES) return "DiVerSeS";
    if(alg == DIV_IHS) return "DiVerSeS+IHS";
    return "unknown algorithm";
}


int ExpConfig::scaleIters(int N) {
    if (max_new_cycles_iter_scale == "linear_h") return N/2;
    if (max_new_cycles_iter_scale == "linear") return N;
    if (max_new_cycles_iter_scale == "log") return N * ceil(log2(N));
    if (max_new_cycles_iter_scale == "sqrt") return N * sqrt(N);
    if ( max_new_cycles_iter_scale == "bounded" ) return max_new_cycles_per_iter;
    if ( max_new_cycles_iter_scale == "unbounded" ) return inf; // take all
    assert(false && "invalid max_new_cycles_per_iter_scale");
    return -1;
}

vector<pair<string, string>> ExpConfig::getConfigEntries() {
    vector<pair<string, string>> entries;
    entries.emplace_back("metadata_filepath",metadata_filepath);
    entries.emplace_back("alg",parseAlgorithm(alg));
    entries.emplace_back("threads",to_string(threads));
    entries.emplace_back("max_time_sec",to_string(max_time_sec));
    entries.emplace_back("ihs_single_iteration_sec",to_string(ihs_single_iteration_sec));
    entries.emplace_back("unhit_cycle_enumeration_type",to_string(unhit_cycle_enumeration_type));
    entries.emplace_back("find_optimal_result",to_string(find_optimal_result));
    entries.emplace_back("use_only_cpsat_lns",to_string(use_only_cpsat_lns));
    entries.emplace_back("ihs_iterations_in_mtz",to_string(ihs_iterations_in_mtz));
    entries.emplace_back("max_cycles_for_hs",to_string(max_cycles_for_hs));
    entries.emplace_back("log_cpsat_search_progress",to_string(log_cpsat_search_progress));
    entries.emplace_back("ihs_max_iterations",to_string(ihs_max_iterations));
    entries.emplace_back("next_sol_max_dst_from_init_sol",to_string(next_sol_max_dst_from_init_sol));
    entries.emplace_back("use_init_sol_as_hint_mode",to_string(use_init_sol_as_hint_mode));
    entries.emplace_back("init_L_for_all_constraints",to_string(init_L_for_all_constraints));
    entries.emplace_back("max_new_cycles_iter_scale",max_new_cycles_iter_scale);
    entries.emplace_back("max_new_cycles_per_iter",to_string(max_new_cycles_per_iter));
    entries.emplace_back("pi_arcs_perc_to_add",to_string(pi_arcs_perc_to_add));
    entries.emplace_back("fill_partial_result_using_greedy_fvs",to_string(fill_partial_result_using_greedy_fvs));
    entries.emplace_back("focus_mostly_onh_heuristics",to_string(focus_mostly_onh_heuristics));
    entries.emplace_back("ihs_init_sol_creation_mode",to_string(ihs_init_sol_creation_mode));
    entries.emplace_back("mtz_auxiliary_cycles_mode",to_string(mtz_auxiliary_cycles_mode));
    entries.emplace_back("use_cycle_trimming",to_string(use_cycle_trimming));
    entries.emplace_back("cycle_trimming_probab",to_string(cycle_trimming_probab));
    entries.emplace_back("cycle_trimming_min_nodes_in_hs",to_string(cycle_trimming_min_nodes_in_hs));

    return entries;
}

void ExpConfig::writeConfig() {
    auto entries = getConfigEntries();
    clog << "Config: " << endl;
    for ( auto & [k,v] : entries ) clog << "\t " << k << " = " << v << endl;
}
