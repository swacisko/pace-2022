//
// Created by sylwe on 08/04/2026.
//

#include "ExpConfig.h"

void ExpConfig::setParametersForAlgorithm() {
    if ( alg == HS1 ) {
        ihs_single_iteration_sec = 10;
        max_cycles_for_hs = 1e8;
        unhit_cycle_enumeration_type = 1;
    }
    if (alg == IHS1 || alg == IHS2) {
        ihs_single_iteration_sec = 3;
        unhit_cycle_enumeration_type = 1;
    }
    if (alg == MTZ1 || alg == MTZ2) {
        max_time_fraction_for_ihs_cycles_in_mtz = 0.2;
        ihs_single_iteration_sec = 2;
        max_cycles_for_hs = 1e8;
        unhit_cycle_enumeration_type = 1;
    }
}

string parseAlgorithm(int alg)  {
    using ::Algorithm;
    if(alg == HS1) return "HS1";
    if(alg == IHS1) return "IHS-1";
    if(alg == IHS2) return "IHS-2";
    if(alg == MTZ1) return "MTZ-1";
    if(alg == MTZ2) return "MTZ-2";
    if(alg == DIVERSES) return "DiVerSeS";
    return "unknown algorithm";
}


vector<pair<string, string>> ExpConfig::getConfigEntries() {
    vector<pair<string, string>> entries;
    entries.emplace_back("metadata_filepath",metadata_filepath);
    entries.emplace_back("threads",to_string(threads));
    entries.emplace_back("max_time_sec",to_string(max_time_sec));
    entries.emplace_back("ihs_single_iteration_sec",to_string(ihs_single_iteration_sec));
    entries.emplace_back("unhit_cycle_enumeration_type",to_string(unhit_cycle_enumeration_type));
    entries.emplace_back("find_optimal_result",to_string(find_optimal_result));
    entries.emplace_back("use_only_cpsat_lns",to_string(use_only_cpsat_lns));
    entries.emplace_back("max_time_fraction_for_ihs_cycles_in_mtz",to_string(max_time_fraction_for_ihs_cycles_in_mtz));
    entries.emplace_back("max_cycles_for_hs",to_string(max_cycles_for_hs));
    entries.emplace_back("log_cpsat_search_progress",to_string(log_cpsat_search_progress));

    return entries;
}
