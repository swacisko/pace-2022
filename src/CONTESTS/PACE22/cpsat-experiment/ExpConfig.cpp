//
// Created by sylwe on 08/04/2026.
//

#include "ExpConfig.h"

void ExpConfig::setParametersForAlgorithm(Algorithm alg) {
}

vector<pair<string, string>> ExpConfig::getConfigEntries() {
    vector<pair<string, string>> entries;
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
