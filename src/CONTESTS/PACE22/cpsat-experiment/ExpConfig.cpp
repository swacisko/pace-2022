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
    entries.emplace_back("ihs_iterations_in_mtz",to_string(ihs_iterations_in_mtz));
    entries.emplace_back("max_cycles_for_hs",to_string(max_cycles_for_hs));
    entries.emplace_back("log_cpsat_search_progress",to_string(log_cpsat_search_progress));
    entries.emplace_back("ihs_max_iterations",to_string(ihs_max_iterations));

    return entries;
}

void ExpConfig::writeConfig() {
    auto entries = getConfigEntries();
    clog << "Config: " << endl;
    for ( auto & [k,v] : entries ) clog << "\t " << k << " = " << v << endl;
}
