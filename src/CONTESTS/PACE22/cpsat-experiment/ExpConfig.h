//
// Created by sylwe on 08/04/2026.
//

#ifndef DIVERSES_EXPCONFIG_H
#define DIVERSES_EXPCONFIG_H

#include "Makros.h"

enum Algorithm {
    HS1 = 0,
    IHS1,
    IHS2,
    MTZ1,
    MTZ2,
    DIVERSES,
};

string parseAlgorithm(int alg) {
    if(alg == HS1) return "HS1";
    if(alg == IHS1) return "IHS-1";
    if(alg == IHS2) return "IHS-2";
    if(alg == MTZ1) return "MTZ-1";
    if(alg == MTZ2) return "MTZ-2";
    if(alg == DIVERSES) return "DiVerSeS";
}

class ExpConfig{
public:

    Algorithm alg = IHS1;
    /**
     * Sets proper values of different parameters to run the given algorithms.
     * This uses the standard configuration for some parameters, their values might be changed later if needed.
     */
    void setParametersForAlgorithm(Algorithm alg);

    int max_time_sec = 10 * 60; // ten minutes default runtime
    int ihs_single_iteration_sec = 2;

    /**
     * Value can be 1 or 2, since there are only two options.
     * 1 - standard version used in DiVerSeS
     * 2 - takes random directed tree, then for each arc that would close a cycles, considers this cycle and makes it
     * chordless if neccessary. This should be (in theory) better for graphs that do not contain short cycles or
     * contain very few of them and the rest is long, since enumeration of long induced cycles using the standard
     * method might be slow...
    */
    int unhit_cycle_enumeration_type = 1;

    /**
     * If true, then optimal result will be found. In such a case all running time bounds set by parameters are ignored,
     * perhaps some other parameters as well.
     */
    bool find_optimal_result = false;

    int threads = 4;

    /**
     * if true, then only LNS will be used from cpsat #CAUTION! This might not work as intended, unless carefully taken care of
     */
    bool use_only_cpsat_lns = false;

    string metadata_filepatg = "";

    /**
     * The fraction of [max_time_sec] which will be spent on finding chordless cycles that will be used to augment
     * the MTZ formulation.
     * #CAUTION! Parameters like ihs_single_iteration_sec must be set manually to make it work reasonably.
     */
    double max_time_fraction_for_ihs_cycles_in_mtz = 0.2;

    /**
     *  Maximum number of cycles to create in HS method, so that we can terminate without exceeding memory limit.
     */
    int max_cycles_for_hs = 1e8;



    vector<pair<string, string>> getConfigEntries() {
        vector<pair<string, string>> entries;
    }



};

#endif //DIVERSES_EXPCONFIG_H