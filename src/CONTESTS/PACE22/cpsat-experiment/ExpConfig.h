//
// Created by sylwe on 08/04/2026.
//

#ifndef DIVERSES_EXPCONFIG_H
#define DIVERSES_EXPCONFIG_H

#include "Makros.h"

enum Algorithm {
    HS = 0,
    IHS = 1, // IHS with cycle_enumeration_type=1
    MTZ = 2,
    DIVERSES = 3,
};

string parseAlgorithm(int alg);
constexpr int inf = 1e9+1;

class ExpConfig{
public:

    Algorithm alg = IHS;
    /**
     * Sets proper values of different parameters to run the given algorithms.
     * This uses the standard configuration for some parameters, their values might be changed later if needed.
     */
    void setParametersForAlgorithm();

    int max_time_sec = 10 * 60; // ten minutes default runtime
    int ihs_single_iteration_sec = 3;

    /**
     * This number of iterations will be done in IHS, until it terminates.
     * Valid result might not be found if this value is too small. This parameter enables finding some subset of cycles.
     */
    int ihs_max_iterations = inf;

    /**
     * If set to some value ofther than inf, then an additional constraint will be added to the solvers.
     * For an initial solution, the next, improved solution can differ at most by [next_sol_max_dst_from_init_sol]
     * from the initial one. This might help speed up search by making it more local.
     */
    int next_sol_max_dst_from_init_sol = inf;

    /**
     * 0 - hint values for all N variables
     * 1 - hint values only for variables of nodes that are in init_sol
     * 2 - hint values only for variables of nodes that are in prev_res
     */
    int use_init_sol_as_hint_mode = 2;

    /**
     * Value can be 1 or 2, since there are only two options.
     * 1 - standard version used in DiVerSeS
     * 2 - takes random directed tree, then for each arc that would close a cycles, considers this cycle and makes it
     * chordless if neccessary. This should be (in theory) better for graphs that do not contain short cycles or
     * contain very few of them and the rest is long, since enumeration of long induced cycles using the standard
     * method might be slow...
     * 3 - both 1 and 2
    */
    int unhit_cycle_enumeration_type = 1;

    /**
     * If true, then optimal result will be found. In such a case all running time bounds set by parameters are ignored,
     * perhaps some other parameters as well.
     */
    bool find_optimal_result = false;

    int threads = 4;

    bool log_cpsat_search_progress = false;

    /**
     * if true, then only LNS will be used from cpsat #CAUTION! This might not work as intended, unless carefully taken care of
     */
    bool use_only_cpsat_lns = false;

    string metadata_filepath;

    /**
     * This number of iterations will be done in IHS to find cycles.
     * The expected time is roughly  ihs_single_iteration_sec * ihs_iterations_in_mtz   plus time required to list cycles.
     */
    int ihs_iterations_in_mtz = inf;

    /**
     * Option used to create auxiliary cycles for MTZ formulation.
     * 0 - no auxiliary cycles
     * 1 - all chordless cycles with length up to  [init_L_for_all_constraints]
     * 2 - cycles found by IHS method run for at most 20% of total time and at most [ihs_max_iterations] iterations.
     */
    int mtz_auxiliary_cycles_mode = 2;

    /**
     *  Maximum number of cycles to create in HS method, so that we can terminate without exceeding memory limit.
     */
    int max_cycles_for_hs = 1e7;

    /**
     * This is the value, for which all cycles of at most this length will be created as constraints.
     */
    int init_L_for_all_constraints = 4;

    string max_new_cycles_iter_scale = "log";
    int max_new_cycles_per_iter = inf;
    int scaleIters(int N);


    vector<pair<string, string>> getConfigEntries();

    void writeConfig();

};

#endif //DIVERSES_EXPCONFIG_H