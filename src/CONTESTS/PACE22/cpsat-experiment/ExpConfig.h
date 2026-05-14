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
    DIV_IHS = 4,
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

    int max_time_sec = 15 * 60; // 15 minutes default runtime
    int ihs_single_iteration_sec = 1;

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
     * 0 - hint values for all N variables in init_sol
     * 1 - hint values only for variables of nodes that are set in init_sol
     * 2 - hint values only for variables of nodes that are set in prev_res
     * 3 - hint values for variables that are in best_fvs - this will reduce variability,
     * but should enable to find better HS in iterations
     */
    int use_init_sol_as_hint_mode = 3;

    /**
     * Value can be 1 or 2, since there are only two options.
     * 1 - standard version used in DiVerSeS
     * 2 - takes random directed tree, then for each arc that would close a cycles, considers this cycle and makes it
     * chordless if neccessary. This should be (in theory) better for graphs that do not contain short cycles or
     * contain very few of them and the rest is long, since enumeration of long induced cycles using the standard
     * method might be slow...
     * 3 - both 1 and 2
    */
    int unhit_cycle_enumeration_type = 3;

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

    /**
     * If true, some CPSAT's parameters will be set to run it much more heurtistically. This will severely impact
     * proving optimality, but should increase the quality/speed of used heuristics.
     */
    bool focus_mostly_onh_heuristics = true;

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
     * 2 - cycles found by IHS method run for at most 20% of total time and at most [ihs_max_iterations] iterations, or
     * for 30 seconds if find_optimal_result is set.
     * Remember that the bound on the number of iterations (at most ihs_iterations_in_mtz) still holds.
     */
    int mtz_auxiliary_cycles_mode = 2;

    /**
     *  Maximum number of cycles to create in HS method, so that we can terminate without exceeding memory limit.
     */
    int max_cycles_for_hs = 2e6;

    /**
     * If true, then cycle trimming will be used in IHS.
     */
    bool use_cycle_trimming = true;
    int cycle_trimming_freq = inf;
    double cycle_trimming_probab = 0.5;
    string cycle_trimming_max_cycles_to_trim_scale = "linear_h";
    int cycle_trimming_min_nodes_in_hs = 3;

    /**
     * This is the value, for which all cycles of at most this length will be created as constraints.
     */
    int init_L_for_all_constraints = 4;

    string max_new_cycles_iter_scale = "linear"; // options: linear_h, linear, log, sqrt, bounded, unbounded
    int max_new_cycles_per_iter = inf;
    static int scaleIters(int N, string scale, int bound = inf);

    /**
     * If true, then any time the CPSAT improves a solution, it is checked, whether it is a valid FVS of the graph V.
     * If so, it is returned. This way it might be possible to check many hitting sets, not just the best one
     * returned at the end.
     */
    bool check_incumbent_cpsat_solutions = true;

    /**
     * This number of pi-arcs will be added.
     * For |A| * pi_arcs_perc_to_add randomly selected arcs, their reverse counterparts will be added to the graph.
     */
    double pi_arcs_perc_to_add = 0.0;

    /**
     * If true, then a hitting set S for a set of cycles that is not a FVS of the graph V will be filled to a valid FVS
     * by adding to it a FVS of the graph V \ S, found by the Diverses solver, using basic reductions and
     * sinkhorn-knopp algorithm.
     */
    static constexpr bool fill_partial_result_using_greedy_fvs = true;


    /**
     * This steers how the initial solution is created in IHS.
     * 0 - no initial solution
     * 1 - uses getUnhitGraphGreedyFVS(V, emptyset) function to create the solution
     * 2 - uses diverses solver - this can be used to check how much the solution returned by diverses can be improved
     */
    int ihs_init_sol_creation_mode = 1;

    vector<pair<string, string>> getConfigEntries();

    void writeConfig();

};

#endif //DIVERSES_EXPCONFIG_H