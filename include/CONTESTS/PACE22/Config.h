//
// Created by sylwester on 12/20/21.
//

#ifndef ALGORITHMSPROJECT_CONFIG_H
#define ALGORITHMSPROJECT_CONFIG_H

#include <utils/Stopwatch.h>
#include <csignal>

class Config{
public:

    static const int agent_flow_remove_largest_flow_node = 0;
    static const int agent_flow_merge_smallest_flow_node = 1;

    int agent_flow_node_selection_type = agent_flow_merge_smallest_flow_node;

    bool agent_flow_alternate_selection_type = false;
    bool solver_improve_alternate_selection_type = false;

    static const int agent_flow_continuous = 0;
    static const int agent_flow_tokens = 1;
    static const int agent_flow_sinkhorn = 2;

    static const int agent_flow_methods_cnt = 3;

    bool agent_flow_alternate_flow_method = false;

    int agent_flow_method = agent_flow_sinkhorn;

    int agent_flow_node_update_frequency = 1;

    int agent_flow_min_distance = 2;
    int agent_flow_max_distance_from_best = 1000; // this should be set to e.g. sqrt(N)


    //************************************************************************************ REDUCER


    bool reducer_use_pie = false;
    bool reducer_use_core = false;
    bool reducer_use_strongly_connected = false;
    bool reducer_use_dome = false;
    bool reducer_use_twins_merge = false;
    bool reducer_use_inoutclique = false;
    bool reducer_use_folding = false;
    bool reducer_use_general_folding = false;
    bool reducer_use_folding_twins = false;
    bool reducer_use_full_bipartite_blocker = false;
    bool reducer_use_edge_neighborhood_blocker = false;
    bool reducer_use_desk = false;
    bool reducer_use_unconfined = false;
    bool reducer_use_nonsimple_cycle_arcs = false;
    bool reducer_use_nonsimple_cycle_arcs_full = false;
    bool reducer_use_domination = false;
    bool reducer_use_domination_3 = false;
    bool reducer_use_domination_4 = false;
    bool reducer_use_domination_5 = false;
    bool reducer_use_domination_6 = false;
    bool reducer_use_domination_6inserter = false;
    bool reducer_use_reverse_triangle_gadgets = false;
    bool reducer_use_mixed_domination = false;
    bool reducer_use_mixed_domination_full = false;
    bool reducer_use_funnel = false;
    bool reducer_use_cycle_folding = false;
    bool reducer_use_spiderweb_gadgets = false;
    bool reducer_use_bottleneck = false;
    bool reducer_use_bottleneck2 = false;
    bool reducer_use_recursive_reducer = false;

    //******************************************************************** ED
    bool reducer_use_ed = false;

    /**
     * 0 - only at the very end
     * 1 - apply only once after reducible and before liftable rules, then at the very end.
     * 2 - apply exhaustively after all reducible reductions are done and before foldable reductions are done,
     *      and at the end
     *
     * When ed is applied at the end and no reducible node is found, edges are added according to type1-constraints
     */
    int ed_application_mode = 0;

    /**
     * Used for order in which nodes are considered by ED rule
     * 0 - order in which the nodes just are, nothing is done
     * 1 - sort nodes considered in ED by degrees, largest to smallest
     * 2 - sort nodes considered in ED by degrees, smallest to largest
     */
    int ed_node_sorting_mode = 1;

    /**
     * Order in which nodes in
     * 0 - order in which the nodes just are, nothing is done
     * 1 - sort nodes in U1 based on their |N(u) \setminus W|, largest to smallest
     * 2 - sort nodes in U1 based on their |N(u) \setminus W|, smallest to largest
     */
    int ed_U_nodes_sorting_mode = 0;

    /**
     * If true, then all detected nodes are moved to the set U simultaneously, in the same iteration, otherise
     * they are moved one in each iteration.
     * Similary for the set S.
     */
    bool ed_move_nodes_to_U_simultaneously = true;
    bool ed_move_nodes_to_S_simultaneously = false;

    /**
     * If true, then nodes are considered to be moved to the set U, even if they do not belong to N(W).
     * This might considerably speed up the process - and not only speed up, but also increase the areas of
     * subgraphs searched by the rule.
     * This should be set to true for best performance (both quality and efficiency).
     */
    bool ed_consider_nodes_to_move_outside_NW = false;

    /**
     * This is the standard concept that must be used.
     */
    constexpr static bool use_ed_domination = true;

    /**
     * If true, then the ``same neighborhood'' appraoch will be used to identify nodes to move to U.
     * This is just a special case of `deficit1 domination'' and works in the same complexity,
     * but should be slighly faster in practice (but possibly much weaker).
     */
    bool ed_use_same_neigh_domination = false;

    /**
     * If true, then the ``deficit1'' approach will be used to find nodes to move to U.
     * This is a GENERALIZATION of the ``same neighborhood'' rule.
     */
    bool ed_use_deficit1_domination = false;

    /**
     * If true, then the ``biset'' approach will be used to find nodes to move to U.
     * CAUTION! This rules is slower than other rules used to determine nodes to move to U.
     * It needs to be checked if in practice it is efficient enough, and perhaps limit it to only some special
     * (smaller) substructures.
     */
    bool ed_use_double_ed_checks = false;
    int ed_double_ed_max_candidates = 10;

    /**
     * If true, then apart from nodes to add to the kernel, edges will be checked for possible removal.
    For each edge {a,b}, consider({a,b}) will be run, and the edge will be removed if it returns true.
     * #CAUTION! This requires careful handling for solution lifting...
     */
    bool ed_use_edge_removal = false;

    /**
     * If true, then type-1 cnstraints will be added on the fly, as the reduction rule executes.
     */
    bool ed_apply_type1_constraints_on_the_fly = false;

    /**
     * If true, then initial sets S of the form {v} will be checked for each node in the graph.
     * This is the default check and should be set to true, unless you really need to disable that for some reason.
     */
    bool ed_use_node_removal = true;

    /**
     * Only nodes u \in U1 with |N(u) \setminus W| <= ed_ext_dom_max_node_neigh will be considered
     */
    int ed_ext_dom_max_node_neigh = 10;

    /**
     * If true, then in the Reducer there will be at the very end considered adding edges to the graph on the fly
     */
    bool ed_use_edge_insertion = true;

    /**
     * If true, then extensive edge insertion will be used.
     * Many candidate nodes x will be considered for each node v, and pairs {v,x} will be checked using consider({v,x}).
     * If it returns true, then edge {v,x} will be added to the graph.
     */
    bool ed_use_extended_edges_insertion = false;

    //******************************************************************** ED

    int reducer_simple_cycle_max_branch_depth = 50; // set this to 1e9 to make full search for simple cycles

    int reducer_max_time_millis = 60'000; // one minute max reduction time


    int reducer_nonsimple_cycle_arcs_full_max_time_millis_per_arc = 100;
    int reducer_nonsimple_cycle_arcs_full_max_time_millis_total = 7'000;

    int reducer_domination4_max_time_millis_per_node = 100;
    int reducer_domination_3_4_max_time_millis_total = 5'000;

    int reducer_mixed_domination_full_max_time_millis_per_node = 100;
    int reducer_mixed_domination_full_max_time_millis_total = 5'000;

    int reducer_domination5_max_time_millis_per_node = 100;
    int reducer_domination_5_max_time_millis_total = 5'000;

    int reducer_domination6_max_neigh_size = 4;
    int reducer_domination6inserter_max_neigh_size = 7;
    int reducer_domination6inserter_distance = 3; // greater values than 6 have no sense (smaller can be a bit just faster)

    int reducer_max_twin_merge_neighborhood_size = 16; // for smaller graphs it may be larger, e.g. 24 seems to be good
    int reducer_max_bottleneck2_neighborhood_size = 10;
    int reducer_max_general_folding_neighborhood_size = 10;
    int reducer_max_general_folding_antiedges = 3; // original value 1e9

    int reducer_max_component_size_for_spiderweb_gadgets = 10;

    void enableAllReductions(){
        reducer_use_core = true;
        reducer_use_dome = true;
        reducer_use_pie = true;
        reducer_use_inoutclique = true;
        reducer_use_twins_merge = true;

        reducer_use_nonsimple_cycle_arcs = true;
        reducer_use_nonsimple_cycle_arcs_full = true;
        reducer_use_folding = true;
        reducer_use_general_folding = true;
        reducer_use_folding_twins = true;
        reducer_use_full_bipartite_blocker = true;
        reducer_use_edge_neighborhood_blocker = true;
        reducer_use_desk = true;
        reducer_use_unconfined = true;
        reducer_use_funnel = true;
        reducer_use_domination = true;
        reducer_use_domination_3 = true;
        reducer_use_domination_4 = true;
        reducer_use_domination_6 = true;
        reducer_use_reverse_triangle_gadgets = true;
        reducer_use_mixed_domination = true;
        reducer_use_mixed_domination_full = true;
        reducer_use_bottleneck = true;
        reducer_use_bottleneck2 = true;
        reducer_use_cycle_folding = true;

        // reducer_use_ed = true;
    }

    void disableAllNonbasicReductions(){
        reducer_use_core = false;
        reducer_use_dome = false;
        reducer_use_pie = false;
        reducer_use_inoutclique = false;
        reducer_use_twins_merge = false;
        reducer_use_strongly_connected = false;

        reducer_use_nonsimple_cycle_arcs = false;
        reducer_use_nonsimple_cycle_arcs_full = false;
        reducer_use_folding = false;
        reducer_use_general_folding = false;
        reducer_use_folding_twins = false;
        reducer_use_full_bipartite_blocker = false;
        reducer_use_edge_neighborhood_blocker = false;
        reducer_use_desk = false;
        reducer_use_unconfined = false;
        reducer_use_funnel = false;
        reducer_use_domination = false;
        reducer_use_domination_3 = false;
        reducer_use_domination_4 = false;
        reducer_use_domination_5 = false;
        reducer_use_domination_6 = false;
        reducer_use_domination_6inserter = false;
        reducer_use_reverse_triangle_gadgets = false;
        reducer_use_mixed_domination = false;
        reducer_use_mixed_domination_full = false;
        reducer_use_bottleneck = false;
        reducer_use_bottleneck2 = false;
        reducer_use_cycle_folding = false;
        reducer_use_spiderweb_gadgets = false;
        reducer_use_recursive_reducer = false;

        reducer_use_ed = false;
    }

    void disableAllConditionalReductions(){
        reducer_use_funnel = false;
        reducer_use_folding = false;
        reducer_use_general_folding = false;
        reducer_use_folding_twins = false;
        reducer_use_full_bipartite_blocker = false;
        reducer_use_desk = false;
        reducer_use_cycle_folding = false;
        reducer_use_spiderweb_gadgets = false;
        reducer_use_reverse_triangle_gadgets = false;
    }

    void disableAllRecursiveReductions(){
        reducer_use_twins_merge = false;
        reducer_use_bottleneck = false;
        reducer_use_bottleneck2 = false;
        reducer_use_general_folding = false;
        reducer_use_recursive_reducer = false;

        reducer_use_spiderweb_gadgets = false;
    }

//    ****************************************************************************************** DFVSSolverH
    bool solverh_use_superpi_vc_ub = true;
    int solverh_improvement_iterations = 10; // this number of improvement iterations will be done in DFVSolverH
    bool solverh_use_sals = true; // this number of improvement iterations will be done in DFVSolverH
    bool solverh_use_sals_improver = false;
    int solverh_sals_improver_iterations = 3; // original value 3

    bool solverh_use_conditional_sals_improver = false;
    int solverh_conditional_sals_improver_min_size = 300;

    double min_density_for_sals3 = 0.067;

    bool solverh_use_reductions_for_each_scc = false;

    bool solverh_use_reductions_initial = false;

    bool solverh_use_reductions_AF = false;

    int solverh_min_graph_size_for_improvements = 30;


    //****************************************************************************************** IHS
    int solverh_ihs_init_max_cycles = 300; // default value 300 seems to work well

    int solverh_ihs_secondary_max_cycles = 1e9;

    int solverh_ihs_max_rescaling_times = 10;

    bool solverh_use_hsls_after_each_node_addition = false;

    int ihs_init_cycle_length = 3; // default value 3 seem to work well

    int ihs_hsls_perm_deviation_frequency = 800; // default value 800 seems to work well

    bool hsls_use_continuous_perm_deviation = true;
    //******************************************************************************************




    //******************************************************************************** VCImprover
    long long vc_improver_milliseconds = 200;

    //********************************************************************************* DFVSImprover
    int dfvsimprover_local_optimum_violation_frequency = 1e9;
    double dfvsimprover_alpha = 0.6;

    double dfvsimprover_local_opt_max_deviation_from_best_relative = 0.003;
    int dfvsimprover_local_opt_max_deviation_from_best_absolute_addition = 3;
    int dfvsimprover_max_iters_without_improvement = 15;


    //*********************************************************************************************
    Stopwatch sw;

    bool tle(){ return sw.tle("main") || sigterm_received; }

    volatile static sig_atomic_t sigterm_received;
    static void terminate(int signum) { sigterm_received = 1; }
    static void addSigtermCheck(){
        struct sigaction action;
        memset(&action, 0, sizeof(struct sigaction));
        action.sa_handler = Config::terminate;
        sigaction(SIGTERM, &action, NULL);
    }

    bool write_logs = true;

};

#endif //ALGORITHMSPROJECT_CONFIG_H
