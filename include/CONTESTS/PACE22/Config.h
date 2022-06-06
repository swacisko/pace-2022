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
        reducer_use_domination_6 = true; // do not enable this by default - very time consuming
        reducer_use_reverse_triangle_gadgets = true;
        reducer_use_mixed_domination = true;
        reducer_use_mixed_domination_full = true;
        reducer_use_bottleneck = true;
        reducer_use_bottleneck2 = true;
        reducer_use_cycle_folding = true;
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
