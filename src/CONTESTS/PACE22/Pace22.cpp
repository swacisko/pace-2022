//
// Created by sylwester on 12/20/21.
//

#include <CONTESTS/PACE22/Reducer.h>
#include <CONTESTS/PACE22/heur/VCImprover.h>
#include <graphs/scc/StronglyConnectedComponents.h>
#include <graphs/GraphInducer.h>
#include <CONTESTS/PACE22/exact/DFVSSolverE.h>
#include <utils/TimeMeasurer.h>
#include <graphs/VertexCover/approximation/LibMVC/fastvc.h>
#include <filesystem>
#include "CONTESTS/PACE22/heur/DFVSSolverH.h"
#include "MemoryUtils.h"
#include "getopt.h"

constexpr bool USE_ONLY_AF = false;

constexpr bool write_solution = true; // if true, then DFVS solution will be written to stdout

bool use_heuristic_solver = true;
bool use_exact_solver = !use_heuristic_solver;

constexpr bool lite_track = false;
bool MUTE_MODE = false;
string input_filepath = "";

static int time_limit_millis = 590'000;

vector<DFVSReduction*> initialReductions(VVI &V, Config cnf, bool heuristic_track){


    constexpr bool use_graph_sparsification = true;

    if(heuristic_track && use_graph_sparsification
        &&  GraphUtils::countEdges(V,true) > 1'500'000
        ){
        // Some unhit cycles will be hit in emergencyExit
        double pi_perc = Utils::getPieEdgesPercentage(V);
        double threshold = 0.6;
        if( pi_perc > threshold ) {
            DEBUG(GraphUtils::countEdges(V,true));
            int length = 4;
            if(pi_perc > 0.7) length = 5;
            if(pi_perc > 0.8) length = 6;

            set<PII> arcs;
            { // this is much faster to get all arcs in some 'short' induced cycle
                vector<Triple<int>> min_lengths = Utils::getLengthOfShortestInducedCycleWithArc(V, length);
                for( auto & tr : min_lengths ) if( tr.rd <= length ) arcs.insert(PII(tr.st, tr.nd));
            }

            for (int i = 0; i < V.size(); i++) V[i].clear();
            VB helper(V.size(), false);
            VPII to_add(ALL(arcs));
            Utils::addEdges(V, to_add, helper);
            DEBUG(GraphUtils::countEdges(V, true));
            clog << "Main elapsed time: " << (cnf.sw.getTime() / 1000) << " sec." << endl;
        }else{
            clog << "Graph to sparse for 'edge-sparsification" << endl;
        }
    }

    cnf.write_logs = false;

    Reducer red(V,cnf);

    {
        red.cnf.enableAllReductions();
        red.cnf.reducer_use_strongly_connected = false;
        red.cnf.reducer_use_spiderweb_gadgets = false;

        red.cnf.reducer_use_recursive_reducer = false; // this is time-consuming!!!

//        red.cnf.reducer_use_domination_6inserter = true; // extremely-time-consuming
//        red.cnf.reducer_domination6inserter_max_neigh_size = 5;
//        red.cnf.reducer_domination6inserter_distance = 3;

//        red.cnf.reducer_use_bottleneck = false;
//        red.cnf.reducer_use_bottleneck2 = false;

        if(heuristic_track){
            red.cnf.reducer_use_inoutclique = false;
            red.cnf.reducer_use_recursive_reducer = false;

            { // disable rarely used reductions
                red.cnf.reducer_use_cycle_folding = false;
                red.cnf.reducer_use_desk = false;
                red.cnf.reducer_use_general_folding = false;
                red.cnf.reducer_use_edge_neighborhood_blocker = false;
                red.cnf.reducer_use_reverse_triangle_gadgets = false;
            }

        }
        else{
            red.cnf.reducer_domination6_max_neigh_size = 20;
        }


        { // time-consuming reductions
            if( GraphUtils::countEdges(V,true) > 1e6 && !Utils::isPIGraph(V) ){
                red.cnf.reducer_use_domination_6 = false; // effective, but time-consuming for large graphs
                red.cnf.reducer_use_domination_6inserter = false; // effective, but time-consuming for large graphs
            }
            if( heuristic_track ){
                red.cnf.reducer_use_bottleneck = false; // only for very small graphs
                red.cnf.reducer_use_bottleneck2 = false; // only for very small graphs
                red.cnf.reducer_use_domination_6inserter = false;
            }
        }
    }


    red.cnf.reducer_max_component_size_for_spiderweb_gadgets = 5;


    red.cnf.reducer_max_twin_merge_neighborhood_size = 16;
    red.cnf.reducer_simple_cycle_max_branch_depth = 50;

    red.cnf.reducer_max_time_millis = 120 * time_limit_millis / 590;
    if(!heuristic_track) red.cnf.reducer_max_time_millis = 500'000;

    int MAX_MILLIS_PER_REDUCTION = 10'000;
    if(heuristic_track) MAX_MILLIS_PER_REDUCTION = 5'000;

    red.cnf.reducer_nonsimple_cycle_arcs_full_max_time_millis_total = MAX_MILLIS_PER_REDUCTION;

    red.cnf.reducer_domination4_max_time_millis_per_node = 100;
    red.cnf.reducer_domination_3_4_max_time_millis_total = MAX_MILLIS_PER_REDUCTION;

    red.cnf.reducer_domination5_max_time_millis_per_node = 100;
    red.cnf.reducer_domination_5_max_time_millis_total = MAX_MILLIS_PER_REDUCTION;

    red.cnf.reducer_mixed_domination_full_max_time_millis_per_node = 100;
    red.cnf.reducer_mixed_domination_full_max_time_millis_total = MAX_MILLIS_PER_REDUCTION;

    auto reductions = red.reduce();
    int red_size_ub = Reducer::getReductionsSizeDiff(reductions);
    DEBUG(red_size_ub);
    red.writeTotals();

    V = red.V;
    {
        DEBUG(GraphUtils::countEdges(V, true));
        VVI underPI = Utils::getUnderlyingPIGraph(V);
        DEBUG( GraphUtils::countEdges( underPI, true ) );
        DEBUG(1.0 * GraphUtils::countEdges(underPI, true) / GraphUtils::countEdges(V, true));

        auto nonpiV = Utils::getNonPIGraph(V);
        int nonpi_arcs = GraphUtils::countEdges(nonpiV, true);
        DEBUG(nonpi_arcs);

        {
            StronglyConnectedComponents scc(nonpiV);
            scc.createStronglyConnectedComponents();
            auto comps = scc.getComponents();
            comps.resize( remove_if( ALL(comps), [](auto & v){ return v.size() <= 2; } ) - comps.begin());

            VI sizes(V.size(),0);
            for( auto& v : comps ) sizes[v.size()]++;

            VI cycles_sizes(V.size(),0);
            for( auto & v : comps ){
                InducedGraph g = GraphInducer::induce(nonpiV, v);
                if( GraphUtils::countEdges(g.V,true) == g.V.size() ) cycles_sizes[v.size()]++;
            }

            clog << "There are " << comps.size() << " connected_components of size >= 3 in nonpiV" << endl;
            for( int i=3; i<V.size(); i++ ){
                if(sizes[i] > 0) clog << "There are " << sizes[i] << " components in nonpiV of size " << i
                     << ", from which " << cycles_sizes[i] << " are 'cycles'" << endl;
            }
        }

        ENDL(2);
    }

    TimeMeasurer::stop("main");
    TimeMeasurer::writeAllMeasurements();
    ENDL(10);

    int total_reduction_diff = red_size_ub;
    ENDL(1);
    DEBUG(total_reduction_diff);
    ENDL(1);

    TimeMeasurer::start("main");
    return reductions;
}



void initializeParams(int argc, char **argv) {
    string time_limit = "time-limit";
    string track = "track";
    string quiet = "quiet";
    string file = "file";

    static struct option long_options[] = {
            {time_limit.c_str(), required_argument, 0, 0},
            {track.c_str(), required_argument, 0, 0},
            {quiet.c_str(), required_argument, 0, 0},
            {file.c_str(), required_argument, 0, 0},
            {0, 0,                                           0, 0}
    };

    while (1) {
        int option_index = 0;
        int c;
        string option, option_name;

        c = getopt_long(argc, argv, "l:",
                        long_options, &option_index);
        if (c == -1) break;
        switch (c) {
            case 0:
                option = string(optarg);
                option_name = string(long_options[option_index].name);

                if (option_name == time_limit) {
                    time_limit_millis = stoi(option);
                }
                if(option_name == track){
                    if( option == "exact" ){
                        use_heuristic_solver = false;
                        use_exact_solver = true;
                    }
                }
                if(option_name == quiet){
                    if( option == "true" ) MUTE_MODE = true;
                }
                if(option_name == file){
                    input_filepath = option;
                }

                break;
            case '?':
                break;
            default:
                printf("?? getopt returned character code 0%o ??\n", c);
        }
    }

    if( input_filepath != "" && !filesystem::exists( input_filepath ) ){
        cerr << "File " << input_filepath << ", provided as input file, does not exist" << endl;
        exit(1);
    }

}


/**
 * MAIN ALGORITHM MAY NOT RUN DETERMINISTICALLY!
 * That is because some algorithms, like NuMVC, are run for e.g. 1'000 milliseconds. They may not find
 * the same solutions if run twice with the same time limit.
 */
int main(int argc, char** argv){
    MemoryUtils::increaseStack();
    Config::addSigtermCheck();

    initializeParams(argc, argv);

    auto old_clog_buf = clog.rdbuf();
    if(MUTE_MODE){
        clog << "MUTE MODE" << endl;
        clog.rdbuf( nullptr );
    }

    Config main_cnf;

    if( use_heuristic_solver ){
        if(lite_track) main_cnf.sw.setLimit("main", 295'000); // heuristic track
//        else main_cnf.sw.setLimit("main", 590'000); // heuristic track - wait for SIGTERM or 590 seconds
        else main_cnf.sw.setLimit("main", time_limit_millis); // heuristic track - wait for SIGTERM or 590 seconds
    }
    else{
        time_limit_millis = 300'000'000;
        main_cnf.sw.setLimit("main", time_limit_millis); // exact track - time limit set to 'almost infinity'
    }

    clog << "Using " << (use_heuristic_solver ? "heuristic" : "exact") << " version of DiVerSeS solver" << endl;
    clog << "Running DiVerSeS for " << (1.0 * time_limit_millis / 1000) << " seconds" << endl;

    main_cnf.sw.start("main");

    TimeMeasurer::start("main");


    VVI V;

    string test_path = ""; // for tests
//    test_path = "pace22_exact/e_127";
//    test_path = "pace22_heur/h_039";

    if(test_path != ""){
        ifstream str(test_path.c_str());
        V = Utils::readGraph(str);
        str.close();
    }else if( input_filepath != "" ){
        clog << "Reading input from file: " << input_filepath << endl;
        ifstream str(input_filepath.c_str());
        V = Utils::readGraph(str);
        str.close();
    }
    else{
        clog << "Reading input from standard input" << endl;
        V = Utils::readGraph(cin);
    }

    DEBUG(V.size());
    DEBUG(GraphUtils::countEdges(V,true));
    assert(GraphUtils::isSimple(V));


    //***************************************************** INITIAL REDUCTIONS AND REMAPPING
    VVI origV = V;
    VI origV_dfvs;
    vector<DFVSReduction*> reductions;
    if(!USE_ONLY_AF) {
    // uncomment following two line not to use the whole set of reductions (WGYC will be used only)
        if (use_exact_solver && Utils::isPIGraph(V) && filesystem::exists("vc_solver") ) {}
        else reductions = initialReductions(V, main_cnf, use_heuristic_solver);

//        reductions = initialReductions(V, main_cnf, use_heuristic_solver);
    }

    InducedGraph g = GraphInducer::induceByNonisolatedNodes(V);
    V = g.V;
    assert(GraphUtils::isSimple(V));
    DEBUG(g.V.size());
    DEBUG(GraphUtils::countEdges(g.V,true));
    DEBUG(GraphUtils::density(g.V,true));


    auto liftSolution = [&]( VI & dfvs ){
        set<int> zb(ALL(dfvs));
        dfvs = VI(ALL(zb));

        ENDL(10);
        clog << "LIFTING SOLUTION!" << endl;
        DEBUG(dfvs.size());
        DEBUG(reductions.size());
        for (int &d : dfvs) d = g.nodes[d];
        int red_dfvs_size_ub = Reducer::getReductionsSizeDiff(reductions);
        int init_dfvs_size = dfvs.size();
        DEBUG(red_dfvs_size_ub);

        int N_ub = origV.size() + V.size();
        Reducer::liftSolution(N_ub, dfvs, reductions);

        clog << "After lifting solution, dfvs.size(): " << dfvs.size() << endl;
        if(! (init_dfvs_size + red_dfvs_size_ub >= dfvs.size() ) ){
            clog << "Condition init_dfvs_size + red_dfvs_size_ub >= dfvs.size() in liftSolution() "
                    "does not hold!" << endl;
        }
        Utils::emergencyExit(origV, dfvs);
    };
    //********************************************************************************

    clog << "Main elapsed seconds: " << main_cnf.sw.getTime() / 1'000 << endl;


    if(V.empty()){
        use_exact_solver = use_heuristic_solver = false;
        VI dfvs = {};
        liftSolution(dfvs);
        origV_dfvs = dfvs;
    }


    if(use_heuristic_solver){ // heuristic solution
        TimeMeasurer::start("DFVSSolverH");
        ENDL(5);
        clog << "USING HEURISTIC SOLVER" << endl;


        DFVSSolverH solver(main_cnf);
        solver.cnf.reducer_max_time_millis = 30'000;
        solver.cnf.reducer_nonsimple_cycle_arcs_full_max_time_millis_total = 2'000;
        solver.cnf.reducer_domination_5_max_time_millis_total = 2'000;
        solver.cnf.reducer_mixed_domination_full_max_time_millis_total = 2'000;

        {
            solver.cnf.enableAllReductions();
//            solver.cnf.disableAllConditionalReductions();
            solver.cnf.disableAllRecursiveReductions();
//            solver.cnf.reducer_use_domination_6 = false;
            solver.cnf.reducer_use_domination_6inserter = false;
            solver.cnf.reducer_use_nonsimple_cycle_arcs_full = false;
            solver.cnf.reducer_use_bottleneck = false;
            solver.cnf.reducer_use_spiderweb_gadgets = false;

            solver.cnf.solverh_use_reductions_for_each_scc = true;
        }


        solver.cnf.agent_flow_max_distance_from_best = sqrt(V.size());
        // this should be default - for larger graphs we cannot merge nodes
        solver.cnf.agent_flow_node_selection_type = Config::agent_flow_remove_largest_flow_node;
        if(V.size() < 10'000) solver.cnf.agent_flow_node_selection_type = Config::agent_flow_merge_smallest_flow_node;

//        solver.cnf.solverh_use_sals_improver = true; // originally false
        solver.cnf.agent_flow_min_distance = 2;

        {
            double pie_edges_perc = Utils::getPieEdgesPercentage(V);
            if (V.size() < 30'000 && pie_edges_perc < 0.3 && GraphUtils::countEdges(V, true) < 150'000) {
                solver.cnf.agent_flow_max_distance_from_best = 10;
            }
        }

        solver.cnf.vc_improver_milliseconds = 2'000;
        if(Utils::isPIGraph(V)) solver.cnf.vc_improver_milliseconds = 10'000;

        solver.cnf.dfvsimprover_max_iters_without_improvement = 50;
        solver.cnf.solverh_improvement_iterations = 5;

        solver.cnf.solverh_min_graph_size_for_improvements = 5;

        solver.cnf.solverh_use_conditional_sals_improver = true;
        solver.cnf.solverh_conditional_sals_improver_min_size = 350;

        VI dfvs(V.size(),0); iota(ALL(dfvs),0);

        if(USE_ONLY_AF){
            solver.cnf.disableAllNonbasicReductions();
            solver.cnf.solverh_use_reductions_initial = true;
            solver.cnf.solverh_use_reductions_for_each_scc = false;
            solver.cnf.solverh_min_graph_size_for_improvements = 1e9; // disabling improvements
        }

        VVI all_cycles;
        dfvs = solver.solveForGraph(V);
        assert(Utils::isFVS(V,dfvs));


        constexpr bool use_ihs = (!USE_ONLY_AF);

        if(use_ihs){

            if(V.size() < 3'000 && dfvs.size() < 300 ){ // only for small graphs
                solver.cnf.ihs_hsls_perm_deviation_frequency = 800;
                solver.cnf.solverh_ihs_init_max_cycles = 100;
                solver.cnf.ihs_init_cycle_length = 3;
                solver.cnf.hsls_use_continuous_perm_deviation = true;
                solver.cnf.solverh_ihs_max_rescaling_times = 15;
                dfvs = solver.solveForStronglyConnectedIterativeHittingSet2(V, dfvs, all_cycles, 0.15);
                assert(Utils::isFVS(V,dfvs));
            }

            solver.cnf.ihs_hsls_perm_deviation_frequency = 800;
            solver.cnf.solverh_ihs_init_max_cycles = 300;
            solver.cnf.ihs_init_cycle_length = 3;
            solver.cnf.hsls_use_continuous_perm_deviation = true;
            dfvs = solver.solveForStronglyConnectedIterativeHittingSet2(V, dfvs, all_cycles, 0.1);
            assert(Utils::isFVS(V,dfvs));

            solver.cnf.ihs_hsls_perm_deviation_frequency = 800;
            solver.cnf.solverh_ihs_init_max_cycles = 150;
            solver.cnf.ihs_init_cycle_length = 3;
            solver.cnf.hsls_use_continuous_perm_deviation = true;
            dfvs = solver.solveForStronglyConnectedIterativeHittingSet2(V, dfvs, all_cycles, 0.15);
            assert(Utils::isFVS(V,dfvs));

            solver.cnf.ihs_hsls_perm_deviation_frequency = 400;
            solver.cnf.solverh_ihs_init_max_cycles = 100;
            solver.cnf.ihs_init_cycle_length = 3;
            solver.cnf.hsls_use_continuous_perm_deviation = true;
            dfvs = solver.solveForStronglyConnectedIterativeHittingSet2(V, dfvs, all_cycles, 0.25);
            assert(Utils::isFVS(V,dfvs));

            solver.cnf.ihs_hsls_perm_deviation_frequency = 600;
            solver.cnf.solverh_ihs_init_max_cycles = 50;
            solver.cnf.ihs_init_cycle_length = 3;
            solver.cnf.hsls_use_continuous_perm_deviation = true;
            dfvs = solver.solveForStronglyConnectedIterativeHittingSet2(V, dfvs, all_cycles, 0.5);
            assert(Utils::isFVS(V,dfvs));

            solver.cnf.ihs_hsls_perm_deviation_frequency = 300;
            solver.cnf.solverh_ihs_init_max_cycles = 200;
            solver.cnf.ihs_init_cycle_length = 4;
            solver.cnf.hsls_use_continuous_perm_deviation = true;
            dfvs = solver.solveForStronglyConnectedIterativeHittingSet2(V, dfvs, all_cycles, 0.5);
            assert(Utils::isFVS(V,dfvs));
        }

        Utils::emergencyExit(V, dfvs);

        clog << "solver.solveForGraph(V).size(): " << dfvs.size() << endl;
        TimeMeasurer::stop("DFVSSolverH");

        liftSolution(dfvs);
        origV_dfvs = dfvs;
    }



    if(use_exact_solver){

        clog << "USING EXACT SOLVER" << endl;
        Config cnf = main_cnf;

        cnf.enableAllReductions();
        cnf.disableAllRecursiveReductions();
        cnf.reducer_use_general_folding = true;
        cnf.reducer_use_domination_6 = false;
        cnf.reducer_use_domination_6inserter = false;
        cnf.reducer_use_nonsimple_cycle_arcs_full = false;
        cnf.reducer_use_spiderweb_gadgets = false;


        cnf.vc_improver_milliseconds = 2'000;
        cnf.dfvsimprover_max_iters_without_improvement = 10;
        cnf.solverh_improvement_iterations = 5;

        cnf.solverh_use_superpi_vc_ub = false;

        // the following parameters will be used to heuristically find cycles using IHS
        cnf.ihs_hsls_perm_deviation_frequency = 400;
        cnf.hsls_use_continuous_perm_deviation = false;
        cnf.ihs_init_cycle_length = 3;
        cnf.solverh_ihs_init_max_cycles = 10;

        DFVSSolverE solverE(&V, cnf);
        solverE.write_logs = true;
        VI exact_dfvs = solverE.solveForInputGraph(V);
//        VI exact_dfvs = solverE.solveForStronglyConnectedIterativeHittingSet2(V, VI(V.size(),0), true); // for IHS tests only

        DEBUG(exact_dfvs.size());
        assert( Utils::isFVS(V, exact_dfvs) );

        liftSolution(exact_dfvs);
        origV_dfvs = exact_dfvs;
    }

    ENDL(10);

//    if(MUTE_MODE) clog.rdbuf(old_clog_buf);

    DEBUG(origV_dfvs.size());

    TimeMeasurer::stop("main");
    TimeMeasurer::write();
    clog << "Main elapsed time (real): " << main_cnf.sw.getTime("main") / 1'000 << endl;


    if( !Utils::isFVS(origV, origV_dfvs) ){
        Utils::emergencyExit( origV, origV_dfvs );
        clog << "After emergency exit, origV_dfvs.size(): " << origV_dfvs.size() << endl;
    }

    {
        set<int> zb(ALL(origV_dfvs));
        origV_dfvs = VI(ALL(zb));
        assert( set<int>(ALL(origV_dfvs)).size() == origV_dfvs.size() );
    }

    DEBUG(origV_dfvs.size());

    if(write_solution){
        for(int d : origV_dfvs) cout << d+1 << "\n";
        cout << flush;
    }

    DEBUG(origV_dfvs.size());

    if(MUTE_MODE) clog.rdbuf(old_clog_buf);

    return 0;
}

