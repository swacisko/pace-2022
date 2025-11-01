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
#include <GraphReader.h>

#include "CONTESTS/PACE22/heur/DFVSSolverH.h"
#include "MemoryUtils.h"
#include "getopt.h"
#include "Makros.h"

#include <omp.h>

constexpr bool USE_ONLY_AF = false;

constexpr bool write_solution = true; // if true, then DFVS solution will be written to stdout

bool use_heuristic_solver = true;
bool use_exact_solver = !use_heuristic_solver;

constexpr bool lite_track = false;
bool MUTE_MODE = false;
string input_filepath = "";

static int time_limit_millis = 600'000;

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


struct ExpData {
    int N, M, pi_arcs, nonpi_arcs;
    double pi_arcs_fraction;
    double nonpi_arcs_fraction;

    static vector<string> getHeader() {
        // vector<string> fields{"N", "M", "pi_arcs", "npi_arcs", "paf"};
        // vector<string> fields{"N", "M", "paf", "time"};
        vector<string> fields{"N", "M", "time"};
        vector<string> res;

        // writeData(cur_data_os, initV,0);
        // writeData(cur_data_os, basic_data, time_basic);
        // writeData(cur_data_os, known, time_known);
        // writeData(cur_data_os, all_dom, time_all_dom);
        // writeData(cur_data_os, only_fast_dom, time_only_fast_dom);
        //
        // writeData(cur_data_os, dom12, time_dom12);
        // writeData(cur_data_os, dom3, time_dom3);
        // writeData(cur_data_os, dom4, time_dom4);
        // writeData(cur_data_os, dom5, time_dom5);
        // writeData(cur_data_os, dom_nonsimple_arcs, time_dom_nonsimple_arcs);
        // writeData(cur_data_os, dom_nonsimple_arcs_full, time_dom_nonsimple_arcs_full);
        // writeData(cur_data_os, dom_mixed, time_dom_mixex);
        // writeData(cur_data_os, dom_mixed_full, time_dom_mixed_full);


        for( string s : { "orig", "basic", "known", "all_dom", "fast_dom", "dom12", "dom3",
            "dom4", "dom5", "dom_nca", "dom_nca_full", "dom_mixed", "dom_mixed_full"} ) {
            for(const auto & f : fields) res.push_back(s + "-" + f);
            // res.push_back(s + "_N / bN" );
            // res.push_back(s + "_M / bM" );
            // res.push_back(s + "_N / kN" );
            // res.push_back(s + "_M / kM" );
        }


        return res;
    }

    void createData(VVI V) {
        N = V.size();
        M = GraphUtils::countEdges(V,true);
        pi_arcs = Utils::countPiEdges(V);
        nonpi_arcs = M - pi_arcs;
        if(M) pi_arcs_fraction = (double)pi_arcs / M;
        else pi_arcs_fraction = 0;

        if(M) nonpi_arcs_fraction = (double)nonpi_arcs / M;
        else nonpi_arcs_fraction = 0;
    }
};

using PEI = pair<ExpData,int>;


tuple<ExpData,
    PEI, PEI, PEI, PEI, PEI, PEI, PEI, PEI, PEI, PEI, PEI, PEI>
    createDataForInstance(VVI V, bool include_single_reduction) {

    auto induceByNonisolated = [&]() {
        InducedGraph g = GraphInducer::induceByNonisolatedNodes(V);
        V = g.V;
    };


    auto assignTimeLimits = [&](auto & red) {
        red.cnf.reducer_max_time_millis = 3'600'000; // 1 hour max time
        int MAX_MILLIS_PER_REDUCTION = 100'000;
        red.cnf.reducer_nonsimple_cycle_arcs_full_max_time_millis_total = MAX_MILLIS_PER_REDUCTION;
        red.cnf.reducer_domination4_max_time_millis_per_node = 1'000;
        red.cnf.reducer_domination_3_4_max_time_millis_total = MAX_MILLIS_PER_REDUCTION;
        red.cnf.reducer_domination5_max_time_millis_per_node = 1'000;
        red.cnf.reducer_domination_5_max_time_millis_total = MAX_MILLIS_PER_REDUCTION;
        red.cnf.reducer_mixed_domination_full_max_time_millis_per_node = 1'000;
        red.cnf.reducer_mixed_domination_full_max_time_millis_total = MAX_MILLIS_PER_REDUCTION;
    };

    ExpData initV_data, basic_data, known_red_data;
    ExpData all_dom, only_fast_dom, dom12, dom3, dom4, dom5,
        dom_nonsimple_arcs, dom_nonsimple_arcs_full, dom_mixed, dom_mixed_full;

    int time_basic, time_known, time_dom12, time_dom3, time_dom4, time_dom5, time_dom_nonsimple_arcs, time_dom_nonsimple_arcs_full,
        time_dom_mixed, time_dom_mixed_full, time_only_fast_dom, time_all_dom;

    induceByNonisolated();
    initV_data.createData(V);

    auto initV = V;

    Stopwatch sw;
    sw.setLimit("red", time_limit_millis);
    sw.start("red");

    { // basic rules
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        Reducer red(V, main_cnf);
        red.disableAllNonbasicReductions();
        assignTimeLimits(red);

        sw.restart("red");
        auto reductions = red.reduce();
        for(auto * r : reductions) delete r;
        V = red.V;
        induceByNonisolated();
        basic_data.createData(V);

        time_basic = sw.getTime("red");
    }

    VVI basicV = V;

    { // known rules
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        // Reducer red(initV, main_cnf);
        red.disableAllNonbasicReductions();
        red.cnf.reducer_use_core = true;
        red.cnf.reducer_use_dome = true;
        red.cnf.reducer_use_pie = true;
        red.cnf.reducer_use_inoutclique = true;
        // red.cnf.reducer_use_nonsimple_cycle_arcs = true;
        // red.cnf.reducer_use_nonsimple_cycle_arcs_full = true;
        assignTimeLimits(red);

        sw.restart("red");
        auto reductions = red.reduce();
        for(auto * r : reductions) delete r;
        V = red.V;
        induceByNonisolated();
        known_red_data.createData(V);

        time_known = sw.getTime("red");
    }




    auto setAllForRed = [&](auto & red,
        const bool use_dom4 = false, const bool use_dom5 = false,
        const bool use_nonsimple_cycle_arcs_full = false,
        bool use_mixed_domination_full = false) {

        red.disableAllNonbasicReductions();

        red.cnf.reducer_use_core = true;
        red.cnf.reducer_use_dome = true;
        red.cnf.reducer_use_pie = true;
        red.cnf.reducer_use_inoutclique = true;

        red.cnf.reducer_use_core = !include_single_reduction;
        red.cnf.reducer_use_dome = !include_single_reduction;
        red.cnf.reducer_use_pie = !include_single_reduction;
        red.cnf.reducer_use_inoutclique = !include_single_reduction;
        red.cnf.reducer_use_nonsimple_cycle_arcs = !include_single_reduction;

        red.cnf.reducer_use_nonsimple_cycle_arcs_full = use_nonsimple_cycle_arcs_full;

        red.cnf.reducer_use_domination = !include_single_reduction;
        red.cnf.reducer_use_domination_3 = !include_single_reduction;

        red.cnf.reducer_use_domination_4 = use_dom4;
        red.cnf.reducer_use_domination_5 = use_dom5;

        red.cnf.reducer_use_mixed_domination = !include_single_reduction;

        red.cnf.reducer_use_mixed_domination_full = use_mixed_domination_full;
    };

    auto createForRed = [&](auto & red, auto & data) {
        assignTimeLimits(red);
        auto reductions = red.reduce();
        for(auto * r : reductions) delete r;
        V = red.V;
        induceByNonisolated();
        data.createData(V);
    };




    if(true){ // all_dom
        auto b = include_single_reduction;
        include_single_reduction = 0;

        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red, true, true, true, true);

        sw.restart("red");
        createForRed(red, all_dom);
        time_all_dom = sw.getTime("red");

        include_single_reduction = b;
    }

    if(true){ // fast_dom
        auto b = include_single_reduction;
        include_single_reduction = 0;

        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red, false,false, false, false);

        sw.restart("red");
        createForRed(red, only_fast_dom);
        time_only_fast_dom = sw.getTime("red");

        include_single_reduction = b;
    }


    if(true){ // reducer_use_domination
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red,!include_single_reduction);
        red.cnf.reducer_use_domination = include_single_reduction;

        sw.restart("red");
        createForRed(red, dom12);
        time_dom12 = sw.getTime("red");
    }

    if(true){ // reducer_use_domination_3
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red,!include_single_reduction);
        red.cnf.reducer_use_domination_3 = include_single_reduction;

        sw.restart("red");
        createForRed(red, dom3);
        time_dom3 = sw.getTime("red");
    }

    if(true){ // reducer_use_domination_4
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single_reduction);
        red.cnf.reducer_use_domination_4 = include_single_reduction;

        sw.restart("red");
        createForRed(red, dom4);
        time_dom4 = sw.getTime("red");
    }

    if(true){ // reducer_use_domination_5
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single_reduction);
        red.cnf.reducer_use_domination_5 = include_single_reduction;

        sw.restart("red");
        createForRed(red, dom5);
        time_dom5 = sw.getTime("red");
    }

    if(true){ // nonsimple arcs
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single_reduction);
        red.cnf.reducer_use_nonsimple_cycle_arcs = include_single_reduction;

        sw.restart("red");
        createForRed(red, dom_nonsimple_arcs);
        time_dom_nonsimple_arcs = sw.getTime("red");
    }

    if(true){ // nonsimple arcs full
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single_reduction);
        red.cnf.reducer_use_nonsimple_cycle_arcs_full = include_single_reduction;

        sw.restart("red");
        createForRed(red, dom_nonsimple_arcs_full);
        time_dom_nonsimple_arcs_full = sw.getTime("red");
    }

    if(true){ // mixed domination
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single_reduction);
        red.cnf.reducer_use_mixed_domination = include_single_reduction;

        sw.restart("red");
        createForRed(red, dom_mixed);
        time_dom_mixed = sw.getTime("red");
    }

    if(true){ // mixed full
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single_reduction);
        red.cnf.reducer_use_mixed_domination_full = include_single_reduction;

        sw.restart("red");
        createForRed(red, dom_mixed_full);
        time_dom_mixed_full = sw.getTime("red");
    }


    return make_tuple(initV_data,
        make_pair(basic_data,time_basic),
        make_pair(known_red_data,time_known),
        make_pair(all_dom, time_all_dom),
        make_pair(only_fast_dom, time_only_fast_dom),
        make_pair(dom12, time_dom12),
        make_pair(dom3, time_dom3),
        make_pair(dom4, time_dom4),
        make_pair(dom5, time_dom5),
        make_pair(dom_nonsimple_arcs, time_dom_nonsimple_arcs),
        make_pair(dom_nonsimple_arcs_full, time_dom_nonsimple_arcs_full),
        make_pair(dom_mixed, time_dom_mixed),
        make_pair(dom_mixed_full, time_dom_mixed_full)
    );
}

int main(int argc, char** argv) {
    MemoryUtils::increaseStack();
    ios_base::sync_with_stdio(0);
    cin.tie(0);


    vector<string> instances_e;
    vector<string> instances_h;
    for( int i=1; i<=200; i++ ) {
        string id = "";
        if( i < 10 ) id += "0";
        if( i < 100 ) id += "0";
        id += to_string(i);

        string se = "e_" + id;
        string sh = "h_" + id;

        instances_e.push_back(se);
        instances_h.push_back(sh);
    }

    // if (false)
    { // #TEST - for running tests
        instances_h.clear();
        // instances_e.erase(instances_e.begin(), instances_e.begin() + 100);
        // instances_e.erase(instances_e.begin()+100, instances_e.end());

        // instances_h = {"e_012", "e_013"};
        // instances_e.resize(140);

        // instances_e.clear();
    }


    // vector<string> instances_all = instances_e;
    vector<string> instances_all = instances_e + instances_h;

    auto header = ExpData::getHeader();




    for ( int include_single_reduction = 0; include_single_reduction <= 1; include_single_reduction++ ) {

        DEBUG(include_single_reduction);
        string suf =  (include_single_reduction ? "_single_included" : "_single_excluded");

        constexpr bool compute = false;
        constexpr bool retain_only_absent_csvs = true;
        constexpr bool create_results_all = true;


        if(compute) {
            reverse(ALL(instances_all));


            if(retain_only_absent_csvs){
                auto fun = [&]( string f ) {
                    string s = "datasets/" + f + ".csv";
                    return !filesystem::exists(s);
                };
                auto it = stable_partition( ALL(instances_all), fun );

                instances_all.resize(it - instances_all.begin());
                DEBUG(instances_all);
            }



            omp_set_num_threads( min( (int)instances_all.size(), 10 ) );

            clog << "There are " << instances_all.size() << " instances to consider" << endl << endl;

            constexpr int chunk = 1;
            #pragma omp parallel for schedule(dynamic,chunk)
            for( int ind = 0; ind < instances_all.size(); ind++ ) {
                string s = instances_all[ind];

                string msg = "Considering instance " + s;
                msg += " in thread id: " + to_string(omp_get_thread_num() );
                clog << msg << endl;
                s = "datasets/" + s;

                if ( !filesystem::exists(s) ) {
                    clog << "Dataset " << s << " DOES NOT EXIST! skipping this instance...." << endl;
                    continue;
                }

                ifstream cur_str(s);
                auto V = Utils::readGraph(cur_str);


                auto [initV,
                    basic, knownd, alld, only_fast,
                    d12, d3, d4, d5,
                    nca, nca_full,
                    mixed, mixed_full]
                = createDataForInstance(V, include_single_reduction);

                auto [basic_data,time_basic] = basic;
                auto [known,time_known] = knownd;
                auto [all_dom, time_all_dom] = alld;
                auto [only_fast_dom, time_only_fast_dom] = only_fast;
                auto [dom12, time_dom12] = d12;
                auto [dom3, time_dom3] = d3;
                auto [dom4, time_dom4] = d4;
                auto [dom5, time_dom5] = d5;
                auto [dom_nonsimple_arcs, time_dom_nonsimple_arcs] = nca;
                auto [dom_nonsimple_arcs_full, time_dom_nonsimple_arcs_full] = nca_full;
                auto [dom_mixed, time_dom_mixex] =  mixed;
                auto [dom_mixed_full, time_dom_mixed_full] = mixed_full;

                ofstream cur_data_os(s + suf + ".csv");
                cur_data_os.precision(2);
                cur_data_os << fixed;

                {
                    int cnt = 0;
                    for( auto h : header ) cur_data_os << (cnt++ ? ", " : "") << h;
                    cur_data_os << endl;
                }

                auto writeData = [&](auto & str, auto & data, int time, bool end_of_line = false) {
                    const auto & r = data;

                    double dN_basicN = ( basic_data.N > 0 ?  data.N / basic_data.N : -1);
                    double dM_basicM = ( basic_data.M > 0 ?  data.M / basic_data.M : -1);

                    double dN_knownN = ( known.N > 0 ?  data.N / known.N : -1);
                    double dM_knownM = ( known.M > 0 ?  data.M / known.M : -1);

                    str << r.N << ", " << r.M << ", " << 1.0 * time / 1000;
                    // str << r.N << ", " << r.M << ", " << r.pi_arcs_fraction << ", " << time << ", ";
                    // str << dN_basicN << ", " << dM_basicM << ", " << dN_knownN << ", " << dM_knownM;
                    if (end_of_line) str << endl;
                    else str << ", ";
                };

                writeData(cur_data_os, initV,0);
                writeData(cur_data_os, basic_data, time_basic);
                writeData(cur_data_os, known, time_known);
                writeData(cur_data_os, all_dom, time_all_dom);
                writeData(cur_data_os, only_fast_dom, time_only_fast_dom);
                writeData(cur_data_os, dom12, time_dom12);
                writeData(cur_data_os, dom3, time_dom3);
                writeData(cur_data_os, dom4, time_dom4);
                writeData(cur_data_os, dom5, time_dom5);
                writeData(cur_data_os, dom_nonsimple_arcs, time_dom_nonsimple_arcs);
                writeData(cur_data_os, dom_nonsimple_arcs_full, time_dom_nonsimple_arcs_full);
                writeData(cur_data_os, dom_mixed, time_dom_mixex);
                writeData(cur_data_os, dom_mixed_full, time_dom_mixed_full, true);

                clog << "\tFinished processing instance " << s << endl;
            }
        }


        if(create_results_all) {

            sort(ALL(instances_all));

            auto writeForStream = [&]( auto & str, string skip_trivial ) {
                str.precision(3);
                str << fixed;

                { // writing header
                    str << "id";
                    int cnt = 1;
                    for( auto h : header ) str << (cnt++ ? ", " : "") << h;
                    str << endl;

                }

                for(auto s : instances_all) {
                    string instance_name = s;
                    s = "datasets/" + s + suf + ".csv";
                    if(!filesystem::exists(s)){
                        // str << instance_name << endl;
                        continue;
                    }

                    string header_line, data_line;
                    {
                        ifstream cur_str(s);
                        getline(cur_str, header_line);
                        getline(cur_str, data_line);

                        // DEBUG(header_line);
                        // DEBUG(data_line);
                    }

                    auto trim = [&](string & e) {
                        int p = 0, q = (int)e.size()-1;
                        while(p < e.size()) {
                            if(e[p] != ' ') break;
                            p++;
                        }
                        while(q >= 0) {
                            if( e[q] != ' ' ) break;
                            q--;
                        }
                        e = e.substr( p, q-p+1 );
                    };

                    auto entries = StandardUtils::split(data_line, ",");
                    for(auto & e : entries) trim(e);

                    int basic_N = stoi(entries[4]);
                    int known_N = stoi(entries[7]);

                    if ( skip_trivial == "basic" && basic_N == 0 ) continue;
                    if ( skip_trivial == "known" && known_N == 0 ) continue;

                    str << instance_name;
                    for(auto e : entries) str << ", " << e;
                    str << endl;
                }
            };

            ofstream res_all_str("res_all" + suf + ".csv");
            writeForStream(res_all_str, "all");

            ofstream res_nontrivial_basic( "res_all_nontrivial_basic" + suf + ".csv" );
            writeForStream(res_nontrivial_basic, "basic" );

            ofstream res_nontrivial_known( "res_all_nontrivial_known" + suf + ".csv" );
            writeForStream(res_nontrivial_known, "known");
        }
    }

    return 0;
}

