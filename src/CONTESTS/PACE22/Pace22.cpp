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
        vector<string> fields{"N", "M", "pi_arcs", "npi_arcs", "paf"};
        vector<string> res;

        // writeData(cur_data_os, dom_and_liftables_no_funnel);
        // writeData(cur_data_os, dom_and_liftables_no_desk);
        // writeData(cur_data_os, dom_and_liftables_no_folding);
        // writeData(cur_data_os, dom_and_liftables_no_funnel_no_folding_twins);
        // writeData(cur_data_os, dom_and_liftables_no_funnel_no_general_folding);
        // writeData(cur_data_os, dom_and_liftables_no_funnel_no_cycle_folding);
        // writeData(cur_data_os, dom_and_liftables_no_funnel_no_rev_triangles);
        // writeData(cur_data_os, dom_and_liftables_no_funnel_no_edge_blocker);
        // writeData(cur_data_os, dom_and_liftables_no_funnel_no_bipartite_blocker);

        for( string s : { "orig", "basic", "known", "dom", "dom+lift", "all-non-dom6",
            "all-no-funnel", "no-desk", "no-folding",
            "no-folding-twins", "no-general-folding", "no-cycle-folding",
            "no-rev-triangles", "no-edge-blocker", "no-bipartite-blocker"} ) {
            for(auto f : fields) {
                res.push_back(s + "-" + f);
            }
        }

        res.push_back( "time-dom-sec" );
        res.push_back("time-all-sec");

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
    PEI, PEI, PEI, PEI, PEI, PEI, PEI, PEI, PEI, PEI, PEI>
    createDataForInstance(VVI V) {

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



    bool include_single = false; // if true, then only single rule will be included, otherwise excluded
    auto setAllForRed = [&](auto & red,
        const bool use_dom4 = false, const bool use_dom5 = false,
        const bool use_nonsimple_cycle_arcs_full = false,
        bool use_mixed_domination_full = false) {

        red.disableAllNonbasicReductions();

        red.cnf.reducer_use_core = !include_single;
        red.cnf.reducer_use_dome = !include_single;
        red.cnf.reducer_use_pie = !include_single;
        red.cnf.reducer_use_inoutclique = !include_single;
        red.cnf.reducer_use_nonsimple_cycle_arcs = !include_single;

        red.cnf.reducer_use_nonsimple_cycle_arcs_full = use_nonsimple_cycle_arcs_full;

        red.cnf.reducer_use_domination = !include_single;
        red.cnf.reducer_use_domination_3 = !include_single;

        red.cnf.reducer_use_domination_4 = use_dom4;
        red.cnf.reducer_use_domination_5 = use_dom5;

        red.cnf.reducer_use_mixed_domination = !include_single;

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


    VVI domV;

    if(true){ // all_dom
        auto b = include_single;
        include_single = false;
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = domV;
        Reducer red(V, main_cnf);
        setAllForRed(red, true, true, true, true);

        sw.restart("red");
        createForRed(red, all_dom);
        time_all_dom = sw.getTime("red");
        include_single = b;
    }

    if(true){ // fast_dom
        auto b = include_single;
        include_single = false;

        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = domV;
        Reducer red(V, main_cnf);
        setAllForRed(red, false,false, false, false);

        sw.restart("red");
        createForRed(red, all_dom);
        time_only_fast_dom = sw.getTime("red");

        include_single = b;
    }


    if(true){ // reducer_use_domination
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = domV;
        Reducer red(V, main_cnf);
        setAllForRed(red,!include_single);
        red.cnf.reducer_use_domination = include_single;

        sw.restart("red");
        createForRed(red, dom12);
        time_dom12 = sw.getTime("red");
    }

    if(true){ // reducer_use_domination_3
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = domV;
        Reducer red(V, main_cnf);
        setAllForRed(red,!include_single);
        red.cnf.reducer_use_domination_3 = include_single;

        sw.restart("red");
        createForRed(red, dom3);
        time_dom3 = sw.getTime("red");
    }

    if(true){ // reducer_use_domination_4
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = domV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single);
        red.cnf.reducer_use_domination_4 = include_single;

        sw.restart("red");
        createForRed(red, dom4);
        time_dom4 = sw.getTime("red");
    }

    if(true){ // reducer_use_domination_5
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = domV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single);
        red.cnf.reducer_use_domination_5 = include_single;

        sw.restart("red");
        createForRed(red, dom5);
        time_dom5 = sw.getTime("red");
    }

    if(true){ // nonsimple arcs
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = domV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single);
        red.cnf.reducer_use_nonsimple_cycle_arcs = include_single;

        sw.restart("red");
        createForRed(red, dom_nonsimple_arcs);
        time_dom_nonsimple_arcs = sw.getTime("red");
    }

    if(true){ // nonsimple arcs full
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = domV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single);
        red.cnf.reducer_use_nonsimple_cycle_arcs_full = include_single;

        sw.restart("red");
        createForRed(red, dom_nonsimple_arcs_full);
        time_dom_nonsimple_arcs_full = sw.getTime("red");
    }

    if(true){ // mixed domination
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = domV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single);
        red.cnf.reducer_use_mixed_domination = include_single;

        sw.restart("red");
        createForRed(red, dom_mixed);
        time_dom_mixed = sw.getTime("red");
    }

    if(true){ // mixed full
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = domV;
        Reducer red(V, main_cnf);
        setAllForRed(red, !include_single);
        red.cnf.reducer_use_mixed_domination_full = include_single;

        sw.restart("red");
        createForRed(red, dom_mixed_full);
        time_dom_mixed_full = sw.getTime("red");
    }


    return make_tuple(initV_data,
        make_pair(basic_data,time_basic),
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

    // {
    //     instances_h.clear();
    //     // instances_e.erase(instances_e.begin(), instances_e.begin() + 100);
    //     // instances_e.erase(instances_e.begin()+100, instances_e.end());
    //
    //     // instances_h = { "h_164", "h_188", "h_189"};
    //     // instances_h = { "h_177", "h_175", "h_176"};
    //     // instances_h = {"h_164", "h_177", "h_175", "h_176"};
    //     // instances_e = {"e_152", "e_154", "e_184", "e_155", "e_157"};
    //
    //     // instances_e.clear();
    // }


    // vector<string> instances_all = instances_e;
    vector<string> instances_all = instances_e + instances_h;

    // clog << "All instances: " << endl;
    // for( auto s : instances_all ) DEBUG(s);

    auto header = ExpData::getHeader();

    // instances_all.resize(250);
    // instances_all = StandardUtils::slice(instances_all, 200 + 197, 400);


    // omp_set_num_threads(1);
    // DEBUG(omp_get_max_threads());
    // exit(1);

    constexpr bool compute = true;
    constexpr bool retain_only_absent_csvs = false;

    constexpr bool create_results_all = true;



    if(compute) {
        reverse(ALL(instances_all));


        if(retain_only_absent_csvs){
            auto fun = [&]( string f ) {
                string s = "datasets/" + f + ".csv";
                return !filesystem::exists(s);
            };
            auto it = partition( ALL(instances_all), fun );

            // for(auto s : instances_all){
            //     if(!filesystem::exists( "datasets/" + s + ".csv" )){
            //         DEBUG(s);
            //     }
            // }

            // DEBUG(instances_all);
            // DEBUG( (int)(it - instances_all.begin()) );
            instances_all.resize(it - instances_all.begin());
            DEBUG(instances_all);


            // DEBUG(filesystem::exists("datasets/h_175.csv"));
            // DEBUG(filesystem::exists("datasets/h_176.csv"));
            // exit(2);
        }

        omp_set_num_threads( min( (int)instances_all.size(), 24 ) );

        clog << "There are " << instances_all.size() << " instances to consider" << endl << endl;

        constexpr int chunk = 1;
        // #pragma omp parallel for schedule(dynamic,chunk) num_threads(14)
        // #pragma omp parallel for num_threads(24) schedule(dynamic,chunk)
        #pragma omp parallel for schedule(dynamic,chunk)
        // for( auto s : instances_all ) {
        for( int ind = 0; ind < instances_all.size(); ind++ ) {
            string s = instances_all[ind];

            string msg = "Considering instance " + s;
            msg += " in thread id: " + to_string(omp_get_thread_num() );
            clog << msg << endl;
            s = "datasets/" + s;

            ifstream cur_str(s);
            auto V = Utils::readGraph(cur_str);

        //     initV_data,
        // make_pair(basic_data,time_basic),
        // make_pair(all_dom, time_all_dom),
        // make_pair(only_fast_dom, time_only_fast_dom),
        // make_pair(dom12, time_dom12),
        // make_pair(dom3, time_dom3),
        // make_pair(dom4, time_dom4),
        // make_pair(dom5, time_dom5),
        // make_pair(dom_nonsimple_arcs, time_dom_nonsimple_arcs),
        // make_pair(dom_nonsimple_arcs_full, time_dom_nonsimple_arcs_full),
        // make_pair(dom_mixed, time_dom_mixed),
        // make_pair(dom_mixed_full, time_dom_mixed_full)

            auto [initV,
                basic, alld, only_fast,
                d12, d3, d4, d5,
                nca, nca_full,
                mixed, mixed_full]
            = createDataForInstance(V);

            auto [basic_data,time_basic] = basic;
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

            ofstream cur_data_os(s + ".csv");
            cur_data_os.precision(3);
            cur_data_os << fixed;

            {
                int cnt = 0;
                for( auto h : header ) cur_data_os << (cnt++ ? ", " : "") << h;
                cur_data_os << endl;
            }

            auto writeData = [&](auto & str, auto & data, int time, bool end_of_line = false) {
                const auto & r = data;
                str << r.N << ", " << r.M << ", " << r.pi_arcs << ", " << r.nonpi_arcs << ", " << r.pi_arcs_fraction;
                str << r.N << ", " << r.M << ", " << time;
                if (end_of_line) str << endl;
                else str << ", ";
            };

            writeData(cur_data_os, initV,0);
            writeData(cur_data_os, basic_data, time_basic);
            writeData(cur_data_os, all_dom, time_all_dom);
            writeData(cur_data_os, only_fast_dom, time_only_fast_dom);

            writeData(cur_data_os, dom12, time_dom12);
            writeData(cur_data_os, dom3, time_dom3);
            writeData(cur_data_os, dom4, time_dom4);
            writeData(cur_data_os, dom5, time_dom5);
            writeData(cur_data_os, dom_nonsimple_arcs, time_dom_nonsimple_arcs);
            writeData(cur_data_os, dom_nonsimple_arcs_full, time_dom_nonsimple_arcs_full);
            writeData(cur_data_os, dom_mixed, time_dom_mixex);
            writeData(cur_data_os, dom_mixed_full, time_dom_mixed_full);

            clog << "\tFinished processing instance " << s << endl;
        }
    }

    exit(1); // convert below !!

    if(create_results_all) {
        ofstream res_all_str("res_all.csv");
        res_all_str.precision(3);
        res_all_str << fixed;
        { // writing header
            res_all_str << "id";
            int cnt = 1;
            for( auto h : header ) res_all_str << (cnt++ ? ", " : "") << h;
            res_all_str << ", dN / bN, dE / bE, dN / kN, dE / kE";
            res_all_str << ", dliftN / kN, dliftE / kE, allN / kN, allE / kE";
            res_all_str << endl;
        }

        for(auto s : instances_all) {
            string instance_name = s;
            s = "datasets/" + s + ".csv";
            if(!filesystem::exists(s)){
                res_all_str << instance_name << endl;
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
            // DEBUG(entries);
            // for(auto e : entries) clog << e << endl;
            for(auto & e : entries) trim(e);
            // DEBUG(entries);
            // for(auto e : entries) clog << e << endl;
            // exit(1);

            res_all_str << instance_name;
            for(auto e : entries) res_all_str << ", " << e;


            int basicN = stoi(entries[5]);
            int knownN = stoi(entries[10]);
            int domN = stoi(entries[15]);
            int dliftN = stoi(entries[20]);
            int allN = stoi(entries[25]);
            double domN_basicN = (basicN ?  1.0*domN / basicN : 0 );
            double domN_knownN = (knownN ? 1.0*domN / knownN : 0 );
            double dliftN_knownN = (knownN ? 1.0*dliftN / knownN : 0 );
            double allN_knownN = (knownN ? 1.0*allN / knownN : 0 );

            int basicE = stoi(entries[6]);
            int knownE = stoi(entries[11]);
            int domE = stoi(entries[16]);
            int dliftE = stoi(entries[21]);
            int allE = stoi(entries[26]);
            double domE_basicE = (basicE ?  1.0*domE / basicE : 0 );
            double domE_knownE = (knownE ? 1.0*domE / knownE : 0 );
            double dliftE_knownE = (knownE ? 1.0*dliftE / knownE : 0 );
            double allE_knownE = (knownE ? 1.0*allE / knownE : 0 );

            res_all_str << ", " << domN_basicN << ", " << domE_basicE << ", " << domN_knownN << ", " << domE_knownE;
            res_all_str << ", " << dliftN_knownN << ", " << dliftE_knownE << ", " << allN_knownN << ", " << allE_knownE;

            res_all_str << endl;
        }
    }

    return 0;
}



/**
 * MAIN ALGORITHM MAY NOT RUN DETERMINISTICALLY!
 * That is because some algorithms, like NuMVC, are run for e.g. 1'000 milliseconds. They may not find
 * the same solutions if run twice with the same time limit.
 */
int main3(int argc, char** argv){
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

