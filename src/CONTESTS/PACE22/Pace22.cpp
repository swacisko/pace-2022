//
// Created by sylwester on 12/20/21.
//

#include <CONTESTS/PACE22/Reducer.h>
#include <graphs/GraphInducer.h>
#include <graphs/VertexCover/approximation/LibMVC/fastvc.h>
#include <filesystem>
#include "CONTESTS/PACE22/heur/DFVSSolverH.h"
#include "MemoryUtils.h"
#include "getopt.h"
#include "Makros.h"
#include <omp.h>

static int time_limit_millis = 600'000;




struct ExpData {
    int N, M, pi_arcs, nonpi_arcs;
    double pi_arcs_fraction;
    double nonpi_arcs_fraction;

    static vector<string> getHeader() {
        vector<string> fields{"N", "E", "time"};
        vector<string> res;

        for( string s : { "orig", "basic", "known", "all_dom", "part-dom-1", "part-dom-2", "part-dom-3"} ) {
            for(const auto & f : fields) res.push_back(s + "-" + f);
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
    PEI, PEI, PEI, PEI, PEI, PEI>
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
    ExpData all_dom, mixed_dom_1, mixed_dom_2, mixed_domination_3, dom4, dom5;

    int time_basic, time_known, time_part_dom1, time_part_dom2, time_part_dom_3, time_dom4, time_dom5, time_all_dom;

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


    auto setKnownRules = [&]( Reducer & red) {
        red.cnf.reducer_use_core = true;
        red.cnf.reducer_use_dome = true;
        red.cnf.reducer_use_pie = true;
        red.cnf.reducer_use_inoutclique = true;

        red.cnf.reducer_use_domination_0_pinodes = true;
        red.cnf.reducer_use_folding = true;
        red.cnf.reducer_use_folding_only_for_pi_nodes = true;
        red.cnf.reducer_use_funnel = true;
        red.cnf.reducer_use_desk = true;
        red.cnf.reducer_use_twins_merge = true;
        red.cnf.reducer_use_folding_twins = true;
    };

    { // known rules
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        red.disableAllNonbasicReductions();
        setKnownRules(red);
        assignTimeLimits(red);

        sw.restart("red");
        auto reductions = red.reduce();
        for(auto * r : reductions) delete r;
        V = red.V;
        induceByNonisolated();
        known_red_data.createData(V);

        time_known = sw.getTime("red");
    }




    auto setAllForRed = [&](Reducer & red) {

        red.disableAllNonbasicReductions();
        setKnownRules(red);

        red.cnf.reducer_use_mixed_domination = !include_single_reduction;
        red.cnf.reducer_use_mixed_domination_full = !include_single_reduction;
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
        setAllForRed(red);

        sw.restart("red");
        createForRed(red, all_dom);
        time_all_dom = sw.getTime("red");

        include_single_reduction = b;
    }


    if(true){ // mixed_domination_1
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red);
        red.cnf.reducer_use_mixed_domination = include_single_reduction;

        sw.restart("red");
        createForRed(red, mixed_dom_1);
        time_part_dom1 = sw.getTime("red");
    }

    if(true){ // mixed_domination_2 - version with dominators
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red);
        red.cnf.reducer_use_mixed_domination2 = include_single_reduction;

        sw.restart("red");
        createForRed(red, mixed_dom_2);
        time_part_dom2 = sw.getTime("red");
    }

    if(true){ // // mixed_domination_3 - full DFS+backtracking based search
        Config main_cnf;
        main_cnf.sw.setLimit("main", time_limit_millis);
        main_cnf.sw.start("main");

        V = basicV;
        Reducer red(V, main_cnf);
        setAllForRed(red);
        red.cnf.reducer_use_domination_3 = include_single_reduction;

        sw.restart("red");
        createForRed(red, mixed_domination_3);
        time_part_dom_3 = sw.getTime("red");
    }



    return make_tuple(initV_data,
        make_pair(basic_data,time_basic),
        make_pair(known_red_data,time_known),
        make_pair(all_dom, time_all_dom),
        make_pair(mixed_dom_1, time_part_dom1),
        make_pair(mixed_dom_2, time_part_dom2),
        make_pair(mixed_domination_3, time_part_dom_3)
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

    if (false)
    { // #TEST - for running tests
        instances_h.clear();
        // instances_e.erase(instances_e.begin(), instances_e.begin() + 100);
        // instances_e.erase(instances_e.begin()+100, instances_e.end());

        // instances_h = {"e_012", "e_013"};
        // instances_e.resize(140);

        // instances_e.clear();
        instances_e.resize( min(instances_e.size() / 2.0, 100.0) );
    }


    // vector<string> instances_all = instances_e;
    vector<string> instances_all = instances_e + instances_h;
    vector<string> in_all_cp = instances_all;

    auto header = ExpData::getHeader();

    const string input_files_path = "pace-2022-datasets";

    // #TEST - now considering only addition of the single reduction rule
    for ( int include_single_reduction = 1; include_single_reduction <= 1; include_single_reduction++ ) {

        DEBUG(include_single_reduction);
        string suf =  (include_single_reduction ? "_single_included" : "_single_excluded");

        constexpr bool compute = true;
        constexpr bool retain_only_absent_csvs = false;
        constexpr bool create_results_all = true;


        if(compute) {
            reverse(ALL(instances_all));


            if(retain_only_absent_csvs){
                clog << "Removing all instances for which result files already exist..." << endl;

                auto fun = [&]( string f ) {
                    // string s = "datasets/" + f + ".csv";
                    string s = input_files_path + "/" + f + suf + ".csv";
                    return !filesystem::exists(s);
                };
                auto it = stable_partition( ALL(instances_all), fun );

                // DEBUG((int)(it - instances_all.begin()));

                instances_all.resize(it - instances_all.begin());
                DEBUG(instances_all);
            }



            omp_set_num_threads( min( (int)instances_all.size(), 8 ) );

            clog << "There are " << instances_all.size() << " instances to consider" << endl << endl;

            constexpr int chunk = 1;
            #pragma omp parallel for schedule(dynamic,chunk)
            for( int ind = 0; ind < instances_all.size(); ind++ ) {
                string s = instances_all[ind];

                string msg = "Considering instance " + s;
                msg += " in thread id: " + to_string(omp_get_thread_num() );
                clog << msg << endl;
                // s = "datasets/" + s;
                s = input_files_path + "/" + s;

                if ( !filesystem::exists(s) ) {
                    clog << "Dataset " << s << " DOES NOT EXIST! skipping this instance...." << endl;
                    continue;
                }

                ifstream cur_str(s);
                auto V = Utils::readGraph(cur_str);


                auto [initV, basic, knownd, alld, pd1,pd2,pd3]
                = createDataForInstance(V, include_single_reduction);

                auto [basic_data,time_basic] = basic;
                auto [known,time_known] = knownd;
                auto [all_dom, time_all_dom] = alld;
                auto [part_dom1, time_part_dom1] = pd1;
                auto [part_dom2, time_part_dom2] = pd2;
                auto [part_dom3, time_part_dom3] = pd3;

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

                    str << r.N << ", " << r.M << ", " << 1.0 * time / 1000;
                    if (end_of_line) str << endl;
                    else str << ", ";
                };

                writeData(cur_data_os, initV,0);
                writeData(cur_data_os, basic_data, time_basic);
                writeData(cur_data_os, known, time_known);
                writeData(cur_data_os, all_dom, time_all_dom);
                writeData(cur_data_os, part_dom1, time_part_dom1);
                writeData(cur_data_os, part_dom2, time_part_dom2);
                writeData(cur_data_os, part_dom3, time_part_dom3);

                clog << "\tFinished processing instance " << s << endl;
            }
        }


        if(create_results_all) {

            instances_all = in_all_cp;

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
                    // s = "datasets/" + s + suf + ".csv";
                    s = input_files_path + "/" + s + suf + ".csv";
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
                    int basic_M = stoi(entries[5]);
                    int known_N = stoi(entries[7]);
                    int known_M = stoi(entries[8]);

                    if ( skip_trivial == "basic" && basic_N == 0 ) continue;
                    if ( skip_trivial == "known" && known_N == 0 ) continue;
                    if ( skip_trivial == "known_improved_basic"
                        // && (basic_N > 0 || ( basic_N == 0 && known_N == 0 ) )
                        && known_N == basic_N && known_M == basic_M ) continue;


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

            ofstream res_known_improved_basic( "res_all_known_improved_basic" + suf + ".csv" );
            writeForStream(res_known_improved_basic, "known_improved_basic");
        }
    }

    return 0;
}

