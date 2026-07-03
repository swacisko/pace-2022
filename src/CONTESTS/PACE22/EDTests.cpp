//
// Created by sylwe on 01/07/2026.
//


#include "EDReducer.h"
#include "GraphInducer.h"
#include "GraphReader.h"
#include "GraphUtils.h"
#include "CONTESTS/PACE22/Config.h"
#include "CONTESTS/PACE22/Reducer.h"
#include "VertexCover/kernelization/KernelizerVC.h"
#include <ranges>

struct ExpData {
    int N0,M0,N1,M1,N2,M2;
    int type1_constraints;
    int type2_constraints;

    int red_init_time_millis;
    int red_noned_time_millis;
    int red_ed_time_millis;

    int solver_max_time_sec;
    int solver_time_granularity;
    int solver_repeats;
    VD noned_results;
    VD ed_results;

    int ed_nodes_reduced = 0;
    int ed_edges_removed = 0;
    int ed_t1_inference_rules_added = 0;
    int ed_total_t2_inference_rules_created = 0;
    int ed_t2_inference_rules_added = 0;

    map<string,string> getEntries() {
        map<string,string> res;
        res["N0"] = to_string(N0);
        res["M0"] = to_string(M0);
        res["N1"] = to_string(N1);
        res["M1"] = to_string(M1);
        res["N2"] = to_string(N2);
        res["M2"] = to_string(M2);

        res["type1_constraints"] = to_string(type1_constraints);
        res["type2_constraints"] = to_string(type2_constraints);

        res["red_init_time_millis"] = to_string(red_init_time_millis);
        res["red_noned_time_millis"] = to_string(red_noned_time_millis);
        res["red_ed_time_millis"] = to_string(red_ed_time_millis);

        res["solver_max_time_sec"] = to_string(solver_max_time_sec);
        res["solver_time_granularity"] = to_string(solver_time_granularity);
        res["solver_repeats"] = to_string(solver_repeats);

        stringstream str;
        for (int d : noned_results) str << d << " ";
        res["noned_results"] = str.str();
        str.clear(); str.str("");

        for (int d : ed_results) str << d << " ";
        res["ed_results"] = str.str();
        str.clear(); str.str("");

        return res;
    }


    void writeToFile(ostream & str) {
        auto mapa = getEntries();
        vector<string> header = { "N0", "M0", "N1", "M1", "N2", "M2",
            "type1_constraints", "type2_constraints",
        "red_init_time_millis", "red_noned_time_millis", "red_ed_time_millis",
        "solver_max_time_sec", "solver_time_granularity", "solver_repeats",
        "noned_results", "ed_results"
        };
        for (auto k : views::keys(mapa)) header.push_back(k);

        auto writeLine = [&](vector<string> & l) {
            int cnt = 0;
            for ( auto s : l ) {
                if (cnt++) str << " ";
                str << s;
            }
            str << "\n";
        };

        writeLine(header); // write header
        vector<string> line;
        for (auto k : header) line.push_back(mapa[k]); // create entries in the order of the header
        writeLine(line); // write entries
    }
};

pair<VI,VD> solveInstanceUsingSatBasedSolver(VPII constraints, int max_time_sec, int granularity, int repeats, string alg = "cpsat") {
    // times[i] is the result found by the solver after time (i+1) * granularity seconds
    const int I = ceil(max_time_sec / granularity);
    VVI iteration_res_times(I);
    VD res_times(I, -1);
    VI res; // valid solution for provided constraints

    // TODO: implement solving constraitns using CPSAT or EvalMaxSAT or whatever other SAT-based solver will be used...
    for (int rep=0; rep<repeats; rep++) {

        auto solve = [&]() {
            if (alg == "cpsat") {

            }

            if (alg == "evalmaxsat") {

            }

            if (alg == "some_anytime_maxsat_solver...") {

            }

            return VI(I,-1);
        };

        VI iter_res_times = solve();
        for ( int i=0; i<iter_res_times.size(); i++ ) iteration_res_times[i].push_back(iter_res_times[i]);
    }

    for (int i=0; i<iteration_res_times.size(); i++ ) {
        const auto & irt = iteration_res_times[i];
        if ( !irt.empty() ) {
            res_times[i] = accumulate(ALL(irt),0.0) / irt.size();
        }
    }

    return make_pair(res,res_times);
}

pair<VPII,ExpData> runVCTestforGraph(VVI V, int solver_max_time_sec, int solver_time_granularity, int solver_repeats, string alg) {
    int N = V.size();

    ExpData exp_data;
    exp_data.solver_time_granularity = solver_time_granularity;
    exp_data.solver_max_time_sec = solver_max_time_sec;
    exp_data.solver_repeats = solver_repeats;

    exp_data.N0 = N;
    exp_data.M0 = GraphUtils::countEdges(V,true);

    DEBUG(PII(exp_data.N0,exp_data.M0));

    if (false){
        Stopwatch sw;
        sw.start("main");
        KernelizerVC kern;
        auto [kern_nodes, edges_removed] = kern.initialKernelization(V);
        VB helper(N);

        auto revV = GraphUtils::reverseGraph(V);
        Utils::removeNodes(V, revV, kern_nodes, helper);

        V = GraphInducer::induceByNonisolatedNodes(V).V;
        N = V.size();
        sw.stop("main");

        exp_data.red_init_time_millis = sw.getTime("main");
    }

    exp_data.N1 = N;
    exp_data.M1 = GraphUtils::countEdges(V,true);

    DEBUG(PII(exp_data.N1,exp_data.M1));

    { // measuring just the VC reduction time WITHOUT ED rule, and the solver results for the non-ed reduced graph
        Stopwatch sw;
        sw.start("main");
        Config cnf;
        cnf.disableAllNonbasicReductions();
        cnf.reducer_use_folding = cnf.reducer_use_folding_twins = cnf.reducer_use_funnel = cnf.reducer_use_desk = true;
        cnf.reducer_use_unconfined =  cnf.reducer_use_twins_merge = cnf.reducer_use_domination = true;
        cnf.reducer_use_general_folding = true; cnf.reducer_max_general_folding_antiedges = 2; cnf.reducer_max_general_folding_neighborhood_size = 10;
        Reducer red(V,cnf);
        auto to_lift = red.reduce();
        sw.stop("main");

        int solution_lift_overhead = Reducer::getReductionsSizeDiff(to_lift);

        exp_data.red_noned_time_millis = exp_data.red_init_time_millis + sw.getTime("main");

        VPII constraints = GraphUtils::getGraphEdges(V);
        auto [res,times] = solveInstanceUsingSatBasedSolver(constraints, solver_max_time_sec, solver_time_granularity, solver_repeats, alg);
        for (auto& d : times) d += solution_lift_overhead;
        exp_data.noned_results = times;
    }


    // now measuring VC reduction time WITH ED rule
    Stopwatch sw;
    sw.start("main");
    Config cnf;
    cnf.disableAllNonbasicReductions();
    cnf.reducer_use_folding = cnf.reducer_use_folding_twins = cnf.reducer_use_funnel = cnf.reducer_use_desk = true;
    cnf.reducer_use_unconfined =  cnf.reducer_use_twins_merge = cnf.reducer_use_domination = true;
    cnf.reducer_use_general_folding = true; cnf.reducer_max_general_folding_antiedges = 2; cnf.reducer_max_general_folding_neighborhood_size = 10;
    cnf.reducer_use_ed = true;
    Reducer red(V,cnf);
    auto to_lift = red.reduce();
    sw.stop("main");

    exp_data.red_ed_time_millis = exp_data.red_init_time_millis + sw.getTime("main");

    exp_data.ed_nodes_reduced = red.ed_nodes_reduced;
    exp_data.ed_edges_removed = red.ed_edges_removed; assert(red.ed_edges_removed == 0);
    exp_data.ed_t1_inference_rules_added = red.ed_t1_inference_rules_added;
    exp_data.ed_t2_inference_rules_added = red.ed_t2_inference_rules_added;
    exp_data.ed_total_t2_inference_rules_created = red.ed_total_t2_inference_rules_created;

    auto newV = red.V;
    auto indg = GraphInducer::induceByNonisolatedNodes(newV);

    V = indg.V;
    N = V.size();

    exp_data.N2 = N;
    exp_data.M2 = GraphUtils::countEdges(V,true);

    DEBUG(PII(exp_data.N2,exp_data.M2));



    // now solve the reduced problem using constraints...


    VPII constraints;
    // TODO: create constraints for the ED-reduced graph.
    // Take into account type2 constraints as well, type1 constrains were already added to the graph

    int solution_lift_overhead = Reducer::getReductionsSizeDiff(to_lift);
    auto [res,times] = solveInstanceUsingSatBasedSolver(constraints,solver_max_time_sec,solver_time_granularity, solver_repeats);
    for (auto& d : times) d += solution_lift_overhead;
    exp_data.ed_results = times;



    return {constraints,exp_data};
}

VVI getTestV1() {
    /**
     v a b x c d p q k y r  s  l
     0 1 2 3 4 5 6 7 8 9 10 11 12
     *
     **/
    VPII edges;
    {
        int v = 0, a = 1, b = 2, x = 3, c = 4, d = 5, p = 6, q = 7, k = 8, y = 9, r = 10, s = 11, l = 12;
        edges = {
            {v,a}, {v,b}, {v,c}, {v,d},
            {a,p}, {a,q}, {a,b}, {b,q},
            {c,x}, {c,y}, {c,d},
            {d,r}, {d,s},
            {p,x}, {p,l},
            {q,x}, {q,l},
            {x,k}, {x,y},
            {k,l},
            {y,r},
            {r,l}, {r,s},
            {s,l}
        };
    }

    return GraphUtils::getGraphForEdges(edges);
}



int main() {
    ios_base::sync_with_stdio(0);

    // VVI V = GraphReader::readGraphStandardEdges(cin);
    VVI V = getTestV1();


    int solver_max_time_sec = 300;
    int solver_time_granularity = 5;
    int solver_repeats = 5;
    string alg = "cpsat";
    auto [constraints,exp_data] = runVCTestforGraph(V,  solver_max_time_sec, solver_time_granularity,solver_repeats, alg);




    exp_data.writeToFile(cout);


    return 0;
}
