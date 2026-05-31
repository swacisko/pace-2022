//
// Created by sylwe on 08/04/2026.
//

#include "CombinatoricUtils.h"
#include "CpsatExp1.h"
#include "ExpConfig.h"
#include "GraphReader.h"
#include "GraphUtils.h"
#include "MemoryUtils.h"
#include "StandardUtils.h"
#include "Stopwatch.h"
#include "CONTESTS/PACE22/Utils.h"
#include "CONTESTS/PACE22/heur/DFVSSolverH.h"
#include "dreyfvs/dfvs.h"
#include "dreyfvs/graph.h"

VVI generateSets(int N, int M, int D) {
    VVI res;
    IntGenerator rnd;
    for ( int i=0; i<M; i++ ) {
        int d = 2 + rnd.nextInt(D-1);
        VI S = CombinatoricUtils::getRandomSubset(N-1, d);
        res.push_back(S);
    }
    return res;
}

static VI solveDreyFVS(VVI V, ExpConfig cnf) {
    Stopwatch sw;
    sw.start("dreyfvs");

    stringstream str;
    int N = V.size();
    int M = GraphUtils::countEdges(V,true);
    str << N << " " << M << " 0" << "\n";
    for (int i=0; i<N; i++) {
        int cnt = 0;
        for (int d : V[i]) str << (cnt++ ? " " : "") << d+1;
        str << "\n";
    }

    clog << "Starting DreyFVS" << endl;

    ExpData exp_data;
    Graph g = Graph::from_istream(str);
    auto sol = computeDFVS(g, cnf.max_time_sec, exp_data, true);
    assert(Utils::isFVS(V,sol));
    sw.stop("dreyfvs");

    return sol;
}


VI solveDiVerSeS(VVI V, ExpConfig cnf) {

    Stopwatch sw;
    sw.setLimit("diverses", cnf.max_time_sec * 1000);
    sw.start("diverses");

    int best_res_size = 0;
    VI best_res;

    int iter_id = 0;
    while (!sw.tle("diverses")) {
        Config diverses_cnf;
        diverses_cnf.write_logs = false;
        int time_left_millis = ( sw.getLimit("diverses") - sw.getTime("diverses") );
        diverses_cnf.sw.setLimit("main", time_left_millis);
        diverses_cnf.sw.start("main");
        diverses_cnf.disableAllNonbasicReductions();
        diverses_cnf.solverh_use_reductions_AF = false;

        Stopwatch sw2;
        sw2.start("iteration");
        DFVSSolverH sh(diverses_cnf);
        auto res = sh.solveForGraph(V);
        sw2.stop("iteration");

        if (best_res_size == 0 || res.size() < best_res_size) {
            best_res_size = res.size();
            best_res = res;
        }

        iter_id++;
        clog << "After iteration " << iter_id << ", res.size(): " << res.size() << ", best_res_size: " << best_res_size << endl;
    }

    DEBUG(best_res_size);

    return best_res;
}

VI testForTopologicalOrder(VI order, VVI sets, string alg = "diverses", int time_sec = 3) {
    int N = order.size();
    VI in_order(N,0);
    for (auto [i,d] : views::enumerate(order)) in_order[d] = i;

    set<PII> arcs;

    for ( auto & cyc : sets ) {
        assert(cyc.size() >= 2);
        sort(ALL(cyc), [&](int a, int b){ return in_order[a] < in_order[b]; });
        for (int i=0; i<(int)cyc.size()-1; i++) arcs.insert( {cyc[i],cyc[i+1]});
        arcs.insert( {cyc.back(),cyc.front()} );
    }

    VPII A(ALL(arcs));
    VVI V = GraphUtils::getGraphForEdges(A,true);

    clog << "\t created graph has " << A.size() << " arcs" << endl;

    VI sol;
    ExpConfig cnf;
    cnf.max_time_sec = time_sec;

    if (alg == "dreyfvs") {
        clog << "\t Looking for a FVS of the digraph using DreyFVS for " << cnf.max_time_sec << " seconds" << endl;
        sol = solveDreyFVS(V, cnf);
    }

    if (alg == "diverses") {
        clog << "\t Looking for a FVS of the digraph using DiVerSeS for " << cnf.max_time_sec << " seconds" << endl;
        sol = solveDiVerSeS(V, cnf);
    }

    clog << "\t found FVS of size " << sol.size() << endl;

    return sol;
}

void testAlgorithms() {

    const int MAXD = 20;
    const int N = 1e3;
    const int MAXM = 1e6;

    const int MAIN_REPS = 20;

    for (int r = 0; r < MAIN_REPS; r++) {
        ENDL(5); ENDLS(10," -->MAIN REPS<-- "); ENDL(5);

        for (int M = 10*N; M <= MAXM; M *= 2) {
            ENDL(5); ENDLS(50,"*"); ENDL(5);
            clog << "Creating instance for M = " << M << " sets" << endl;
            VVI sets = generateSets(N,M,MAXD);

            {
                map<int,int> sets_sizes;
                for (auto & cyc : sets) sets_sizes[cyc.size()]++;
                DEBUG(sets_sizes);
                int total_elements = accumulate(ALL(sets_sizes),0,[](int a, auto & b){ return a + b.first*b.second; });
                DEBUG(total_elements);
            }

            const int SEC = 10;
            clog << "Running CP-SAT for " << SEC << " seconds to find HS" << endl;
            VI hs;
            {
                SatParameters params;
                params.set_max_time_in_seconds(SEC);
                params.set_num_search_workers(6);

                Model solver_model;
                solver_model.Add(NewSatParameters(params));

                CpModelBuilder model;
                vector<BoolVar> nodes;
                for (int i=0; i<N; i++) nodes.push_back(model.NewBoolVar());
                for ( auto & c : sets ) {
                    vector<BoolVar> cyc_vars;
                    for ( int d : c ) cyc_vars.push_back(nodes[d]);
                    model.AddAtLeastOne(cyc_vars);
                }

                auto objective = LinearExpr::Sum(nodes);
                model.Minimize(objective);
                auto model_proto = model.Build();
                CpSolverResponse response = SolveCpModel(model_proto, &solver_model);
                for (int i = 0; i < nodes.size(); ++i) if (SolutionBooleanValue(response, nodes[i])) hs.push_back(i);
            }

            clog << "\t found HS of size " << hs.size() << endl;

            VI ord = CombinatoricUtils::getRandomPermutation(N);
            VI cnt(N,0);
            for (auto & cyc : sets) for (auto & d : cyc) cnt[d]++;

            {
                clog << "\tCreating order by setting elements in HS at the end" << endl;
                auto in_hs = StandardUtils::toVB(N,hs);
                sort(ALL(ord),[&](int a, int b) { return in_hs[a] && !in_hs[b]; });
                auto fvs = testForTopologicalOrder(ord, sets);
                clog << "\t Found FVS of size " << fvs.size() << " for a graph created for HS of size " << hs.size() << endl;
            }


            if (true){
                VI best_res;

                IntGenerator rnd;
                const int T = 100;
                clog << "\t Running " << T << " iterations with HS order rearrangement" << endl;
                sort(ALL(ord),[&](int a, int b) { return cnt[a] > cnt[b]; });

                for ( int t=0; t<T; t++ ) {
                    auto sol = testForTopologicalOrder(ord, sets, "diverses", 3 + (7.0 * (t+1) / T));
                    // auto sol = testForTopologicalOrder(ord, sets, "dreyfvs", 2 + (3.0 * (t+1) / T));

                    if (best_res.empty() || sol.size() <= best_res.size()) best_res = sol;

                    clog << "\tAfter iteration " << t+1 << " found FVS of size: " << sol.size() << endl;

                    auto in_best_res = StandardUtils::toVB(N,best_res);
                    StandardUtils::shuffle(ord,rnd);
                    // if (rnd.nextInt(4) == 0) sort(ALL(ord), [&](int a, int b){ return cnt[a] < cnt[b]; });
                    // if (rnd.nextInt(4) == 1) sort(ALL(ord), [&](int a, int b){ return cnt[a] > cnt[b]; });
                    stable_partition(ALL(ord), [&](int d){ return in_best_res[d]; });


                    // if (t & 1)
                    { // moving to the front elements in solution for which for some set only they hit it
                        VB sole_element(N,false);
                        for ( auto & cyc : sets ) {
                            int c = 0;
                            for (int d : cyc) c += in_best_res[d];
                            assert(c>=1);
                            int el = -1;
                            for (int d : cyc) if ( in_best_res[d] ) el=d;
                            assert(el != -1);
                            if (c == 1) sole_element[el] = true;
                        }
                        // stable_partition(ALL(ord), [&](int d){ return sole_element[d]; });
                        stable_sort(ALL(ord), [&](int a, int b) {
                            if ( in_best_res[a] != in_best_res[b] ) return (bool)in_best_res[a];
                            else {
                                if (t&1) return (bool)sole_element[b];
                                else return !sole_element[b];
                            }
                        });

                        assert(is_partitioned(ALL(ord),[&](int d){ return in_best_res[d]; }));
                    }

                }
            }


            if (false){
                IntGenerator rnd;
                const int T = 50;
                clog << "\t Running " << T << " iterations with HS order rearrangement using reinforcement learning" << endl;
                sort(ALL(ord),[&](int a, int b) { return cnt[a] > cnt[b]; });
                int prev_res = inf;
                VI occ(N,0);

                for ( int t=0; t<T; t++ ) {
                    auto sol = testForTopologicalOrder(ord, sets);
                    clog << "\tAfter iteration " << t+1 << " found FVS of size: " << sol.size() << endl;

                    auto in_sol = StandardUtils::toVB(N,sol);
                    VI tab = occ;
                    if (t & 1) for (int & d : tab) d += rnd.nextInt(5);
                    sort(ALL(ord), [&](int a, int b){ return tab[a] < tab[b]; });

                    prev_res = sol.size();
                    for (int d : sol) occ[d]++;
                }
            }

        }
    }

}


int main(int argc, char** argv){
    MemoryUtils::increaseStack();

    testAlgorithms();


    return 0;
}
