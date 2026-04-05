//
// Created by sylwester on 12/20/21.
//

#include <filesystem>
#include <utility>
#include "MemoryUtils.h"
#include "getopt.h"
#include "GraphReader.h"
#include "GraphUtils.h"
#include "Stopwatch.h"
#include "StandardUtils.h"
#include "CONTESTS/PACE22/Utils.h"

#include "ortools/sat/cp_model.h"

using namespace operations_research;
using namespace operations_research::sat;

constexpr int inf = 1e9+1;



constexpr bool lite_track = false;
bool MUTE_MODE = false;

int threads = 6;
static int time_limit_millis = 15'000;


void initializeParams(int argc, char **argv) {
    string time_limit = "time";
    string threads = "threads";
    string quiet = "quiet";

    static struct option long_options[] = {
            {time_limit.c_str(), required_argument, 0, 0},
            {threads.c_str(), required_argument, 0, 0},
            {quiet.c_str(), required_argument, 0, 0},
            {0, 0,                                           0, 0}
    };

    while (1) {
        int option_index = 0;
        int c;
        string option, option_name;

        c = getopt_long(argc, argv, "l:", long_options, &option_index);
        if (c == -1) break;
        switch (c) {
            case 0:
                option = string(optarg);
                option_name = string(long_options[option_index].name);

                if(option_name == time_limit) time_limit_millis = stoi(option);
                if(option_name == quiet) if( option == "true" ) MUTE_MODE = true;
                if(option_name == threads) if( option == "true" ) ::threads = stoi(option);

                break;
            case '?':
                break;
            default:
                printf("?? getopt returned character code 0%o ??\n", c);
        }
    }
}


/**
 * Iterative Hitting-Set approach - the most straightforward type.
 * For each cycle length L, starting from 1, considers all cycles of length <= L, then finds HS of those cycles.
 * If the found HS is not a FVS of V, then increases L and repeats.
 */
pair<string,VI> solveCPSAT1(VVI V, int max_l, int time_limit_millis, bool find_optimal = false) {
    clog << "Solving using CPSAT, model v1" << endl;

    int N = V.size();

    VI prev_res;

    auto solveHS = [&](auto & cycles) -> VI {
        SatParameters params;
        if (!find_optimal) params.set_max_time_in_seconds(time_limit_millis / 1000.0);
        params.set_num_search_workers(threads);   // deterministic runs
        params.set_log_search_progress(false);
        Model solver_model;
        solver_model.Add(NewSatParameters(params));

        CpModelBuilder model;
        vector<BoolVar> nodes;
        for (int i=0; i<N; i++) nodes.push_back(model.NewBoolVar());
        for ( auto & c : cycles ) {
            vector<BoolVar> cyc_vars;
            for ( int d : c ) cyc_vars.push_back(nodes[d]);
            model.AddAtLeastOne(cyc_vars);
        }

        // VB in_res = StandardUtils::toVB(N,prev_res);
        // for(int i=0; i<N; i++) model.AddHint(nodes[i],in_res[i]);
        for(int d : prev_res) model.AddHint(nodes[d],1);

        if ( find_optimal ) model.AddGreaterOrEqual(LinearExpr::Sum(nodes), (int)prev_res.size());

        model.Minimize(LinearExpr::Sum(nodes));
        CpSolverResponse response = SolveCpModel(model.Build(), &solver_model);

        VI res;
        for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
        return res;
    };

    Stopwatch timer;
    timer.setLimit("all_iters", max_l * time_limit_millis);
    timer.start("all_iters");

    if(find_optimal) max_l = inf;

    VI res;
    for (int L=3; L <= max_l ; L++) {
        if ( !find_optimal && timer.tle("all_iters")) {
            clog << "Solver did not find optimal value for given set of cycles in admissible time" << endl;
            break;
        }


        clog << endl << "Considering cycles of length <= " << L << endl;
        auto cycles = Utils::getAllSimpleCycles3(V,L );
        clog << "\t there are " << cycles.size() << " such cycles" << endl;

        VI sol = solveHS(cycles);
        if ( Utils::isFVS(V,sol) ){ res = sol; break; }

        prev_res = sol;

        clog << "\t found sol.size(): " << sol.size() << ", but it is not a FVS, increasing cycle length" << endl;
    }

    string status;
    if( Utils::isFVS(V,res) ) {
        if( find_optimal ) status = "OPTIMAL";
        else status = "FEASIBLE";
    }
    else status = "INCORRECT";

    return {status, res};
}

/**
 * Uses the MTZ formulation augmented with all cycles of length <= 3 to speed up propagation.
 */
pair<string,VI> solveCPSAT2(VVI V, int MAX_RANK_VALUE = 1e9, VI init_sol = {}, bool use_only_lns = false) {
    clog << "Solving using CPSAT, model v2" << endl;

    int N = V.size();

    VB in_init_sol = StandardUtils::toVB(N,init_sol);
    if(!init_sol.empty()) MAX_RANK_VALUE = N;

    CpModelBuilder model;

    vector<IntVar> ranks;
    vector<BoolVar> nodes;
    for (int i = 0; i < N; ++i) {
        ranks.push_back(model.NewIntVar({0,MAX_RANK_VALUE}));
        nodes.push_back(model.NewBoolVar());
    }

    auto arcs = GraphUtils::getGraphEdges(V,true);
    for (auto [a,b] : arcs) {
        auto cond_var = model.NewBoolVar();
        model.AddLessThan(ranks[a], ranks[b]).OnlyEnforceIf(cond_var);
        model.AddBoolOr( {nodes[a], nodes[b], cond_var} );
    }

    { // here we add all pi-edges or triangles to make the propagation faster
        int L0 = 3;
        auto cyc = Utils::getAllSimpleCycles3(V,L0);
        clog << "\t adding " << cyc.size() << " constraints for all simple cycles of length <= " << L0 << endl;
        for(auto & v : cyc) {
            if( v.size() == 2 ) model.AddBoolOr({nodes[v[0]], nodes[v[1]]});
            if( v.size() == 3 ) model.AddBoolOr({nodes[v[0]], nodes[v[1]], nodes[v[2]]});
        }
    }

    if(!init_sol.empty()) {
        // for( int i=0; i<N; i++ ) model.AddHint(nodes[i], in_init_sol[i]);
        for( int d : init_sol ) model.AddHint(nodes[d], 1);

        // now toposort to create initial ranks values
        VI topo, deg(N,0);
        for( int i=0; i<N; i++ ) if(!in_init_sol[i]) for(int d : V[i]) if(!in_init_sol[d]){ deg[d]++; }
        for( int i=0; i<N; i++ ) if( !in_init_sol[i] && deg[i] == 0 ) topo.push_back(i);
        for(int i=0; i<topo.size(); i++) {
            int v = topo[i];
            for(int d : V[v]) if(!in_init_sol[d]) {
                deg[d]--;
                if(deg[d] == 0) topo.push_back(d);
            }
        }
        for(int i=0; i<topo.size(); i++) model.AddHint( ranks[topo[i]], i );
    }

    SatParameters params;
    params.set_max_time_in_seconds(time_limit_millis / 1000.0);
    params.set_num_search_workers(threads);   // deterministic runs
    if(use_only_lns) params.set_use_lns_only(true);
    params.set_log_search_progress(false);
    Model solver_model;
    solver_model.Add(NewSatParameters(params));

    model.Minimize(LinearExpr::Sum(nodes));
    CpSolverResponse response = SolveCpModel(model.Build(), &solver_model);

    string status;
    if ( response.status() == OPTIMAL ) status = "OPTIMAL";
    if ( response.status() == UNKNOWN ) status = "UNKNOWN";
    if ( response.status() == FEASIBLE ) status = "FEASIBLE";
    if ( response.status() == INFEASIBLE ) {
        status = "INFEASIBLE";
        assert(false && "status cannot be infeasible, unless model is incorrect");
    }

    VI res;
    if ( response.status() == OPTIMAL || response.status() == FEASIBLE ) {
        for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
    }



    return {status,res};
}

/**
 * Iterative HS.
 * Starting with L = 2 and res = {}, finds all cycles in the graph G[V \ res] of length <= L.
 * If there is just a small number of such cycles, increases L and repeats.
 * Then finds HS of the set of all cycles found so far.
 */
pair<string,VI> solveCPSAT3(VVI V, int total_time_seconds, int max_time_seconds_per_iter, bool find_optimal) {
    clog << "Solving using CPSAT, model v3" << endl;

    int N = V.size();

    VVI cycles; // = Utils::getAllSimpleCycles3(V,2);
    VI prev_res;

    VI best_fvs;

    auto updateCyclesAndHit = [&](VVI new_cycles) -> VI {
        cycles += std::move(new_cycles);

        SatParameters params;
        if (!find_optimal) {
            params.set_max_time_in_seconds(max_time_seconds_per_iter);
            clog << "\t running cpsat solver with time limit of " << max_time_seconds_per_iter << " sec." << endl;
        }
        params.set_num_search_workers(threads);   // deterministic runs
        // if( !find_optimal && !prev_res.empty() ) params.set_use_lns_only(true);
        params.set_log_search_progress(false);
        Model solver_model;
        solver_model.Add(NewSatParameters(params));

        CpModelBuilder model;
        vector<BoolVar> nodes;
        for (int i=0; i<N; i++) nodes.push_back(model.NewBoolVar());
        for ( auto & c : cycles ) {
            vector<BoolVar> cyc_vars;
            for ( int d : c ) cyc_vars.push_back(nodes[d]);
            model.AddAtLeastOne(cyc_vars);
        }

        // VB in_res = StandardUtils::toVB(N,prev_res);
        // for(int i=0; i<N; i++) model.AddHint(nodes[i],in_res[i]);
        for(int d : prev_res) model.AddHint(nodes[d], 1);

        if(find_optimal && !prev_res.empty()) model.AddGreaterOrEqual(LinearExpr::Sum(nodes), (int)prev_res.size());

        model.Minimize(LinearExpr::Sum(nodes));
        CpSolverResponse response = SolveCpModel(model.Build(), &solver_model);
        int F = 2;
        while(response.status() == UNKNOWN) {
            clog << "Response status unknown, increasing max_time_per_iter to " << F*max_time_seconds_per_iter << endl;

            params.set_max_time_in_seconds(F*max_time_seconds_per_iter);
            Model solver_model;
            solver_model.Add(NewSatParameters(params));
            F *= 2;
            response = SolveCpModel(model.Build(), &solver_model);
        }


        VI res;
        for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
        return res;
    };


    // VI res;
    Stopwatch timer;
    timer.setLimit("all_iters", find_optimal ? inf : total_time_seconds * 1000);
    timer.start("all_iters");


    int L = 2;
    while(true) {
        if(timer.tle("all_iters")) break;

        VVI H = V;
        VVI revH = GraphUtils::reverseGraph(H);
        VB helper(N);
        Utils::removeNodes(H, revH, prev_res, helper);
        VVI nonpiH = Utils::getNonPIGraph(H);
        VVI revnonpiH = GraphUtils::reverseGraph(nonpiH);

        clog << endl << "Looking for new cycles, cycles.size(): " << cycles.size() << ", prev_res.size(): "
             << prev_res.size() << ", time: " << (int)timer.getTime("all_iters") / 1000 << endl;

        int all_arcs = GraphUtils::countEdges(V,true);
        set<PII> zb;
        for( auto & cyc : cycles ) for( int j=cyc.size()-1, i=0; i < cyc.size(); j = i++ ) zb.insert( {cyc[j], cyc[i]} );
        int arcs = zb.size();
        clog << "\t arcs in constraints: " << arcs << " / " << all_arcs << endl;

        VVI new_cycles = Utils::getAllSimpleCycles3(H, L);

        // we cannot add more than MAX_NEW_CYCLES_PER_ITERATION_PERC * cycles.size() new cycles in each iteration
        // this is here, because when increasing the length size, we might get an awful lot of new cycles of
        // that length, we do not want that, we want to keep number of cycles used for constraints as small as possible
        constexpr double MAX_NEW_CYCLES_PER_ITERATION_PERC = 0.1;
        if( L >= 5 && cycles.size() >= 100 && new_cycles.size() > MAX_NEW_CYCLES_PER_ITERATION_PERC * cycles.size() ) {
            StandardUtils::shuffle(new_cycles);
            new_cycles.resize( MAX_NEW_CYCLES_PER_ITERATION_PERC * cycles.size() );
        }

        if( new_cycles.empty() && !Utils::isFVS(V,prev_res) ) {
            L++;
            clog << endl << "---> INCREASING LENGTH, L: " << L << endl << endl;
            continue;
        }

        clog << "\t there are " << new_cycles.size() << " new cycles found, all cycles: " << cycles.size() << endl;


        VI sol = updateCyclesAndHit(new_cycles);
        if ( Utils::isFVS(V,sol) ) {
            if( best_fvs.empty() || sol.size() < best_fvs.size() ) best_fvs = sol;

            if(find_optimal) {
                // res = sol;
                assert( best_fvs.size() == sol.size() );
                break;
            }else {
                max_time_seconds_per_iter++;
                clog << endl << "--> Found a valid FVS, increasing max_time_seconds_per_iter to "
                     << max_time_seconds_per_iter << " secc." << endl << endl;
            }
        }

        prev_res = sol;

        clog << "\t found sol.size(): " << sol.size() << ", best_fvs.size(): " << best_fvs.size() << endl;
    }

    timer.stop("all_iters");
    timer.write("all_iters");

    string status;
    if( Utils::isFVS(V,best_fvs) ) {
        if( find_optimal ) status = "OPTIMAL";
        else status = "FEASIBLE";
    }
    else status = "INCORRECT";

    return {status, best_fvs};
}


VVI readDirectedExample() {
    int N,M,c;
    cin >> N >> M >> c;
    cin.ignore();
    // cin.ignore();

    VVI V(N);
    for (int i=0; i<N; i++) {
        string s;
        getline(cin,s);
        auto l = StandardUtils::split(s, " ");
        if (!l.empty()) {
            for ( auto & v : l ) if (!v.empty() && v != " ") {
                V[i].push_back(stoi(v)-1);
                assert(V[i].back() >= 0);
                assert(V[i].back() < N);
            }
        }
    }

    return V;
}


int main(int argc, char** argv){
    MemoryUtils::increaseStack();

    initializeParams(argc, argv);

    auto old_clog_buf = clog.rdbuf();
    if(MUTE_MODE){
        clog << "MUTE MODE" << endl;
        clog.rdbuf( nullptr );
    }

    VVI V = readDirectedExample();
    int N = V.size();


    DEBUG(V.size());
    DEBUG(GraphUtils::countEdges(V,true));


    constexpr bool check_heuristic_algorithms = true;
    constexpr bool check_exact_algorithms = true;

    if(check_heuristic_algorithms) {
        Stopwatch sw;

        //******************************

        sw.start("cpsat-3");
        // auto[status0,res0] = solveCPSAT3(V, inf, inf, true); // exact solution
        auto[status0,res0] = solveCPSAT3(V, 30, 1, false); // heuristic approach
        sw.stop("cpsat-3");

        DEBUG(status0);
        DEBUG(res0.size());

        //******************************


        sw.start("cpsat-2-N");
        auto[status2,res2] = solveCPSAT2(V, N);
        sw.stop("cpsat-2-N");

        DEBUG(status2);
        DEBUG(res2.size());
        // DEBUG(res2);


        //******************************

        sw.start("cpsat-2-inf");
        auto[status3,res3] = solveCPSAT2(V, inf);
        sw.stop("cpsat-2-inf");

        DEBUG(status3);
        DEBUG(res3.size());
        // DEBUG(res3);

        //******************************

        sw.start("cpsat-2-N/10");
        auto[status4,res4] = solveCPSAT2(V, V.size()/10);
        sw.stop("cpsat-2-N/10");

        DEBUG(status4);
        DEBUG(res4.size());
        // DEBUG(res4);

        //******************************

        sw.start("cpsat-2-only-lns-from-res4");
        auto[status5,res5] = solveCPSAT2(V, N, res4, true);
        sw.stop("cpsat-2-only-lns-from-res4");

        DEBUG(status5);
        DEBUG(res5.size());
        // DEBUG(res4);

        //******************************

        sw.start("cpsat-1");
        auto[status1,res1] = solveCPSAT1(V,25,time_limit_millis);
        sw.stop("cpsat-1");

        DEBUG(status1);
        DEBUG(res1.size());
        // DEBUG(res1);


        if ( status1 == "OPTIMAL" && status2 == "OPTIMAL" ) assert( res1.size() == res2.size() );

        ENDL(3);
        DEBUG(res1.size());
        DEBUG(res2.size());
        DEBUG(res3.size());
        DEBUG(res4.size());
        DEBUG(res5.size());

        if( !res1.empty() && status1 != "INCORRECT" ) assert(Utils::isFVS(V,res1));
        if( !res2.empty() && status2 != "INCORRECT" ) assert(Utils::isFVS(V,res2));
        if( !res3.empty() && status3 != "INCORRECT" ) assert(Utils::isFVS(V,res3));
        if( !res4.empty() && status4 != "INCORRECT" ) assert(Utils::isFVS(V,res4));
        if( !res5.empty() && status5 != "INCORRECT" ) assert(Utils::isFVS(V,res5));
        if( !res0.empty() && status0 != "INCORRECT" ) assert(Utils::isFVS(V,res0));

        sw.write("cpsat-1");
        sw.write("cpsat-2-N");
        sw.write("cpsat-2-inf");
        sw.write("cpsat-2-N/10");
        sw.write("cpsat-2-only-lns-from-res4");
        sw.write("cpsat-3");
    }



    if(check_exact_algorithms) {
        time_limit_millis = inf;

        Stopwatch sw;

        //******************************

        sw.start("cpsat-3");
        auto[status0,res0] = solveCPSAT3(V, inf, inf, true);
        sw.stop("cpsat-3");

        DEBUG(status0);
        DEBUG(res0.size());

        //******************************

        // we set time_limit_millis = inf, so this will find optimal result
        sw.start("cpsat-2-N");
        auto[status2,res2] = solveCPSAT2(V, N);
        sw.stop("cpsat-2-N");

        DEBUG(status2);
        DEBUG(res2.size());

        //******************************

        sw.start("cpsat-1");
        auto[status1,res1] = solveCPSAT1(V,25,inf, false);
        sw.stop("cpsat-1");

        DEBUG(status1);
        DEBUG(res1.size());

        ENDL(3);
        DEBUG(res0.size());
        DEBUG(res1.size());
        DEBUG(res2.size());

        assert(Utils::isFVS(V,res0));
        assert(Utils::isFVS(V,res1));
        assert(Utils::isFVS(V,res2));

        sw.write("cpsat-1");
        sw.write("cpsat-2-N");
        sw.write("cpsat-3");

        if( !res1.empty() && status1 != "INCORRECT" ) assert(Utils::isFVS(V,res1));
        if( !res2.empty() && status2 != "INCORRECT") assert(Utils::isFVS(V,res2));
        if( !res0.empty() && status0 != "INCORRECT") assert(Utils::isFVS(V,res0));
    }



    return 0;
}

