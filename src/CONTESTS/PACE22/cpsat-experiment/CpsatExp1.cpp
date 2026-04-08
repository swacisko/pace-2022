//
// Created by sylwe on 08/04/2026.
//

#include "CpsatExp1.h"

#include "GraphUtils.h"
#include "StandardUtils.h"
#include "Stopwatch.h"
#include "CONTESTS/PACE22/Utils.h"
#include "ortools/sat/cp_model.h"
using namespace operations_research::sat;


inline ExpData CpsatExp1::solveHS1(VVI V, ExpConfig cnf) {
     clog << "Solving using CPSAT, model v1" << endl;

    int N = V.size();

    VI prev_res;

    auto solveHS = [&](auto & cycles) -> VI {
        SatParameters params;
        if (!cnf.find_optimal_result) params.set_max_time_in_seconds(time_limit_millis / 1000.0);
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

        if ( cnf.find_optimal_result ) model.AddGreaterOrEqual(LinearExpr::Sum(nodes), (int)prev_res.size());

        model.Minimize(LinearExpr::Sum(nodes));
        CpSolverResponse response = SolveCpModel(model.Build(), &solver_model);

        VI res;
        for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
        return res;
    };

    Stopwatch timer;
    timer.setLimit("all_iters", max_l * time_limit_millis);
    timer.start("all_iters");

    if(cnf.find_optimal_result) max_l = inf;

    VI res;
    for (int L=3; L <= max_l ; L++) {
        if ( !cnf.find_optimal_result && timer.tle("all_iters")) {
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
        if( cnf.find_optimal_result ) status = "OPTIMAL";
        else status = "FEASIBLE";
    }
    else status = "INCORRECT";

    return {status, res};
}

inline ExpData CpsatExp1::solveIHS(VVI V, ExpConfig cnf, int cycle_enumeration_type) {
    clog << "Solving using CPSAT, model v3" << endl;

    int N = V.size();

    VVI cycles; // = Utils::getAllSimpleCycles3(V,2);
    VI prev_res;

    VI best_fvs;

    auto updateCyclesAndHit = [&](VVI new_cycles) -> VI {
        cycles += std::move(new_cycles);

        SatParameters params;
        if (!cnf.find_optimal_result) {
            params.set_max_time_in_seconds(max_time_seconds_per_iter);
            clog << "\t running cpsat solver with time limit of " << max_time_seconds_per_iter << " sec." << endl;
        }
        params.set_num_search_workers(cnf.threads);   // deterministic runs
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

        if(cnf.find_optimal_result && !prev_res.empty()) model.AddGreaterOrEqual(LinearExpr::Sum(nodes), (int)prev_res.size());

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
    timer.setLimit("all_iters", cnf.find_optimal_result ? inf : total_time_seconds * 1000);
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

            if(cnf.find_optimal_result) {
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
        if( cnf.find_optimal_result ) status = "OPTIMAL";
        else status = "FEASIBLE";
    }
    else status = "INCORRECT";

    return {status, best_fvs};
}

inline ExpData CpsatExp1::solveMTZ(VVI V, ExpConfig cnf, int auxiliary_cycles_mode) {
    clog << "Solving using CPSAT, model v2" << endl;

    int N = V.size();
    int MAX_RANK_VALUE = N;

    VB in_init_sol = StandardUtils::toVB(N,init_sol);

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
    params.set_num_search_workers(cnf.threads);   // deterministic runs
    if(cnf.use_only_cpsat_lns) params.set_use_lns_only(true);
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

inline ExpData CpsatExp1::solveDiVerSeS(VVI V, ExpConfig cnf) {
}