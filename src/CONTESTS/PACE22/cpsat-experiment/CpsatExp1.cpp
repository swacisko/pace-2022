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

SatParameters getDefaultSatParameters(ExpConfig cnf) {
    SatParameters params;
    if (!cnf.find_optimal_result) params.set_max_time_in_seconds(cnf.max_time_sec);
    if(cnf.use_only_cpsat_lns) params.set_use_lns_only(true);
    params.set_num_search_workers(cnf.threads);
    params.set_log_search_progress(cnf.log_cpsat_search_progress);
    return params;
}

VVI CpsatExp1::getUnhitChordlessCycles(VVI &V, VI &S, int max_l, int enumeration_option) {
    if (enumeration_option == 1) {
        VVI H = V;
        VVI revH = GraphUtils::reverseGraph(H);
        VB helper(V.size());
        Utils::removeNodes(H, revH, S, helper);
        return Utils::getAllSimpleCycles3(H,max_l);
    }else {
        assert(false && "not implemented yet");
    }
}

ExpData CpsatExp1::solveHS1(VVI V, ExpConfig cnf) {
     clog << "Solving using CPSAT, model v1" << endl;

    ExpData exp_data;
    int N = V.size();

    VI prev_res;

    auto solveHS = [&](auto & cycles) -> VI {
        SatParameters params = getDefaultSatParameters(cnf);
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

        VB in_res = StandardUtils::toVB(N,prev_res);
        for(int i=0; i<N; i++) model.AddHint(nodes[i],in_res[i]);
        // for(int d : prev_res) model.AddHint(nodes[d],1);

        if ( cnf.find_optimal_result ) model.AddGreaterOrEqual(LinearExpr::Sum(nodes), (int)prev_res.size());

        model.Minimize(LinearExpr::Sum(nodes));
        CpSolverResponse response = SolveCpModel(model.Build(), &solver_model);

        VI res;
        for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
        return res;
    };

    constexpr int max_l = inf;

    Stopwatch timer;
    timer.setLimit("all_iters", max_l * cnf.max_time_sec*1'000 );
    timer.start("all_iters");

    VI res;
    for (int L=3; L <= max_l ; L++) {
        if ( !cnf.find_optimal_result && timer.tle("all_iters")) {
            clog << "Solver did not find optimal value for given set of cycles in admissible time" << endl;
            break;
        }


        clog << endl << "Considering cycles of length <= " << L << endl;
        auto cycles = Utils::getAllSimpleCycles3(V,L );
        clog << "\t there are " << cycles.size() << " such cycles" << endl;
        if (cycles.size() > cnf.max_cycles_for_hs) break;

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

    // return {status, res};
    return exp_data;
}

ExpData CpsatExp1::solveIHS(VVI V, ExpConfig cnf) {
    VVI cycles;
    VI res;
    return solveIHS(V,cnf,cycles,res);
}

ExpData CpsatExp1::solveIHS(VVI V, ExpConfig cnf, VVI & cycles, VI & res) {
    clog << "Solving using CPSAT, model v3" << endl;

    int cycle_enumeration_type = cnf.unhit_cycle_enumeration_type;
    ExpData exp_data;

    int N = V.size();

    cycles.clear();
    VI prev_res;

    VI best_fvs;

    auto updateCyclesAndHit = [&](VVI new_cycles) -> VI {
        cycles += std::move(new_cycles);

        SatParameters params = getDefaultSatParameters(cnf);
        if (!cnf.find_optimal_result) clog << "\t running cpsat solver with time limit of " << cnf.ihs_single_iteration_sec << " sec." << endl;
        else clog << "\t running cpsat solver without time limit, looking for optimal result" << endl;

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

        VB in_res = StandardUtils::toVB(N,prev_res);
        for(int i=0; i<N; i++) model.AddHint(nodes[i],in_res[i]);
        // for(int d : prev_res) model.AddHint(nodes[d], 1);

        if(cnf.find_optimal_result && !prev_res.empty()) model.AddGreaterOrEqual(LinearExpr::Sum(nodes), (int)prev_res.size());

        model.Minimize(LinearExpr::Sum(nodes));
        CpSolverResponse response = SolveCpModel(model.Build(), &solver_model);
        int F = 2;
        while(response.status() == UNKNOWN) {
            clog << "Response status unknown, increasing max_time_per_iter to " << F * cnf.ihs_single_iteration_sec << endl;

            params.set_max_time_in_seconds(F * cnf.ihs_single_iteration_sec);
            Model solv_model;
            solv_model.Add(NewSatParameters(params));
            F *= 2;
            response = SolveCpModel(model.Build(), &solv_model);
        }


        VI res;
        for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
        return res;
    };


    Stopwatch timer;
    timer.setLimit("all_iters", cnf.find_optimal_result ? inf : cnf.max_time_sec * 1000);
    timer.start("all_iters");


    int L = 2;
    while(true) {
        if(timer.tle("all_iters")) break;

        // VVI H = V;
        // VVI revH = GraphUtils::reverseGraph(H);
        // VB helper(N);
        // Utils::removeNodes(H, revH, prev_res, helper);
        // VVI nonpiH = Utils::getNonPIGraph(H);
        // VVI revnonpiH = GraphUtils::reverseGraph(nonpiH);

        clog << endl << "Looking for new cycles, cycles.size(): " << cycles.size() << ", prev_res.size(): "
             << prev_res.size() << ", time: " << (int)timer.getTime("all_iters") / 1000 << endl;

        // VVI new_cycles = Utils::getAllSimpleCycles3(H, L);
        VVI new_cycles = getUnhitChordlessCycles(V,prev_res,L,cycle_enumeration_type);

        int all_arcs = GraphUtils::countEdges(V,true);
        set<PII> zb;
        for( auto & cyc : cycles ) for( int j=cyc.size()-1, i=0; i < cyc.size(); j = i++ ) zb.insert( {cyc[j], cyc[i]} );
        int arcs = zb.size();
        clog << "\t arcs in constraints: " << arcs << " / " << all_arcs << endl;


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

        int max_time_seconds_per_iter = cnf.ihs_single_iteration_sec;
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

    // return {status, best_fvs};
    if (!best_fvs.empty()) res = best_fvs;
    else res = prev_res;

    return exp_data;
}

ExpData CpsatExp1::solveMTZ(VVI V, ExpConfig cnf, int auxiliary_cycles_mode) {
    clog << "Solving using CPSAT, model v2" << endl;

    ExpData exp_data;

    int N = V.size();
    int MAX_RANK_VALUE = min(1ll*inf,1ll*(N+2)*(int)sqrt(N));

    VI init_sol = {};
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

    if (auxiliary_cycles_mode == 1){ // here we add all pi-edges or triangles to make the propagation faster
        int L0 = 3;
        auto cyc = Utils::getAllSimpleCycles3(V,L0);
        clog << "\t adding " << cyc.size() << " constraints for all simple cycles of length <= " << L0 << endl;
        for(auto & v : cyc) {
            if( v.size() == 2 ) model.AddBoolOr({nodes[v[0]], nodes[v[1]]});
            if( v.size() == 3 ) model.AddBoolOr({nodes[v[0]], nodes[v[1]], nodes[v[2]]});
        }
    }else if (auxiliary_cycles_mode == 2) {
        auto new_cnf = cnf;
        new_cnf.max_time_sec *= new_cnf.max_time_fraction_for_ihs_cycles_in_mtz;
        if ( new_cnf.ihs_single_iteration_sec > 2 ) {
            new_cnf.ihs_single_iteration_sec *= new_cnf.max_time_fraction_for_ihs_cycles_in_mtz;
            new_cnf.ihs_single_iteration_sec = max(new_cnf.ihs_single_iteration_sec,2);
        }
        VVI cycles;
        VI res;
        auto r = solveIHS(V,new_cnf,cycles, res);
        for(auto & v : cycles) {
            vector<BoolVar> cyc_vars;
            for (int d : v) cyc_vars.push_back(nodes[d]);
            model.AddBoolOr(cyc_vars);
        }
        init_sol = res;
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
        for(int i=0; i<topo.size(); i++) model.AddHint( ranks[topo[i]], (i+1)*sqrt(N) );
    }

    SatParameters params = getDefaultSatParameters(cnf);
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



    // return {status,res};
    return exp_data;
}

ExpData CpsatExp1::solveDiVerSeS(VVI V, ExpConfig cnf) {
    ExpData exp_data;
    assert(false && "not implemented");

    return exp_data;
}

ExpData CpsatExp1::solve(VVI V, ExpConfig cnf) {
    auto alg = cnf.alg;
    if (alg == Algorithm::HS) return solveHS1(V,cnf);
    if (alg == Algorithm::IHS) return solveIHS(V,cnf);
    if (alg == Algorithm::MTZ) return solveMTZ(V,cnf);
    if (alg == Algorithm::DIVERSES) return solveDiVerSeS(V,cnf);

    return ExpData{};
}
