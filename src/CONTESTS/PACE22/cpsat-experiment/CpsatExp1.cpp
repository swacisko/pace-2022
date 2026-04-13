//
// Created by sylwe on 08/04/2026.
//

#include "CpsatExp1.h"

#include "GraphUtils.h"
#include "StandardUtils.h"
#include "Stopwatch.h"
#include "CONTESTS/PACE22/Utils.h"
#include "ortools/sat/cp_model.h"
#include "scc/StronglyConnectedComponents.h"
using namespace operations_research::sat;

SatParameters getDefaultSatParameters(ExpConfig cnf, int time_in_sec) {
    SatParameters params;
    if (!cnf.find_optimal_result) params.set_max_time_in_seconds(time_in_sec);
    if(cnf.use_only_cpsat_lns) params.set_use_lns_only(true);
    params.set_num_search_workers(cnf.threads);
    params.set_log_search_progress(cnf.log_cpsat_search_progress);
    return params;
}


VVI CpsatExp1::getUnhitChordlessCycles(VVI &V, VI &S, int max_l, int max_millis, int enumeration_option) {
    if (enumeration_option == 1) {
        VVI H = V;
        VVI revH = GraphUtils::reverseGraph(H);
        VB helper(V.size());
        Utils::removeNodes(H, revH, S, helper);
        clog << "\t Looking for cycles for at most " << max_millis << " millis" << endl;
        return Utils::getAllSimpleCycles3(H,max_l, max_millis);
    }else {
        auto indg = getUnhitGraph(V,S);
        VVI H = indg.V;
        int N = H.size();
        vector<LL> hashes(N);
        IntGenerator rnd;
        for (int i=0; i<N; i++) hashes[i] = rnd.rand();
        unordered_set<LL> cycle_hashes;

        Stopwatch sw;
        sw.setLimit("cycles",max_millis);
        sw.start("cycles");

        auto getCycles = [&](int v) -> VVI {
            VB was(N);
            VI par(N,-1);
            VB on_path(N);
            for (auto & vec : H) StandardUtils::shuffle(vec,rnd);
            VVI cycles;

            function<void(int)> dfs = [&](int num) {
                was[num] = true;
                on_path[num] = true;
                for ( int d : H[num] ) {
                    if ( was[d] && on_path[d] ) { // create a cycle
                        VI cyc(1,num);
                        int p = num;
                        while (p != d) cyc.push_back(p = par[p]);
                        reverse(ALL(cyc));
                        cycles.push_back(cyc);
                    }else if ( !was[d] ) {
                        par[d] = num;
                        dfs(d);
                    }
                }
                on_path[num] = false;
            };
            dfs(v);

            return cycles;
        };

        VI ind_on_cycle(N,-1);
        VB temp(N);

        auto makeChordless = [&](VI cyc)-> VI {
            for (int i=0; i<cyc.size(); i++) ind_on_cycle[cyc[i]] = i;
            int v = cyc[0];
            VI res(1,v);
            temp[v] = true;

            while (true) {
                int u = v;
                for ( int d : H[u] ) if (ind_on_cycle[d] != -1 && ind_on_cycle[d] > ind_on_cycle[v] ) v = d;
                res.push_back(v);
                temp[v] = true;
                int max_ind = -1;
                u = -1;
                for ( int d : H[v] ) if ( ind_on_cycle[d] != -1 && temp[d] && ind_on_cycle[d] < ind_on_cycle[v] && ind_on_cycle[d] > max_ind ) {
                    max_ind = ind_on_cycle[d];
                    u = d;
                }
                if (u != -1) {
                    reverse(ALL(res));
                    while (res.back() != u) res.pop_back();
                    reverse(ALL(res));
                    break;
                }
            }

            for (int d : cyc) temp[d] = false;
            for (int d : cyc) ind_on_cycle[d] = -1;

            return res;
        };

        auto getHash = [&](VI cyc)-> LL {
            LL h = 0;
            for (int d : cyc) h ^= hashes[d];
            return h;
        };

        auto isChordless = [&](VI cyc)-> bool {
            bool chordless = true;
            for (int i=0; i<cyc.size()-1; i++) ind_on_cycle[cyc[i]] = i;
            for ( int v : cyc ) for (int d : H[v]) if (ind_on_cycle[d] != -1) {
                if ( ind_on_cycle[d] > ind_on_cycle[v]+1 ) chordless = false;
                if ( ind_on_cycle[d] < ind_on_cycle[v]-1 && (v != cyc.back() || ind_on_cycle[d] != 0) ) chordless = false;
                if ( ind_on_cycle[d] == ind_on_cycle[v] - 1 && cyc.size() != 2 ) chordless = false;
            }
            for (int i=0; i<cyc.size()-1; i++) ind_on_cycle[cyc[i]] = -1;
            return chordless;
        };

        VVI cycles;
        VI perm(N); iota(ALL(perm),0);
        StandardUtils::shuffle(perm,rnd);

        for (int v : perm) {
            if ( sw.tle("cycles") ) break;
            clog << "Getting cycles for node v: " << v << endl;
            auto new_cycles = getCycles(v);
            int cycles_added = 0;
            int new_chordless_cycles = 0;
            double suml1 = 0, suml2 = 0;
            for ( auto & cyc : new_cycles ) {
                suml1 += cyc.size();
                int l = cyc.size();
                cyc = makeChordless(cyc);
                new_chordless_cycles += ( l == cyc.size() );
                assert(isChordless(cyc));
                auto h = getHash(cyc);
                if ( !cycle_hashes.contains(h) ) {
                    cycle_hashes.insert(h);
                    cycles.push_back(cyc);
                    cycles_added++;
                    suml2 += cyc.size();
                }
            }
            if (!new_cycles.empty()) {
                clog << "\t found " << new_cycles.size() << " new cycles of avg_length: " << suml1 / new_cycles.size()
                     << " from which " << new_chordless_cycles << " were already chordless --> only added "
                     << cycles_added << " new chordless cycles, with avg_length: " << suml2 / cycles_added << endl;
            }
        }

        for (auto & cyc : cycles) indg.remapNodes(cyc);
        return cycles;
    }
}

InducedGraph CpsatExp1::getUnhitGraph(VVI &V, VI &S) {
    int N = V.size();
    VVI H = V;
    VVI revH = GraphUtils::reverseGraph(H);
    VB helper(N);
    Utils::removeNodes(H, revH, S, helper);

    StronglyConnectedComponents scc(V);
    scc.createStronglyConnectedComponents();
    auto comps = scc.getComponents();
    VI in_comp = StandardUtils::layersToPartition(comps);
    VPII arcs_to_remove;
    VPII all_arcs = GraphUtils::getGraphEdges(V, true);
    for( auto & [a,b] : all_arcs ) if( in_comp[a] != in_comp[b] ) arcs_to_remove.emplace_back(a,b);

    fill(ALL(helper),false);
    Utils::removeEdges( H, arcs_to_remove, helper );

    return GraphInducer::induceByNonisolatedNodes(H);
}

PII CpsatExp1::getUnhitGraphSizes(VVI &V, VI &S) {
    auto H = getUnhitGraph(V,S);
    return {H.V.size(), GraphUtils::countEdges(H.V,true)};
}

void CpsatExp1::addCycleConstraints(CpModelBuilder &model, VVI & cycles, vector<BoolVar> & nodes) {
    for ( auto & c : cycles ) {
        vector<BoolVar> cyc_vars;
        for ( int d : c ) cyc_vars.push_back(nodes[d]);
        model.AddAtLeastOne(cyc_vars);
    }
}

void CpsatExp1::addInitialSolutionHint(CpModelBuilder &model, vector<BoolVar> &nodes, VI &init_sol, VI &prev_res, ExpConfig cnf) {
    int N = nodes.size();
    VB in_init_sol = StandardUtils::toVB(N,init_sol);
    if (cnf.use_init_sol_as_hint_mode == 0) for(int i=0; i<N; i++) model.AddHint(nodes[i],in_init_sol[i]);
    if (cnf.use_init_sol_as_hint_mode == 1) for(int i=0; i<N; i++) if (in_init_sol[i]) model.AddHint(nodes[i],1);
    if (cnf.use_init_sol_as_hint_mode == 2) for(int d : prev_res) model.AddHint(nodes[d],1);
}

void CpsatExp1::addMaxHammingDstConstraint(CpModelBuilder &model, vector<BoolVar> &nodes, VI &init_sol, ExpConfig cnf) {
    int N = nodes.size();
    VB in_init_sol = StandardUtils::toVB(N,init_sol);
    if (cnf.find_optimal_result) {
        clog << endl << "#CAUTION!!! Using next_sol_max_dst_from_init_sol, but trying to find optimal solution!"
             << " Optimality might not be preserved!" << endl;
    }
    LinearExpr hamming_dist;
    for (int i = 0; i < N; ++i) {
        hamming_dist += in_init_sol[i] ? (1 - nodes[i]) : nodes[i];
    }
    model.AddLessOrEqual(hamming_dist, cnf.next_sol_max_dst_from_init_sol);
}

VI CpsatExp1::rerunModelUntilFeasibleOrTle(CpModelProto &model_proto, vector<BoolVar> &nodes,
    CpSolverResponse & response, Stopwatch & timer, string timer_option, int init_time, ExpConfig cnf) {

    while(response.status() == CpSolverStatus::UNKNOWN && !timer.tle(timer_option)) {
        init_time = ceil(min( 3000.0 * init_time, timer.getLimit(timer_option) - timer.getTime(timer_option) ) / 1000);
        clog << "Response status unknown, increasing max_time_per_iter to " << init_time << endl;

        SatParameters params = getDefaultSatParameters(cnf, init_time);
        Model solv_model;
        solv_model.Add(NewSatParameters(params));
        response = SolveCpModel(model_proto, &solv_model);
    }

    int N = nodes.size();
    if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {
        VI res;
        for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
        return res;
    }

    return {};
}

VI CpsatExp1::solveCpsatForCycles(VVI &V, VVI &cycles, VI &prev_res, VI &init_sol, Stopwatch &timer,
    string timer_option, ExpConfig cnf) {

    int N = V.size();
    int F = ( Utils::isFVS(V,prev_res) ? 5 : 1 );
    int time = ceil(min( 1000.0 * F * cnf.ihs_single_iteration_sec, timer.getLimit(timer_option) - timer.getTime(timer_option) ) / 1000);
    time = max(time,1);
    SatParameters params = getDefaultSatParameters(cnf, time);
    if (cnf.find_optimal_result) clog << "\t running cpsat solver without time limit, looking for optimal result" << endl;
    else clog << "\t running cpsat solver with time limit of " << time << " sec." << endl;

    Model solver_model;
    solver_model.Add(NewSatParameters(params));

    CpModelBuilder model;
    vector<BoolVar> nodes;
    for (int i=0; i<N; i++) nodes.push_back(model.NewBoolVar());
    addCycleConstraints(model,cycles,nodes);

    addInitialSolutionHint(model,nodes,init_sol,prev_res, cnf);
    assert( isHS(cycles, init_sol) );

    if (cnf.next_sol_max_dst_from_init_sol != inf && !prev_res.empty()) addMaxHammingDstConstraint(model,nodes,init_sol,cnf);

    model.Minimize(LinearExpr::Sum(nodes));
    auto model_proto = model.Build();
    CpSolverResponse response = SolveCpModel(model_proto, &solver_model);

    return rerunModelUntilFeasibleOrTle( model_proto,nodes, response,timer, timer_option,time,cnf );
}

VI CpsatExp1::getUnhitCyclesHSGreedy(VVI &cycles) {
    // clog << "Looking for greedy HS" << endl;
    int N = 0, M = cycles.size();
    for (auto & c : cycles) for (int d : c) N = max(N,d);
    N++;
    VVI A(N), B(M);
    for ( int i=0; i<M; i++ ) for ( int d : cycles[i] ) {
        A[d].push_back(i);
        B[i].push_back(d);
    }

    VI deg(N,0);
    for (int i=0; i<N; i++) deg[i] = A[i].size();

    VI res;
    VB hit(M,false);

    priority_queue<PII> zb;
    for (int i=0; i<N; i++) zb.emplace( deg[i],i );
    while (!zb.empty()) {
        auto [d,v] = zb.top();
        zb.pop();
        if ( d != deg[v] ) continue;
        assert(deg[v] >= 0);
        if ( deg[v] <= 0 ) continue;

        res.push_back(v);
        for (int d : A[v]) if (!hit[d]) {
            hit[d] = true;
            for (int dd : cycles[d]) {
                deg[dd]--;
                zb.emplace(deg[dd],dd);
            }
        }
    }
    // clog << "\tfound greedy HS" << endl;
    return res;
}

bool CpsatExp1::isHS(VVI &cycles, VI &S) {
    if (cycles.empty()) return true;

    int N = 0;
    for (auto & c : cycles) for (int d : c) N = max(N,d);
    for (int d : S) N = max(N,d);
    N++;
    VB was(N);

    for ( int d : S ) was[d] = true;
    for ( auto & c : cycles ) {
        bool hit = false;
        for ( int d : c ) hit |= was[d];
        if (!hit) return false;
    }
    return true;
}

ExpData CpsatExp1::solveHS(VVI V, ExpConfig cnf) {
     clog << "Solving using CpsatExp1::solveHS" << endl;
    cnf.writeConfig();

    ExpData exp_data;
    int N = V.size();

    VI prev_res;

    Stopwatch timer;
    string timer_option = "solveHS";
    timer.setLimit(timer_option, 1000 * cnf.max_time_sec );
    timer.start(timer_option);

    // auto solveHS = [&](auto & cycles, VI unhit_cycles_hs) -> VI {
    //     int F = ( Utils::isFVS(V,prev_res) ? 5 : 1 );
    //     int time = ceil(min( 1000.0 * F * cnf.ihs_single_iteration_sec, timer.getLimit(timer_option) - timer.getTime(timer_option) ) / 1000);
    //     time = max(time,1);
    //     SatParameters params = getDefaultSatParameters(cnf, time);
    //     if (!cnf.find_optimal_result) clog << "\t running cpsat solver with time limit of " << time << " sec." << endl;
    //     else clog << "\t running cpsat solver without time limit, looking for optimal result" << endl;
    //
    //     Model solver_model;
    //     solver_model.Add(NewSatParameters(params));
    //
    //     CpModelBuilder model;
    //     vector<BoolVar> nodes;
    //     for (int i=0; i<N; i++) nodes.push_back(model.NewBoolVar());
    //     addCycleConstraints(model,cycles,nodes);
    //
    //     VI init_sol = prev_res +  unhit_cycles_hs;
    //     addInitialSolutionHint(model,nodes,init_sol,prev_res, cnf);
    //     assert( isHS(cycles, init_sol) );
    //
    //     if (cnf.next_sol_max_dst_from_init_sol != inf && !prev_res.empty()) addMaxHammingDstConstraint(model,nodes,init_sol,cnf);
    //
    //     model.Minimize(LinearExpr::Sum(nodes));
    //     auto model_proto = model.Build();
    //     // CpSolverResponse response = SolveCpModel(model.Build(), &solver_model);
    //     CpSolverResponse response = SolveCpModel(model_proto, &solver_model);
    //
    //     // while(response.status() == CpSolverStatus::UNKNOWN && !timer.tle(timer_option)) {
    //     //     time = ceil(min( 3000.0 * time, timer.getLimit(timer_option) - timer.getTime(timer_option) ) / 1000);
    //     //     clog << "Response status unknown, increasing max_time_per_iter to " << time << endl;
    //     //
    //     //     params.set_max_time_in_seconds(time);
    //     //     Model solv_model;
    //     //     solv_model.Add(NewSatParameters(params));
    //     //     response = SolveCpModel(model_proto, &solv_model);
    //     // }
    //     //
    //     // if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {
    //     //     VI res;
    //     //     for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
    //     //     return res;
    //     // }
    //     // return {};
    //
    //     return rerunModelUntilFeasibleOrTle( model_proto,nodes, response,timer, timer_option,time,cnf );
    // };

    auto solveHSLocal = [&](auto & cycles, VI unhit_cycles_hs) -> VI {
        VI init_sol = prev_res +  unhit_cycles_hs;
        return solveCpsatForCycles(V,cycles,prev_res,init_sol,timer, timer_option, cnf);
    };

    constexpr int max_l = inf;


    for (int L=3; L <= max_l ; L++) {
        if ( !cnf.find_optimal_result && timer.tle(timer_option)) {
            clog << "Solver did not find optimal value for given set of cycles in admissible time" << endl;
            break;
        }

        exp_data.iterations.emplace_back();
        exp_data.iterations.back().res_size_before_impr = prev_res.size();

        clog << endl << "Considering cycles of length <= " << L << ", time: " << timer.getTime(timer_option) / 1000 << endl;

        Stopwatch s;

        s.start("cycles");
        int millis = max(10.0,timer.getLimit(timer_option) - timer.getTime(timer_option));
        auto cycles = Utils::getAllSimpleCycles3(V,L, millis );
        sort(ALL(cycles),[&](auto & c1, auto & c2){ return c1.size() < c2.size(); });
        s.stop("cycles");
        clog << "\t there are " << cycles.size() << " such cycles, found in time: " << s.getTime("cycles") / 1000 << endl;
        if (cycles.size() > cnf.max_cycles_for_hs) break;

        {
            int all_arcs = GraphUtils::countEdges(V,true);
            set<PII> zb;
            for( auto & cyc : cycles ) for( int j=cyc.size()-1, i=0; i < cyc.size(); j = i++ ) zb.insert( {cyc[j], cyc[i]} );
            int arcs = zb.size();
            clog << "\t arcs in constraints: " << arcs << " / " << all_arcs << endl;
        }

        s.start("hs_greedy");
        VVI new_cycles;
        for (auto & cyc : cycles) if ( L <= 3 || cyc.size() == L) new_cycles.push_back(cyc);
        VI unhit_cycles_hs = getUnhitCyclesHSGreedy(new_cycles);
        assert( isHS(new_cycles, unhit_cycles_hs) );
        s.stop("hs_greedy");
        clog << "\t found HS for " << new_cycles.size() << " unhit cycles of size " << unhit_cycles_hs.size() << endl;


        VI sol = solveHSLocal(cycles, unhit_cycles_hs);
        clog << "\t found hs of size sol.size(): " << sol.size() << endl;
        clog << "\t "; DEBUG(Utils::isFVS(V,sol));

        while ( !cnf.find_optimal_result && Utils::isFVS(V,sol) && !timer.tle(timer_option) ) {
            cnf.ihs_single_iteration_sec *= 3;
            sol = solveHSLocal(cycles, unhit_cycles_hs);
            clog << "\t found hs of size sol.size(): " << sol.size() << endl;
            clog << "\t "; DEBUG(Utils::isFVS(V,sol));
        }

        {
            set<PII> arcs;
            for (auto & cyc : cycles) for ( int i=0, j=(int)cyc.size()-1; i < cyc.size(); j = i++ ) arcs.insert( {cyc[j],cyc[i]});

            map<int,int> cycles_of_length;
            for (auto & cyc : cycles) cycles_of_length[cyc.size()]++;

            // gather statistics
            exp_data.iterations.back().res_size_after_impr = sol.size();
            exp_data.iterations.back().unhit_cycle_enumeration_time_millis = s.getTime("cycles");
            exp_data.iterations.back().hs_greedy_time = s.getTime("hs_greedy");
            exp_data.iterations.back().distinct_arcs_in_all_cycles = arcs.size();
            exp_data.iterations.back().res_valid = Utils::isFVS(V,sol);
            exp_data.iterations.back().res_optimal = cnf.find_optimal_result;
            exp_data.iterations.back().time_since_start_millis = timer.getTime(timer_option);
            exp_data.iterations.back().total_cycles = cycles.size();
            exp_data.iterations.back().unhit_graph_sizes = getUnhitGraphSizes(V,sol);
            exp_data.iterations.back().cycles_of_length = cycles_of_length;
            exp_data.iterations.back().new_cycles_added = new_cycles.size();
            exp_data.iterations.back().new_cycles_found = new_cycles.size();
        }

        if ( cnf.find_optimal_result && Utils::isFVS(V,sol) ) break;

        prev_res = sol;

        clog << "\t found sol.size(): " << sol.size() << ", but it is not a FVS, increasing cycle length" << endl;
    }

    clog << "\t final solution size: " << exp_data.iterations.back().res_size_after_impr << endl;
    timer.stop(timer_option);
    timer.write(timer_option);

    return exp_data;
}

ExpData CpsatExp1::solveIHS(VVI V, ExpConfig cnf) {
    VVI cycles;
    VI res;
    return solveIHS(V,cnf,cycles,res);
}

ExpData CpsatExp1::solveIHS(VVI V, ExpConfig cnf, VVI & cycles, VI & res) {
    clog << "Solving using CpsatExp1::solveIHS" << endl;
    cnf.writeConfig();

    int cycle_enumeration_type = cnf.unhit_cycle_enumeration_type;
    ExpData exp_data;

    int N = V.size();

    cycles.clear();
    VI prev_res;

    VI best_fvs;

    Stopwatch timer;
    string timer_option = "solveIHS";
    timer.setLimit(timer_option, cnf.find_optimal_result ? inf : cnf.max_time_sec * 1000);
    timer.start(timer_option);

    // auto updateCyclesAndHit = [&](VVI new_cycles, VI unhit_cycles_hs) -> VI {
    //     cycles += new_cycles;
    //
    //     int F = ( Utils::isFVS(V,prev_res) ? 5 : 1 );
    //     int time = ceil(min( 1000.0 * F * cnf.ihs_single_iteration_sec, timer.getLimit(timer_option) - timer.getTime(timer_option) ) / 1000);
    //     time = max(time,1);
    //     SatParameters params = getDefaultSatParameters(cnf, time);
    //     if (!cnf.find_optimal_result) clog << "\t running cpsat solver with time limit of " << time << " sec." << endl;
    //     else clog << "\t running cpsat solver without time limit, looking for optimal result" << endl;
    //
    //     Model solver_model;
    //     solver_model.Add(NewSatParameters(params));
    //
    //     CpModelBuilder model;
    //     vector<BoolVar> nodes;
    //     for (int i=0; i<N; i++) nodes.push_back(model.NewBoolVar());
    //     for ( auto & c : cycles ) {
    //         vector<BoolVar> cyc_vars;
    //         for ( int d : c ) cyc_vars.push_back(nodes[d]);
    //         model.AddAtLeastOne(cyc_vars);
    //     }
    //
    //     // this theoretically is bounding the solution from below, but it apparently can hurt performance...
    //     // if(cnf.find_optimal_result && !prev_res.empty()) model.AddGreaterOrEqual(LinearExpr::Sum(nodes), (int)prev_res.size());
    //
    //     VI init_sol = prev_res +  unhit_cycles_hs;
    //     addInitialSolutionHint(model,nodes,init_sol,prev_res, cnf);
    //     assert( isHS(cycles, init_sol) );
    //
    //     if (cnf.next_sol_max_dst_from_init_sol != inf && !prev_res.empty()) addMaxHammingDstConstraint(model,nodes,init_sol,cnf);
    //
    //     model.Minimize(LinearExpr::Sum(nodes));
    //     auto model_proto = model.Build();
    //     CpSolverResponse response = SolveCpModel(model_proto, &solver_model);
    //
    //     // while(response.status() == UNKNOWN && !timer.tle(timer_option)) {
    //     //     time = ceil(min( 3000.0 * time, timer.getLimit(timer_option) - timer.getTime(timer_option) ) / 1000);
    //     //     clog << "Response status unknown, increasing max_time_per_iter to " << time << endl;
    //     //
    //     //     params.set_max_time_in_seconds(time);
    //     //     Model solv_model;
    //     //     solv_model.Add(NewSatParameters(params));
    //     //     response = SolveCpModel(model_proto, &solv_model);
    //     // }
    //     //
    //     // if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {
    //     //     VI res;
    //     //     for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
    //     //     return res;
    //     // }
    //     //
    //     // return {};
    //
    //     return rerunModelUntilFeasibleOrTle( model_proto,nodes, response,timer, timer_option,time,cnf );
    // };

    auto solveHSLocal = [&](auto & cycles, VI unhit_cycles_hs) -> VI {
        VI init_sol = prev_res +  unhit_cycles_hs;
        return solveCpsatForCycles(V,cycles,prev_res,init_sol,timer, timer_option, cnf);
    };



    int L = 2;
    int iters_done = 0;
    int old_cycles = 0;

    while(true) {
        if(timer.tle(timer_option)) break;
        if (iters_done++ > cnf.ihs_max_iterations) break;

        clog << endl << "Looking for new cycles, iter #" << iters_done << ", cycles.size(): " << cycles.size() << ", prev_res.size(): "
             << prev_res.size() << ", time: " << (int)timer.getTime(timer_option) / 1000 << endl;

        exp_data.iterations.emplace_back();
        exp_data.iterations.back().res_size_before_impr = prev_res.size();

        Stopwatch s;

        s.start("cycles");
        VVI new_cycles = getUnhitChordlessCycles(V,prev_res,L,200*cnf.ihs_single_iteration_sec, cycle_enumeration_type );
        s.stop("cycles");


        clog << "\t there are " << new_cycles.size() << " new cycles found, found in time "
             << s.getTime("cycles") / 1000 << endl;

        {
            int all_arcs = GraphUtils::countEdges(V,true);
            set<PII> zb;
            for( auto & cyc : cycles ) for( int j=cyc.size()-1, i=0; i < cyc.size(); j = i++ ) zb.insert( {cyc[j], cyc[i]} );
            int arcs = zb.size();
            clog << "\t arcs in constraints: " << arcs << " / " << all_arcs << endl;
        }


        // we cannot add more than MAX_NEW_CYCLES_PER_ITERATION_PERC * cycles.size() new cycles in each iteration
        // this is here, because when increasing the length size, we might get an awful lot of new cycles of
        // that length, we do not want that, we want to keep number of cycles used for constraints as small as possible
        // constexpr int MIN_NEW_CYCLES = 100;
        int MIN_NEW_CYCLES = sqrt(GraphUtils::countEdges(V,true));
        int I = exp_data.iterations.size();
        if ( new_cycles.size() < MIN_NEW_CYCLES &&
            ( exp_data.iterations.size() <= 3 || exp_data.iterations.back().new_cycles_found != exp_data.iterations[I-3].new_cycles_found )
            ) new_cycles.clear();

        constexpr double MAX_NEW_CYCLES_PER_ITERATION_PERC = 1;
        const int MAX_NEW_CYCLES = N * log(N);

        bool cond1 = ( new_cycles.size() > MAX_NEW_CYCLES );
        // bool cond2 = ( cycles.size() >= N && new_cycles.size() > MAX_NEW_CYCLES_PER_ITERATION_PERC * cycles.size());
        if( cond1 ) {
            StandardUtils::shuffle(new_cycles);
            new_cycles.resize( MAX_NEW_CYCLES );
        }

        if( new_cycles.empty() && !Utils::isFVS(V,prev_res) ) {
            L++;
            clog << endl << "---> INCREASING LENGTH, L: " << L << endl << endl;
            continue;
        }

        clog << "\t adding " << new_cycles.size() << " new cycles to cycles, cycles.size(): " << cycles.size() << endl;

        s.start("hs_greedy");
        VI unhit_cycles_hs = getUnhitCyclesHSGreedy(new_cycles);
        assert( isHS(new_cycles, unhit_cycles_hs) );
        s.stop("hs_greedy");
        clog << "\t found HS for " << new_cycles.size() << " unhit cycles of size " << unhit_cycles_hs.size() << endl;

        // VI sol = updateCyclesAndHit(new_cycles, unhit_cycles_hs);
        cycles += new_cycles;
        VI sol = solveHSLocal(cycles, unhit_cycles_hs);
        clog << "\t found hs of size sol.size(): " << sol.size() << endl;
        clog << "\t "; DEBUG(Utils::isFVS(V,sol));


        if ( Utils::isFVS(V,sol) ) {
            if( best_fvs.empty() || sol.size() < best_fvs.size() ) best_fvs = sol;

            if(cnf.find_optimal_result) {
                assert( best_fvs.size() == sol.size() );
                break;
            }else {
                cnf.ihs_single_iteration_sec *= 3;
                clog << endl << "--> Found a valid FVS, increasing max_time_seconds_per_iter to " << cnf.ihs_single_iteration_sec << " secc." << endl << endl;
            }
        }

        {
            set<PII> arcs;
            for (auto & cyc : cycles) for ( int i=0, j=(int)cyc.size()-1; i < cyc.size(); j = i++ ) arcs.insert( {cyc[j],cyc[i]});

            map<int,int> cycles_of_length;
            for (auto & cyc : cycles) cycles_of_length[cyc.size()]++;

            // gather statistics
            exp_data.iterations.back().res_size_after_impr = sol.size();
            exp_data.iterations.back().unhit_cycle_enumeration_time_millis = s.getTime("cycles");
            exp_data.iterations.back().hs_greedy_time = s.getTime("hs_greedy");
            exp_data.iterations.back().distinct_arcs_in_all_cycles = arcs.size();
            exp_data.iterations.back().res_valid = Utils::isFVS(V,sol);
            exp_data.iterations.back().res_optimal = cnf.find_optimal_result;
            exp_data.iterations.back().time_since_start_millis = timer.getTime(timer_option);
            exp_data.iterations.back().total_cycles = cycles.size();
            exp_data.iterations.back().unhit_graph_sizes = getUnhitGraphSizes(V,sol);
            exp_data.iterations.back().cycles_of_length = cycles_of_length;
            exp_data.iterations.back().new_cycles_added = cycles.size() - old_cycles;
            exp_data.iterations.back().new_cycles_found = new_cycles.size();
        }

        prev_res = sol;
        old_cycles = cycles.size();

        clog << "\t found sol.size(): " << sol.size() << ", best_fvs.size(): " << best_fvs.size()
             << ", while best_fvs.size(): " << best_fvs.size() << endl;
    }

    timer.stop(timer_option);
    timer.write(timer_option);


    if (!best_fvs.empty()) res = best_fvs;
    else res = prev_res;

    clog << "\t final solution size: " << exp_data.iterations.back().res_size_after_impr << endl;

    return exp_data;
}

ExpData CpsatExp1::solveMTZ(VVI V, ExpConfig cnf, int auxiliary_cycles_mode) {
    clog << "Solving using CpsatExp1::solveMTZ" << endl;
    cnf.writeConfig();

    ExpData exp_data;

    Stopwatch timer;
    timer.start("solveMTZ");

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
        sort(ALL(cyc),[&](auto & c1, auto & c2){ return c1.size() < c2.size(); });
        clog << "\t adding " << cyc.size() << " constraints for all simple cycles of length <= " << L0 << endl;
        addCycleConstraints(model,cyc,nodes);
    }else if (auxiliary_cycles_mode == 2) {
        auto new_cnf = cnf;
        new_cnf.ihs_max_iterations = cnf.ihs_iterations_in_mtz;
        VVI cycles;
        VI res;
        auto r = solveIHS(V,new_cnf,cycles, res);
        addCycleConstraints(model,cycles,nodes);
        init_sol = res;
    }

    if(!init_sol.empty()) {
        in_init_sol = StandardUtils::toVB(N,init_sol);
        addInitialSolutionHint(model,nodes,init_sol,init_sol,cnf); // here we intentionally add init_sol instead of prev_res

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

    SatParameters params = getDefaultSatParameters(cnf, cnf.max_time_sec);
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

    clog << "\t final solution size: " << exp_data.iterations.back().res_size_after_impr << endl;
    timer.stop("solveMTZ");
    timer.write("solveMTZ");

    return exp_data;
}

ExpData CpsatExp1::solveDiVerSeS(VVI V, ExpConfig cnf) {
    ExpData exp_data;
    assert(false && "not implemented");

    return exp_data;
}

ExpData CpsatExp1::solve(VVI V, ExpConfig cnf) {
    auto alg = cnf.alg;
    if (alg == Algorithm::HS) return solveHS(V,cnf);
    if (alg == Algorithm::IHS) return solveIHS(V,cnf);
    if (alg == Algorithm::MTZ) return solveMTZ(V,cnf);
    if (alg == Algorithm::DIVERSES) return solveDiVerSeS(V,cnf);

    return ExpData{};
}


