//
// Created by sylwe on 08/04/2026.
//

#include "CpsatExp1.h"

#include "GraphUtils.h"
#include "StandardUtils.h"
#include "Stopwatch.h"
#include "CONTESTS/PACE22/Reducer.h"
#include "CONTESTS/PACE22/Utils.h"
#include "CONTESTS/PACE22/heur/DFVSSolverH.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"
#include "scc/StronglyConnectedComponents.h"
#include "mutex"
using namespace operations_research::sat;

SatParameters getDefaultSatParameters(ExpConfig cnf, int time_in_sec) {
    SatParameters params;
    if (!cnf.find_optimal_result) params.set_max_time_in_seconds(time_in_sec);
    if(cnf.use_only_cpsat_lns) params.set_use_lns_only(true);
    params.set_num_search_workers(cnf.threads);
    params.set_log_search_progress(cnf.log_cpsat_search_progress);

    if(cnf.focus_mostly_onh_heuristics) {
        IntGenerator rnd;
        params.set_randomize_search(true);
        params.set_random_seed(rnd.nextInt(inf));

        params.set_use_lns(true);
        params.set_symmetry_level(0); // disables detecting symmetries - often costs more than it helps, especially in hard instances
        params.set_use_sat_inprocessing(false);

        params.set_linearization_level(0);
        params.set_cut_level(0);
        params.set_add_lp_constraints_lazily(false);

        // params.set_optimize_with_core(false); // without this it seems to work better
        params.set_use_exact_lp_reason(false);

        // The [set_use_feasibility_jump] enables a local search / repair heuristic inside CP-SAT. Instead of only doing: systematic branching +
        // propagation the solver also does: start from a (possibly infeasible) assignment iteratively “repair” it by flipping
        // variables try to reach feasibility and then improve objective
        params.set_use_feasibility_jump(true);
        params.set_random_branches_ratio(0.1); // controls how often CP-SAT ignores its heuristics and makes a random branching decision.
        params.set_random_polarity_ratio(0.015); // controls how often CP-SAT assigns a variable, not based on heuristics but randomly to true or false

        params.set_max_number_of_conflicts(100'000);
        // params.set_max_deterministic_time(0.1);
    }

    return params;
}


bool ExpData::foundValidHSResult(vector<IterationEntry> &entries) {
    return ranges::any_of(entries, [&](auto & ie) {
        return ie.hs_valid_fvs;
    });
}

// #define ADD_ENTRY(x) entries.back().emplace(#x,to_string(x))
// #define ADD_ENTRY(x) entries.back().emplace(#x,to_string("I".x))
#define ADD_ENTRY(x,y) entries.back().emplace(#y,to_string(x.y))

vector<map<string, string>> ExpData::getIterationEntries() {
    vector<map<string,string>> entries;

    for (int i=0; i<iterations.size(); i++) {
        entries.emplace_back();
        const auto & I = iterations[i];

        entries.back().emplace("iteration_id",to_string(i+1));
        ADD_ENTRY(I,hs_size_before_impr);
        ADD_ENTRY(I,hs_size_after_impr);
        ADD_ENTRY(I,hs_valid_fvs);
        ADD_ENTRY(I,improved_best_res);
        ADD_ENTRY(I,res_optimal);
        ADD_ENTRY(I,res_lower_bound);
        ADD_ENTRY(I,distinct_arcs_in_all_cycles);
        ADD_ENTRY(I,unhit_cycle_enumeration_time_millis);
        ADD_ENTRY(I,hs_greedy_time);
        entries.back().emplace("unhit_graph_sizes", to_string(I.unhit_graph_sizes.first) + " " + to_string(I.unhit_graph_sizes.second));
        ADD_ENTRY(I,unhit_graph_greedy_dfvs_size);
        ADD_ENTRY(I,unhit_graph_greedy_dfvs_time);
        ADD_ENTRY(I,total_cycles);
        ADD_ENTRY(I,new_cycles_added);
        ADD_ENTRY(I,new_cycles_found);
        ADD_ENTRY(I,time_since_start_millis);
        ADD_ENTRY(I,iteration_time);
        ADD_ENTRY(I,max_cycle_length);
        ADD_ENTRY(I,best_result_so_far);
        // ADD_ENTRY(I,cycle_hs_valid_dfvs);

        string str;
        for (auto [k,v] : I.cycles_of_length) {
            str += "(" + to_string(k) + "->" + to_string(v) + ") ";
        }
        entries.back().emplace("cycles_of_length", str);
    }

    return entries;
}

void ExpData::updateBestResultSoFar() {
    // iterations[0].best_result_so_far = (iterations[0].hs_valid_fvs ? iterations[0].hs_size_after_impr : 0);
    // for (int i=1; i < iterations.size(); i++) iterations[i].best_result_so_far = inf;
    // for ( int i=1; i<iterations.size(); i++ ) {
    //     auto & I0 = iterations[i-1];
    //     auto & I1 = iterations[i];
    //     if ( I1.hs_valid_fvs ) I1.best_result_so_far = I1.hs_size_after_impr;
    //     I1.best_result_so_far = min(I1.best_result_so_far, I0.best_result_so_far);
    // }
    // for ( auto & I : iterations ) if ( I.best_result_so_far == inf ) I.best_result_so_far = 0;

    iterations[0].best_result_so_far = iterations[0].full_sol_size;
    for ( int i=1; i<iterations.size(); i++ ) {
        auto & I0 = iterations[i-1];
        auto & I1 = iterations[i];
        I1.best_result_so_far = min(I1.full_sol_size, I0.best_result_so_far);
        I1.improved_best_res = (I1.best_result_so_far < I0.best_result_so_far);
    }

}

void ExpData::writeToFile(ExpConfig cnf) {
    ofstream str(cnf.metadata_filepath);

    auto entries = getIterationEntries();
    using VS = vector<string>;
    VS header;
    for ( auto s : views::keys(entries[0]) ) header.push_back(s);
    for ( int i=0; i<header.size(); i++ ) if ( header[i] == "iteration_id" ) {
        swap(header[0], header[i]);
        sort(header.begin()+1, header.end());
        break;
    }

    map<string,string> config_entries;
    vector<string> config_header;
    {
        auto cnf_entries = cnf.getConfigEntries();
        for(auto [a,b] : cnf_entries) config_entries.emplace(a,b);
        for(const auto & s : views::keys(cnf_entries)) config_header.push_back(s);
        // DEBUG(cnf_entries); DEBUG(config_header); DEBUG(config_entries);
    }

    for (auto [i,s] : views::enumerate(header)) str << (i > 0 ? "," : "") << s;
    for (const auto& s : config_header ) str << "," << s;
    str << "\n";

    for ( auto & I : entries ) {
        for (const auto & [i,key] : header | views::enumerate) {
            if (i > 0) str << ",";
            str << I[key];
        }
        for (const auto & key : config_header) str << "," << config_entries[key];
        str << "\n";
    }

    str.close();
}

VVI CpsatExp1::getUnhitChordlessCycles(VVI &V, VI &S, int max_l, int max_millis, int enumeration_option) {
    clog << "\t Looking for cycles for at most " << max_millis << " millis" << endl;
    VVI all_cycles;

    if (enumeration_option == 3) max_millis /= 2;

    if (enumeration_option == 1 || enumeration_option == 3) {
        VVI H = V;
        VVI revH = GraphUtils::reverseGraph(H);
        VB helper(V.size());
        Utils::removeNodes(H, revH, S, helper);
        auto cycles = Utils::getAllSimpleCycles3(H,max_l, max_millis);
        if (enumeration_option == 1) return cycles;
        all_cycles += cycles;
    }

    int N = V.size();
    vector<LL> hashes(N);
    IntGenerator rnd;
    for (int i=0; i<N; i++) hashes[i] = rnd.rand();
    auto getHash = [&](VI cyc)-> LL {
        LL h = 0;
        for (int d : cyc) h ^= hashes[d];
        return h;
    };

    if (enumeration_option == 2 || enumeration_option == 3){
        auto indg = getUnhitGraph(V,S);
        VVI H = indg.V;
        int N = H.size();

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
            // clog << "Getting cycles for node v: " << v << endl;
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
            // if (!new_cycles.empty()) {
            //     clog << "\t found " << new_cycles.size() << " new cycles of avg_length: " << suml1 / new_cycles.size()
            //          << " from which " << new_chordless_cycles << " were already chordless --> only added "
            //          << cycles_added << " new chordless cycles, with avg_length: " << suml2 / cycles_added << endl;
            // }
        }

        for (auto & cyc : cycles) indg.remapNodes(cyc);
        if ( enumeration_option == 2 ) return cycles;
        all_cycles += cycles;
    }

    map<LL,VI> unique_cycles;
    for ( auto cyc : all_cycles ) unique_cycles[getHash(cyc)] = cyc;
    VVI res;
    for( auto & cyc : views::values(unique_cycles) ) res.push_back(cyc);
    clog << "\t\t there are " << res.size() << " unique cycles, out of " << all_cycles.size() << " all cycles" << endl;
    return res;
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

VI CpsatExp1::getUnhitGraphGreedyFVS(VVI &V, VI &S) {
    auto indg = getUnhitGraph(V,S);

    Config diverses_cnf;
    diverses_cnf.write_logs = false;
    diverses_cnf.agent_flow_min_distance = 3;
    // diverses_cnf.agent_flow_node_update_frequency = 1;
    // diverses_cnf.agent_flow_max_distance_from_best = 1;
    diverses_cnf.agent_flow_max_distance_from_best = 1 + sqrt(V.size());
    diverses_cnf.agent_flow_method = Config::agent_flow_sinkhorn;
    diverses_cnf.agent_flow_node_selection_type = Config::agent_flow_remove_largest_flow_node;
    diverses_cnf.agent_flow_alternate_selection_type = diverses_cnf.solver_improve_alternate_selection_type = true;

    diverses_cnf.disableAllNonbasicReductions();
    // diverses_cnf.reducer_use_pie = diverses_cnf.reducer_use_dome = diverses_cnf.reducer_use_folding = true;
    diverses_cnf.reducer_use_pie = diverses_cnf.reducer_use_dome = true;

    Reducer red(indg.V, diverses_cnf);
    auto reductions = red.reduce();
    auto newV = red.V;
    assert(newV.size() == indg.V.size());

    DFVSSolverH sh(diverses_cnf);
    auto res = sh.solveByAgentFlow(newV);
    assert(Utils::isFVS(newV,res));

    Reducer::liftSolution(newV.size(), res, reductions);

    assert(Utils::isFVS(indg.V,res));
    indg.remapNodes(res);

    return res;
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
    if (cnf.use_init_sol_as_hint_mode == 1 || cnf.use_init_sol_as_hint_mode == 3) for(int i=0; i<N; i++) if (in_init_sol[i]) model.AddHint(nodes[i],1);
    if (cnf.use_init_sol_as_hint_mode == 2) for(int d : prev_res) model.AddHint(nodes[d],1);
}

void CpsatExp1::addInitialSolutionSizeConstraint(CpModelBuilder &model, vector<BoolVar> &nodes, VI &init_sol) {
    model.AddLessOrEqual(LinearExpr::Sum(nodes), init_sol.size());
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

tuple<VI,CpSolverStatus,VI> CpsatExp1::rerunModelUntilFeasibleOrTle(VVI & V, CpModelProto &model_proto, vector<BoolVar> &nodes,
    CpSolverResponse & response, Stopwatch & timer, string timer_option, VI & init_sol, int init_time, ExpConfig& cnf) {

    VI inter_fvs;

    while(response.status() == CpSolverStatus::UNKNOWN && !timer.tle(timer_option)) {
        if (response.status() == CpSolverStatus::UNKNOWN) cnf.ihs_single_iteration_sec++;

        init_time = ceil(min( 3000.0 * init_time, timer.getLimit(timer_option) - timer.getTime(timer_option) ) / 1000);
        clog << "Response status unknown, increasing max_time_per_iter to " << init_time << endl;

        SatParameters params = getDefaultSatParameters(cnf, init_time);
        Model solv_model;
        solv_model.Add(NewSatParameters(params));
        if ( cnf.check_incumbent_cpsat_solutions ) {
            mutex log_mutex;
            solv_model.Add(NewFeasibleSolutionObserver(
                [&](const CpSolverResponse& r) {
                    int cnt = 0;
                    for (int i = 0; i < nodes.size(); ++i) cnt += SolutionBooleanValue(r, nodes[i]);
                    if ( !init_sol.empty() && cnt >= init_sol.size() ) return;

                    VI temp; temp.reserve(nodes.size());
                    for (int i = 0; i < nodes.size(); ++i) if (SolutionBooleanValue(r, nodes[i])) temp.push_back(i);
                    if ( Utils::isFVS(V,temp) ) {
                        std::lock_guard<std::mutex> lk(log_mutex);
                        inter_fvs = temp;
                        clog << "\t\t\t found new inter_fvs of size: " << inter_fvs.size() << endl;
                    }
                    // else { // uncomment this to check the a 'fixed' solution every time a new hitting set is found
                    //     std::lock_guard<std::mutex> lk(log_mutex);
                    //     auto unhit_graph_dfvs = getUnhitGraphGreedyFVS(V,temp);
                    //     if (temp.size() + unhit_graph_dfvs.size() < inter_fvs.size() ||
                    //         (inter_fvs.empty() && temp.size() + unhit_graph_dfvs.size() < init_sol.size())
                    //         ) {
                    //         inter_fvs = temp + unhit_graph_dfvs;
                    //         assert(Utils::isFVS(V,inter_fvs));
                    //         clog << "\t\t\t found SUPPL. inter_fvs of size: " << inter_fvs.size() << endl;
                    //     }
                    // }
                }
            ));
        }
        response = SolveCpModel(model_proto, &solv_model);
    }

    int N = nodes.size();
    if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {
        VI res;
        for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
        return {res,response.status(), inter_fvs};
    }

    return { {}, CpSolverStatus::UNKNOWN, VI{} };
}

tuple<VI,CpSolverStatus, VI> CpsatExp1::solveCpsatForCycles(VVI &V, VVI &cycles, VI &prev_res, VI &init_sol, Stopwatch &timer,
    string timer_option, ExpConfig& cnf) {

    int N = V.size();
    int F = ( Utils::isFVS(V,prev_res) ? 2 : 1 );
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
    // addInitialSolutionSizeConstraint(model, nodes, init_sol); // #TEST
    assert( isHS(cycles, init_sol) );

    if (cnf.next_sol_max_dst_from_init_sol != inf && !prev_res.empty()) addMaxHammingDstConstraint(model,nodes,init_sol,cnf);

    VI inter_fvs;
    if ( cnf.check_incumbent_cpsat_solutions ) {
        mutex log_mutex;
        solver_model.Add(NewFeasibleSolutionObserver(
            [&](const CpSolverResponse& r) {
                int cnt = 0;
                for (int i = 0; i < nodes.size(); ++i) cnt += SolutionBooleanValue(r, nodes[i]);
                if ( !init_sol.empty() && cnt >= init_sol.size() ) return;

                VI temp; temp.reserve(cnt);
                for (int i = 0; i < nodes.size(); ++i) if (SolutionBooleanValue(r, nodes[i])) temp.push_back(i);
                if ( Utils::isFVS(V,temp) ) {
                    std::lock_guard<std::mutex> lk(log_mutex);
                    inter_fvs = temp;
                    clog << "\t\t\t found new inter_fvs of size: " << inter_fvs.size() << endl;
                }
                // else { // uncomment this to check the a 'fixed' solution every time a new hitting set is found
                //     std::lock_guard<std::mutex> lk(log_mutex);
                //     auto unhit_graph_dfvs = getUnhitGraphGreedyFVS(V,temp);
                //     if (temp.size() + unhit_graph_dfvs.size() < inter_fvs.size() ||
                //         (inter_fvs.empty() && temp.size() + unhit_graph_dfvs.size() < init_sol.size())
                //         ) {
                //         inter_fvs = temp + unhit_graph_dfvs;
                //         assert(Utils::isFVS(V,inter_fvs));
                //         clog << "\t\t\t found SUPPL. inter_fvs of size: " << inter_fvs.size() << endl;
                //     }
                // }
            }
        ));
    }

    model.Minimize(LinearExpr::Sum(nodes));
    auto model_proto = model.Build();
    CpSolverResponse response = SolveCpModel(model_proto, &solver_model);
    // if (response.status() == CpSolverStatus::UNKNOWN) cnf.ihs_single_iteration_sec++;

    // return rerunModelUntilFeasibleOrTle( V, model_proto,nodes, response,timer, timer_option,time,cnf );
    auto [sol,resp,fvs] = rerunModelUntilFeasibleOrTle( V, model_proto,nodes,
        response,timer, timer_option, init_sol, time,cnf );
    if ( !inter_fvs.empty() && (fvs.empty() || inter_fvs.size() < fvs.size()) ) fvs = inter_fvs;
    return {sol,resp,fvs};
}

VI CpsatExp1::getUnhitCyclesHSGreedy(VVI &cycles) {
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


    auto solveHSLocal = [&](auto & cycles, VI unhit_cycles_hs) {
        VI init_sol = prev_res +  unhit_cycles_hs;
        return solveCpsatForCycles(V,cycles,prev_res,init_sol,timer, timer_option, cnf);
    };

    constexpr int max_l = inf;


    for (int L=cnf.init_L_for_all_constraints; L <= max_l ; L++) {
        if ( !cnf.find_optimal_result && timer.tle(timer_option)) {
            clog << "Solver did not find optimal value for given set of cycles in admissible time" << endl;
            break;
        }


        clog << endl << "Considering cycles of length <= " << L << ", time: " << timer.getTime(timer_option) / 1000 << endl;

        Stopwatch s;
        s.start("iteration");

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
        for (auto & cyc : cycles) if ( L <= cnf.init_L_for_all_constraints || cyc.size() == L) new_cycles.push_back(cyc);
        VI unhit_cycles_hs = getUnhitCyclesHSGreedy(new_cycles);
        assert( isHS(new_cycles, unhit_cycles_hs) );
        s.stop("hs_greedy");
        clog << "\t found HS for " << new_cycles.size() << " unhit cycles of size " << unhit_cycles_hs.size() << endl;

        if (timer.tle(timer_option)) break;

        auto [sol,response_status, inter_fvs] = solveHSLocal(cycles, unhit_cycles_hs);
        clog << "\t found hs of size sol.size(): " << sol.size() << endl;
        clog << "\t "; DEBUG(Utils::isFVS(V,sol));

        if (cnf.find_optimal_result) assert(response_status == CpSolverStatus::OPTIMAL);


        if ( !cnf.find_optimal_result && Utils::isFVS(V,sol) && !timer.tle(timer_option) ) {
            cnf.ihs_single_iteration_sec *= 3;
            tie(sol,response_status, inter_fvs) = solveHSLocal(cycles, unhit_cycles_hs);
            clog << "\t found hs of size sol.size(): " << sol.size() << endl;
            clog << "\t "; DEBUG(Utils::isFVS(V,sol));
        }

        s.stop("iteration");

        VI unhit_graph_dfvs;
        VI full_sol;
        if ( cnf.fill_partial_result_using_greedy_fvs ) {
            s.start("unhit_graph_greedy_dfvs_time");
            unhit_graph_dfvs = getUnhitGraphGreedyFVS(V,sol);
            full_sol = sol + unhit_graph_dfvs;
            assert(Utils::isFVS(V,full_sol));
            s.stop("unhit_graph_greedy_dfvs_time");
            DEBUG(unhit_graph_dfvs.size());
            DEBUG(full_sol.size());
            DEBUG(s.getTime("unhit_graph_greedy_dfvs_time"));
        }

        if ( !inter_fvs.empty() && inter_fvs.size() < full_sol.size() ) full_sol = inter_fvs;

        {
            exp_data.iterations.emplace_back();
            exp_data.iterations.back().hs_size_before_impr = prev_res.size();
            exp_data.iterations.back().max_cycle_length = L;

            set<PII> arcs;
            for (auto & cyc : cycles) for ( int i=0, j=(int)cyc.size()-1; i < cyc.size(); j = i++ ) arcs.insert( {cyc[j],cyc[i]});

            map<int,int> cycles_of_length;
            for (auto & cyc : cycles) cycles_of_length[cyc.size()]++;

            // gather statistics
            exp_data.iterations.back().full_sol_size = full_sol.size();
            exp_data.iterations.back().hs_size_after_impr = sol.size();
            exp_data.iterations.back().unhit_cycle_enumeration_time_millis = s.getTime("cycles");
            exp_data.iterations.back().hs_greedy_time = s.getTime("hs_greedy");
            exp_data.iterations.back().distinct_arcs_in_all_cycles = arcs.size();
            exp_data.iterations.back().hs_valid_fvs = Utils::isFVS(V,sol);
            exp_data.iterations.back().res_optimal = (response_status == CpSolverStatus::OPTIMAL && sol.size() == full_sol.size());
            exp_data.iterations.back().time_since_start_millis = timer.getTime(timer_option);
            exp_data.iterations.back().total_cycles = cycles.size();
            exp_data.iterations.back().unhit_graph_sizes = getUnhitGraphSizes(V,sol);
            exp_data.iterations.back().cycles_of_length = cycles_of_length;
            exp_data.iterations.back().new_cycles_added = new_cycles.size();
            exp_data.iterations.back().new_cycles_found = new_cycles.size();
            exp_data.iterations.back().iteration_time = s.getTime("iteration");
            if ( !Utils::isFVS(V,sol) && cnf.fill_partial_result_using_greedy_fvs ) {
                exp_data.iterations.back().unhit_graph_greedy_dfvs_size = unhit_graph_dfvs.size();
                exp_data.iterations.back().unhit_graph_greedy_dfvs_time = s.getTime("unhit_graph_greedy_dfvs_time");
            }
        }

        if ( response_status == CpSolverStatus::OPTIMAL && Utils::isFVS(V,sol) ) break;

        prev_res = sol;

        clog << "\t found sol.size(): " << sol.size() << ", but it is not a FVS, increasing cycle length" << endl;
    }

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


    ExpData exp_data;

    int N = V.size();

    cycles.clear();
    VI prev_res;


    Stopwatch timer;
    string timer_option = "solveIHS";
    timer.setLimit(timer_option, cnf.find_optimal_result ? inf : cnf.max_time_sec * 1000);
    timer.start(timer_option);


    VI best_fvs;
    if(cnf.ihs_init_sol_creation_mode != 0) best_fvs = createInitialSolution(V,cnf);



    int L = cnf.init_L_for_all_constraints;
    int iters_done = 0;
    int old_cycles = 0;
    int iters_without_new_cycles = 0;
    int iters_streak_with_hs_valid_res = 0;

    auto solveHSLocal = [&](auto & cycles, VI unhit_cycles_hs) {
        VI init_sol = prev_res +  unhit_cycles_hs;
        if (cnf.use_init_sol_as_hint_mode == 3 && !best_fvs.empty() && best_fvs.size() < init_sol.size() ) init_sol = best_fvs;
        auto new_cnf = cnf;
        // if (cnf.use_init_sol_as_hint_mode == 3) new_cnf.use_init_sol_as_hint_mode = 1; // hint only variable set to 1
        if (cnf.use_init_sol_as_hint_mode == 3) new_cnf.use_init_sol_as_hint_mode = 0; // hint all variables to 0 or 1
        auto r = solveCpsatForCycles(V,cycles,prev_res,init_sol,timer, timer_option, new_cnf);
        cnf.ihs_single_iteration_sec = new_cnf.ihs_single_iteration_sec;
        return r;
    };

    const bool USE_CYCLE_TRIMMING = cnf.use_cycle_trimming;
    const int CYCLE_TRIMMING_FREQ = cnf.cycle_trimming_freq;
    const double CYCLE_TRIMMING_PROBAB = cnf.cycle_trimming_probab;
    const int MIN_NODES_IN_HS = cnf.cycle_trimming_min_nodes_in_hs;
    auto trimCycles = [&](VI & hs, VVI & cycles, int max_to_trim, int min_nodes_in_hs) {
        // clog << "\t Trimming cycles" << endl;
        int trimmed = 0;

        IntGenerator rnd;
        StandardUtils::shuffle(cycles,rnd);
        sort(ALL(cycles),[&](auto & c1, auto & c2){ return c1.size() < c2.size(); });
        VB in_hs = StandardUtils::toVB(N,hs);
        for ( int i=(int)cycles.size()-1; i>=0; i-- ) {
            int cnt = accumulate(ALL(cycles[i]),0, [&](int a, int b){ return a + in_hs[b]; });
            if (cnt >= min_nodes_in_hs) {
                swap(cycles[i],cycles.back());
                cycles.pop_back();
                max_to_trim--;
                trimmed++;
                if (max_to_trim == 0) break;
            }
        }
        clog << "\t trimmed " << trimmed << " cycles with at least " << min_nodes_in_hs << " nodes in hs" << endl;
    };

    IntGenerator rnd;
    while(true) {
        if(timer.tle(timer_option)) break;
        if (iters_done++ > cnf.ihs_max_iterations) break;

        if (cycles.size() > cnf.max_cycles_for_hs) break;

        clog << endl << "Looking for new cycles, iter #" << iters_done << ", L: " << L << ", cycles.size(): "
             << cycles.size()
             << ", prev_res.size(): " << prev_res.size() << ", best_fvs.size(): " << best_fvs.size()
             << ", time: " << (int)timer.getTime(timer_option) / 1000 << endl;


        if ( USE_CYCLE_TRIMMING && (iters_done % CYCLE_TRIMMING_FREQ == 0 || rnd.nextInt(1000) < 1000*CYCLE_TRIMMING_PROBAB ) ) {
            auto hs = best_fvs;

            // int R = (CYCLE_TRIMMING_FREQ == 1 ? 1 : 4);
            int R = 4;
            for (int i=0; i<R; i++) {
                trimCycles(hs,cycles,cnf.scaleIters(N) / (i+1), MIN_NODES_IN_HS + R-1-i); // orig

                // if (CYCLE_TRIMMING_FREQ > 1) trimCycles(hs,cycles,cnf.scaleIters(N) / (i+1), MIN_NODES_IN_HS + R-1-i); // orig
                // else trimCycles(hs,cycles,cnf.scaleIters(N) / (i+2), MIN_NODES_IN_HS + R-1-i);
            }
        }

        exp_data.iterations.emplace_back();
        // exp_data.iterations.back().hs_size_before_impr = prev_res.size();
        exp_data.iterations.back().max_cycle_length = L;

        Stopwatch s;
        s.start("iteration");

        s.start("cycles");
        int cycle_enumeration_type = cnf.unhit_cycle_enumeration_type;
        if ( cycle_enumeration_type == 3 && (iters_done % 2 == 0) ) cycle_enumeration_type = 1; // take every second iteration, just to be able to increase L after some time #original
        VVI new_cycles = getUnhitChordlessCycles(V,prev_res,L,200*cnf.ihs_single_iteration_sec, cycle_enumeration_type );


        {
            int all_arcs = GraphUtils::countEdges(V,true);
            set<PII> zb;
            for( auto & cyc : cycles ) for( int j=cyc.size()-1, i=0; i < cyc.size(); j = i++ ) zb.insert( {cyc[j], cyc[i]} );
            int arcs = zb.size();
            clog << "\t arcs in constraints: " << arcs << " / " << all_arcs << endl;
        }
        int new_cycles_found = new_cycles.size();
        s.stop("cycles");
        clog << "\t there are " << new_cycles.size() << " new cycles found, found in time "
             << s.getTime("cycles") / 1000 << endl;

        // we cannot add more than MAX_NEW_CYCLES_PER_ITERATION_PERC * cycles.size() new cycles in each iteration
        // this is here, because when increasing the length size, we might get an awful lot of new cycles of
        // that length, we do not want that, we want to keep number of cycles used for constraints as small as possible
        // constexpr int MIN_NEW_CYCLES = 100;
        int MIN_NEW_CYCLES = 2*sqrt(GraphUtils::countEdges(V,true));
        if (new_cycles.size() < MIN_NEW_CYCLES && iters_without_new_cycles <= 2 ) {
            new_cycles.clear();
            iters_without_new_cycles++;
        }


        const int MAX_NEW_CYCLES = cnf.scaleIters(N);

        bool cond1 = ( new_cycles.size() > MAX_NEW_CYCLES );
        if( ( cnf.unhit_cycle_enumeration_type >= 2 || L > cnf.init_L_for_all_constraints) &&
            cond1
            ) {
            StandardUtils::shuffle(new_cycles);
            sort(ALL(new_cycles),[&](auto & c1, auto & c2){ return c1.size() < c2.size(); });
            if ( L == cnf.init_L_for_all_constraints ) {
                int p = 0;
                while ( p < new_cycles.size() && new_cycles[p].size() <= L ) p++;
                new_cycles.resize( max(p,MAX_NEW_CYCLES) );
            }else {
                new_cycles.resize( MAX_NEW_CYCLES );
            }
        }

        if( new_cycles.empty() && !Utils::isFVS(V,prev_res) ) {
            L++;
            clog << endl << "---> INCREASING LENGTH, L: " << L << endl << endl;
            exp_data.iterations.pop_back();
            continue;
        }

        iters_without_new_cycles = 0;

        clog << "\t adding " << new_cycles.size() << " new cycles to cycles, cycles.size(): " << cycles.size() << endl;

        s.start("hs_greedy");
        VI unhit_cycles_hs = getUnhitCyclesHSGreedy(new_cycles);
        assert( isHS(new_cycles, unhit_cycles_hs) );
        s.stop("hs_greedy");
        clog << "\t found HS for " << new_cycles.size() << " unhit cycles of size " << unhit_cycles_hs.size() << endl;

        if (timer.tle(timer_option)) {
            exp_data.iterations.pop_back();
            break;
        }

        exp_data.iterations.back().hs_size_before_impr = prev_res.size() + unhit_cycles_hs.size();

        cycles += new_cycles;
        auto [sol,response_status, inter_fvs] = solveHSLocal(cycles, unhit_cycles_hs);
        clog << "\t found hs of size sol.size(): " << sol.size() << endl;
        clog << "\t "; DEBUG(Utils::isFVS(V,sol));
        if (cnf.find_optimal_result) assert(response_status == CpSolverStatus::OPTIMAL);

        VI full_sol;
        VI unhit_graph_dfvs;

        if ( Utils::isFVS(V,sol) ) {
            iters_streak_with_hs_valid_res++;

            if( best_fvs.empty() || sol.size() < best_fvs.size() ) best_fvs = sol;
            if(cnf.find_optimal_result) {
                assert( best_fvs.size() == sol.size() );
                break;
            }else {
                // cnf.ihs_single_iteration_sec++;
                cnf.ihs_single_iteration_sec += iters_streak_with_hs_valid_res;
                clog << endl << "--> Found a valid FVS, increasing max_time_seconds_per_iter to " << cnf.ihs_single_iteration_sec << " sec." << endl << endl;
            }

        }else iters_streak_with_hs_valid_res = 0;

        if ( cnf.fill_partial_result_using_greedy_fvs ) {
            s.start("unhit_graph_greedy_dfvs_time");
            unhit_graph_dfvs = getUnhitGraphGreedyFVS(V,sol);
            full_sol = sol + unhit_graph_dfvs;
            assert(Utils::isFVS(V,full_sol));
            s.stop("unhit_graph_greedy_dfvs_time");
            // DEBUG(sol.size());
            // DEBUG(unhit_graph_dfvs.size());
            // DEBUG(full_sol.size());
            // DEBUG(s.getTime("unhit_graph_greedy_dfvs_time"));
        }

        if ( !inter_fvs.empty() && inter_fvs.size() < full_sol.size() ) {
            clog << "\t\t\t inter_fvs improves full_sol!" << endl;
            full_sol = inter_fvs;
        }

        s.stop("iteration");

        auto addStats = [&]() {
            set<PII> arcs;
            for (auto & cyc : cycles) for ( int i=0, j=(int)cyc.size()-1; i < cyc.size(); j = i++ ) arcs.insert( {cyc[j],cyc[i]});

            map<int,int> cycles_of_length;
            for (auto & cyc : cycles) cycles_of_length[cyc.size()]++;

            // gather statistics
            exp_data.iterations.back().full_sol_size = full_sol.size();
            exp_data.iterations.back().hs_size_after_impr = sol.size();
            exp_data.iterations.back().unhit_cycle_enumeration_time_millis = s.getTime("cycles");
            exp_data.iterations.back().hs_greedy_time = s.getTime("hs_greedy");
            exp_data.iterations.back().distinct_arcs_in_all_cycles = arcs.size();
            exp_data.iterations.back().hs_valid_fvs = Utils::isFVS(V,sol);
            exp_data.iterations.back().res_optimal = (response_status == CpSolverStatus::OPTIMAL && sol.size() == full_sol.size());
            exp_data.iterations.back().time_since_start_millis = timer.getTime(timer_option);
            exp_data.iterations.back().total_cycles = cycles.size();
            exp_data.iterations.back().unhit_graph_sizes = getUnhitGraphSizes(V,sol);
            exp_data.iterations.back().cycles_of_length = cycles_of_length;
            exp_data.iterations.back().new_cycles_added = new_cycles_found;
            exp_data.iterations.back().new_cycles_found = new_cycles.size();
            exp_data.iterations.back().iteration_time = s.getTime("iteration");
            if ( !Utils::isFVS(V,sol) && cnf.fill_partial_result_using_greedy_fvs ) {
                exp_data.iterations.back().unhit_graph_greedy_dfvs_size = unhit_graph_dfvs.size();
                exp_data.iterations.back().unhit_graph_greedy_dfvs_time = s.getTime("unhit_graph_greedy_dfvs_time");
            }
        };
        addStats();

        prev_res = sol;
        old_cycles = cycles.size();
        if (full_sol.size() <= best_fvs.size() || best_fvs.empty()) {
            best_fvs = full_sol;
        }

        clog << "\t found sol.size(): " << sol.size() << ", unhit_graph_fvs_size: " << unhit_graph_dfvs.size()
             << ", full_sol.size(): " << full_sol.size() << ", best_fvs.size(): " << best_fvs.size() << endl;

        if ( response_status == CpSolverStatus::OPTIMAL && Utils::isFVS(V,sol) ) break;

        // constexpr int B = 50'000;
        // cnf.ihs_single_iteration_sec = max( cnf.ihs_single_iteration_sec, 1 + (int)cycles.size() / B );
    }

    timer.stop(timer_option);
    timer.write(timer_option);


    if (!best_fvs.empty()) res = best_fvs;
    else res = prev_res;

    return exp_data;
}

ExpData CpsatExp1::solveMTZ(VVI V, ExpConfig cnf, int auxiliary_cycles_mode) {
    clog << "Solving using CpsatExp1::solveMTZ" << endl;
    cnf.writeConfig();

    ExpData exp_data;

    Stopwatch timer;
    string timer_option = "solveMTZ";
    timer.setLimit(timer_option, cnf.find_optimal_result ? inf : cnf.max_time_sec * 1000);
    timer.start(timer_option);

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

    VVI augmented_cycles;

    if (auxiliary_cycles_mode == 1){ // here we add all pi-edges or triangles to make the propagation faster
        int L0 = cnf.init_L_for_all_constraints;
        auto cyc = Utils::getAllSimpleCycles3(V,L0);
        sort(ALL(cyc),[&](auto & c1, auto & c2){ return c1.size() < c2.size(); });
        clog << "\t adding " << cyc.size() << " constraints for all simple cycles of length <= " << L0 << endl;
        addCycleConstraints(model,cyc,nodes);
        augmented_cycles = cyc;
        clog << "\t augmenting MTZ model with " << augmented_cycles.size() << " cycle constraints" << endl;
    }else if (auxiliary_cycles_mode == 2) {
        auto new_cnf = cnf;
        new_cnf.ihs_max_iterations = cnf.ihs_iterations_in_mtz;
        new_cnf.max_time_sec = 0.2 * ( timer.getLimit(timer_option) - timer.getTime(timer_option) ) / 1000;
        if (cnf.find_optimal_result ) cnf.max_time_sec = 30;
        VVI cycles;
        VI res;
        auto r = solveIHS(V,new_cnf,cycles, res);
        addCycleConstraints(model,cycles,nodes);
        init_sol = res;
        augmented_cycles = cycles;
        clog << "\t augmenting MTZ model with " << augmented_cycles.size() << " cycle constraints obtained using IHS" << endl;
    }

    if( auxiliary_cycles_mode != 2 && cnf.ihs_init_sol_creation_mode != 0 ) init_sol = createInitialSolution(V,cnf);

    clog << "In MTZ, init_sol.size(): " << init_sol.size() << endl;

    if(!init_sol.empty()) {
        in_init_sol = StandardUtils::toVB(N,init_sol);
        auto new_cnf = cnf;
        if(cnf.use_init_sol_as_hint_mode == 3) new_cnf.use_init_sol_as_hint_mode = 0;
        addInitialSolutionHint(model,nodes,init_sol,init_sol,new_cnf); // here we intentionally add init_sol instead of prev_res

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

    clog << "\t Starting to solve MTZ model" << endl;

    int time_sec_left = ( timer.getLimit(timer_option) - timer.getTime(timer_option) ) / 1000;
    SatParameters params = getDefaultSatParameters(cnf, time_sec_left);
    Model solver_model;
    solver_model.Add(NewSatParameters(params));

    {
        auto addStats = [&](int objective, int objective_bound, int time_sec) {
            exp_data.iterations.emplace_back();
            exp_data.iterations.back().full_sol_size = objective;
            exp_data.iterations.back().res_optimal = (objective == objective_bound);
            exp_data.iterations.back().time_since_start_millis = time_sec;
            exp_data.iterations.back().res_lower_bound = objective_bound;
        };

        mutex log_mutex;
        solver_model.Add(NewFeasibleSolutionObserver(
            [&](const CpSolverResponse& r) {
                std::lock_guard<std::mutex> lk(log_mutex);
                // std::cout << "obj=" << r.objective_value() << " bound=" << r.best_objective_bound() << " time=" << r.wall_time() << '\n';
                addStats(r.objective_value(), r.best_objective_bound(), (int)r.wall_time() * 1000);
            }
        ));
    }

    model.Minimize(LinearExpr::Sum(nodes));
    CpSolverResponse response = SolveCpModel(model.Build(), &solver_model);

    if ( response.status() == INFEASIBLE ) assert(false && "status cannot be infeasible, unless model is incorrect");

    VI res;
    if ( response.status() == OPTIMAL || response.status() == FEASIBLE ) {
        for (int i=0; i<N; i++) if ( SolutionBooleanValue(response,nodes[i]) ) res.push_back(i);
    }

    {
        map<int,int> cycles_of_length;
        for (auto & cyc : augmented_cycles) cycles_of_length[cyc.size()]++;

        exp_data.iterations.emplace_back();
        exp_data.iterations.back().full_sol_size = res.size();
        // exp_data.iterations.back().hs_size_after_impr = res.size();
        // exp_data.iterations.back().hs_valid_fvs = Utils::isFVS(V,res);
        exp_data.iterations.back().res_optimal = (response.status() == OPTIMAL);
        exp_data.iterations.back().time_since_start_millis = timer.getTime(timer_option);
        exp_data.iterations.back().total_cycles = augmented_cycles.size();
        exp_data.iterations.back().cycles_of_length = cycles_of_length;
    }

    timer.stop(timer_option);
    timer.write(timer_option);

    return exp_data;
}

ExpData CpsatExp1::solveDiVerSeS(VVI V, ExpConfig cnf) {
    ExpData exp_data;

    Stopwatch sw;
    sw.setLimit("diverses", cnf.max_time_sec * 1000);
    sw.start("diverses");


    int best_res_size = 0;
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

        if (best_res_size == 0 || res.size() < best_res_size) best_res_size = res.size();

        exp_data.iterations.emplace_back();
        exp_data.iterations.back().full_sol_size = res.size();
        exp_data.iterations.back().hs_size_after_impr = res.size();
        exp_data.iterations.back().hs_valid_fvs = Utils::isFVS(V,res);
        exp_data.iterations.back().iteration_time = sw2.getTime("iteration");

        iter_id++;
        clog << "After iteration " << iter_id << ", res.size(): " << res.size() << ", best_res_size: " << best_res_size << endl;
    }

    DEBUG(best_res_size);


    return exp_data;
}

VI CpsatExp1::createInitialSolution(VVI &V, ExpConfig cnf) {
    VI init_sol;
    if(cnf.ihs_init_sol_creation_mode == 1) {
        clog << "Looking for initial solution using Agent-Flow approach" << endl;
        init_sol = getUnhitGraphGreedyFVS(V,init_sol);
        clog << "Found initial solution of size " << init_sol.size() << endl;
    }
    if(cnf.ihs_init_sol_creation_mode == 2) {
        Config diverses_cnf;
        diverses_cnf.write_logs = false;
        diverses_cnf.sw.setLimit("main", 1000 * cnf.max_time_sec);
        diverses_cnf.sw.start("main");
        diverses_cnf.disableAllNonbasicReductions();
        diverses_cnf.solverh_use_reductions_AF = false;
        int iters = 1;
        clog << "Looking for initial solution using " << iters << " iteration(s) of DiVerSeS" << endl;
        while(iters--) {
            DFVSSolverH sh(diverses_cnf);
            auto r = sh.solveForGraph(V);
            assert(Utils::isFVS(V,r));
            if( init_sol.empty() || r.size() < init_sol.size() ) init_sol = r;
        }
        clog << "Found initial solution of size " << init_sol.size() << endl;
    }
    return init_sol;
}

ExpData CpsatExp1::solve(VVI V, ExpConfig cnf) {
    auto alg = cnf.alg;
    if (alg == Algorithm::HS) return solveHS(V,cnf);
    if (alg == Algorithm::IHS) return solveIHS(V,cnf);
    if (alg == Algorithm::MTZ) return solveMTZ(V,cnf, cnf.mtz_auxiliary_cycles_mode);
    if (alg == Algorithm::DIVERSES) return solveDiVerSeS(V,cnf);
    if (alg == Algorithm::DIV_IHS) {
        cnf.ihs_init_sol_creation_mode = 2;
        return solveIHS(V,cnf);
    }

    return ExpData{};
}


