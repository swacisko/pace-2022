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
#include "IntGenerator.h"
#include "components/ConnectedComponents.h"
#include "VertexCover/VCUtils.h"

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"

using namespace operations_research::sat;


enum Alg {
    CPSAT_SAT = 0,
    CPSAT_LP = 1,
    CPSAT_DEF = 2,
    NUMVC = 3,
    HIGHS = 4,
};

string parseAlgorithm(Alg alg) {
    if (alg == CPSAT_SAT) return "cpsat-sat";
    if (alg == CPSAT_LP) return "cpsat-lp";
    if (alg == CPSAT_DEF) return "cpsat-def";
    if (alg == NUMVC) return "numvc";
    if (alg == HIGHS) return "highs";

    assert(false && "incorrect algorithm");
}

struct ExpData {
    int N0=-1, M0=-1, N1=-1, M1=-1, N2=-1, M2=-1, N3=-1, M3=-1, N4=-1, M4=-1;
    int red1_offset = 0, red_noned_offset = 0, red_ed_offset = 0, red_ed2_offset = 0;
    // int type1_constraints = -1;
    // int type2_constraints = -1;

    int red_init_time_millis = -1;
    int red_noned_time_millis = -1;
    int red_ed_time_millis = -1;
    int red_ed2_time_millis = -1;

    int solver_max_time_sec = -1;
    int solver_time_granularity = -1;
    int solver_repeats = -1;
    VD noned_results;
    VD ed_results;
    VD ed2_results;

    int ed_nodes_reduced = -1;
    int ed_edges_removed = -1;
    int ed2_nodes_reduced = -1;
    int ed2_edges_removed = -1;
    int ed_t1_inference_rules_added = -1;
    int ed2_t1_inference_rules_added = -1;
    int ed_total_t2_inference_rules_created = -1;
    int ed2_total_t2_inference_rules_created = -1;
    int ed_t2_inference_rules_added = -1;
    int ed2_t2_inference_rules_added = -1;

    string metadata_filepath = "";

    int reducer_max_time_millis = 900 * 1000; // 15 minutes

    int noned_folds = 0, noned_funnels = 0, noned_unconfined = 0, noned_ext_dom = 0, noned_desks = 0, noned_twins = 0, noned_dominations = 0;
    int ed_folds = 0, ed_funnels = 0, ed_unconfined = 0, ed_ext_dom = 0, ed_desks = 0, ed_twins = 0, ed_dominations = 0;

    Alg alg = CPSAT_SAT;

    bool run_noned = true;
    bool run_ed = true;
    bool run_ed2 = true;

    bool use_def1_dom = true;
    bool ed_use_edge_removal = false;
    bool ed_use_double_ed_checks = false;

    map<string,string> getEntries() {
        map<string,string> res;
        res["N0"] = to_string(N0); res["M0"] = to_string(M0);
        res["N1"] = to_string(N1); res["M1"] = to_string(M1);
        res["N2"] = to_string(N2); res["M2"] = to_string(M2);
        res["N3"] = to_string(N3); res["M3"] = to_string(M3);
        res["N4"] = to_string(N4); res["M4"] = to_string(M4);

        res["red1_offset"] = to_string(red1_offset);
        res["red_noned_offset"] = to_string(red_noned_offset);
        res["red_ed_offset"] = to_string(red_ed_offset);
        res["red_ed2_offset"] = to_string(red_ed2_offset);

        res["ed_nodes_reduced"] = to_string(ed_nodes_reduced);
        res["ed2_nodes_reduced"] = to_string(ed2_nodes_reduced);
        res["ed_edges_removed"] = to_string(ed_edges_removed);
        res["ed2_edges_removed"] = to_string(ed2_edges_removed);
        res["ed_t1_inference_rules_added"] = to_string(ed_t1_inference_rules_added);
        res["ed2_t1_inference_rules_added"] = to_string(ed2_t1_inference_rules_added);
        res["ed_total_t2_inference_rules_created"] = to_string(ed_total_t2_inference_rules_created);
        res["ed2_total_t2_inference_rules_created"] = to_string(ed2_total_t2_inference_rules_created);
        res["ed_t2_inference_rules_added"] = to_string(ed_t2_inference_rules_added);
        res["ed2_t2_inference_rules_added"] = to_string(ed2_t2_inference_rules_added);

        res["red_init_time_millis"] = to_string(red_init_time_millis);
        res["red_noned_time_millis"] = to_string(red_noned_time_millis);
        res["red_ed_time_millis"] = to_string(red_ed_time_millis);
        res["red_ed2_time_millis"] = to_string(red_ed2_time_millis);

        res["solver_max_time_sec"] = to_string(solver_max_time_sec);
        res["solver_time_granularity"] = to_string(solver_time_granularity);
        res["solver_repeats"] = to_string(solver_repeats);
        res["metadata_filepath"] = metadata_filepath;
        res["run_noned"] = to_string(run_noned);
        res["run_ed"] = to_string(run_ed);
        res["run_ed2"] = to_string(run_ed2);
        res["use_def1_dom"] = to_string(use_def1_dom);
        res["ed_use_edge_removal"] = to_string(ed_use_edge_removal);
        res["ed_use_double_ed_checks"] = to_string(ed_use_double_ed_checks);

        res["algorithm"] = parseAlgorithm(alg);

        // ed_folds = 0, ed_funnels = 0, ed_unconfined = 0, ed_ext_dom = 0, ed_desks = 0, ed_twins
        res["ed_folds"] = to_string(ed_folds); res["ed_funnels"] = to_string(ed_funnels);
        res["ed_unconfined"] = to_string(ed_unconfined); res["ed_ext_dom"] = to_string(ed_ext_dom);
        res["ed_desks"] = to_string(ed_desks); res["ed_twins"] = to_string(ed_twins);
        res["ed_dominations"] = to_string(ed_dominations);

        res["noned_folds"] = to_string(noned_folds); res["noned_funnels"] = to_string(noned_funnels);
        res["noned_unconfined"] = to_string(noned_unconfined); res["noned_ext_dom"] = to_string(noned_ext_dom);
        res["noned_desks"] = to_string(noned_desks); res["noned_twins"] = to_string(noned_twins);
        res["noned_dominations"] = to_string(noned_dominations);

        stringstream str;
        for (auto d : noned_results) str << d << " ";
        res["noned_results"] = str.str();
        str.clear(); str.str("");

        for (auto d : ed_results) str << d << " ";
        res["ed_results"] = str.str();
        str.clear(); str.str("");

        for (auto d : ed2_results) str << d << " ";
        res["ed2_results"] = str.str();
        str.clear(); str.str("");

        return res;
    }


    void writeToFile(ostream & str, bool debug_entries = false, const bool debug_only = false) {
        str.precision(2);

        auto mapa = getEntries();
        vector<string> header = { "N0", "M0", "N1", "M1", "N2", "M2", "N3", "M3", "N4", "M4",
            "red1_offset", "red_noned_offset", "red_ed_offset", "red_ed2_offset",
            "ed_nodes_reduced", "ed2_nodes_reduced",
            "ed_edges_removed", "ed2_edges_removed",
            "ed_t1_inference_rules_added", "ed2_t1_inference_rules_added",
            "ed_total_t2_inference_rules_created", "ed2_total_t2_inference_rules_created",
            "ed_t2_inference_rules_added", "ed2_t2_inference_rules_added",
            "red_init_time_millis", "red_noned_time_millis",
            "red_ed_time_millis", "red_ed2_time_millis",
            "solver_max_time_sec", "solver_time_granularity", "solver_repeats", "algorithm",
            "noned_results", "ed_results", "ed2_results",
            "noned_folds", "noned_funnels", "noned_unconfined", "noned_ext_dom", "noned_desks", "noned_twins", "noned_dominations",
            "ed_folds", "ed_funnels", "ed_unconfined", "ed_ext_dom", "ed_desks", "ed_twins", "ed_dominations",
            "metadata_filepath", "run_noned", "run_ed", "run_ed2",
            "use_def1_dom", "ed_use_edge_removal", "ed_use_double_ed_checks"
        };

        auto writeLine = [&](vector<string> & l) {
            int cnt = 0;
            for ( auto s : l ) {
                if (cnt++) str << ",";
                str << s;
            }
            str << "\n";
        };

        writeLine(header); // write header
        vector<string> line;
        for (auto k : header) line.push_back(mapa[k]); // create entries in the order of the header
        if (!debug_only) writeLine(line); // write entries

        if (debug_entries || debug_only) for (auto k : header) clog << k << ": " << mapa[k] << endl;
    }
};


static void updateResTimes(auto & res_times) {
    for ( int i=0; i+1<res_times.size(); i++ ) if ( res_times[i] != -1 ) {
        if ( res_times[i+1] == -1 ) res_times[i+1] = res_times[i];
        else res_times[i+1] = min( res_times[i+1], res_times[i] );
    }
}


pair<VI,VI> solveByFastVC(VPII & constraints, ExpData & exp_data ) {
    int max_sec = exp_data.solver_max_time_sec;
    int res_measure_freq_sec = exp_data.solver_time_granularity;

    clog << "Looking for solution using NuMVC/FastVC for " << exp_data.solver_max_time_sec << " seconds" << endl;

    VVI V = GraphUtils::getGraphForEdges(constraints);
    auto vc = VCUtils::getMinCVUsingFastVC(V, exp_data.solver_max_time_sec * 1000);
    for (int & d : vc) d--;

    VI res = vc;
    VI res_times(ceil(1.0*max_sec/res_measure_freq_sec)+1, vc.size());

    return {res,res_times};
}

pair<VI,VI> solveByHIGHS(VPII & constraints,ExpData & exp_data ) {
    int max_sec = exp_data.solver_max_time_sec;
    int res_measure_freq_sec = exp_data.solver_time_granularity;

    clog << "Looking for solution using HiGHS for " << exp_data.solver_max_time_sec << " seconds" << endl;

    VI res, res_times(ceil(1.0*max_sec/res_measure_freq_sec)+1, -1);
    return {};
}

// pair<VI,VI> solveByCPSAT(VPII & constraints, int max_sec, int res_measure_freq_sec ) {
pair<VI,VI> solveByCPSAT(VPII & constraints, ExpData &exp_data ) {
    int max_sec = exp_data.solver_max_time_sec;
    int res_measure_freq_sec = exp_data.solver_time_granularity;

    VI res, res_times(ceil(1.0*max_sec/res_measure_freq_sec)+1, -1 );

    clog << "Looking for solution using CPSAT (" << parseAlgorithm(exp_data.alg) << ") for "
         << exp_data.solver_max_time_sec << " seconds" << endl;

    int N = 0;
    for ( auto [a,b] : constraints ) N = max( N, 1 + max( abs(a), abs(b) ) );
    Stopwatch sw;
    sw.start("cpsat");

    SatParameters params;
    params.set_max_time_in_seconds(max_sec);

    bool use_sat_heavy_computation = (exp_data.alg == CPSAT_SAT);
    if (use_sat_heavy_computation) {
        clog << "Running CPSAT with SAT orientation" << endl;
        constexpr int workers = 2;
        params.set_num_workers(workers);  // Preferred over deprecated num_search_workers.

        params.clear_subsolvers();
        params.set_num_full_subsolvers(workers);

        params.add_subsolvers("quick_restart_no_lp");
        params.add_subsolvers("no_lp");

        // Explicitly disable incomplete primal heuristics.
        params.set_use_lns(false);
        params.set_use_rins_lns(false);
        params.set_use_feasibility_pump(false);
        params.set_use_feasibility_jump(false);
        params.set_num_violation_ls(0);

        // params.set_log_search_progress(true);
        // params.set_log_subsolver_statistics(true);
    }

    bool use_lp_heavy_compuation =  (exp_data.alg == CPSAT_LP);
    if (use_lp_heavy_compuation) {
        clog << "Running CPSAT with LP orientation" << endl;
        constexpr int workers = 2;
        params.set_num_workers(workers);

        params.clear_subsolvers();
        params.set_num_full_subsolvers(workers);

        // Make both workers construct the LP constraints eagerly.
        params.set_add_lp_constraints_lazily(false);

        params.add_subsolvers("quick_restart_max_lp");
        params.add_subsolvers("max_lp");

        // Disable incomplete primal heuristics.
        params.set_use_lns(false);
        params.set_use_rins_lns(false);
        params.set_use_feasibility_pump(false);
        params.set_use_feasibility_jump(false);
        params.set_num_violation_ls(0);

        // params.set_log_search_progress(true);
        // params.set_log_subsolver_statistics(true);
    }

    bool use_default_computation = (exp_data.alg == CPSAT_DEF);
    if (use_default_computation) {
        clog << "Running CPSAT with DEFAULT orientation" << endl;
        constexpr int workers = 2;
        params.set_num_workers(workers);
    }



    Model solver_model;
    CpModelBuilder model;
    auto sat_params = NewSatParameters(params);
    solver_model.Add(sat_params);

    vector<BoolVar> nodes;
    for (int i=0; i<N+1; i++) nodes.push_back(model.NewBoolVar());
    model.AddEquality(nodes[0], 0);

    for ( auto [a,b] : constraints ) {
        vector<BoolVar> cnstr; cnstr.reserve(2);

        if (a > 0) cnstr.push_back(nodes[a]);
        else cnstr.push_back(nodes[-a].Not());

        if (b > 0) cnstr.push_back(nodes[b]);
        else cnstr.push_back(nodes[-b].Not());

        model.AddBoolOr(cnstr );
    }

    mutex log_mutex;
    solver_model.Add(NewFeasibleSolutionObserver(
        [&](const CpSolverResponse& r) {
            int cnt = 0;
            log_mutex.lock();
            for (int i = 0; i < nodes.size(); ++i) cnt += SolutionBooleanValue(r, nodes[i]);

            const double t = r.wall_time();              // seconds
            // clog << "Found new results of size " << cnt << " at time " << t << " seconds " << endl;

            int ind = ceil(1.0 * t / res_measure_freq_sec);
            if ( ind < res_times.size() ) res_times[ind] = cnt;
            log_mutex.unlock();
        }
    ));

    model.Minimize( LinearExpr::Sum(nodes) );

    const CpSolverResponse response = SolveCpModel(model.Build(), &solver_model);
    updateResTimes(res_times);
    sw.stop("cpsat");

    if (response.status() == CpSolverStatus::OPTIMAL ||
        response.status() == CpSolverStatus::FEASIBLE) {
        for (int i=1; i<N; i++) if (SolutionBooleanValue(response,nodes[i])) res.push_back(i-1);
        if (response.status() == CpSolverStatus::OPTIMAL) clog << "Found OPTIMAL result of size " << res.size();
        if (response.status() == CpSolverStatus::FEASIBLE) clog << "Found feasible result of size " << res.size();
        clog << " in time " << response.wall_time() << " seconds" << endl;
    }
    else {
        if ( response.status() == CpSolverStatus::INFEASIBLE ) clog << "CPSAT solver status INFEASIBLE" << endl;
        if ( response.status() == CpSolverStatus::MODEL_INVALID ) clog << "CPSAT solver status MODEL_INVALID" << endl;
        if ( response.status() == CpSolverStatus::UNKNOWN ) clog << "CPSAT solver status UNKNOWN" << endl;
    }

    int ind = ceil(1.0 * sw.getTime("cpsat") / (1000*res_measure_freq_sec));
    if ( ind < res_times.size() ) res_times[ind] = res.size();

    clog << "Returning results found by cpsat, res.size() " << res.size() << ", res_times: " << res_times << endl;

    return {res,res_times};
}


pair<VI,VD> solveInstanceUsingSolver(VPII constraints, ExpData & exp_data) {

    int max_time_sec = exp_data.solver_max_time_sec;
    int res_measure_freq_sec = exp_data.solver_time_granularity;
    Alg alg = exp_data.alg;
    int repeats = exp_data.solver_repeats;

    // times[i] is the result found by the solver after time (i+1) * granularity seconds
    const int I = ceil(1.0 * max_time_sec / res_measure_freq_sec) + 1;
    VVI iteration_res_times(I, VI());
    VD res_times(I, -1);
    VI res; // valid solution for provided constraints

    for (int rep=0; rep<repeats; rep++) {

        auto solve = [&]() {
            if (alg == CPSAT_SAT || alg == CPSAT_LP || alg == CPSAT_DEF) return solveByCPSAT(constraints, exp_data);
            // if (alg == "evalmaxsat") return solveByEvalMaxSAT(constraints, max_time_sec, res_measure_freq_sec);
            if (alg == HIGHS) return solveByHIGHS(constraints, exp_data);
            if (alg == NUMVC) return solveByFastVC( constraints, exp_data);
            return pair<VI,VI>{};
        };

        auto [iter_res, iter_res_times] = solve();
        if (res.empty() || iter_res.size() < res.size()) res = iter_res;
        assert(iter_res_times.size() == I);

        for ( int i=0; i<I; i++ ) {
            assert(i < iter_res_times.size());
            assert(i < iteration_res_times.size());
            iteration_res_times[i].push_back(iter_res_times[i]);
        }
    }

    for (int i=0; i<iteration_res_times.size(); i++ ) {
        const auto & irt = iteration_res_times[i];
        assert(!irt.empty());
        if ( !irt.empty() ) res_times[i] = accumulate(ALL(irt),0.0) / irt.size();

    }

    return make_pair(res,res_times);
}

static void runVCTestforGraph(VVI V, ExpData & exp_data) {
    int N = V.size();
    auto initV = V;

    auto exp_data_cnf = exp_data;

    exp_data.N0 = N;
    exp_data.M0 = GraphUtils::countEdges(V);

    DEBUG(PII(exp_data.N0,exp_data.M0));



    auto writeConnCompInfo = [&](VVI & V, string msg = "") {
        return;

        if (msg != "") clog << msg << endl;
        auto [cmp_cnt, sizes] = GraphUtils::getConnectedcomponentsInfo(V);
        clog << "There are " << cmp_cnt << " nontrivial connected components, with the following size distribution" << endl;
        int cnt = 0;
        for (auto [k,v] : sizes) {
            clog << k << ": " << v << "  |  ";
            if (++cnt % 10 == 0) clog << "\n";
        }
        ENDL(1);
    };


    constexpr bool run_ed_test = false;
    if (run_ed_test){
        Config cnf;
        cnf.disableAllNonbasicReductions();
        cnf.reducer_use_ed = true;
        EDReducer edred(N,cnf);
        auto res = edred.reduce(V);
        exit(4);
    }


    bool run_initial_vc_kernel = true;
    if (run_initial_vc_kernel){
        Stopwatch sw;
        sw.start("main");
        KernelizerVC kern;
        kern.use_crown_and_lp_checks = false;
        auto [kern_nodes, edges_removed] = kern.initialKernelization(V);
        VB helper(N);

        auto revV = GraphUtils::reverseGraph(V);
        GraphUtils::removeNodes(V, kern_nodes, helper);

        V = GraphInducer::induceByNonisolatedNodes(V).V;
        N = V.size();
        sw.stop("main");

        exp_data.red_init_time_millis = sw.getTime("main");

        exp_data.N1 = N;
        exp_data.M1 = GraphUtils::countEdges(V);
        exp_data.red1_offset = kern_nodes.size();

        writeConnCompInfo(V, "Connected components after init-kernelization");

        DEBUG(exp_data.red_init_time_millis);
    }


    DEBUG(PII(exp_data.N1,exp_data.M1));
    DEBUG(exp_data.red1_offset);
    ENDL(3);

    double numvc_time_check_sec = 4;

    auto checkByNuMVC = [&](VVI & V, int additional_offset = 0) {
        if ( GraphUtils::countEdges(V) == 0 ) return VI{};

        double T = numvc_time_check_sec;
        clog << endl << "********************** NUMVC CHECK" << endl;
        clog << "Running NuMVC/FastVC for " << T << " seconds to check size" << endl;
        auto vc = VCUtils::getMinCVUsingFastVC(V, T * 1000);
        DEBUG(vc.size());
        clog << "Solution size with additional offset: " << vc.size() + additional_offset << endl;
        clog << "********************** NUMVC CHECK" << endl << endl;

        return vc;
    };

    constexpr bool use_numvc_for_testing = false;

    if constexpr(use_numvc_for_testing){ // numvc/fastvc testing
        int t = numvc_time_check_sec;
        numvc_time_check_sec *= 1;
        auto fastvc_sol = checkByNuMVC(V, 0);
        assert(VCUtils::isVertexCover( V, fastvc_sol ));
        numvc_time_check_sec = t;
    }


    const bool test_noned_vc_rules = exp_data.run_noned;
    if (test_noned_vc_rules){ // measuring just the VC reduction time WITHOUT ED rule, and the solver results for the non-ed reduced graph
        clog << endl << "***************** CHECKING FULL NON-ED RULES" << endl;

        Stopwatch sw;
        sw.start("main");
        Config cnf;
        cnf.disableAllNonbasicReductions();
        cnf.reducer_use_folding = cnf.reducer_use_funnel = cnf.reducer_use_desk = true;
        cnf.reducer_use_unconfined = true;
        cnf.reducer_use_twins = true;
        cnf.reducer_use_general_folding = true; cnf.reducer_max_general_folding_antiedges = 1; cnf.reducer_max_general_folding_neighborhood_size = 5;
        auto edges = GraphUtils::getGraphEdges(V);
        cnf.reducer_max_time_millis = exp_data_cnf.reducer_max_time_millis;

        Reducer red(edges,cnf);
        auto reduced_instance = red.reduce();

        sw.stop("main");
        red.writeTotals();
        clog << "Full NON-ED graph reduction took " << sw.getTime("main") / 1000 << " seconds" << endl;

        exp_data.red_noned_time_millis = exp_data.red_init_time_millis + sw.getTime("main");
        VVI & coreV = reduced_instance.getCoreV();

        exp_data.noned_dominations = red.total_dominations_done;
        exp_data.noned_desks = red.total_desks_done;
        exp_data.noned_folds = red.total_folds_done;
        exp_data.noned_twins = red.total_twins_done;
        exp_data.noned_ext_dom = red.ed_nodes_reduced;
        exp_data.noned_funnels = red.total_funnels_done;
        exp_data.noned_unconfined = red.total_unconfined_nodes;

        exp_data.N2 = coreV.size();
        exp_data.M2 = GraphUtils::countEdges(coreV);
        DEBUG(PII(exp_data.N2,exp_data.M2));

        int noned_solution_lift_overhead = reduced_instance.getReductionsOffset();
        DEBUG(noned_solution_lift_overhead);
        exp_data.red_noned_offset = noned_solution_lift_overhead;

        if constexpr(use_numvc_for_testing){ // numvc/fastvc testing
            auto fastvc_sol = checkByNuMVC(coreV, exp_data.red_noned_offset);
            assert(VCUtils::isVertexCover( coreV, fastvc_sol ));
            fastvc_sol = reduced_instance.liftSolution(fastvc_sol);
            assert(VCUtils::isVertexCover( V, fastvc_sol ));
        }

        VPII constraints = GraphUtils::getGraphEdges(coreV);
        for (auto & [a,b] : constraints){a++; b++;}
        auto [solver_vc,times] = solveInstanceUsingSolver(constraints, exp_data_cnf);
        DEBUG(times);
        for (auto& d : times) if (d != -1) d += noned_solution_lift_overhead;
        exp_data.noned_results = times;
        DEBUG(exp_data.noned_results);
        DEBUG(solver_vc.size());

        if ( !solver_vc.empty() || coreV.empty() ) {
            solver_vc = reduced_instance.liftSolution(solver_vc);
            assert(VCUtils::isVertexCover( V, solver_vc ));
        }

        clog << "After lifting, solver_vc.size(): " << solver_vc.size() << endl;

        writeConnCompInfo(coreV, "Connected components after full NON-ED reduction");
    }

    auto testEDRules = [&](ExpData & exp_data, bool add_constraints = false) {
        constexpr bool test_ed_vc_rules = true;

        if (test_ed_vc_rules) {
            clog << endl << "***************** CHECKING FULL ED RULES" << endl;

            // now measuring VC reduction time WITH ED rule
            Stopwatch sw;
            sw.start("main");
            Config cnf;
            cnf.disableAllNonbasicReductions();
            cnf.reducer_use_folding = cnf.reducer_use_funnel = cnf.reducer_use_desk = true;
            cnf.reducer_use_unconfined = true;
            cnf.reducer_use_twins = true;
            // cnf.reducer_use_general_folding = true; cnf.reducer_max_general_folding_antiedges = 2; cnf.reducer_max_general_folding_neighborhood_size = 10; // original
            cnf.reducer_use_general_folding = true; cnf.reducer_max_general_folding_antiedges = 1; cnf.reducer_max_general_folding_neighborhood_size = 5;
            cnf.reducer_use_ed = true;
            cnf.ed_consider_nodes_to_move_outside_NW = true;
            cnf.ed_use_same_neigh_domination = true;
            cnf.ed_use_deficit1_domination = exp_data.use_def1_dom;
            cnf.ed_use_edge_removal = exp_data.ed_use_edge_removal;

            cnf.ed_use_double_ed_checks = exp_data.ed_use_double_ed_checks; // time-consuming, especially for denser graphs... use for sparse graphs only
            // cnf.ed_use_edge_removal = true;
            cnf.ed_apply_type1_constraints_on_the_fly = add_constraints;

            Reducer red(GraphUtils::getGraphEdges(V),cnf);
            auto reduced_instance = red.reduce();
            sw.stop("main");
            red.writeTotals();
            clog << "Full ED graph reduction took " << sw.getTime("main") / 1000 << " seconds" << endl;

            exp_data.red_ed_time_millis = exp_data.red_init_time_millis + sw.getTime("main");

            exp_data.ed_nodes_reduced = red.ed_nodes_reduced;
            exp_data.ed_edges_removed = red.ed_edges_removed; assert(red.ed_edges_removed == 0);
            exp_data.ed_t1_inference_rules_added = red.ed_t1_inference_rules_added;
            exp_data.ed_t2_inference_rules_added = red.ed_t2_inference_rules_added;
            exp_data.ed_total_t2_inference_rules_created = red.ed_total_t2_inference_rules_created;

            exp_data.ed_dominations = red.total_dominations_done;
            exp_data.ed_desks = red.total_desks_done;
            exp_data.ed_folds = red.total_folds_done;
            exp_data.ed_twins = red.total_twins_done;
            exp_data.ed_ext_dom = red.ed_nodes_reduced;
            exp_data.ed_funnels = red.total_funnels_done;
            exp_data.ed_unconfined = red.total_unconfined_nodes;

            VVI & coreV = reduced_instance.getCoreV();
            exp_data.N3 = coreV.size();
            exp_data.M3 = GraphUtils::countEdges(coreV);
            DEBUG(PII(exp_data.N3,exp_data.M3));
            writeConnCompInfo(coreV, "Connected components after full NON-ED reduction");

            int ed_solution_lift_overhead = reduced_instance.getReductionsOffset();
            DEBUG(ed_solution_lift_overhead);
            exp_data.red_ed_offset = ed_solution_lift_overhead;

            if constexpr(use_numvc_for_testing){ // numvc/fastvc testing
                auto fastvc_sol = checkByNuMVC(coreV, exp_data.red_ed_offset);
                assert(VCUtils::isVertexCover( coreV, fastvc_sol ));
                fastvc_sol = reduced_instance.liftSolution(fastvc_sol);
                assert(VCUtils::isVertexCover( V, fastvc_sol ));
            }


            // now solve the reduced problem using constraints...


            VPII constraints;
            auto edges = GraphUtils::getGraphEdges(coreV);
            for ( auto [a,b] : edges ) constraints.emplace_back(a+1,b+1);

            bool add_t2_constraints = (add_constraints && exp_data.alg != NUMVC); // for numvc we do not want to add t2-constraints
            // bool add_t2_constraints = true;
            if (add_t2_constraints) {

                clog << "Preparing to create type-2 constraints" << endl;
                EDReducer edred(coreV.size(),cnf);
                edred.resetAllUsedTechniques();
                edred.cnf.ed_use_node_removal = true;
                edred.cnf.gather_t2_inf_rules = true;
                auto removed_nodes = edred.reduce(coreV);

                // assert(removed_nodes.empty() && " this should hold if the preprocessing finished and did not terminate due to timeout");
                DEBUG(removed_nodes.size());
                if (!removed_nodes.empty()) clog << "Found reducible nodes... why?!" << endl;

                if (removed_nodes.empty()) {
                    clog << "Creating and adding type-2 constraints" << endl;
                    auto rules = edred.all_inf_rules_2_found;
                    clog << "\t found " << rules.size() << " rules" << endl;

                    sort(ALL(rules),[&](auto & a, auto & b){ return a.second.size() > b.second.size(); });

                    VB was(N);
                    // int cnstr_added_test = 0;
                    for ( auto & [v,r] : rules ) {
                        // DEBUG(v); DEBUG(r); ENDL(1);

                        // check if the neighborhood is blocked
                        bool blocked = was[v];
                        for (int d : V[v]) blocked |= was[d];
                        if (blocked) continue;
                        for (int d : r) blocked |= was[d];
                        if (blocked) continue;
                        for (int d : r) for (int dd : V[d]) blocked |= was[dd];
                        if (blocked) continue;

                        // mark neighborhoods as blocked
                        was[v] = true;
                        for (int d : V[v]) was[d] = true;
                        for (int d : r) was[d] = true;
                        for (int d : r) for (int dd : V[d]) was[dd] = true;

                        for (int d : r) {
                            // clog << "\tadding t2-inf-rule, v: " << v << ", d: " << d << endl;
                            constraints.emplace_back( -(d+1), v+1 ); // if d belongs, then v also belongs
                            exp_data.ed_t2_inference_rules_added++;
                            // break;
                        }

                        // if (++cnstr_added_test > 0) break;
                    }

                    clog << "Created and added " << exp_data.ed_t2_inference_rules_added << " type-2 constraints" << endl;
                }
            }


            auto [solver_vc,times] = solveInstanceUsingSolver(constraints, exp_data_cnf);
            DEBUG(times);
            for (auto& d : times) if (d != -1) d += ed_solution_lift_overhead;
            exp_data.ed_results = times;

            DEBUG(solver_vc.size());
            DEBUG(exp_data.ed_results);

            if (!solver_vc.empty() || coreV.empty()) {
                solver_vc = reduced_instance.liftSolution(solver_vc);
                assert(VCUtils::isVertexCover( V, solver_vc ));
            }
            clog << "After lifting, solver_vc.size(): " << solver_vc.size() << endl;
        }
    };

    // now checking the impact without adding t1-constraints on the fly
    if (exp_data.run_ed) testEDRules(exp_data, false);


    if (exp_data.run_ed2) {
        // now checking the impact with adding t1-constraints on the fly and t2 constraints on the used solver
        ExpData dummy_exp_data = exp_data;
        testEDRules( dummy_exp_data, true);
        exp_data.N4 = dummy_exp_data.N3;
        exp_data.M4 = dummy_exp_data.M3;
        exp_data.red_ed2_offset = dummy_exp_data.red_ed_offset;
        exp_data.ed2_results = dummy_exp_data.ed_results;
        exp_data.red_ed2_time_millis = dummy_exp_data.red_ed_time_millis;
        exp_data.ed2_t1_inference_rules_added = dummy_exp_data.ed_t1_inference_rules_added;
        exp_data.ed2_t2_inference_rules_added = dummy_exp_data.ed_t2_inference_rules_added;
        exp_data.ed2_total_t2_inference_rules_created = dummy_exp_data.ed_total_t2_inference_rules_created;
        exp_data.ed2_nodes_reduced = dummy_exp_data.ed_nodes_reduced;
        exp_data.ed2_edges_removed = dummy_exp_data.ed_edges_removed;
    }

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
            {a,p}, {a,q}, {a,b},
            {b,q}, {b,x},
            {c,x}, {c,y}, {c,d},
            {d,r}, {d,s},
            {p,x}, {p,l},
            {q,x}, {q,l},
            {x,k}, {x,y},
            {k,l}, {k,y},
            {y,r},
            {r,l}, {r,s},
            {s,l}
        };
    }

    return GraphUtils::getGraphForEdges(edges,false);
}


VVI getRandomGraph( int N, int M ) {
    set<PII> zb;
    IntGenerator rnd;
    while (zb.size() < M) {
        int a = rnd.nextInt(N);
        int b = rnd.nextInt(N);
        if (a > b) swap(a,b);
        if (a != b) zb.insert(PII(a,b));
    }

    VPII edges(ALL(zb));
    VVI V = GraphUtils::getGraphForEdges(edges,false);
    return V;
}

void trimAllDegreeOneNodes(VVI & V) {
    int N = V.size();
    VI deg(N,0);
    for (int i=0; i<N; i++) deg[i] = V[i].size();

    deque<int> q;
    for (int i=0; i<N; i++) if (deg[i] == 1) {
        int b = V[i][0];
        deg[b]--;
        if (deg[b] == 1) q.push_back(b);
        GraphUtils::removeNodeFromGraph(V,i);
        deg[i] = 0;
        assert(deg[i] == V[i].size());

        while (!q.empty()) {
            int a = q.back();
            q.pop_back();
            if ( deg[a] != 1 ) {
                if (deg[a] != 0) DEBUG(PII(a,deg[a]));
                assert(deg[a] == 0);
                continue;
            }
            assert(deg[a] == 1);

            int b = V[a][0];
            deg[b]--;
            if (deg[b] == 1) q.push_back(b);
            GraphUtils::removeNodeFromGraph(V,a);
            deg[a] = 0;

            if (deg[b] != V[b].size()) {
                DEBUG(b);
                DEBUG(PII(deg[b], V[b].size()));
                assert(deg[b] == V[b].size());
            }
        }
    }
}


ExpData parseArguments(int argc, char ** argv) {
    ExpData cnf{};

    class ArgParser {
    public:
        // Store flags: --verbose, --help
        unordered_map<string, bool> flags;

        // Store options with values: --input=..., --threads=...
        unordered_map<string, string> options;

        // Which names are flags/options
        unordered_set<string> flag_names;
        unordered_set<string> option_names;
        unordered_set<string> required_options;

        void addFlag(const string &name) {
            flag_names.insert(name);
            flags[name] = false;
        }

        void addOption(const string &name, bool required) {
            option_names.insert(name);
            options[name] = "";
            if (required) required_options.insert(name);
        }

        using VARIANT = variant<bool*,int*,double*,string*>;
        void findAndAssign(string name, string type, VARIANT data) {
            if (!hasProvidedOption(name)) return;


            if (type == "bool") {
                bool* ptr = get<bool*>(data);
                auto isTrue = [&](string s) { return s == "True" || s == "true" || s == "1"; };
                *ptr = isTrue(getOption(name));
            }
            else if (type == "int") {
                int* ptr = get<int*>(data);
                *ptr = stoi(getOption(name));
            }else if (type == "double") {
                double* ptr = get<double*>(data);
                *ptr = stod(getOption(name));
            }else if (type == "string") {
                string* ptr = get<string*>(data);
                *ptr = getOption(name);
            }
        }

        void parse(int argc, char **argv) {
            for (int i = 1; i < argc; i++) {
                string arg = argv[i];

                if (!startsWithDoubleDash(arg)) throw runtime_error("Unknown positional or malformed argument: " + arg);

                string inner = arg.substr(2); // strip "--"
                size_t eq = inner.find('='); // Split on '='
                string name, value;

                if (eq == string::npos) { // No '=' → must be a flag (e.g., --verbose)
                    name = inner;

                    if (flag_names.contains(name)) {
                        flags[name] = true;
                    } else if (option_names.contains(name)) {
                        throw runtime_error("Missing '=value' for option --" + name + " (expected --" + name + "=VALUE)");
                    } else {
                        throw runtime_error("Unknown argument: --" + name);
                    }
                } else {
                    // Has '=' → must be an option: --name=value
                    name = inner.substr(0, eq);
                    value = inner.substr(eq + 1);

                    if (flag_names.contains(name)) {
                        throw runtime_error("Flag --" + name + " does not take a value (remove '=...').");
                    } else if (option_names.contains(name)) {
                        if (value.empty()) {
                            throw runtime_error("Missing value for option --" + name + " (use --" + name + "=VALUE)");
                        }
                        options[name] = value;
                    } else {
                        throw runtime_error("Unknown argument: --" + name);
                    }
                }
            }
        }

        bool getFlag(const string &name) const {
            auto it = flags.find(name);
            if (it == flags.end()) throw runtime_error("Flag not registered: " + name);
            return it->second;
        }

        bool hasProvidedOption(const string &name) const { return options.find(name)->second != ""; }

        string getOption(const string &name) const {
            auto it = options.find(name);
            if (it == options.end()) throw runtime_error("Option not registered: " + name);
            return it->second;
        }

        void printHelp(const string &progName) const {
            cout << "Usage: " << progName << " [options]\n\n";
            cout << "Options:\n";
            for (auto &f : flag_names) cout << "  --" << f << "\n";
            for (auto &o : option_names) cout << "  --" << o << "=<value>\n";
            cout << "\n";
        }

    private:
        static bool startsWithDoubleDash(const string &s) {
            return s.size() >= 2 && s[0] == '-' && s[1] == '-';
        }
    };


    ArgParser ap;
    ap.addOption("alg", false);
    ap.addOption("mtd", true);
    ap.addOption("time", false); // max time for solvers in seconds
    ap.addOption("rep", false); // max time for solvers in seconds
    ap.addOption("gran", false); // granularity
    ap.addOption("run_noned", false); // noned tests
    ap.addOption("run_ed", false); // ed tests
    ap.addOption("run_ed2", false); // ed2 tests
    ap.addOption("use_def1_dom", false); // ed2 tests
    ap.addOption("ed_use_edge_removal", false); // ed2 tests
    ap.addOption("ed_use_double_ed_checks", false); // ed2 tests

    ap.parse(argc, argv);
    for ( const string& opt : ap.required_options ) if( !ap.hasProvidedOption(opt) ) {
        clog << "Option " << opt << " is not provided, but is mandatory!" << endl;
    }
    for ( const string& opt : ap.required_options ) assert( ap.hasProvidedOption(opt) );

    string alg;
    ap.findAndAssign("alg", "string", &alg);
    std::transform(alg.begin(), alg.end(), alg.begin(), [](unsigned char c){ return std::tolower(c); });
    DEBUG(alg);
    if ( alg == "cpsat-sat" ) cnf.alg = CPSAT_SAT;
    if ( alg == "cpsat-def" ) cnf.alg = CPSAT_DEF;
    if ( alg == "cpsat-lp" ) cnf.alg = CPSAT_LP;
    if ( alg == "highs" ) cnf.alg = HIGHS;
    if ( alg == "numvc" ) cnf.alg = NUMVC;

    ap.findAndAssign("mtd", "string", &cnf.metadata_filepath);
    ap.findAndAssign("time", "int", &cnf.solver_max_time_sec);
    ap.findAndAssign("gran", "int", &cnf.solver_time_granularity);
    ap.findAndAssign("rep", "int", &cnf.solver_repeats);
    ap.findAndAssign("run_noned", "bool", &cnf.run_noned);
    ap.findAndAssign("run_ed", "bool", &cnf.run_ed);
    ap.findAndAssign("run_ed2", "bool", &cnf.run_ed2);
    ap.findAndAssign("use_def1_dom", "bool", &cnf.use_def1_dom);
    ap.findAndAssign("ed_use_edge_removal", "bool", &cnf.ed_use_edge_removal);
    ap.findAndAssign("ed_use_double_ed_checks", "bool", &cnf.ed_use_double_ed_checks);


    return cnf;
}

void trimDeg1Nodes(VVI & V) {
    V = GraphInducer::induceByNonisolatedNodes(V).V;
    int N = V.size(), M = GraphUtils::countEdges(V);
    clog << "Before trimming degree 1 nodes, N:" << N << ", M: " << M << endl;

    trimAllDegreeOneNodes(V);
    V = GraphInducer::induceByNonisolatedNodes(V).V;
    N = V.size(), M = GraphUtils::countEdges(V);
    clog << "After trimming degree 1 nodes, N:" << N << ", M: " << M << endl;

    int c = 0;
    for (int i=0; i<N; i++) c += (V[i].size() == 1);
    assert(c == 0);
}

int main(int argc, char** argv) {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout << fixed;
    clog << fixed;



    // int N0 = 3'000, M0 = 4'700;
    // double C = 3;
    // N0 *= C; M0 *= C;
    // VVI V = getRandomGraph(N0, M0);

    VVI V = GraphReader::readGraphStandardEdges(cin);
    // VVI V = GraphReader::readGraphDIMACSWunweighed(cin,true);
    // VVI V = getTestV1();

    V = GraphUtils::makeSimple(V); // making the graph simple, as at the input it might not...

    constexpr bool trim_deg1_nodes = false; // can be used to obtain a slightly more difficult instance
    if (trim_deg1_nodes) trimDeg1Nodes(V);

    assert(GraphUtils::isSimple(V));


    ExpData exp_data = parseArguments(argc, argv);
    exp_data.writeToFile(clog, true, true); // just log the parameters from the input

    runVCTestforGraph(V,  exp_data);

    ENDL(5);
    clog << "FINISHED TESTS!" << endl;

    ofstream mtd_str(exp_data.metadata_filepath);
    exp_data.writeToFile(mtd_str, true);

    return 0;
}
