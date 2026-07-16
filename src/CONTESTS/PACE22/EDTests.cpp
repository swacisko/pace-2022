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

struct ExpData {
    int N0=-1, M0=-1, N1=-1, M1=-1, N2=-1, M2=-1, N3=-1, M3=-1, N4=-1, M4=-1;
    int red1_offset, red_noned_offset = 0, red_ed_offset = 0, red_ed2_offset = 0;
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

        stringstream str;
        for (int d : noned_results) str << d << " ";
        res["noned_results"] = str.str();
        str.clear(); str.str("");

        for (int d : ed_results) str << d << " ";
        res["ed_results"] = str.str();
        str.clear(); str.str("");

        for (int d : ed2_results) str << d << " ";
        res["ed2_results"] = str.str();
        str.clear(); str.str("");

        return res;
    }


    void writeToFile(ostream & str, bool debug_entries = false) {
        str.precision(2);

        auto mapa = getEntries();
        vector<string> header = { "N0", "M0", "N1", "M1", "N2", "M2", "N3", "M3", "N4", "M4",
            "red1_offset", "red_noned_offset", "red_ed_offset", "red_ed2_offset",
            "ed_total_t2_inference_rules_created",
            "ed_t2_inference_rules_added",
            "ed_nodes_reduced", "ed2_nodes_reduced",
            "ed_edges_removed", "ed2_edges_removed",
            "ed_t1_inference_rules_added", "ed2_t1_inference_rules_added",
            "red_init_time_millis",
            "red_noned_time_millis",
            "red_ed_time_millis", "red_ed2_time_millis",
            "solver_max_time_sec",
            "solver_time_granularity",
            "solver_repeats",
            "noned_results",
            "ed_results", "ed2_results"
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
        writeLine(line); // write entries

        if (debug_entries) for (auto k : header) clog << k << ": " << mapa[k] << endl;
    }
};


static void updateResTimes(auto & res_times) {
    for ( int i=0; i+1<res_times.size(); i++ ) if ( res_times[i] != -1 ) {
        if ( res_times[i+1] == -1 ) res_times[i+1] = res_times[i];
        else res_times[i+1] = min( res_times[i+1], res_times[i] );
    }
}


pair<VI,VI> solveByFastVC(VPII & constraints, int max_sec, int res_measure_freq_sec ) {
    VI res, res_times(ceil(1.0*max_sec/res_measure_freq_sec)+1, -1);
    return {};
}

pair<VI,VI> solveByHIGHS(VPII & constraints, int max_sec, int res_measure_freq_sec ) {
    VI res, res_times(ceil(1.0*max_sec/res_measure_freq_sec)+1, -1);
    return {};
}

pair<VI,VI> solveByCPSAT(VPII & constraints, int max_sec, int res_measure_freq_sec ) {
    VI res, res_times(ceil(1.0*max_sec/res_measure_freq_sec)+1, -1 );

    clog << "Running cpsat for at most " << max_sec << " seconds" << endl;

    int N = 0;
    for ( auto [a,b] : constraints ) N = max( N, 1 + max( abs(a), abs(b) ) );
    Stopwatch sw;
    sw.start("cpsat");

    SatParameters params;
    params.set_num_search_workers(6);
    params.set_max_time_in_seconds(max_sec);

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
        else cnstr.push_back(nodes[a].Not());

        if (b > 0) cnstr.push_back(nodes[b]);
        else cnstr.push_back(nodes[b].Not());

        model.AddBoolOr(cnstr );
    }

    mutex log_mutex;
    solver_model.Add(NewFeasibleSolutionObserver(
        [&](const CpSolverResponse& r) {
            int cnt = 0;
            log_mutex.lock();
            for (int i = 0; i < nodes.size(); ++i) cnt += SolutionBooleanValue(r, nodes[i]);

            const double t = r.wall_time();              // seconds
            const double dt = r.deterministic_time();    // deterministic solver time

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

pair<VI,VI> solveByEvalMaxSAT(VPII & constraints, int max_sec, int res_measure_freq_sec ) {
    VI res, res_times(ceil(1.0*max_sec/res_measure_freq_sec), -1);
    return {};
}

pair<VI,VD> solveInstanceUsingSolver(VPII constraints, int max_time_sec, int res_measure_freq_sec, int repeats, string alg = "cpsat") {
    // times[i] is the result found by the solver after time (i+1) * granularity seconds
    const int I = ceil(1.0 * max_time_sec / res_measure_freq_sec) + 1;
    VVI iteration_res_times(I, VI());
    VD res_times(I, -1);
    VI res; // valid solution for provided constraints

    for (int rep=0; rep<repeats; rep++) {

        auto solve = [&]() {
            if (alg == "cpsat") return solveByCPSAT(constraints, max_time_sec, res_measure_freq_sec);
            if (alg == "evalmaxsat") return solveByEvalMaxSAT(constraints, max_time_sec, res_measure_freq_sec);
            if (alg == "highs") return solveByHIGHS(constraints, max_time_sec, res_measure_freq_sec);
            if (alg == "numvc" || alg == "fastvc") return solveByFastVC( constraints, max_time_sec, res_measure_freq_sec);
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

static ExpData runVCTestforGraph(VVI V, int solver_max_time_sec, int solver_time_granularity, int solver_repeats, string alg) {
    int N = V.size();
    auto initV = V;

    ExpData exp_data;
    exp_data.solver_time_granularity = solver_time_granularity;
    exp_data.solver_max_time_sec = solver_max_time_sec;
    exp_data.solver_repeats = solver_repeats;

    exp_data.N0 = N;
    exp_data.M0 = GraphUtils::countEdges(V);

    DEBUG(PII(exp_data.N0,exp_data.M0));



    auto writeConnCompInfo = [&](VVI & V, string msg = "") {
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

    { // numvc/fastvc testing
        int t = numvc_time_check_sec;
        numvc_time_check_sec *= 2;
        auto fastvc_sol = checkByNuMVC(V, 0);
        assert(VCUtils::isVertexCover( V, fastvc_sol ));
        numvc_time_check_sec = t;
    }


    constexpr bool test_noned_vc_rules = true;
    if (test_noned_vc_rules){ // measuring just the VC reduction time WITHOUT ED rule, and the solver results for the non-ed reduced graph
        clog << endl << "***************** CHECKING FULL NON-ED RULES" << endl;

        Stopwatch sw;
        sw.start("main");
        Config cnf;
        cnf.disableAllNonbasicReductions();
        cnf.reducer_use_folding = cnf.reducer_use_folding_twins = cnf.reducer_use_funnel = cnf.reducer_use_desk = true;
        cnf.reducer_use_unconfined = true;
        cnf.reducer_use_twins = true;
        cnf.reducer_use_general_folding = true; cnf.reducer_max_general_folding_antiedges = 1; cnf.reducer_max_general_folding_neighborhood_size = 5;
        auto edges = GraphUtils::getGraphEdges(V);
        Reducer red(edges,cnf);
        auto reduced_instance = red.reduce();

        sw.stop("main");
        red.writeTotals();
        clog << "Full NON-ED graph reduction took " << sw.getTime("main") / 1000 << " seconds" << endl;

        exp_data.red_noned_time_millis = exp_data.red_init_time_millis + sw.getTime("main");
        VVI & coreV = reduced_instance.getCoreV();

        exp_data.N2 = coreV.size();
        exp_data.M2 = GraphUtils::countEdges(coreV);
        DEBUG(PII(exp_data.N2,exp_data.M2));

        int noned_solution_lift_overhead = reduced_instance.getReductionsOffset();
        DEBUG(noned_solution_lift_overhead);
        exp_data.red_noned_offset = noned_solution_lift_overhead;

        { // numvc/fastvc testing
            auto fastvc_sol = checkByNuMVC(coreV, exp_data.red_noned_offset);
            assert(VCUtils::isVertexCover( coreV, fastvc_sol ));
            fastvc_sol = reduced_instance.liftSolution(fastvc_sol);
            assert(VCUtils::isVertexCover( V, fastvc_sol ));
        }

        VPII constraints = GraphUtils::getGraphEdges(coreV);
        for (auto & [a,b] : constraints){a++; b++;}
        auto [solver_vc,times] = solveInstanceUsingSolver(constraints, solver_max_time_sec, solver_time_granularity, solver_repeats, alg);
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

    auto testEDRules = [&](auto & exp_data, bool add_constraints = false) {
        constexpr bool test_ed_vc_rules = true;
        if (test_ed_vc_rules) {
            clog << endl << "***************** CHECKING FULL ED RULES" << endl;

            // now measuring VC reduction time WITH ED rule
            Stopwatch sw;
            sw.start("main");
            Config cnf;
            cnf.disableAllNonbasicReductions();
            cnf.reducer_use_folding = cnf.reducer_use_folding_twins = cnf.reducer_use_funnel = cnf.reducer_use_desk = true;
            cnf.reducer_use_unconfined = true;
            cnf.reducer_use_twins = true;
            // cnf.reducer_use_general_folding = true; cnf.reducer_max_general_folding_antiedges = 2; cnf.reducer_max_general_folding_neighborhood_size = 10; // original
            cnf.reducer_use_general_folding = true; cnf.reducer_max_general_folding_antiedges = 1; cnf.reducer_max_general_folding_neighborhood_size = 5;
            cnf.reducer_use_ed = true;
            cnf.ed_consider_nodes_to_move_outside_NW = true;
            cnf.ed_use_same_neigh_domination = true;
            cnf.ed_use_deficit1_domination = true;

            cnf.ed_use_double_ed_checks = false; // time-consuming, especially for denser graphs... use for sparse graphs only
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

            VVI & coreV = reduced_instance.getCoreV();
            exp_data.N3 = coreV.size();
            exp_data.M3 = GraphUtils::countEdges(coreV);
            DEBUG(PII(exp_data.N3,exp_data.M3));
            writeConnCompInfo(coreV, "Connected components after full NON-ED reduction");

            int ed_solution_lift_overhead = reduced_instance.getReductionsOffset();
            DEBUG(ed_solution_lift_overhead);
            exp_data.red_ed_offset = ed_solution_lift_overhead;

            { // numvc/fastvc testing
                auto fastvc_sol = checkByNuMVC(coreV, exp_data.red_ed_offset);
                assert(VCUtils::isVertexCover( coreV, fastvc_sol ));
                fastvc_sol = reduced_instance.liftSolution(fastvc_sol);
                assert(VCUtils::isVertexCover( V, fastvc_sol ));
            }


            // now solve the reduced problem using constraints...


            VPII constraints;
            auto edges = GraphUtils::getGraphEdges(coreV);
            for ( auto [a,b] : edges ) constraints.emplace_back(a+1,b+1);

            if (add_constraints) {

            }


            auto [solver_vc,times] = solveInstanceUsingSolver(constraints,solver_max_time_sec,solver_time_granularity, solver_repeats);
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
    testEDRules(exp_data, false);


    // now checking the impact with adding t1-constraints on the fly
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



    return exp_data;
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


int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout << fixed;
    clog << fixed;


    int N0 = 3'000, M0 = 4'700;
    double C = 3;
    N0 *= C; M0 *= C;

    VVI V = GraphReader::readGraphStandardEdges(cin);
    // VVI V = GraphReader::readGraphDIMACSWunweighed(cin,true);
    // VVI V = getTestV1();
    // VVI V = getRandomGraph(N0, M0);

    V = GraphUtils::makeSimple(V); // making the graph simple, as at the input it might not...

    assert(GraphUtils::isSimple(V));

    // int solver_max_time_sec = 300;
    // int solver_time_granularity = 10;
    int solver_max_time_sec = 10;
    int solver_time_granularity = 1;
    int solver_repeats = 1;
    string alg = "cpsat";

    auto exp_data = runVCTestforGraph(V,  solver_max_time_sec, solver_time_granularity,solver_repeats, alg);


    ENDL(5);
    clog << "FINISHED TESTS!" << endl;


    exp_data.writeToFile(cout, true);


    return 0;
}
