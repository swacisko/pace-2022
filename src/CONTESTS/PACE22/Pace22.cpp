//
// Created by sylwester on 12/20/21.
//

#include <filesystem>
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




pair<string,VI> solveCPSAT1(VVI V, int max_l, bool find_optimal = false) {
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

        for (int d : prev_res) model.AddHint(nodes[d],1);
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
        for( int i=0; i<N; i++ ) model.AddHint(nodes[i], in_init_sol[i]);

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

    Stopwatch sw;


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
    auto[status1,res1] = solveCPSAT1(V,25);
    sw.stop("cpsat-1");

    DEBUG(status1);
    DEBUG(res1.size());
    // DEBUG(res1);

    //******************************

    if ( status1 == "OPTIMAL" && status2 == "OPTIMAL" ) assert( res1.size() == res2.size() );

    ENDL(3);
    DEBUG(res1.size());
    DEBUG(res2.size());
    DEBUG(res3.size());
    DEBUG(res4.size());
    DEBUG(res5.size());

    assert(Utils::isFVS(V,res1));
    assert(Utils::isFVS(V,res2));
    assert(Utils::isFVS(V,res3));
    assert(Utils::isFVS(V,res4));
    assert(Utils::isFVS(V,res5));

    sw.write("cpsat-1");
    sw.write("cpsat-2-N");
    sw.write("cpsat-2-inf");
    sw.write("cpsat-2-N/10");
    sw.write("cpsat-2-only-lns-from-res4");


    return 0;
}

