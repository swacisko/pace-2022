//
// Created by sylwester on 12/20/21.
//

// #include <graphs/GraphWriter.h>
// #include <CONTESTS/PACE22/Reducer.h>
#include "CONTESTS/PACE22/Utils.h"

#include "GraphUtils.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"

using namespace operations_research::sat;

#include "IntGenerator.h"


namespace Utils{


    void writeNeighborhood(VVI & V, int v) {
        DEBUG(v);
        DEBUG(V[v]);
        for (int u : V[v]) clog << "V[" << u << "]: " << V[u] << endl;
    }

    void writeNeighborhood(VVI & V, VVI & revV, int v){
        DEBUG(v);
        DEBUG(V[v]);
        DEBUG(revV[v]);
        ENDL(1);

        set<int> zb(ALL(V[v])); zb.insert(ALL(revV[v]));
        for(int d : zb){
            DEBUG(d);
            DEBUG(V[d]);
            DEBUG(revV[d]);
            ENDL(1);
        }
        ENDL(1);
    }



    bool hasLoop( VVI & V, int a ){
        for( int d : V[a] ) if(d == a) return true;
        return false;
    }

    void writeRemainingGraph(VVI &V) {
        clog << "Graph V:" << endl;
        for(int i=0; i<V.size(); i++){
            if(!V[i].empty()){
                clog << i << " --> ";
                for(int d : V[i]) clog << d << " ";
                clog << endl;
            }
        }
    }




    LL getSetHash(int N, VI &s, int seed) {
        VLL hashes(N);
        IntGenerator rnd(seed);
        for( int i=0; i<N; i++ ) hashes[i] = rnd.rand();
        return accumulate(ALL(s), 0ll, [&](LL h, int b){ return h ^ hashes[b]; });
    }

    VI getMinVcCPSAT(VVI &V, int thread_workers) {
        SatParameters params;
        params.set_num_search_workers(thread_workers);

        CpModelBuilder model;

        int N = V.size();
        vector<BoolVar> nodes;
        for (int i=0; i<N; i++) nodes.push_back(model.NewBoolVar());

        VPII edges = GraphUtils::getGraphEdges(V);

        for ( auto [a,b] : edges ) {
            vector<BoolVar> cnstr; cnstr.reserve(2);
            cnstr.push_back(nodes[a]);
            cnstr.push_back(nodes[b]);
            model.AddBoolOr(cnstr );
        }

        model.Minimize( LinearExpr::Sum(nodes) );

        Model solver_model;
        auto sat_params = NewSatParameters(params);
        solver_model.Add(sat_params);

        const CpSolverResponse response = SolveCpModel(model.Build(), &solver_model);
        assert (response.status() == CpSolverStatus::OPTIMAL);

        VI vc;
        vc.reserve(N);
        for (int i=1; i<N; i++) if (SolutionBooleanValue(response,nodes[i])) vc.push_back(i);

        return vc;
    }



}
