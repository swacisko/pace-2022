//
// Created by sylwester on 12/20/21.
//

#ifndef ALGORITHMSPROJECT_DFVSSOLVERH_H
#define ALGORITHMSPROJECT_DFVSSOLVERH_H

#include <graphs/GraphUtils.h>
#include "CONTESTS/PACE22/Utils.h"
#include "CONTESTS/PACE22/Config.h"

class DFVSSolverH{
public:

    DFVSSolverH(Config c) : cnf(c) {}

    VI solveForGraph(VVI V);

    VI solveForGraph2(VVI V);

    VI solveForBiconnectedGraph(VVI V, bool allow_improvements_here = true, const bool use_vc_improvement = false);

    VI solveByAgentFlow(VVI V);

    VI solveByAgentFlowFast(VVI V);

    VI solveByAgentFlowAllWithinDistance(VVI V, const int minimize = 2);

    VI solveByLowerBounding(VVI V);

    VI solveByVCForAllPIGraph(VVI V, int milliseconds);

    VI improve1(VVI V, VI dfvs, int iterations, double alpha);

    template<class _V>
    VD sinkhorn(VVI & V, VVI & revV, _V & values) {
        int N = V.size();

        for( int i=0; i<N; i++ ){
            values[i][i] = 1;
            for(int d : V[i]) values[i][d] = 1;
        }

        auto normalizeRows = [&](){
            for( int r=0; r<N; r++ ){
                double s = values[r][r];
                for( int d : V[r] ) if(d!=r) s += values[r][d];

                values[r][r] /= s;
                for( int d : V[r] ) if(d!=r) values[r][d] /= s;
            }
        };

        auto normalizeColumns = [&](){
            for(int c=0; c<N; c++){
                double s = values[c][c];
                for( int d : revV[c] ) if(d!=c) s += values[d][c];

                values[c][c] /= s;
                for( int d : revV[c] ) if(d!=c) values[d][c] /= s;
            }
        };

        int reps = 10 + log(N);
        while(reps--){
            normalizeRows();
            normalizeColumns();
        }

        VD res(N);
        for(int i=0; i<N; i++) res[i] = values[i][i];
        return res;
    }

    static VD sinkhorn(VVI & V);

    VI solveForStronglyConnectedIterativeHittingSet2(VVI V, VI heur_dfvs, VVI & all_cycles, double F = 1.0);

    Config cnf;

    int origV_edges = -1;
    int origV_N = -1;
};

#endif //ALGORITHMSPROJECT_DFVSSOLVERH_H
