//
// Created by sylwester on 12/30/20.
//

#include <graphs/GraphUtils.h>
#include "graphs/cycles/CycleCounter.h"

namespace CycleCounter{

    VVI getAllTrianglesInDirectedGraph( VVI V ){
        int N = (int)V.size();

        VVI revV = GraphUtils::reverseGraph(V);

        VI order(N);
        iota(ALL(order),0);
        sort(ALL(order), [&]( int a, int b ){ return V[a].size() > V[b].size(); } );
        VI inOrder(N);
        for(int i=0; i<N; i++) inOrder[order[i]] = i;

        VVI triangles;
        triangles.reserve(3*N);
        VB was(N,false);

        for(int a : order){
            for( int d : revV[a] ) was[d] = true;
            for( int b : V[a] ){
                if( inOrder[b] < inOrder[a] ) continue;
                for(int c : V[b]) if(was[c] && inOrder[c] > inOrder[a]) triangles.push_back( VI({a,b,c}) );
            }
            for( int d : revV[a] ) was[d] = false;
        }

        return triangles;
    }

}