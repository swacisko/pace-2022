//
// Created by sylwester on 8/13/19.
//

#include <Constants.h>
#include "graphs/toposort/TopoSort.h"
#include "StandardUtils.h"

TopoSort::TopoSort(VVI &structure) {
    V = structure;
    N = V.size();
}


void TopoSort::DFS( int num ){
    kol.PB( num );
    was[num] = 1;

    REP(i,SIZE( V[num] ) ){
        deg[ V[num][i] ]--;
        if( !deg[V[num][i]] ) DFS( V[num][i] );
    }
}

void TopoSort::topoSort(){

    VI perm(V.size());
    iota(ALL(perm),0);

    UniformIntGenerator gen( 0, 1e9 );
    StandardUtils::shuffle( perm, gen.getRNG() );

    for( VI & v : V ){
        StandardUtils::shuffle( v, gen.getRNG() );
    }



    REP(i,SIZE(V)){
        int d = perm[i];
        if( !deg[d] && !was[d] ) DFS(d);
    }
}

VI TopoSort::sortTopologically(){
    was = VB( SIZE(V),false );
    deg = VI( SIZE(V) );
    REP(i,SIZE(V)) REP(k,SIZE(V[i])) deg[ V[i][k] ]++;
    topoSort();

    return kol;
}


