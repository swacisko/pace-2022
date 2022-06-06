//
// Created by sylwester on 8/13/19.
//

#include <graphs/scc/StronglyConnectedComponents.h>
#include <Constants.h>

StronglyConnectedComponents::StronglyConnectedComponents(VVI &structure) : V(structure) {
    N = V.size();
}

StronglyConnectedComponents::StronglyConnectedComponents(VVI &V, VVI &revV) {
    this->V = V;
    this->revV = revV;
    N = V.size();
}


void StronglyConnectedComponents::createStronglyConnectedComponents(){
    was = VB(N,0);
    compParent = VI(N,Constants::INF);
    PostOrder = VI(N);
    int ponum = 0;
    REP( i,N ) if( !was[i] ) PO_DFS( i,ponum );

    VI sorted( N );
    REP( i,N ) sorted[ PostOrder[i] ] = i;

    if(revV.empty()) revV = transpose( V );

    swap(V, revV);

    was = VB(N,0);
    FORD( i,N-1,0 ){
        temp_comp.clear();
        if( !was[ sorted[i] ] ) Add_PO_DFS( sorted[i] );
        if( SIZE(temp_comp) ) Comps.PB( temp_comp );
    }

    swap(V,revV);

    for( int i=0; i<Comps.size(); i++ ){
        for(int d : Comps[i]){
            compParent[d] = i;
        }
    }
}


void StronglyConnectedComponents::PO_DFS( int num, int & ponum ){
    was[ num ] = true;
    REP( i,SIZE(V[num]) ) if( !was[ V[num][i] ] ) PO_DFS( V[num][i], ponum );
    PostOrder[ num ] = ponum++;
}

void StronglyConnectedComponents::Add_PO_DFS( int num ){
    was[num] = true;
    temp_comp.PB( num );
    REP(i, SIZE(V[num]) ) if( !was[ V[num][i] ] ) Add_PO_DFS( V[num][i] );
}

VVI StronglyConnectedComponents::transpose( VVI & v ){
    VVI g( SIZE(v) );
    REP( i,SIZE(v) )	REP( k,SIZE(v[i]) )	g[ v[i][k] ].PB(i);
    return g;
}

