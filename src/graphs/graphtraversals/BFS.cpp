//
// Created by sylwester on 8/26/19.
//

#include <graphs/graphtraversals/BFS.h>

#include "graphs/graphtraversals/BFS.h"

namespace BFS{

    VVI getBfsLayers(VVI &V, VI sources) {
        int N = SIZE(V);

        VB was(N,false);
        for( int beg : sources ) was[beg] = true;
        VI layer = sources;
        VI temp;

        VVI res;

        do{
            res.PB( layer );

            REP( i, SIZE(layer) ){
                int n = layer[i];


                REP( k, SIZE( V[n] ) ){
                    int d = V[n][k];
                    if( !was[d] ){
                        was[d] = true;
                        temp.PB(d);
                    }
                }
            }

            swap( layer,temp );
            temp.clear();

        }while( !layer.empty() );

        return res;
    }

    VVI getBfsLayers(VVI & V, int source){ return getBfsLayers(V, VI({source})); }


}