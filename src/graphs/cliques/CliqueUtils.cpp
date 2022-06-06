//
// Created by sylwester on 3/25/20.
//

#include <graphs/cliques/CliqueUtils.h>

#include "Makros.h"

namespace CliqueUtils{


    bool isClique( VVI& V, VI& clq){
        unordered_set<int> zb(ALL(clq));
        for( int p : clq ){
            int t = 0;
            for(int d : V[p]){
                if( zb.count(d) ) t++;
            }
            if( t != clq.size()-1 ) return false;
        }
        return true;
    }

    bool isClique(VVI &V, VI &clq, VB &helper) {
        if(clq.empty()) return false;

        for( int d : clq ) if( V[d].size()+1 < clq.size() ) return false;

        for (int d : clq) helper[d] = true;
        for (int p : clq) {
            int t = 0;
            for (int d : V[p]) {
                if (helper[d]) t++;
            }
            if (t != clq.size() - 1) {
                for (int d : clq) helper[d] = false;
                return false;
            }
        }
        for (int d : clq) helper[d] = false;
        return true;
    }


}
