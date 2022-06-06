//
// Created by sylwester on 5/19/22.
//

#include <combinatorics/CombinatoricUtils.h>
#include <utils/StandardUtils.h>
#include "CONTESTS/PACE22/exact/CliqueCoverLS.h"

CliqueCoverLS::CliqueCoverLS(VVI &VV, Config &c) {
    V = VV;
    cnf = c;
}

VVI CliqueCoverLS::coverLS(int iterations) {
    int N = V.size();

    VI in_clique(N,0);
    iota(ALL(in_clique),0);

    VI clique_size(N,1);
    VI edges_to_clique(N,0);

    VB was(N,false), was2(N,false);
    VI cand;
    VI temp;

    int total_cliques = N;

    /**
     * Moves node v to cluster cl.
     */
    auto moveNode = [&]( int v, int cl ){
        clique_size[in_clique[v]]--;
        if(clique_size[in_clique[v]] == 0) total_cliques--;
        in_clique[v] = cl;
        clique_size[in_clique[v]]++;
    };

    for( int iter = 1; iter < iterations; iter++ ){

        VI perm = CombinatoricUtils::getRandomPermutation(N);

        if( iter & 1 ) {
            sort(ALL(perm), [&](int a, int b) {
                     return clique_size[in_clique[a]] < clique_size[in_clique[b]];
                 }
            );
        }

        for (int v : perm) { // for each node in order

            { // creating candidates where node v can be moved
                cand.clear();
                temp.clear();
                for (int d : V[v]) {
                    int cld = in_clique[d];
                    edges_to_clique[cld]++;
                    if (edges_to_clique[cld] == clique_size[cld]) cand.push_back(cld);
                    if(!was[cld]){
                        was[cld] = true;
                        temp.push_back(cld);
                    }
                }

                for(int d : temp){
                    was[d] = false;
                    edges_to_clique[d] = 0; // clearing edges_to_clique
                }
            }

            int best_cl = -1;
            if(!cand.empty()){ // selecting best candidate to move node v
                best_cl = *max_element(ALL(cand), [&](int a, int b){
                    return clique_size[a] < clique_size[b]; // selecting largest cluster
                });
            }

            if(best_cl != -1){ // moving node v to selected clique
                moveNode(v, best_cl);
            }
        }

//        { // #TEST - just an assertion - remove after tests
//            VI cl(N,0);
//            for( int i=0; i<N; i++ ) cl[in_clique[i]]++;
//            for( int i=0; i<N; i++ ) assert(cl[i] == clique_size[i]);
//        }

        if(cnf.write_logs){
            clog << "\rAfter iteration " << iter << ", total_cliques: " << total_cliques << flush;
        }
    }

    if(cnf.write_logs){
        ENDL(1);
        DEBUG(total_cliques);
    }

    VVI cliques = StandardUtils::partitionToLayers(in_clique);
    VVI res; res.reserve(total_cliques);
    for( VI & cl : cliques ) if(!cl.empty()) res.push_back(cl);

    if(res.size() != total_cliques){
        DEBUG(res.size());
        DEBUG(total_cliques);
        assert(res.size() == total_cliques);
    }

    return res;
}
