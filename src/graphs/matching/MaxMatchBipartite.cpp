//
// Created by sylwester on 8/8/19.
//

#include <graphs/GraphUtils.h>
#include <graphs/matching/MaxMatchBipartite.h>
#include <Constants.h>

#include "graphs/matching/MaxMatchBipartite.h"




VI MaxMatchBipartite::getMaximalMatchingOfMinimalSizeRandom(VVI &G, int iterations) {
    VPII edges = GraphUtils::getGraphEdges(G);
    VI res;
    int m = G.size();

    for( int ITER = 0; ITER < iterations; ITER++ ){
        VI temp(G.size(),-1);
        random_shuffle( edges.begin(), edges.end() );
        int M = 0;
        for( PII e : edges ){
            int a = e.first;
            int b = e.second;
            if( temp[a] == -1 && temp[b] == -1 ){
                temp[a] = b;
                temp[b] = a;
                M++;
            }
        }

        if( M < m ){
            m = M;
            res = temp;
        }
    }

    return res;
}


VI MaxMatchBipartite::getMaximumMatchingInBipartition(VVI &G, VB &bipartition, bool fastSearch) {
    VI matching(G.size(),-1);

    if( fastSearch ){
        VVI paths = getMaximalSetOfDisjointAugmentingPaths( G,bipartition,matching );

        while( !paths.empty() ){

            for( VI & path : paths ){
                applyAugmentingPath( matching,path );
            }

            paths = getMaximalSetOfDisjointAugmentingPaths( G,bipartition,matching );
        }

    }else{
        VI was(G.size(),-1);
        VI path;
        int currentWas = 0;

        for( int i=0; i<G.size(); i++ ){
            if( was[i] == currentWas || bipartition[i] == true ) continue; // if i already was in that node or it is in the other bipartition, i do not need to check it.
            bool found = findAugmentingPath( G,i, bipartition, matching, path, was, currentWas );

            if( found ){
                applyAugmentingPath( matching, path );
                path.clear();
                currentWas++;
            }
        }
    }

    return matching;
}


VPII MaxMatchBipartite::convertToPairs( VI & matching ){
    VPII res;
    for( int i=0; i<matching.size(); i++ ){
        if( matching[i] > i ) res.push_back( {i,matching[i]} );
    }
    return res;
}




bool MaxMatchBipartite::findAugmentingPath( VVI & G, int beg, VB& bipartition, VI& matching, VI& path, VI& was, int currentWas ) {

    if( beg < 0 ){
        cerr << "beg < 0 in findAugmentingPath" << endl;
        exit(1);
    }
    if( was[beg] == currentWas ) return false; // i do not process the same node more than once
    was[beg] = currentWas;

    path.push_back(beg);

    if( matching[beg] == -1 && path.size() >= 2 ) return true;


    if( path.size() % 2 == 0 ){
        if( findAugmentingPath( G, matching[beg], bipartition, matching, path,was,currentWas ) ) return true;
    }else{
        for(int p : G[beg]){
            if( bipartition[p] == bipartition[beg] ){
                cerr << "That should not happen if the graph structure is bipartite!" << endl;
                exit(1);
            }

            if( was[ p ] == currentWas || bipartition[p] == bipartition[beg] ) continue;

            if( findAugmentingPath( G, p, bipartition, matching, path, was, currentWas ) ) return true;
        }
    }

    path.pop_back();
    return false;
}


void MaxMatchBipartite::applyAugmentingPath( VI& matching, VI& path ){
    if( path.size() % 2 == 1 ){
        cerr << "error! augmenting path cannot have odd length! "; for(int p : path) cerr << p << " "; cerr << endl;
        exit(1);
    }
    for( int i=1; i<path.size(); i+=2 ){
        int a = path[i-1];
        int b = path[i];

        if( a == b ){
            cerr << "WTF in applyAugmentingPath" << endl;
            exit(1);
        }

        matching[a] = b;
        matching[b] = a;
    }
}


VI MaxMatchBipartite::getMaximumHallViolator(VVI &V, VB &bipartition, VI &matching) {
    VI res;

    VB was(V.size(),false);
    VI neigh;
    for( int i=0; i<V.size(); i++ ){
        if( matching[i] == -1 && bipartition[i] == false ){
            neigh.push_back(i);
            was[i] = true;
        }
    }

    for( int i=0; i<neigh.size(); i++ ){
        int p = neigh[i];

        if( bipartition[p] == false ){
            res.push_back(p);

            for( int d : V[p] ){
                if( !was[d] ){
                    neigh.push_back(d);
                    was[d] = true;
                }
            }
        }else{
            if( matching[p] == -1 ){
                cerr << "ERROR, unmatched node " << p << ", augmenting path exists and hall violator may not even exist! " << endl;
                exit(1);
            }else{
                int d = matching[p];

                neigh.push_back(d);
                was[d] = true;
            }
        }
    }

    return res;
}

VVI MaxMatchBipartite::getMaximalSetOfDisjointAugmentingPaths(VVI & G, VB &bipartition, VI &matching) {
    VVI disjointPaths;

    VVI layerG = createLayerGraph( G, bipartition,matching );

    VB was(G.size(),false); // this should not be neccessary

    for( int i=0; i<bipartition.size(); i++){
        if( bipartition[i] == false && matching[i] == -1 ) getMaximalSetOfDisjointAugmentingPaths( layerG,bipartition,was,i,disjointPaths );
    }

    return disjointPaths;
}

VVI MaxMatchBipartite::createLayerGraph(VVI &G, VB &bipartition, VI &matching) {
    VVI layerG = VVI(G.size());

    VI neigh;
    VI dst(G.size(),Constants::INF);
    for( int i=0; i<G.size(); i++ ){
        if(bipartition[i] == false && matching[i] == -1 ){
            neigh.push_back(i);
            dst[i] = 0;
        }
    }

    for(int i=0; i<neigh.size(); i++){
        int p = neigh[i];

        if( bipartition[p] == true ){

            int d = matching[p];
            if( d != -1 ){ // if the node is already matched, i add the edge to the graph
                layerG[p].push_back(d);
                if( dst[d] == Constants::INF ){
                    dst[d] = 1 + dst[p];
                    neigh.push_back( d );
                }
            }
        }else{
            for( int d : G[p] ){
                if( bipartition[d] == bipartition[p] ) continue;

                if( dst[d] == Constants::INF ){
                    layerG[p].push_back(d);
                    dst[d] = 1 + dst[p];
                    neigh.push_back( d );
                }else if( dst[d] == 1 + dst[p] || ( bipartition[d] == true && matching[d] == -1 ) ){
                    layerG[p].push_back(d);
                }
            }
        }
    }

    return layerG;
}

bool MaxMatchBipartite::getMaximalSetOfDisjointAugmentingPaths(VVI &layerG, VB &bipartition, VB &was, int p, VVI &augmentingPaths) {
    if( was[p] ) return false;
    was[p] = true;

    if( bipartition[p] == true && layerG[p].size() == 0 ){ // i found an augmenting path
        augmentingPaths.push_back( VI( 1,p ) );
        return true;
    }

    for( int i=layerG[p].size()-1; i >= 0; i-- ){
        int d = layerG[p][i];
        if( was[d] ){
            layerG[p].pop_back();
            continue;
        }

        if( getMaximalSetOfDisjointAugmentingPaths( layerG, bipartition,was,d,augmentingPaths ) ){ // if a path is found
            augmentingPaths.back().push_back(p);
            layerG[p].clear();
            return true;
        }else{
            layerG[p].pop_back();
        }
    }
    return false;
}

