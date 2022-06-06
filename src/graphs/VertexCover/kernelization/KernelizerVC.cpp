//
// Created by sylwester on 8/28/19.
//

#include <graphs/matching/MaxMatchBipartite.h>
#include <graphs/GraphUtils.h>
#include <graphs/VertexCover/kernelization/KernelizerVC.h>
#include <combinatorics/CombinatoricUtils.h>


pair<VI, VPII> KernelizerVC::loopKernelization(VVI &G) {
    VI nodes;
    for( int i=0; i<G.size(); i++ ){
        for( int p : G[i] ){
            if( p == i ){
                nodes.push_back(i);
                break;
            }
        }
    }

    set<PII> nodeDegrees; // only to call below function
    VPII edges = removeNodesFromGraph(G, nodes, nodeDegrees );

    return { nodes,edges };
}


pair<VI, VPII> KernelizerVC::crownDecomposition(VVI &G) {

    MaxMatchBipartite matcher;
    VI maximalMatching = matcher.getMaximalMatchingOfMinimalSizeRandom(G, 10);

    VB bipartition(G.size(),false);
    for( int i=0; i<maximalMatching.size(); i++ ){
        if( maximalMatching[i] == -1 ) bipartition[i] = true;
    }

    VPII M1 = matcher.convertToPairs( maximalMatching );
    VI maximumMatchingBipartite = matcher.getMaximumMatchingInBipartition( G, bipartition );
    VPII M2 = matcher.convertToPairs( maximumMatchingBipartite );
    VI X = getMinVCForBipartition(G,maximumMatchingBipartite, bipartition); // X is vertex cover of the graph with given bipartition

    VI H;
    for( int p : X ){
        if( bipartition[p] == false && maximalMatching[p] != -1 ) H.push_back(p);
    }

    pair<VI, VPII> res;
    res.first = H;
    set<PII> nodeDegrees; // this is here only to call below function
    res.second = removeNodesFromGraph( G,H, nodeDegrees );

    return res;
}

pair<VI, VPII> KernelizerVC::initialKernelization(VVI &G) {

    pair<VI, VPII> res;
    auto kernel = loopKernelization(G);
    res = kernel;

    int ITERATIONS = 1;
    int iter = 0;

    int crownAndLPChecks = 0;

    do{
        iter++;

        kernel = degree1Kernelization(G);
        res.first.insert( res.first.end(), kernel.first.begin(), kernel.first.end() );
        res.second.insert( res.second.end(), kernel.second.begin(), kernel.second.end() );

        kernel = dominationRule(G,G.size());
        res.first.insert( res.first.end(), kernel.first.begin(), kernel.first.end() );
        res.second.insert( res.second.end(), kernel.second.begin(), kernel.second.end() );

        if( kernel.first.size() > 0 ) continue;

        crownAndLPChecks++;
        if( crownAndLPChecks <= 1 ){

            kernel = crownDecomposition(G); // CROWN
            if( kernel.first.size() > 0 ) iter = 0;
            res.first.insert( res.first.end(), kernel.first.begin(), kernel.first.end() );
            res.second.insert( res.second.end(), kernel.second.begin(), kernel.second.end() );

            kernel = lpDecomposition(G); // LP
            if( kernel.first.size() > 0 ) iter = 0;
            res.first.insert( res.first.end(), kernel.first.begin(), kernel.first.end() );
            res.second.insert( res.second.end(), kernel.second.begin(), kernel.second.end() );

            if( kernel.first.size() > 0 ) continue;

        }

    }while( kernel.first.size() > 0 || iter < ITERATIONS );

    return res;
}


VI KernelizerVC::getMinVCForBipartition(VVI &G, VI &matching, VB &bipartition) {

    VI U;
    for( int i=0; i<G.size(); i++ ){
        if( matching[i] == -1 && bipartition[i] == false ) U.push_back(i);
    }

    VI Z;
    VB was(G.size(),false);
    VI neigh = U; // this if for bfs to mark all alternating paths.
    for( int i=0; i<neigh.size(); i++ ){
        int d = neigh[i];
        if( was[d] ) continue;
        was[d] = true;
        Z.push_back(d);

        for( int p : G[d] ){
            if( was[p] || bipartition[p] == bipartition[d] ) continue;
            Z.push_back(p);
            neigh.push_back( matching[p] );

            was[p] = true;
        }
    }

    VB Z2(G.size(),false);
    for( int p : Z ) Z2[p] = true;

    VI K;

    for( int i=0; i<G.size(); i++ ){
        if( G[i].size() == 0 ) continue;

        if( bipartition[i] == false && !Z2[i] ){
            K.push_back(i);
        }else if( bipartition[i] == true && Z2[i] ){
            K.push_back(i);
        }
    }

    return K;
}

pair<VI, VPII> KernelizerVC::lpDecomposition(VVI &G) {
    VVI layerG( 2*G.size() );

    for( int i=0; i<G.size(); i++ ){
        for( int d : G[i] ){
            GraphUtils::addEdge( layerG, i, G.size()+d );
        }
    }

    VB bipartition(layerG.size(),false);
    for( int i=0; i<G.size(); i++ ) bipartition[i] = true;

    MaxMatchBipartite matcher;
    VI matching = matcher.getMaximumMatchingInBipartition(layerG,bipartition);

    pair<VI, VPII> kernel;
    VI minVC = getMinVCForBipartition( layerG, matching, bipartition );
    VB minVC2( bipartition.size(), false );
    for(int p : minVC) minVC2[p] = true;

    for( int i=0; i<G.size(); i++ ){
        if( minVC2[i] && minVC2[G.size()+i] ) kernel.first.push_back(i);
        else if( !minVC2[i] && !minVC2[ G.size()+i ] ){
            auto v = G[i];
            for(auto d : v){
                GraphUtils::removeEdge(G,i,d);
            }
        }
    }

    set<PII> nodeDegrees;
    kernel.second = removeNodesFromGraph( G, kernel.first, nodeDegrees );

    return kernel;
}

VPII KernelizerVC::removeNodesFromGraph(VVI &G, VI &toRemove, set<PII> &nodeDegrees) {
    VPII res;
    VI nodes;
    for( int p : toRemove ){
        VPII temp = removeNodeFromGraph( G, p );
        res.insert( res.end(), temp.begin(), temp.end() );
        for( auto p : temp ){
            nodes.push_back(p.first);
            nodes.push_back(p.second);
        }
    }

    sort( nodes.begin(), nodes.end() );
    nodes.resize( unique( nodes.begin(), nodes.end() ) - nodes.begin() );
    for( int p : nodes ) if( G[p].size() > 0 ) nodeDegrees.insert( { G[p].size(),p } );

    return res;
}

VPII KernelizerVC::removeNodeFromGraph( VVI & G, int toRemove ){
    VPII res;
    for( int p : G[toRemove] ) res.push_back( { toRemove, p } );
    for( PII p : res ){
        int a = p.first;
        int b = p.second;
        for( int i=G[a].size()-1; i>=0; i-- ){
            if( G[a][i] == b ){
                swap(G[a][i], G[a].back());
                G[a].pop_back();
                break;
            }
        }

        swap(a,b);
        for( int i=G[a].size()-1; i>=0; i-- ){
            if( G[a][i] == b ){
                swap(G[a][i], G[a].back());
                G[a].pop_back();
                break;
            }
        }
    }
    return res;
}


pair<VI, VPII> KernelizerVC::degree1Kernelization(VVI &V) {

    pair<VI,VPII> kernel;
    VI neigh;
    VI perm = CombinatoricUtils::getRandomPermutation(V.size());
    for(int i : perm){
        if( V[i].size() == 1 ){
            neigh.push_back(i);
        }
    }

    for( int i=0; i<neigh.size(); i++ ){
        int p = neigh[i];
        if( V[p].empty() ) continue;

        int d = V[p][0];

        kernel.first.push_back(d);

        for( int x : V[d] ){
            GraphUtils::removeEdge( V, x,d, true );
            kernel.second.push_back( {x,d} );
        }

        for( int x : V[d] ){
            if( V[x].size() == 1 ) neigh.push_back(x);
        }

        V[d].clear();
    }

    return kernel;
}

pair<VI, VPII> KernelizerVC::dominationRule(VVI &V, int maxDegree) {

    pair<VI,VPII> res;
    VB neigh(V.size(),false);
    VI perm = CombinatoricUtils::getRandomPermutation(V.size());

    for( int i : perm){

        if( V[i].size() >= 2 && V[i].size() <= maxDegree ){

            for( int d : V[i] ) neigh[d] = true;
            neigh[i] = true;

            for( int d : V[i] ){

                bool isDominated = true;
                for( int p : V[d] ){
                    if( !neigh[p] ){
                        isDominated = false;
                        break;
                    }
                }

                if( isDominated ){
                    res.first.push_back(i);
                    for( int p : V[i] ){
                        res.second.push_back( {i,p} );
                        neigh[p] = false;
                    }
                    neigh[i] = false;

                    GraphUtils::removeNodeFromGraph(V, i);
                    break;
                }
            }

            for( int t : V[i] ) neigh[t] = false;
            neigh[i] = false;
        }
    }

    return res;
}

