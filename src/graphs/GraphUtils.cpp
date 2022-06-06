//
// Created by sylwester on 8/8/19.
//


#include <graphs/GraphUtils.h>
#include <graphs/components/ConnectedComponents.h>

#include "combinatorics/CombinatoricUtils.h"
#include <CONTESTS/PACE22/Utils.h>
#include "graphs/cliques/CliqueExtension.h"


VI GraphUtils::getComplimentaryNodes( VVI & V, VI & nodes ){
    VB inNodes( V.size(),false );
    for(auto p : nodes) inNodes[p] = true;
    VI res;
    for( int i=0; i<V.size(); i++ ){
        if( !inNodes[i] ) res.push_back(i);
    }
    return res;
}

VPII GraphUtils::getGraphEdges( VVI & V, bool directed ){
    if(directed) return getDirectedGraphEdges(V);

    VPII res;
    res.reserve( countEdges(V) );
    for( int i=0; i<V.size(); i++ ){
        for( int d : V[i] ){
            if(d>i) res.push_back( {i,d} );
        }
    }
    return res;
}

VPII GraphUtils::getDirectedGraphEdges( VVI & V ){
    VPII res;
    res.reserve( countEdges(V) );
    for( int i=0; i<V.size(); i++ ){
        for( int d : V[i] ){
            res.push_back( {i,d} );
        }
    }
    return res;
}


VVI GraphUtils::transposeGraph(VVI &v) {
    VVI g( SIZE(v) );
    REP( i,SIZE(v) )	REP( k,SIZE(v[i]) )	g[ v[i][k] ].PB(i);
    return g;
}

void GraphUtils::addEdge(VVI &V, int a, int b, bool directed) {
    V[a].push_back(b);
    if( !directed ) V[b].push_back(a);
}

void GraphUtils::removeEdge(VVI &V, int a, int b, bool directed) {
    auto rem = [ &V ](int a, int b){
        for( int i=(int)V[a].size()-1; i>=0; i-- ){
            if( V[a][i] == b ){
                swap(V[a][i], V[a].back() );
                V[a].pop_back();
                break;
            }
        }
    };


    rem(a,b);
    if(!directed) rem(b,a);

}

void GraphUtils::removeNodeFromGraph(VVI &V, int a) {
    for( int d : V[a] ) removeEdge( V,d,a,true );
    V[a].clear();
}

int GraphUtils::countEdges(VVI &V, const bool directed) {
    int res = 0;
    for(auto& v : V) res += v.size();
    if(!directed) return res >> 1;
    else return res;
}

VVI GraphUtils::sortNodeNeighborhoods( VVI & V ){
    VVI V2(V.size());
    for(int i=0; i<V.size(); i++) V2[i].reserve(V[i].size());
    for( int i=0; i<V.size(); i++ ){
        for( int d : V[i] ) V2[d].push_back(i);
    }
    return V2;
}


VVI GraphUtils::getComplimentaryGraph(VVI &V) {
    VVI V2(V.size());
    for( int i=0; i<V.size(); i++ ){
        VI tmp = V[i];
        tmp.push_back(i);
        sort(ALL(tmp));
        V2[i] = getComplimentaryNodes( V,tmp );
    }

    bool correct = true;
    for( int i=0; i<V.size(); i++ ){
        VB was(V.size(),false);
        for( int d : V[i] ) was[d] = true;
        for(int d : V2[i]) was[d] = true;
        if( V[i].size() + V2[i].size() != V.size() -1 ) correct = false;
        for(int k=0; k<V.size(); k++) if( k != i && !was[k]) correct = false;
    }

    assert(correct);

    return V2;
}

void GraphUtils::writeGraphHumanReadable(VVI &V) {
    clog << "****" << endl;
    for( int i=0; i<V.size(); i++ ){
        clog << i << ": ";
        VI neigh = V[i];
        sort(ALL(neigh));
        for(int d : neigh) clog << d << "  ";
        clog << endl;
    }
}

void GraphUtils::removeEdges(VVI &V, VPII &edges, bool directed) {

    if( directed == false ){
        int E = edges.size();
        for( int i=0; i<E; i++ ) edges.emplace_back( edges[i].second, edges[i].first ); // adding reverse edges to remove
    }

    sort( ALL(edges) );

    for( int i=0; i<edges.size(); i++ ){
        int p = i;
        unordered_set<int> toRemove;
        while( p < edges.size() && edges[p].first == edges[i].first ){
            toRemove.insert( edges[p].second );
            p++;
        }

        int t = edges[i].first;
        for( int k=(int)V[t].size()-1; k>=0; k-- ){
            int d = V[t][k];
            if( toRemove.count(d) ){
//                cerr << "Removing edge " << t << " -> " << V[t][k] << endl;
                swap( V[t][k], V[t].back() );
                V[t].pop_back();
            }
        }

        i = p-1;
    }

}

VI GraphUtils::getNeighborhoodExclusive(VVI &V, VI &A, VB &helper) {
    VI res;
    for( int a : A ) helper[a] = true;
    for( int a : A ){
        for( int w : V[a] ){
            if(!helper[w]){
                res.push_back(w);
                helper[w] = true;
            }
        }
    }

    for( int a : A ) helper[a] = false;
    for(int d : res) helper[d] = false;
    return res;
}

VVI GraphUtils::getGraphForEdges(VPII edges, bool directed) {
    int N = 0;
    for(auto & [a,b] : edges) N = max(N, max(a,b));
    VVI V(N+1);
    for( auto & [a,b] : edges ) addEdge(V,a,b,directed);
    return V;
}

bool GraphUtils::isSimple(VVI V) {
    int N = V.size();
     VB helper(N,false);

     for( int i=0; i<N; i++ ){
         for( int d : V[i] ){
             if( d == i || helper[d] ) return false;
             helper[d] = true;
         }

         for(int d : V[i]) helper[d] = false;
     }

     return true;
}

double GraphUtils::density(VVI &V, bool directed) {
    LL denom = V.size();
    denom *= denom-1;
    denom >>= 1;

    double density = countEdges(V, directed);
    return density / denom;
}


