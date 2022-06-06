//
// Created by sylwester on 8/8/19.
//

#include <graphs/GraphInducer.h>

#include "graphs/GraphInducer.h"

InducedGraph GraphInducer::induce( VVI & V, VI & nodes ){
    InducedGraph g;
    g.nodes = nodes;
    g.par = &V;
    int N = SIZE(nodes);

    g.perm = unordered_map<int,int>();
    g.perm.reserve( nodes.size() * 2 );
    REP( i, N ) g.perm[ nodes[i] ] = i;

    g.V = VVI( nodes.size() );
    for( int i=0; i<nodes.size(); i++ ){
        for( int d : V[ nodes[i] ] ){
            auto it = g.perm.find(d);
            if( it != g.perm.end() ){
                int indD = it->second;
                g.V[ i ].push_back( indD );
            }
        }
    }
    return g;
}

ostream& operator<<(ostream& str, InducedGraph& g){
    str << "Par: " << *g.par << endl
        << "Nodes: " << g.nodes << endl
        << "Perm: " << g.perm << endl
        << "V: " << g.V << endl;
    return str;
}

InducedGraphPI GraphInducer::induce(VVPII &V, VI &nodes) {
    InducedGraphPI g;
    g.nodes = nodes;
    g.par = &V;
    int N = SIZE(nodes);

    g.perm = unordered_map<int,int>();
    g.perm.reserve( nodes.size() * 2 );
    for(int i=0; i<N; i++) g.perm[ nodes[i] ] = i;

    g.V = VVPII(nodes.size() );
    for( int i=0; i<nodes.size(); i++ ){
        for( auto pr : V[ nodes[i] ] ){
            int d = pr.first;
            int w = pr.second;
            auto it = g.perm.find(d);
            if( it != g.perm.end() ){
                int indD = it->second;
                g.V[ i ].push_back({indD,w} );
            }
        }
    }
    return g;
}


ostream& operator<<(ostream& str, InducedGraphPI& g){
    str << "Par: " << *g.par << endl
        << "Nodes: " << g.nodes << endl
        << "Perm: " << g.perm << endl
        << "V: " << g.V << endl;
    return str;
}


InducedGraph GraphInducer::induce( VVI & V, VPII & edges, bool directed ){
    InducedGraph g;
    g.edges = edges;
    g.par = &V;

    g.perm = unordered_map<int,int>();
    g.perm.reserve( edges.size() / 4 );
    unordered_map<int,int>::iterator it;

    REP( i, SIZE(edges) ){
        int a = edges[i].ST;
        int b = edges[i].ND;

        it = g.perm.find(a);
        if( it == g.perm.end() ){
            g.nodes.PB( a );
            g.perm[a] = SIZE( g.nodes )-1;
            g.V.PB( V[a] );
        }

        it = g.perm.find(b);
        if( it == g.perm.end() ){
            g.nodes.PB( b );
            g.perm[b] = SIZE( g.nodes )-1;
            g.V.PB( V[b] );
        }
    }

    set< PII > edgeSet( ALL(edges) );
    REP( i, SIZE(g.V) ){
        REP( k, SIZE( g.V[i] ) ){

            int a = g.nodes[i];
            int b = g.V[i][k];
            bool isEdge = edgeSet.count( MP(a,b) );
            if( directed == false && isEdge == false ) isEdge = edgeSet.count( MP(b,a) );

            if( !isEdge ){
                swap( g.V[i][k], g.V[i].back() );
                g.V[i].pop_back();
                k--;
            }else{
                g.V[i][k] = g.perm[b];
            }
        }
    }

    return g;
}

InducedGraph GraphInducer::induceByNonisolatedNodes(VVI &V) {
    VI nodes;
    VI degs(V.size(),0);
    for( int i=0; i<V.size(); i++ ){
        degs[i] += V[i].size();
        for(int d : V[i]) degs[d]++;
    }

    for( int i=0; i<V.size(); i++ ) if(degs[i]) nodes.push_back(i);

    return GraphInducer::induce(V,nodes);
}
