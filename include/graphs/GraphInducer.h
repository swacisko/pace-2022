//
// Created by sylwester on 8/8/19.
//

#ifndef ALGORITHMSPROJECT_GRAPHINDUCER_H
#define ALGORITHMSPROJECT_GRAPHINDUCER_H

#include "Makros.h"


struct InducedGraph{
    VVI *par;
    VI nodes;

    unordered_map<int,int> perm;

    VPII edges;
    VVI V;

    friend ostream& operator<<(ostream& str, InducedGraph& g);

    void write(){
        cerr << "Graph induced by: " << flush; WRITE(nodes);
        if( !edges.empty() ){
            cerr << "induced by edges:" << endl;
            REP(i,SIZE(edges)) cerr << WRP(edges[i]) << endl;
        }
        WRITE_ALL( V, "Graph structure",0 );

    }
};

struct InducedGraphPI{
    VVPII *par;
    VI nodes;

    unordered_map<int,int> perm;

    VPII edges;
    VVPII V;

    friend ostream& operator<<(ostream& str, InducedGraphPI& g);
};


class GraphInducer{
public:

    static InducedGraph induce( VVI & V, VI & nodes );

    static InducedGraphPI induce(VVPII & V, VI & nodes );

    static InducedGraph induce( VVI & V, VPII & edges, bool directed = false );

    static InducedGraph induceByNonisolatedNodes( VVI & V );

};

#endif //ALGORITHMSPROJECT_GRAPHINDUCER_H
