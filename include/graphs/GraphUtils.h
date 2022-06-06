//
// Created by sylwester on 8/8/19.
//

#ifndef ALGORITHMSPROJECT_GRAPHUTILS_H
#define ALGORITHMSPROJECT_GRAPHUTILS_H

#include "Makros.h"

class GraphUtils {
public:

    static VI getComplimentaryNodes( VVI & V, VI & nodes );

    static int countEdges(VVI & V, const bool directed = false);

    template< class wtype>
    static int countEdges( vector< vector< pair<int,wtype> > >& V){
        int res = 0;
        for(auto& v : V) res += v.size();
        return res >> 1;
    }

    static VPII getGraphEdges( VVI & V, bool directed = false );
    static VPII getDirectedGraphEdges( VVI & V );

    static VVI getGraphForEdges(VPII edges, bool directed = false);

    template<class wtype>
    static vector<tuple<int,int,wtype>> getGraphEdges( vector< vector< pair<int,wtype> > >& V, bool directed = false ){
        vector<tuple<int,int,wtype>> res;
        res.reserve( countEdges(V) );
        for( int i=0; i<V.size(); i++ ){
            for( auto & [d,w] : V[i] ){
                if(directed || d>i ) res.emplace_back( i,d,w );
            }
        }
        return res;
    }

    static VI getNeighborhoodExclusive( VVI & V, VI & A, VB & helper );

    static VVI transposeGraph(VVI &v);

    static VVI reverseGraph(VVI & V){ return transposeGraph(V); }

    static VVI getComplimentaryGraph( VVI & V );

    static void addEdge( VVI & V, int a, int b, bool directed = false );

    static void removeEdge( VVI & V, int a, int b, bool directed = false );

    static void removeEdges(VVI& V, VPII& edges, bool directed = false);

    static void removeNodeFromGraph(VVI &V, int a);

    static VVI sortNodeNeighborhoods( VVI & V );

    static void writeGraphHumanReadable(VVI& V);

    static bool isSimple(VVI V);

    static double density( VVI & V, bool directed );
};


#endif //ALGORITHMSPROJECT_GRAPHUTILS_H
