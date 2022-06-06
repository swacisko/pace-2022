//
// Created by sylwester on 9/5/19.
//

#include <graphs/GraphWriter.h>

#include <graphs/GraphUtils.h>

namespace GraphWriter {


    void writeGraphStandardEdges(VVI &V, ostream &out, int addToId) {
        out << V.size() << " " << GraphUtils::countEdges(V) << endl;
        auto edges = GraphUtils::getGraphEdges(V);

        random_shuffle(ALL(edges));

        for (auto p : edges) {
            out << p.first + addToId << " " << p.second + addToId << endl;
        }
    }


    void writeGraphDIMACS(VVI &V, ostream &out, bool edgeFoolowE, int addToId ) {
        out << "p td " << V.size() << " " << GraphUtils::countEdges(V) << endl;
        auto edges = GraphUtils::getGraphEdges(V);
        for (auto p : edges) {
            if (edgeFoolowE) out << "e ";

            out << p.first+addToId << " " << p.second+addToId << endl;
        }
    }

}