//
// Created by sylwester on 8/8/19.
//

#ifndef ALGORITHMSPROJECT_GRAPHREADER_H
#define ALGORITHMSPROJECT_GRAPHREADER_H

#include "Makros.h"

namespace GraphReader {

    extern VVI readGraphStandardEdges(istream &cin, bool directed = false);

    extern VVI readGraphDIMACSWunweighed(istream &cin, bool edgeFoolowE = false);
};


#endif //ALGORITHMSPROJECT_GRAPHREADER_H
