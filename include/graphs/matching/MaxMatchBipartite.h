//
// Created by sylwester on 8/8/19.
//

#ifndef ALGORITHMSPROJECT_MAXMATCH_H
#define ALGORITHMSPROJECT_MAXMATCH_H

#include "Makros.h"

class MaxMatchBipartite {

public:

    VI getMaximalMatchingOfMinimalSizeRandom(VVI &G, int iterations);

    VI getMaximumMatchingInBipartition(VVI &G, VB &bipartition, bool fastSearch = true);

    VVI getMaximalSetOfDisjointAugmentingPaths( VVI & B, VB & bipartition, VI & matching );

    VPII convertToPairs( VI & matching );

    VI getMaximumHallViolator( VVI & V, VB & bipartition, VI & matching );

    bool findAugmentingPath( VVI & G, int beg, VB& bipartition, VI& matching, VI& path, VI& was, int currentWas );

    void applyAugmentingPath( VI& matching, VI& path );

private:

    VVI createLayerGraph(VVI &G, VB &bipartition, VI &matching);

    bool getMaximalSetOfDisjointAugmentingPaths(VVI &layerG, VB &bipartition, VB &was, int p, VVI &augmentingPaths);

};


#endif //ALGORITHMSPROJECT_MAXMATCH_H
