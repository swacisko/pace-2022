//
// Created by sylwester on 8/28/19.
//

#ifndef ALGORITHMSPROJECT_KERNELIZERVC_H
#define ALGORITHMSPROJECT_KERNELIZERVC_H

#include "Makros.h"

class KernelizerVC{
public:

    //******************************** BEGINNING OF KERNELIZATION SECTION
    pair<VI, VPII> initialKernelization(VVI &G);

    pair<VI, VPII> loopKernelization(VVI &G);

    pair<VI, VPII> dominationRule( VVI & V, int maxDegree );

    pair<VI, VPII> crownDecomposition(VVI &G);

    pair<VI,VPII> lpDecomposition(VVI& G);

    pair<VI, VPII> degree1Kernelization( VVI & V );

private:

    VPII removeNodesFromGraph(VVI &G, VI &toRemove, set<PII> &nodeDegrees);

    VPII removeNodeFromGraph( VVI & G, int toRemove );

    //******************************** END OF KERNELIZATION SECTION

    VI getMinVCForBipartition(VVI &G, VI &matching, VB &bipartition);

};

#endif //ALGORITHMSPROJECT_KERNELIZERVC_H
