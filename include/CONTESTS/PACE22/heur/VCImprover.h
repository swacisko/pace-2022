//
// Created by sylwester on 12/26/21.
//

#ifndef ALGORITHMSPROJECT_VCIMPROVER_H
#define ALGORITHMSPROJECT_VCIMPROVER_H

#include <CONTESTS/PACE22/Config.h>
#include "Makros.h"

class VCImprover{
public:

    VCImprover(Config c) : cnf(c){}

    VI improveByVertexCoverMinimizeArcs(VVI V, VI dfvs, int update_frequency );

    VI improveByVertexCoverRandom(VVI V, VI dfvs );

    VI createOrderMinimizeArcs( VVI V, VI dfvs, int update_frequency );

    VI createOrderByDFAS(VVI V, VI dfvs );

    VI mergeOrders( VVI V, VI ord1, VI ord2 );

    VI createOrderIterativeMerge( VVI & V, VI init_ord, int iterations );

    VI improveByVertexCoverForSemiTopologicalOrder(VVI V, VI order, VI current_dfvs = {} );

    bool minimize_found_vc = true;

    Config cnf;
};
#endif //ALGORITHMSPROJECT_VCIMPROVER_H
