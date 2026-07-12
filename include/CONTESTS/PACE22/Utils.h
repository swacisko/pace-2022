//
// Created by sylwester on 12/20/21.
//

#ifndef ALGORITHMSPROJECT_UTILS_H
#define ALGORITHMSPROJECT_UTILS_H

#include "Makros.h"
#include "CollectionOperators.h"
#include "Config.h"

namespace Utils{

    extern bool hasLoop( VVI & V,  int a );

    extern void writeRemainingGraph(VVI & V);

    extern void writeNeighborhood(VVI & V, int v);
    extern void writeNeighborhood(VVI & V, VVI & revV, int v);

    LL getSetHash( int N, VI & s, int seed = 8'592'374 );


    VI getMinVcCPSAT( VVI & V, int thread_workers = 1 );

}

#endif //ALGORITHMSPROJECT_UTILS_H
