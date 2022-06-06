//
// Created by sylwester on 8/13/19.
//

#ifndef ALGORITHMSPROJECT_TOPOSORT_H
#define ALGORITHMSPROJECT_TOPOSORT_H

#include "Makros.h"

class TopoSort {
public:

    TopoSort( VVI & structure );

    VI sortTopologically();

private :

    int N;
    VVI V;
    VI deg;
    VI kol;
    VB was;

    void DFS( int num );

    void topoSort();
};


#endif //ALGORITHMSPROJECT_TOPOSORT_H
