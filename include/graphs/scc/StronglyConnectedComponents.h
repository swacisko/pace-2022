//
// Created by sylwester on 8/13/19.
//

#ifndef ALGORITHMSPROJECT_STRONGLYCONNECTEDCOMPONENTS_H
#define ALGORITHMSPROJECT_STRONGLYCONNECTEDCOMPONENTS_H

#include "Makros.h"

class StronglyConnectedComponents {
public:

    StronglyConnectedComponents( VVI & structure );

    StronglyConnectedComponents( VVI & V, VVI & revV );

    void createStronglyConnectedComponents();

    VVI getComponents(){ return Comps; }

    VI PostOrder;
    VVI V;
    VVI revV;

    int N;
    VVI Comps;

    VI compParent;

    VI temp_comp;

    VB was;

    void PO_DFS( int num, int & ponum );

    void Add_PO_DFS( int num );

    VVI transpose( VVI & v );

};


#endif //ALGORITHMSPROJECT_STRONGLYCONNECTEDCOMPONENTS_H
