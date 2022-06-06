//
// Created by sylwester on 3/16/20.
//

#ifndef ALGORITHMSPROJECT_CONNECTEDCOMPONENTS_H
#define ALGORITHMSPROJECT_CONNECTEDCOMPONENTS_H

#include "Makros.h"

namespace ConnectedComponents{

    extern void dfs(VVI& V, int &num, int &par, VB& was, VVI & comps);

    extern VVI getConnectedComponents( VVI & V );

    extern VI getConnectedComponentForNode(VVI &V, int v, VB &was);
}

#endif //ALGORITHMSPROJECT_CONNECTEDCOMPONENTS_H
