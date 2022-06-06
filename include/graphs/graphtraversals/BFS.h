//
// Created by sylwester on 8/26/19.
//

#ifndef ALGORITHMSPROJECT_BFS_H
#define ALGORITHMSPROJECT_BFS_H

#include <Constants.h>
#include "Makros.h"

namespace BFS{

    extern VVI getBfsLayers(VVI & V, VI sources);

    extern VVI getBfsLayers(VVI & V, int source);

}


#endif //ALGORITHMSPROJECT_BFS_H
