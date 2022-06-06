//
// Created by sylwester on 9/5/19.
//

#ifndef ALGORITHMSPROJECT_GRAPHWRITER_H
#define ALGORITHMSPROJECT_GRAPHWRITER_H

#include "Makros.h"

namespace GraphWriter{

    extern void writeGraphDIMACS( VVI & V, ostream & out, bool edgeFoolowE = false, int addToId = 0 );
}

#endif //ALGORITHMSPROJECT_GRAPHWRITER_H
