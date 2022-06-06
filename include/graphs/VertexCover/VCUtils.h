//
// Created by sylwester on 8/7/19.
//

#ifndef ALGORITHMSPROJECT_VCUTILS_H
#define ALGORITHMSPROJECT_VCUTILS_H

#include "Makros.h"

class VCUtils{
public:

    static bool isVertexCover( VVI & V, VI & vc );

    static bool isIndependentSet( VVI & V, VI & S );

    static VI getMinCVUsingFastVC( VVI V, int milliseconds,  VI init_vc = {}  );

};




#endif //ALGORITHMSPROJECT_VCUTILS_H
