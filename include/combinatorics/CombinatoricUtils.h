//
// Created by sylwester on 8/27/19.
//

#ifndef ALGORITHMSPROJECT_COMBINATORICUTILS_H
#define ALGORITHMSPROJECT_COMBINATORICUTILS_H

#include <utils/RandomNumberGenerators.h>
#include "Makros.h"

namespace CombinatoricUtils{

    extern VI getRandomPermutation( int N );

    extern VI getRandomPermutation( int N, unsigned seed );

    extern VI getFullSetDifference( int N, VI  A );

    extern VI getRandomSequence( int U, int N );

    extern VI getRandomSubset( int U, int L );

    extern VI getRandomSubset( int U, int L, unsigned seed );

}


#endif //ALGORITHMSPROJECT_COMBINATORICUTILS_H
