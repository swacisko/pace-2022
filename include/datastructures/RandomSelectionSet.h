//
// Created by sylwester on 3/26/22.
//

#ifndef ALGORITHMSPROJECT_RANDOMSELECTIONSET_H
#define ALGORITHMSPROJECT_RANDOMSELECTIONSET_H

#include <utils/RandomNumberGenerators.h>
#include "Makros.h"

class RandomSelectionSet{
public:

    RandomSelectionSet(int N );

    int getRandomElement();

    void set( int id, LL val );

    int lowerBound( LL val );

    LL get( int id );

    LL getProbab(int id){ return probab[id]; }

private:

    int N;
    UniformIntGenerator rnd;

    VLL probab;

    VLL cover;

    VI cnt_num;

    VLL max_el;

    void add( int id, LL val );

    int L(int x){ return x<<1; }
    int R(int x){ return (x<<1)+1; }
    int par(int x){ return x >> 1; }

};

#endif //ALGORITHMSPROJECT_RANDOMSELECTIONSET_H
