//
// Created by sylwe on 07/08/2025.
//

#ifndef ELDREDUCER_H
#define ELDREDUCER_H

#include "Makros.h"

class ELDReducer {
public:
    ELDReducer( int NN ) : N(NN), V(N), inW0(N), inW1(N), was(N), was2(N), visited(N), helper(N) {}

    VI reduceELD(VVI & V0);

private:

    int N;
    VVI V;

    VB inW0, inW1, visited, was, was2, helper;
    VI temp, temp2;


    void clearAll();
};

#endif //ELDREDUCER_H
