//
// Created by sylwester on 4/8/22.
//

#ifndef ALGORITHMSPROJECT_HITTINGSETLS_H
#define ALGORITHMSPROJECT_HITTINGSETLS_H


#include <CONTESTS/PACE22/Config.h>
#include <utils/RandomNumberGenerators.h>
#include "Makros.h"

class HittingSetLS{
public:

    HittingSetLS( VVI sets, VI hs, Config cnf );

    VI hsImprovementLS2( int iters, int lower_bound = -1, int deviate_perm_frequency = 200,
                         int add_violation_frequency = 500, int remove_violation_frequency = 700 );

    int persistent_search_iterations_left = 0;

private:

    VVI V;

    VI hs;

    VB in_hs;

    Config cnf;

    VVI sets;

    int N;

    int S;

    VI covered_by;

    VLL hashes;

    UniformIntGenerator rnd;

};

#endif //ALGORITHMSPROJECT_HITTINGSETLS_H
