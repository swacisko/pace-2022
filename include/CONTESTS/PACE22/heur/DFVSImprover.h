//
// Created by sylwester on 1/7/22.
//

#ifndef ALGORITHMSPROJECT_DFVSIMPROVER_H
#define ALGORITHMSPROJECT_DFVSIMPROVER_H

#include <CONTESTS/PACE22/Config.h>
#include "Makros.h"
#include "SALS2.h"

class DFVSImprover{
public:

    DFVSImprover(VVI V, Config c);

    VI improve( VI dfvs, const int iterations );

    VI improveBySALSForSparseGraphs( VI dfvs, int iterations, double alpha, double sals_T0,
                                     double sals_alpha, int sals_maxmvt, int sals_maxfail);

    VI improveBySALSForSparseGraphsDensifier( VI dfvs, int iterations, double init_density = 0.01, double step = 0.005,
                                         double sals_T0 = 0.25, double sals_alpha = 0.99, int sals_maxmvt_factor = 100,
                                         int sals_maxfail = 30);

    VI fullImprovement(VI dfvs, int origV_N, int origV_edges, const bool use_vc_improvement = true);

    VI mergeNodesUntilDensityReached( VVI & V, double density_threshold );

    VI getOutOfLocalOptimum(VI dfvs, double alpha);

    Config cnf;
    VVI V;

    VI best_dfvs;
};

#endif //ALGORITHMSPROJECT_DFVSIMPROVER_H
