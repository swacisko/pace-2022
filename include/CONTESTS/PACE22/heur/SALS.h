//
// Created by sylwester on 1/5/22.
//

#ifndef ALGORITHMSPROJECT_SALS_H
#define ALGORITHMSPROJECT_SALS_H

#include <utils/RandomNumberGenerators.h>
#include "Makros.h"

class SALS{
public:
    SALS( VVI V, VI dfvs );

    VI localSearch(double T0, double alpha, int maxMvt, int maxFail, int max_iters = 1e9);

    static double findInitialTemperature( VVI V, VI& dfvs, int max_mvt );

    bool cool_down_quickly_until_improved = true;

    bool use_uniform_node_selection = true;

    bool write_logs = true;

private:
    int N;
    VVI V, revV;

    VI order;
    VI in_order;

    VI conflicts;

    VLL dfvs_pref_sum;
    double SCALE = 1;
    VD probabs;


    double T;
    int nbMvt, nbFail;

    VI dfvs;
    VI best_dfvs;

    int dfvs_check_marker = -1;

    double last_improvement_T = 2.0;

    VPII delta;
    VB is_valid_delta;

    VB helper;

    UniformIntGenerator rnd;
    UniformDoubleGenerator rnd_d;

    int getNodeForMove();

    void updateOrder();
    void updateOrder(int pos, VI & conflicts);

    void updateValidDelta( VI & conflicts );

    void updateNode(int v);

    int evaluateNode(int v, int before);

    void applyMove(int v, int before);

    void createDfvsPrefSum();
};

#endif //ALGORITHMSPROJECT_SALS_H
