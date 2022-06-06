//
// Created by sylwester on 1/15/22.
//

#ifndef ALGORITHMSPROJECT_SALS2_H
#define ALGORITHMSPROJECT_SALS2_H


#include <utils/RandomNumberGenerators.h>
#include <CONTESTS/PACE22/Config.h>
#include "Makros.h"

class SALS2{
public:
    SALS2( VVI V, VI dfvs, Config c, bool use_fast_updates = false);

    void initialize(VI dfvs);

    VI localSearch(double T0, double alpha, int maxMvt, int maxFail, int max_iters = 1e9);


    const bool use_fast_updates = false;
    int max_iters = 1e9;

    bool cool_down_quickly_until_improved = true;

    static double findInitialTemperature( VVI V, VI& dfvs, Config c, int max_mvt, bool use_fast_updates = true );

    bool find_init_temperature_mode = false;

    Config cnf;

private:
    int N;
    VVI V, revV;

    VI in_order;

    VPII temp_order;

    VI node_on_pos;
    const int SPARSITY = 100;

    VI conflicts;

    double T;
    int nbMvt, nbFail;

    int dfvs_size;
    VI best_dfvs;

    double last_improvement_T = 2.0;

    VPII delta;
    VB is_valid_delta;

    VPII insertion_point;

    VB helper;


    UniformIntGenerator rnd;
    UniformDoubleGenerator rnd_d;

    static unsigned long long x, y, z;
    unsigned long long xorshf96();

    VLL random_ints;
    int random_int_index;
    int getRandomInt(int mod);

    vector<double> random_doubles;
    int random_double_index;
    double getRandomDouble();

    void preprocessRandoms();


    VI getCurrentDFVS();

    int getNodeForMove();

    void sparsifyOrder();

    void insertIntoOrder(int v, int pos);

    void shiftRight(int pos);
    void shiftLeft(int pos);

    void removeFromOrder(int v);

    void updateValidDelta( VI & conflicts );

    void updateNode(int v);

    int evaluateNode(int v, int before);

    void applyMove(int v, int before);

};

#endif //ALGORITHMSPROJECT_SALS2_H
