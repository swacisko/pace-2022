//
// Created by sylwester on 3/25/22.
//

#ifndef ALGORITHMSPROJECT_SALS3_H
#define ALGORITHMSPROJECT_SALS3_H



#include <utils/RandomNumberGenerators.h>
#include <datastructures/RandomSelectionSet.h>
#include <CONTESTS/PACE22/Config.h>
#include "Makros.h"

class SALS3{
public:
    SALS3( VVI V, VI dfvs, Config c);

    void initialize(VI dfvs);

    VI localSearch(double T0, double alpha, int maxMvt, int maxFail, int max_iters = 1e9);

    int max_iters = 1e9;

    bool cool_down_quickly_until_improved = true;

    static double findInitialTemperature( VVI V, VI& dfvs, Config cnf, int max_mvt);

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

    int total_backgoing_arcs = 0;

    VI backgoing_arcs_for_node;


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

    void moveRandomNodeToOrder();

    int getNodeToRemoveFromOrder();

    RandomSelectionSet prs_dfvs;

    bool use_metropolis_sa = true;
public:
    void setUseMetropolisSa(bool useMetropolisSa);

private:
    static const LL INVERSE_PRS = 16ll * 9 * 25 * 7 * 11 * 13;
    void updateDFVSNodePrsProbabilities(int v );
    function<LL(LL)> score_fun_dfvs = [&](LL x){ return x*x; };

    RandomSelectionSet prs_order;
    function<LL(LL)> score_fun_order = [&](LL x){ return x * x; };

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

#endif //ALGORITHMSPROJECT_SALS3_H
