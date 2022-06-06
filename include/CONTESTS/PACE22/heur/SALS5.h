//
// Created by sylwester on 4/21/22.
//

#ifndef ALGORITHMSPROJECT_SALS5_H
#define ALGORITHMSPROJECT_SALS5_H


#include <utils/RandomNumberGenerators.h>
#include <datastructures/RandomSelectionSet.h>
#include <CONTESTS/PACE22/Config.h>
#include "Makros.h"

class SALS5{
public:
    SALS5( VVI V, VI dfvs, Config c);

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

    VVI buckets;

    VI node_on_pos;
    const int SPARSITY = 50;

    VI conflicts;

    double T;
    int nbMvt, nbFail;

    int dfvs_size;
    VI best_dfvs;

    double last_improvement_T = 2.0;

    struct INS_PT{
        INS_PT( int uu, int c, bool b ){ u = uu; cnt = c; before = b; }
        int u, cnt;
        bool before;
    };
    vector<vector<INS_PT>> delta;
    VB is_valid_delta;

    VB helper;

    UniformIntGenerator rnd;
    UniformDoubleGenerator rnd_d;

    struct PT{
        PT( int p, int uu, bool c ){ pos = p, u = uu, in_neigh = c; }
        int pos, u;
        bool in_neigh;
    };
    vector<PT> points;

    static unsigned long long x, y, z;
    unsigned long long xorshf96();

    VLL random_ints; // just to spped up random number generation
    int random_int_index;
    int getRandomInt(int mod);

    vector<double> random_doubles; // just to speed up random number generation
    int random_double_index;
    double getRandomDouble();


public:
    void setUseMetropolisSa(bool useMetropolisSa);

private:
    static const LL INVERSE_PRS = 16ll * 9 * 25 * 7 * 11 * 13;


    VI getCurrentDFVS();

    int getNodeForMove();

    void sparsifyOrder();

    void insertIntoOrder(int v, int pos);

    void shiftRight(int pos);
    void shiftLeft(int pos);

    void removeFromOrder(int v);

    void updateNode(int v);

    void applyMove(int v, int pos);
};


#endif //ALGORITHMSPROJECT_SALS5_H
