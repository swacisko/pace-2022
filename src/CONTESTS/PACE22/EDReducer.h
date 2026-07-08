//
// Created by sylwe on 07/08/2025.
//

#ifndef ELDREDUCER_H
#define ELDREDUCER_H

#include "Makros.h"
#include "CONTESTS/PACE22/Config.h"

class EDReducer {
public:
    explicit EDReducer( int NN, Config c )
    : N(NN), inS(N), inU(N), inU1(N), inW(N), was(N), helper(N), marked(N), marked2(N), cnf(c) {
        temp.reserve(N);
        temp2.reserve(N);
        cnt = VI(N,0);
    }

    VI reduce(VVI V0);


    VI inf_rules_1, inf_rules_2;

    /**
     * If true, then edges will be added on the fly...
     */
    bool apply_type1_constraints_on_the_fly = false;

    int inf_rules_1_added = 0;
    int inf_rules_2_created = 0;

// private:

    bool write_logs = false;


    int N;
    VVI V;
    Config cnf;

    VB inS, inU, inU1, inW, was, helper, marked, marked2;
    VI temp, temp2,S,U,U1,W;

    VI cnt;

    /**
     * Checks whether node u can be safely added to the solution.
     * Additionally creates constraints that can be used if it cannot.
     */
    bool considerNode(int v);

    /**
     *  For node u considers all w \in N(u) \setminus W and finds all nodes x
     *  for which N(u) \setminus W \subseteq N(x).
     *
     *  Uses the provided [marked] bitvector to mark nodes. This is provided to enable easy implementation of the
     *  double-ED rule.
     */
    void markDominationNodes(VB& marked, bool check_double_ed = true);

    /**
     * Finds, using approaches marked in the [cnf], all the nodes that can be moved to U.
     * Those nodes must be marked earlier in the [marked] array using the [markDominationNodes] function.
     */
    VI findNodesToMoveToU();

    /**
     * Performs the next step.
     * Returns 1 if an ext-dominator was found.
     * Returns 0 if the ext-dominator was not found, but the algorithm did not terminate
     * Returns -1 if the ext-dominator was not found and the algorithm did terminate.
     */
    int nextStep();

    /**
     * Checks for given marked nodes in the [marked] bitvector, whether an ext-dominator exists.
     */
    bool existsExtDominator(VB & marked);

    /**
     * Moves node u from N(W) to U.
     * Must be u in N(W)
     */
    void moveToU(int u);

    /**
     * Moves node u from N(U) to S.
     * Must be u \in N(W)
     */
    void moveToS(int u);

    /**
     * Clears all arrays to prepare for checking next node.
     */
    void clearAll();


    void checkEmptyArraysAssertions(bool check_marked, bool check_marked2);
};

#endif //ELDREDUCER_H
