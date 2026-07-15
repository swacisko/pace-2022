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

    /**
     * Uses techniques specified in the [cnf] object to reduce/apply changes to the graph.
     *
     * The most standard approach is node removal - checking initS = {v} for all nodes v.
     *
     * Another approach is edge removal - checking initS = {u,v} for edges in the graph. This requires complicated
     * solution lifting though.
     *
     * It is also possible to use edge insertion - for given initial single-nide set S = {v},
     * whenever a node x is moved to the set U, we can add edge {v,x} to the graph.
     *
     * We can also use extended edge insertion - for some set of candidates C (usually N^3(v)),
     * we check whether initS = {v,x} will yield true, for x \in C. If so, then we can add edge {v,x} to the graph.
     */
    VI reduce(VVI V0);

    /**
     * True, if recent call to [reduce] made any changes to the processed graph.
     * @return
     */
    bool madeChangesInLastReduce();

    /**
     * Resets all techniques that might be used in [reduce].
     * After this, used techniques need to be set manually, via this->cnf object.
     * Without setting it manually, the [reduce] function will do nothing, as all techniques are disabled.
     */
    void resetAllUsedTechniques();


    VI inf_rules_1, inf_rules_2;
    VPII all_inf_rules_1_found, all_inf_rules_2_found;

    int last_reduce_inf_rules_1_added = 0;
    int last_reduce_inf_rules_2_created = 0;

    int last_reduce_edges_removed = 0;
    int last_reduce_nodes_removed = 0;

    VVI getV(){return V;}
    Config cnf;

    bool write_logs = false;


private:


    int N;
    VVI V;

    VB inS, inU, inU1, inW, was, helper, marked, marked2;
    VI temp, temp2,S,U,U1,W;
    VI cnt;

    bool check_double_ed = false;


    /**
     * Counts and returns |N(u) \setminus W|
     */
    int getNonWNeighborhoodSize(int u);
    bool hasNonWIntersectionAtMost( int u, int val ){ return getNonWNeighborhoodSize(u) <= val; }

    /**
     * Coutns and returns |N(u) \cap S|
     */
    int getSIntersection(int u);

    /**
     * Checks whether node u can be safely added to the solution.
     * Additionally creates constraints that can be used if it cannot.
     */
    bool consider(VI initS);

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
     * Removes from U1 all nodes that have more than one neighbor in S.
     */
    void updateU1();

    /**
     * Clears all arrays to prepare for checking next node.
     */
    void clearAllForConsider();


    void checkEmptyArraysAssertions(bool check_marked, bool check_marked2);
};

#endif //ELDREDUCER_H
