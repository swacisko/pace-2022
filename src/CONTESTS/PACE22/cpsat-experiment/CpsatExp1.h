//
// Created by sylwe on 08/04/2026.
//

#ifndef DIVERSES_CPSATEXP1_H
#define DIVERSES_CPSATEXP1_H
#include "ExpConfig.h"
#include "Makros.h"

class ExpData {
public:
    class IterationEntry {
    public:
        int res_size;
        bool res_valid;
        bool res_optimal = false;

        int distinct_arcs_in_all_cycles;
        int unhit_cycle_enumeration_time_millis;
        int unhit_graph_size;
        int unhit_graph_dfvs_size;
        map<int,int> cycles_of_length;
        int total_cycles;
    };

    vector<IterationEntry> iterations;
    static bool foundValidResult(vector<IterationEntry> & entries);

};

class CpsatExp1 {
public:
    /**
     * Iterative Hitting-Set approach - the most straightforward type.
     * For each cycle length L, starting from 1, considers all cycles of length <= L, then finds HS of those cycles.
     * If the found HS is not a FVS of V, then increases L and repeats.
     */
    static ExpData solveHS1(VVI V, ExpConfig cnf);

    /**
   * Iterative HS.
   * Starting with L = 2 and res = {}, finds all cycles in the graph G[V \ res] of length <= L.
   * If there is just a small number of such cycles, increases L and repeats.
   * Then finds HS of the set of all cycles found so far.
   */
    static ExpData solveIHS(VVI V, ExpConfig cnf, int cycle_enumeration_type = 0);

    /**
     * Solves the problem using the MTZ formulation. It might additionally find some cycles
     * and add it to speed up propagation.
     * If auxiliary_cycles_mode == 0, then no cycles are added, if 1, then all cycles of length <= 3 are added,
     * if 2, cycles are found using the IHS approach.
     */
    static ExpData solveMTZ(VVI V, ExpConfig cnf, int auxiliary_cycles_mode = 0);

    /**
     * Runs the DiVerSeS solver to find out high-quality baseline solutions.
     */
    static ExpData solveDiVerSeS(VVI V, ExpConfig cnf);


};


#endif //DIVERSES_CPSATEXP1_H
