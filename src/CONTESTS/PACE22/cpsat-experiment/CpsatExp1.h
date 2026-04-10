//
// Created by sylwe on 08/04/2026.
//

#ifndef DIVERSES_CPSATEXP1_H
#define DIVERSES_CPSATEXP1_H
#include <ortools/util/bitset.h>

#include "ExpConfig.h"
#include "GraphInducer.h"
#include "Makros.h"
#include "ortools/sat/cp_model.h"
using namespace operations_research::sat;

class ExpData {
public:
    class IterationEntry {
    public:
        IterationEntry() {}
        // IterationEntry( int _res_size_before_impr, int _res_size_after_impr, bool _res_valid, bool _res_optimal,
        //     int _distinct_arcs_in_all_cycles, int _unhit_cycle_enumeration_time_millis,
        //     PII _unhit_graph_sizes, int _unhit_graph_dfvs_size, map<int,int> _cycles_of_length,
        //     int _total_cycles, int _time_since_start_millis, bool _improved_res, int _new_cycles_added) {
        //     res_size_before_impr = _res_size_before_impr;
        //     res_size_after_impr = _res_size_after_impr;
        //     res_valid = _res_valid;
        //     res_optimal = _res_optimal;
        //     distinct_arcs_in_all_cycles = _distinct_arcs_in_all_cycles;
        //     unhit_cycle_enumeration_time_millis = _unhit_cycle_enumeration_time_millis;
        //     unhit_graph_sizes = _unhit_graph_sizes;
        //     unhit_graph_dfvs_size = _unhit_graph_dfvs_size;
        //     cycles_of_length = _cycles_of_length;
        //     total_cycles = _total_cycles;
        //     time_since_start_millis = _time_since_start_millis;
        //     improved_res = _improved_res;
        //     new_cycles_added = _new_cycles_added;
        // }

        int res_size_before_impr = -1;
        int res_size_after_impr = -1;
        bool res_valid = false;
        bool improved_res = false;
        bool res_optimal = false;

        int distinct_arcs_in_all_cycles = -1;
        int unhit_cycle_enumeration_time_millis = -1;
        int hs_greedy_time = -1;
        PII unhit_graph_sizes = {-1,-1};
        int unhit_graph_dfvs_size = -1;
        map<int,int> cycles_of_length;
        int total_cycles = -1;
        int new_cycles_added = -1;
        int time_since_start_millis = -1;
    };

    int N,M;
    double avg_deg, max_deg, min_deg;
    VI cyc_of_length, ind_cyc_of_length;

    vector<IterationEntry> iterations;
    static bool foundValidResult(vector<IterationEntry> & entries);

};

class CpsatExp1 {
public:

    /**
     *  Finds and returns all unhit chordless cycles of length at most max_l in the graph G[V \ S].
     *  If enumeration_option == 1, then Utils::getAllSimpleCycles3 is used.
     *  If 2, then w new method is used, based on dfs tree traversal.
     */
    static VVI getUnhitChordlessCycles(VVI &V, VI & S, int max_l, int max_millis, int enumeration_option = 1);

    /**
     * Creates a graph H = G[V\S], then removes from it all arcs that belong to different strongly connected components.
     * Then contracts all nodes with in-degree 1 or out-degree 1. Returns with it an associated map,
     * that lets reconstruct the
     */
    static InducedGraph getUnhitGraph(VVI & V, VI & S);
    static PII getUnhitGraphSizes(VVI & V, VI & S);

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
   * All cycles found by this function will be stored in [cycles] vector
   */
    static ExpData solveIHS(VVI V, ExpConfig cnf);
    static ExpData solveIHS(VVI V, ExpConfig cnf, VVI & cycles, VI & res);

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


    static ExpData solve(VVI V, ExpConfig cnf);

    static void addCycleConstraints(CpModelBuilder &model, VVI & cycles, vector<BoolVar> & nodes);
    static VI getUnhitCyclesHSGreedy(VVI & cycles);
    static bool isHS(VVI & cycles, VI & S);

};


#endif //DIVERSES_CPSATEXP1_H
