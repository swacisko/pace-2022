//
// Created by sylwe on 08/04/2026.
//

#ifndef DIVERSES_CPSATEXP1_H
#define DIVERSES_CPSATEXP1_H
#include <ortools/util/bitset.h>

#include "ExpConfig.h"
#include "GraphInducer.h"
#include "Makros.h"
#include "Stopwatch.h"
#include "ortools/sat/cp_model.h"
using namespace operations_research::sat;

class ExpData {
public:
    class IterationEntry {
    public:
        IterationEntry() {}

        int full_sol_size = -1;
        int hs_size_before_impr = -1;
        int hs_size_after_impr = -1;
        bool hs_valid_fvs = false;
        // bool cycle_hs_valid_dfvs = false;
        bool improved_best_res = false;
        bool res_optimal = false;
        int res_lower_bound = -1; // for mtz only

        int distinct_arcs_in_all_cycles = -1;
        int unhit_cycle_enumeration_time_millis = -1;
        int hs_greedy_time = -1;
        PII unhit_graph_sizes = {-1,-1};
        int unhit_graph_greedy_dfvs_size = 0;
        int unhit_graph_greedy_dfvs_time = -1;
        map<int,int> cycles_of_length;
        int total_cycles = -1;
        int new_cycles_added = -1;
        int new_cycles_found = -1;
        int time_since_start_millis = -1;
        int iteration_time = -1;
        int max_cycle_length = -1;
        int best_result_so_far = -1;
    };

    int N,M;
    double avg_deg, max_deg, min_deg;
    VI cyc_of_length, ind_cyc_of_length;

    vector<IterationEntry> iterations;
    static bool foundValidHSResult(vector<IterationEntry> & entries);

    vector<map<string,string>> getIterationEntries();
    void updateBestResultSoFar();

    void writeToFile(string filename);

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
    static VI getUnhitGraphGreedyFVS(VVI & V, VI & S);

    /**
     * Iterative Hitting-Set approach - the most straightforward type.
     * For each cycle length L, starting from 1, considers all cycles of length <= L, then finds HS of those cycles.
     * If the found HS is not a FVS of V, then increases L and repeats.
     */
    static ExpData solveHS(VVI V, ExpConfig cnf);

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
    static void addInitialSolutionHint(CpModelBuilder &model, vector<BoolVar> & nodes, VI & init_sol, VI &prev_res, ExpConfig cnf);
    static void addMaxHammingDstConstraint(CpModelBuilder &model, vector<BoolVar> & nodes, VI & init_sol, ExpConfig cnf);
    static tuple<VI,CpSolverStatus,VI> rerunModelUntilFeasibleOrTle(VVI & V, CpModelProto &proto, vector<BoolVar> & nodes, CpSolverResponse & response,
        Stopwatch & timer, string timer_option, VI & init_sol, int init_time, ExpConfig cnf);
    static tuple<VI,CpSolverStatus, VI> solveCpsatForCycles( VVI & V, VVI & cycles, VI &prev_res, VI & init_sol, Stopwatch & timer, string timer_option, ExpConfig& cnf );
    static VI getUnhitCyclesHSGreedy(VVI & cycles);
    static bool isHS(VVI & cycles, VI & S);

};


#endif //DIVERSES_CPSATEXP1_H
