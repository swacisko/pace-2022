//
// Created by sylwester on 12/20/21.
//

#ifndef ALGORITHMSPROJECT_UTILS_H
#define ALGORITHMSPROJECT_UTILS_H

#include "Makros.h"
#include "CollectionOperators.h"
#include "Config.h"

namespace Utils{

    extern bool addEdge(VVI & V, int a, int b, bool just_push_back = false);

    extern void addEdges(VVI & V, VPII & edges, VB & helper);

    extern void addEdges(VVI & V, VVI & revV, VPII arcs, VB & helper);

    extern void removeEdge( VVI & V, int a, int b );

    extern void removeNode( VVI & V, VVI & revV, int a, VB & helper);

    extern void removeNodes( VVI & V, VVI & revV, VI & nodes, VB & helper);

    void removeEdges( VVI & V, VPII & edges, VB & helper );

    extern void removeEdges( VVI & V, VVI & revV, VPII arcs, VB & helper );

    extern void merge( VVI & V, VVI & revV, int a, VB & helper );

    extern void partialMerge( VVI & V, VVI & revV, VVI & nonpiV, VVI & revnonpiV, VVI & piV, int v, VB & helper );

    extern bool hasLoop( VVI & V,  int a );

    extern VPII getPIEdges( VVI & V, VVI & revV, int a, VB & helper);

    extern VPII getAllPIEdges( VVI & V, VVI & revV, VB & helper);

    extern VI getNonPISuccessors( VVI & V, VVI & revV, int a, VB & helper );

    extern VI getNonPIPredecessors( VVI & V, VVI & revV, int a, VB & helper );

    extern bool isCorresponding( VVI V, VVI revV );

    extern bool hasCycle( VVI & V );

    extern bool hasPathFromTo(VVI & V, VVI & revV, VB & in_V, VI sources, VI ends, VB & helper1, VB & helper2);

    extern bool hasCycle( VVI & V, VI & A, VB& was, VB& on_path, VB & helper );

    extern bool isFVS(VVI & V, VI S);
    extern bool isFVS2(VVI & V, VI & S);

    extern void writeRemainingGraph(VVI & V);

    extern void writeNeighborhood(VVI & V, VVI & revV, int v);

    VI getRedundantNodes(VVI V, VI dfvs, int max_check_length, const bool minimize, const bool return_on_first = false);

    VI removeRedundantNodes( VI & dfvs, VI & redundant, VB & helper );

    VI findAndRemoveRedundantNodes( VVI & V, VI & dfvs, int max_check_length = 1e9, const bool minimize = true );

    tuple<VVI, VPII> transformDFAStoDFVS(VVI & V);

    int countBackEdgesForNodeOrder(VVI & V, VI & ord);

    VVI readGraph( istream & str );

    VVI getCycleAdjacencyGraph( VVI & cycles );

    VVI getSmallLowerBoundCycles(VVI V);

    VI getLowerBoundCycles( VVI cycV, int milliseconds = 1e9 );

    VI getLowerBoundByVCOnPIGraph(VVI &V, VVI &revV);

    VI getUpperBoundByVCOnSuperPIGraph(VVI V, int milliseconds, bool use_reductions = true, bool use_extensive_reductions = false);

    bool isPIGraph(VVI V);

    bool isPIGraph(VVI &V, VVI & revV, VB & helper);

    VVI getUnderlyingPIGraph(VVI V);
    VVI getUnderlyingPIGraph(VVI &V, VVI & revV);

    VVI getSuperPIGraph(VVI V);

    VVI getNonPIGraph(VVI V);
    VVI getNonPIGraph(VVI &V, VVI & revV);

    void contractNodeToNode(VVI & V, VVI & revV, int b, int c, VB & helper );

    bool containsSimpleCycleWithNode(VVI & V, VVI & revV, VVI& nonpiV, VVI & revnonpiV,
                                      VB & in_V, VI & marker, VB& is_end, int v, const int max_time_millis = 100,
                                      int max_branch_depth = 1e9);

    VVI getAllSimpleCycles( VVI & V, VVI & revV, VVI& nonpiV, VVI & revnonpiV, VI A,
                            VB & in_V, VI & marker, VB& is_end, int max_cycle_length = 1e9 );

    VVI getAllSimpleCycles2( VVI & V, VVI & revV, VVI& nonpiV, VVI & revnonpiV, VI A,
                            VB & in_V, VI & marker, VB& is_end, int max_cycle_length = 1e9, int millis = 1e9 );

    VVI getAllSimpleCycles3( VVI V, int max_cycle_length = 1e9, int millis = 1e9, bool random_A = true );

    vector<Triple<int>> getLengthOfShortestInducedCycleWithArc( VVI V, int max_length, int millis = 1e9 );

    VI getIntersectionOfAllInducedNonpiCyclesWithNode( VVI & V, VVI & revV, VVI& nonpiV, VVI & revnonpiV, int v,
                              VI & marker, VB& is_end, VI & cnt_marker, int max_cycle_length = 1e9, int millis = 1e9 );

    VI getAllNoninclusionNodesForNode( VVI & V, VVI & revV, VVI & nonpiV, VVI& revnonpiV,
                                       VB & in_V, VI & marker, VB& is_end, VB & in_Np, int a, VI bs,
                                       int max_time_millis_total, int max_time_millis_per_node,
                                       int max_branch_depth = 1e9);

    bool isPiNode( VVI & V, VVI & revV, VVI & piV, int a );

    VI findMinHittingSet(VVI & sets, VI upper_bound_solution = {}, int lower_bound = 0 );

    VI hsImprovementLS2( VVI sets, VI hs, Config cnf, int iters, int lower_bound = -1, bool deviate_perm = true,
                         int persistent_search_iterations_left = 0);

    double getPieEdgesPercentage(VVI & V);
    double getPieEdgesPercentage(VVI & V, VVI & revV);

    LL countInducedCycles(VVI & V, int max_length, int max_millis = 1e9);

    LL getSetHash( int N, VI & s, int seed = 8'592'374 );

    int countPiEdges(VVI & V);

    VI getMinDFVSDefault( VVI & V );

    void emergencyExit( VVI & V, VI & dfvs, const bool write_communicate = true );

    bool isPartiallyDominated( VVI & V, VVI & revV, VVI & nonpiV, VVI & revnonpiV, int b, int c, VB & helper );

    bool isFullyDominated( VVI & V, VVI & revV, int b, int c, VB & helper1, VB & helper2 );
}

#endif //ALGORITHMSPROJECT_UTILS_H
