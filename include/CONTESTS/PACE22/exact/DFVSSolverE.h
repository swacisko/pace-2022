//
// Created by sylwester on 1/10/22.
//

#ifndef ALGORITHMSPROJECT_DFVSSOLVERE_H
#define ALGORITHMSPROJECT_DFVSSOLVERE_H

#include <utils/RandomNumberGenerators.h>
#include "Makros.h"
#include "CONTESTS/PACE22/Config.h"

class DFVSSolverE{
public:

    DFVSSolverE(VVI * orig, Config c);

    VI solveForInputGraph( VVI V );

    VI solve(VVI V, VVI revV, VI partial_dfvs, int partial_dfvs_size, int rec_depth, int &lb_total, int &ub_total);

    VI solve2(VVI V, VVI revV, int partial_dfvs_size, int rec_depth, int &lb_total, int &ub_total);

    VI solveForStronglyConnected(VVI V, VVI revV, VI partial_dfvs, int partial_dfvs_size, int rec_depth, int& lb_total, int &ub_total);

    VI solveForStronglyConnectedIterativeHittingSet2(VVI V, VI heur_dfvs,  bool use_exact_algorithm = true, bool return_heuristic_solution = false);

    int getBranchingNode(VVI &V, VVI &revV, int rec_depth);

    int getOptimalSolutionSize();
    void setKnownUpperBound( int ub ){ known_upper_bound = ub; }
    void setFindOnlySizeOfOptimalSolution( bool b ){ find_only_size_of_optimal_solution = b; }

    VLL hashes;
    UniformIntGenerator rnd;

    bool write_logs = false;

    Config cnf;

private:

    //***************//***************//***************//***************//***************
    int known_upper_bound = 1e9;

    bool find_only_size_of_optimal_solution = false;
    int optimal_solution_size = -1;

    //***************//***************//***************//***************//***************

    VVI *origV;

    void indent(int d);

    VVI vcs;
    VI vcs_last_rd;

    bool hasRedundantNodes(VI partial_dfvs);

    bool checkLBUB( VVI & V, VVI & revV, int rec_depth, int partial_size, int & ub, VI & vc_dfvs );

    int getLBCliqueCover( VVI & V, int rec_depth);
};

#endif //ALGORITHMSPROJECT_DFVSSOLVERE_H
