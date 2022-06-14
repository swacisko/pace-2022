//
// Created by sylwester on 1/10/22.
//

#include <graphs/scc/StronglyConnectedComponents.h>
#include <graphs/GraphInducer.h>
#include <CONTESTS/PACE22/Reducer.h>
#include <graphs/GraphUtils.h>
#include <CONTESTS/PACE22/heur/DFVSSolverH.h>
#include <utils/TimeMeasurer.h>
#include <graphs/VertexCover/VCUtils.h>
#include <combinatorics/CombinatoricUtils.h>
#include <utils/Stopwatch.h>
#include <CONTESTS/PACE22/heur/DFVSImprover.h>
#include <CONTESTS/PACE22/exact/CliqueCoverLS.h>
#include <filesystem>
#include "CONTESTS/PACE22/exact/DFVSSolverE.h"

DFVSSolverE::DFVSSolverE( VVI* orig, Config c) : rnd(0,1'000'000'000ll * 1'000'000'000) {
    cnf = c;
    origV = orig;

    hashes = VLL(orig->size());
    for(int i=0; i<orig->size(); i++) hashes[i] = rnd.rand();
}


VI DFVSSolverE::solve2(VVI V, VVI revV, int partial_dfvs_size, int rec_depth, int& lb_total, int& ub_total){

    ub_total = min( ub_total, int( partial_dfvs_size + V.size() ) );

    vector<DFVSReduction*> reductions;
    {
        Reducer red(V,cnf);
        red.cnf.write_logs = false;

        reductions = red.reduce(revV);
        if(write_logs){ indent(rec_depth); clog << "Depth: " << rec_depth << ", reductions done" << endl; }
        revV = red.revV;
        V = red.V;
    }

    if( GraphUtils::countEdges(V, true) == 0 ){
        VI dfvs;
        Reducer::liftSolution(V.size(), dfvs, reductions);
        return dfvs;
    }

    int reductions_size_diff = Reducer::getReductionsSizeDiff(reductions);
    partial_dfvs_size += reductions_size_diff;


    VI vc_dfvs;
    bool needs_further_branching = true;

    bool check_lower_bound = true; // here can be a condition, such as "check LB only if rec_depth % 3 == 2"
//    constexpr bool use_clique_cover_lb_instead_of_wgyc = false;
    const bool use_clique_cover_lb_instead_of_wgyc = ( !filesystem::exists("vc_solver") );
    bool push_vc = (!use_clique_cover_lb_instead_of_wgyc);

    int current_lb = 0;

    if( check_lower_bound ){
        if(use_clique_cover_lb_instead_of_wgyc){
            current_lb = getLBCliqueCover(V, rec_depth);
            needs_further_branching = (partial_dfvs_size + current_lb < ub_total);
        }
        else{
            needs_further_branching = checkLBUB(V, revV, rec_depth, partial_dfvs_size, ub_total, vc_dfvs); // #TEST
        }
    }

    if( !needs_further_branching ){
        VI temp = VI(V.size(),0); iota(ALL(temp),0);
        if(write_logs){ indent(rec_depth); clog << "LOWER BOUND WORKS!" << endl; }
        return temp;
    }
    else{
        if( Utils::isFVS(V, vc_dfvs ) ){ // #TEST - only to use if no subgraph inducing was done
            VI dfvs = vc_dfvs;
            Reducer::liftSolution(V.size(), dfvs, reductions);
            return dfvs;
        }
    }

    if(cnf.write_logs){
        indent(rec_depth);
        clog << "Partial_size + current_lb: " << partial_dfvs_size << " + "  << current_lb << " = " <<
             partial_dfvs_size +current_lb << ", UB: " << ub_total << endl;
    }

    // pushing latest vertex cover of piV
    if(push_vc){
        vcs.push_back(vc_dfvs);
        vcs_last_rd.push_back(rec_depth);
    }else{
//        if(cnf.write_logs){ indent(rec_depth); clog << "Partial_size: " << partial_dfvs_size << endl; }
    }

    VI res = solveForStronglyConnected(V, revV, {}, partial_dfvs_size, rec_depth, lb_total, ub_total);

    Reducer::liftSolution(V.size(), res, reductions);


    if(push_vc){
        vcs.pop_back();
        vcs_last_rd.pop_back();
    }
    ub_total = min( ub_total, int( partial_dfvs_size + res.size() ) );

    if(write_logs){ indent(rec_depth); clog << "Upper bound: " << ub_total << endl; }

    return res;
}


VI DFVSSolverE::solve(VVI V, VVI revV, VI partial_dfvs, int partial_dfvs_size, int rec_depth, int& lb_total, int& ub_total){
    return solve2(V, revV, partial_dfvs_size, rec_depth, lb_total, ub_total);

    ub_total = min( ub_total, int( partial_dfvs_size + V.size() ) );

    VI res;
    { // it seems that reductions done here are a time-consuming part done by DFVSSolverE
        Reducer red(V,cnf);
        red.cnf.disableAllConditionalReductions();
        red.cnf.disableAllRecursiveReductions();

        if(write_logs){
            indent(rec_depth);
            clog << "Depth: " << rec_depth << ", using reductions..." << flush;
        }

        auto reductions = red.reduce(revV);
        assert(reductions.size() <= 1); // this will hold if no conditional reductions were done
        VI red_dfvs;
        if(!reductions.empty()){
            KernelizedNodesReduction * knr = (KernelizedNodesReduction*) reductions[0];
            red_dfvs = knr->getKer();
            Reducer::clearReductionObjects(reductions);
        }

        if(write_logs) clog << "reductions done" << endl;
        res = red_dfvs;
        revV = red.revV;
        V = red.V;
    }

    assert(cnf.reducer_use_twins_merge == false);

    if( (rec_depth % 4) == 3) {
        if (hasRedundantNodes(partial_dfvs + res)) {
            clog << "REDUNDANT NODES IN SOLVE!!" << endl;
            if (write_logs) {
                indent(rec_depth);
                clog << "REDUNDANT NODES after reduce()" << endl;
                DEBUG(partial_dfvs);
            }
            VI temp = VI(V.size(), 0);
            iota(ALL(temp), 0);
            return temp;
        }
    }

    VI vc_dfvs;
    bool needs_further_branching = true;

    // using solve() instead of solve2() we need vc_solver
    bool check_lower_bound = true; // here can be a condition, such as "check LB only if rec_depth % 3 == 2"
    bool push_vc = check_lower_bound; // original
    if( check_lower_bound ){
        needs_further_branching = checkLBUB(V, revV, rec_depth, partial_dfvs_size + res.size(), ub_total, vc_dfvs);
    }

    if( !needs_further_branching ){
        VI temp = VI(V.size(),0);
        iota(ALL(temp),0);

        if(write_logs){
            indent(rec_depth);
            clog << "LOWER BOUND WORKS!" << endl;
        }
        return temp;
    }
    else{
        if( Utils::isFVS(V, vc_dfvs ) ){
            return res + vc_dfvs;
        }
    }

    const bool use_disjoint_comps = false;
    VVI comps;

    if(use_disjoint_comps) {
        TimeMeasurer::start("Finding SCC");
        StronglyConnectedComponents scc(V);
        scc.createStronglyConnectedComponents();
        comps = scc.getComponents();
        TimeMeasurer::stop("Finding SCC");
    }else{
        comps = VVI(1, VI(V.size()));
    }

    // pushing latest vertex cover of piV
    if(push_vc){
        vcs.push_back(vc_dfvs);
        vcs_last_rd.push_back(rec_depth);
    }

    if( comps.size() == 1 ){
        res += solveForStronglyConnected(V, revV, partial_dfvs + res, partial_dfvs_size + res.size(), rec_depth, lb_total, ub_total);
    }
    else{
        /**
         * #CAUTION!!
         * If using inducing, then we cannot use vc_dfvs being a DFVS, neither can we use hasRedundantNodes() without
         * properly remapping all nodes.
         */
        VI dfvs;
        for( VI & cmp : comps ){
            if( cmp.size() == 1 ) continue;
            if( cmp.size() == 2 ){
                res.push_back(cmp[0]);
                continue;
            }

            InducedGraph g = GraphInducer::induce(V,cmp);
            VVI grevV = GraphUtils::reverseGraph(g.V);

            VI dfvs_g = solveForStronglyConnected( g.V, grevV, partial_dfvs + res, partial_dfvs_size + res.size() + dfvs.size(), rec_depth, lb_total, ub_total );
            assert(Utils::isFVS(g.V, dfvs_g));

            for(int & d : dfvs_g) d = g.nodes[d];
            dfvs += dfvs_g;
        }
        res += dfvs;
    }

    if(write_logs){
        indent(rec_depth);
        clog << "Upper bound: " << ub_total << endl;
    }

    if(push_vc){
        vcs.pop_back();
        vcs_last_rd.pop_back();
    }
    ub_total = min( ub_total, int( partial_dfvs_size + res.size() ) );

    return res;
}

VI DFVSSolverE::solveForStronglyConnected(VVI V, VVI revV, VI partial_dfvs, int partial_dfvs_size, int rec_depth, int &lb_total, int& ub_total) {

    int branch_v = getBranchingNode(V, revV, rec_depth);

    if( branch_v == -1  ){
        DEBUG(vcs.size());

        VVI piV = Utils::getUnderlyingPIGraph(V);
        assert( Utils::isFVS(piV, vcs.back()) );
        assert( VCUtils::isVertexCover(piV, vcs.back()) );
        assert( Utils::isFVS(V, vcs.back()) );
        // with lower-bounding and taking a VC if it is a DFVS together with partial_dfvs, we should not enter here at all
        assert(false);
    }

    assert(branch_v >=0 && branch_v < V.size() && V[branch_v].size() > 0);

    bool best_initialized = false;
    VI best_dfvs;

    const int MERGE = 0, REMOVE = 1;
    VI order = {MERGE, REMOVE};

    for(int opt : order) {
        if (opt == REMOVE) { // remove branch_v from V
            VVI H = V;
            VB helper(H.size(), false);

            VVI revH = revV;
            Utils::removeNode(H, revH, branch_v, helper);

            if(write_logs){
                indent(rec_depth);
                clog << "Depth: " << rec_depth << ", branching on node " << branch_v << ", remove node" << endl;
            }

            // we only  branch on this node if it will not create any redundant nodes
            VI rem_dfvs = solve(H, revH, partial_dfvs + VI({branch_v}),  partial_dfvs_size + 1, rec_depth + 1, lb_total, ub_total);
            rem_dfvs.push_back(branch_v);
            assert(Utils::isFVS(V, rem_dfvs));

            if( !best_initialized || rem_dfvs.size() < best_dfvs.size() ){
                best_dfvs = rem_dfvs;
                best_initialized = true;
            }
        }

        if(opt == MERGE){ // merge node branch_v
            if(write_logs){
                indent(rec_depth);
                clog << "Depth: " << rec_depth << ", branching on node " << branch_v << ", merge node" << endl;
            }

            VVI H = V;
            VVI revH = revV;
            VB helper(H.size(), false);
            Utils::merge(H, revH, branch_v, helper);

            VI merge_dfvs = solve(H, revH, partial_dfvs, partial_dfvs_size, rec_depth+1, lb_total, ub_total );

            assert(Utils::isFVS(V, merge_dfvs));
            if(!best_initialized || merge_dfvs.size() < best_dfvs.size()){
                best_dfvs = merge_dfvs;
                best_initialized = true;
            }
        }
    }

    assert(Utils::isFVS(V, best_dfvs));

    return best_dfvs;
}


int DFVSSolverE::getBranchingNode(VVI &V, VVI &revV, int rec_depth) {
    TimeMeasurer::start("getBranchingNode()");

    bool use_fast_branching_node = false;

    if(use_fast_branching_node){
        LL best_val = 0, v = -1;
        VB helper(V.size(),false);

        for( int i=0; i<V.size(); i++ ){
            int pie_e = Utils::getPIEdges(V, revV, i, helper).size();
            assert( pie_e <= max(V[i].size(), revV[i].size()) );

            LL score = (V[i].size()-pie_e) * (revV[i].size()-pie_e);

            if( score > best_val ){
                best_val = score;
                v = i;
            }
        }
        assert(v != -1);
        TimeMeasurer::stop("getBranchingNode()");
        return v;
    }

    int N = V.size();
    VVI H;

//    const bool select_branching_node_from_nonpi = true; // original
    const bool select_branching_node_from_nonpi = (!Utils::isPIGraph(V));

    if(select_branching_node_from_nonpi){  // selecting a branching node form nonpiV \ vcs.back()
        H = Utils::getNonPIGraph(V);
        VVI revH = GraphUtils::reverseGraph(H);
        VB helper(N,false);

        const bool remove_vc_from_H = true;
        if( ( !vcs.empty() && vcs_last_rd.back() == rec_depth ) && remove_vc_from_H ){
            Utils::removeNodes(H, revH, vcs.back(), helper);

            if( !Utils::hasCycle(H) ){
                if(write_logs) clog << "Graph H has no cycles after removing vcs.back() !!" << endl;
                // if nonpi graph H has no cycles after removing latest vc, then this vc must be a DFVS of V

                assert(GraphUtils::isSimple(H));
                assert(GraphUtils::isSimple(V));

                assert( Utils::isFVS(H, vcs.back()) );
                assert( Utils::isFVS(V, vcs.back()) );
                assert(false);
            }
        }else{
            // we did not find VC on previous level, so we cannot remove it from H, hence we take H = V
            H = V;
        }
    }else{
        H = V; // #original
    }

    const bool select_node_from_H = true;

    if(select_node_from_H) {

        // since the reductions were done, then H is a union of strongly connected components
        const bool use_reds = true;
        if(use_reds){
            Reducer red(H, cnf);
            red.cnf.disableAllNonbasicReductions();
            red.cnf.disableAllConditionalReductions();
            red.cnf.disableAllRecursiveReductions();

            auto reductions = red.reduce();
            assert(Utils::isCorresponding(red.V, red.revV));
            VI red_dfvs = Reducer::convertKernelizedReductions(reductions);

            int E = GraphUtils::countEdges(red.V, true);
            if (E == 0 && red_dfvs.empty()) {
                TimeMeasurer::stop("getBranchingNode()");
                assert(!Utils::hasCycle(red.V));
                assert(!Utils::hasCycle(H));
                return -1;
            }

            if(E == 0 || !red_dfvs.empty()) {
                VVI revH = GraphUtils::reverseGraph(H);
                return *max_element(ALL(red_dfvs), [&](int a, int b) {
                    int shortcuts_a = H[a].size() * revH[a].size();
                    int shortcuts_b = H[b].size() * revH[b].size();

                    return shortcuts_a < shortcuts_b;
                });
            }

            H = red.V;
        }
    }

    VD tokens = DFVSSolverH::sinkhorn(H);

    double m = 1e9, v = -1;
    for( int i=0; i<H.size(); i++ ){
        if( !H[i].empty() ){
            if( tokens[i] < m ){
                m = tokens[i];
                v = i;
            }
        }
    }

    TimeMeasurer::stop("getBranchingNode()");
    return v; // we return the node that is in maximal number of DCU's
}

void DFVSSolverE::indent(int d) {
    for(int i=0; i<d; i++) clog << "\t";
}

bool DFVSSolverE::hasRedundantNodes(VI partial_dfvs) {
    TimeMeasurer::start("hasRedundantNodes()");
    VI red = Utils::getRedundantNodes(*origV, partial_dfvs, 1e9, true, true);
    TimeMeasurer::stop("hasRedundantNodes()");
    return red.size() > 0;
}

int DFVSSolverE::getLBCliqueCover(VVI &V, int rec_depth) {
    InducedGraph g = GraphInducer::induceByNonisolatedNodes(V);
    VVI pigV = Utils::getUnderlyingPIGraph(g.V);
    CliqueCoverLS clqc(pigV,cnf);
    clqc.cnf.write_logs = false;
    VVI cliques = clqc.coverLS(50);
    int local_lb = pigV.size() - cliques.size();

    return local_lb;
}


bool DFVSSolverE::checkLBUB(VVI &V, VVI & revV, int rec_depth, int partial_size, int &ub_total, VI & vc_dfvs) {

    int local_lb = 0;
    {
        TimeMeasurer::start("Lower bounding");
        vc_dfvs = Utils::getLowerBoundByVCOnPIGraph(V, revV);
        local_lb = vc_dfvs.size();

        const bool use_cycle_lb = false;
        const bool use_hs_lb = false;
        const int max_length = 4;

        if(use_cycle_lb || use_hs_lb){ // uncomment this to use cycle-based lower-bounding. This may be useful only for sparse graphs...
            VVI cycles;

            { // get all induced cycles up to some length
                VVI H = V;
                VVI revH = GraphUtils::reverseGraph(H);

                VVI nonpiH = Utils::getNonPIGraph(H);
                VVI revnonpiH = GraphUtils::reverseGraph(nonpiH);
                VB in_V(V.size(),false);
                VI marker(V.size(),0);
                VB is_end(V.size(),false);

                VI A(V.size()); iota(ALL(A),0);
                cycles = Utils::getAllSimpleCycles2(H, revH, nonpiH, revnonpiH, A, in_V,
                                                   marker, is_end, max_length);
            }

            if(use_cycle_lb) {
                VVI cycV = Utils::getCycleAdjacencyGraph(cycles);
                VI lb_cyc = Utils::getLowerBoundCycles(cycV, 50); // use FastVC/NuMVC 50 millis to find mis
                local_lb = max(local_lb, (int) lb_cyc.size());
            }

            if(use_hs_lb){
                VI hs = Utils::findMinHittingSet(cycles);
                if( hs.size() > local_lb ){
                    local_lb = hs.size();
                }

                vc_dfvs = hs;
            }
        }

        TimeMeasurer::stop("Lower bounding");
    }

    if(write_logs){
        indent(rec_depth);
        DEBUG(partial_size + local_lb);
    }

    if( partial_size + local_lb >= ub_total ){
        if(write_logs){
            indent(rec_depth);
            clog << "******** upper bound exceeded: " << partial_size << " + " << local_lb << " = "
                << partial_size + local_lb << " >= " << ub_total << endl;
        }
        return false; // no need to branch any further
    }

    return true; // we need to branch further
}

VI DFVSSolverE:: solveForInputGraph(VVI V) {

    VI exact_all;

    { // using all reductions in initial graph reduction

        Reducer red(V,cnf);
        red.cnf.enableAllReductions();
//        red.cnf.disableAllNonbasicReductions();
        red.cnf.disableAllConditionalReductions();
        red.cnf.disableAllRecursiveReductions();
        red.cnf.reducer_use_domination_6inserter = false;

        auto reductions = red.reduce();
        VI red_dfvs = Reducer::convertKernelizedReductions(reductions);

        exact_all = red_dfvs;
        V = red.V;
    }

    StronglyConnectedComponents scc(V);
    scc.createStronglyConnectedComponents();
    VVI comps = scc.getComponents();

    int total_millis_per_graph = cnf.vc_improver_milliseconds;
    int E = GraphUtils::countEdges(V,true);

    for( VI & cmp : comps ) {
        if( cmp.size() == 1 ) continue;
        if( cmp.size() == 2){
            exact_all.push_back(cmp[0]);
            continue;
        }

        InducedGraph g = GraphInducer::induce( V, cmp );

        VI heur_dfvs;

        heur_dfvs.resize(g.V.size());
        iota(ALL(heur_dfvs),0);
        //  HEURISTIC UPPER BOUND
        if( !Utils::isPIGraph(g.V) ){
            DFVSSolverH solverH(cnf);
            solverH.cnf.write_logs = write_logs;
            solverH.cnf.disableAllRecursiveReductions();
            solverH.cnf.disableAllConditionalReductions();

            solverH.cnf.vc_improver_milliseconds = 5 +
                    total_millis_per_graph * GraphUtils::countEdges(g.V,true) / E;

            heur_dfvs = solverH.solveForGraph(g.V);
            if(write_logs) DEBUG(heur_dfvs.size());

            assert(Utils::isFVS(g.V, heur_dfvs));
        }else if(V.size() > 100){
            heur_dfvs = Utils::getUpperBoundByVCOnSuperPIGraph(g.V, 1'000);
            if(write_logs) DEBUG(heur_dfvs.size());
        }

        int lb = 0;
        int ub = heur_dfvs.size();


        TimeMeasurer::start("DFVSSolverE");
        DFVSSolverE solverE(&g.V, cnf);
        solverE.cnf.reducer_use_twins_merge = false;

        solverE.write_logs = write_logs;

        VVI revGV = GraphUtils::reverseGraph(g.V);

        VB helper(g.V.size(),false);
        int pi_arcs = Utils::getAllPIEdges(g.V, revGV, helper).size();
        int all_arcs = GraphUtils::countEdges(g.V, true);

        double threshold = 0.3;

        VI exact_partial;
        if( filesystem::exists("findminhs")
            &&    (1.0 * pi_arcs / all_arcs < threshold && g.V.size() > 40) ){ // for small graphs we use branching
            assert( Utils::isFVS(g.V, heur_dfvs) );
            DFVSImprover impr(g.V, cnf);
            if( 1.0 * pi_arcs / all_arcs < 0.1 ) {
                heur_dfvs = impr.improveBySALSForSparseGraphs(heur_dfvs, 5, 0.1, 0.33,
                                                              0.99, // #TEST - uncomment after tests
                                                              150 * g.V.size(), 30);
            }
            int Egv = GraphUtils::countEdges( g.V, true );
            if(E == Egv && find_only_size_of_optimal_solution){
                // here we set known upper bound to g.V
                solverE.find_only_size_of_optimal_solution = find_only_size_of_optimal_solution;
                solverE.known_upper_bound = min(known_upper_bound - exact_all.size(), heur_dfvs.size() );
                exact_partial = solverE.solveForStronglyConnectedIterativeHittingSet2(g.V, heur_dfvs, true);
            }else{
                DEBUG(E); DEBUG(Egv); DEBUG(find_only_size_of_optimal_solution);
                exact_partial = solverE.solveForStronglyConnectedIterativeHittingSet2(g.V, heur_dfvs, false);
            }

            if( E == Egv && solverE.find_only_size_of_optimal_solution ){
                // here exact_all is just red_dfvs
                optimal_solution_size = solverE.getOptimalSolutionSize() + exact_all.size();
                return {};
            }
            assert( Utils::isFVS(g.V, exact_partial) );
        }
        else{
            exact_partial = solverE.solve(g.V, revGV, {}, 0, 0, lb, ub);
        }
        
        if (exact_partial.size() >= heur_dfvs.size()) {
            if(write_logs) clog << "Heuristic solver found optimum for given SCC" << endl;
            exact_partial = heur_dfvs;
        }

        for( int & d : exact_partial ) d = g.nodes[d];

        exact_all += exact_partial;
        TimeMeasurer::stop("DFVSSolverE");
    }

    assert(Utils::isFVS2(V, exact_all));
    assert(Utils::isFVS(V, exact_all));

    optimal_solution_size = exact_all.size();
    return exact_all;
}




int DFVSSolverE::getOptimalSolutionSize() {
    return optimal_solution_size;
}


VI DFVSSolverE::solveForStronglyConnectedIterativeHittingSet2(VVI V, VI heur_dfvs, bool use_exact_algorithm,
                                                              bool return_heuristic_solution) {

    bool switch_to_exact_after_hs_is_dfvs = false;
    constexpr bool find_cycles_heuristically = false; // original
    const bool take_only_smallest_cycles = false; // original

    int upper_bound = heur_dfvs.size();
    bool used_ihs = false;

    auto improveHeurDFVSByIHS = [&](){
        VVI all_cycles;
        DFVSSolverH solver(cnf);
        solver.cnf.solverh_ihs_init_max_cycles = 10;
        solver.cnf.ihs_hsls_perm_deviation_frequency = 400;
        VI t1 = solver.solveForStronglyConnectedIterativeHittingSet2(V, heur_dfvs, all_cycles, 0.1); // for larger graphs
        if(t1.size() < heur_dfvs.size()) heur_dfvs = t1;

        solver.cnf.solverh_ihs_init_max_cycles = 50;
        solver.cnf.ihs_hsls_perm_deviation_frequency = 800;
        t1 = solver.solveForStronglyConnectedIterativeHittingSet2(V, heur_dfvs, all_cycles, 0.2); // for larger graphs
        if(t1.size() < heur_dfvs.size()) heur_dfvs = t1;

        upper_bound = min(upper_bound, (int)heur_dfvs.size());
        used_ihs = true;
    };

    if(false) improveHeurDFVSByIHS();

    int N = V.size();
    VB helper(N,false);

    VI hs;
    int iter = 0;
    int max_cycle_length = 3;
    VVI all_cycles;


    if(find_cycles_heuristically){
        VVI all_cycles1, all_cycles2;
        DFVSSolverH solver(cnf);

        VI t1,t2;
        {
            solver.cnf.solverh_ihs_init_max_cycles = 50;
            solver.cnf.ihs_hsls_perm_deviation_frequency = 800;
            t1 = solver.solveForStronglyConnectedIterativeHittingSet2(V, heur_dfvs, all_cycles1, 0.1);

            solver.cnf.solverh_ihs_init_max_cycles = 10;
            solver.cnf.ihs_hsls_perm_deviation_frequency = 400;
            t2 = solver.solveForStronglyConnectedIterativeHittingSet2(V, t1, all_cycles2, 0.2);
        }

        if(cnf.write_logs){ DEBUG(t1.size());DEBUG(t2.size()); }

        if( t1.size() < heur_dfvs.size() ){
            heur_dfvs = t1;
            hs = t1;
        }
        if( t2.size() < heur_dfvs.size() ){
            heur_dfvs = t2;
            hs = t2;
        }
        upper_bound = min(upper_bound, (int)heur_dfvs.size());
        if(cnf.write_logs) DEBUG(upper_bound);


        if(take_only_smallest_cycles) {
            if(t1.size() != t2.size()){
                if(t1.size() < t2.size()){
                    all_cycles = all_cycles1;
                    hs = t1;
                }
                else{
                    all_cycles = all_cycles2;
                    hs = t2;
                }
            }else{
                if( all_cycles1.size() < all_cycles2.size() ){
                    all_cycles = all_cycles1;
                    hs = t1;
                }
                else{
                    all_cycles = all_cycles2;
                    hs = t2;
                }
            }
        }else{ // take union of all found cycles
            set<LL> present_cycles_hashes;
            VVI all_cycles12 = all_cycles1 + all_cycles2;
            for (VI &C : all_cycles12) {
                LL hash = 0;
                for (int d : C) hash ^= hashes[d];
                if (present_cycles_hashes.count(hash) == 0) {
                    present_cycles_hashes.insert(hash);
                    all_cycles.push_back(C);
                }
            }

            if( t1.size() < t2.size() ) hs = t1;
            else hs = t2;
        }
    }

    bool previous_hs_proved_optimal = false;
    int previous_optimal_hs_size = -1;


    int ITR1 = 1e6; // original
    int ITR2 = 5e4; // original
    int ITR3 = 3e6; // original
    int ITR4 = 3e4; // original
    int ITR5 = 1e5; // original
    int MAX_CYCLES = 10; // original

    auto rescaleITRS = [&](double F) {
        ITR1 *= F;
        ITR2 *= F;
        ITR3 *= F;
        ITR4 *= F;
        ITR5 *= F;
    };

    rescaleITRS(0.2); // original 0.2

    if(cnf.write_logs){
        DEBUG(MAX_CYCLES);
        DEBUG(all_cycles.size());
    }

    int previous_cycle_search_millis = 0;
    const int max_cycle_search_time_millis = 1'000;

    label_loop_beginning:
    do{

        iter++;

        if(cnf.write_logs) {
            clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;
            DEBUG(iter); DEBUG(hs.size());
            DEBUG(all_cycles.size());
            DEBUG(max_cycle_length);
        }

        VB in_V(N,false);
        VB is_end(N,false);
        VI marker(N,0);
        VI A = CombinatoricUtils::getRandomPermutation(N);

        while(true){ // add all unhit induced cycles of length <= max_cycle_length
            // H is a graph nonpiV / vc
            VVI H = V;
            VVI revH = GraphUtils::reverseGraph(H);

            if (iter > 1) Utils::removeNodes(H, revH, hs, helper);

            VVI nonpiH = Utils::getNonPIGraph(H);
            VVI revnonpiH = GraphUtils::reverseGraph(nonpiH);


            Stopwatch sw;
            sw.start("getAllSimpleCycles");
            VVI cycles;
            if(previous_cycle_search_millis <= max_cycle_search_time_millis) {
                cycles = Utils::getAllSimpleCycles2(H, revH, nonpiH, revnonpiH, A, in_V, marker,
                                                    is_end, max_cycle_length);
            }
            else cycles = Utils::getAllSimpleCycles3(H, max_cycle_length); // much faster!

            sw.stop("getAllSimpleCycles");
            sw.write("getAllSimpleCycles");
            if(previous_cycle_search_millis <= max_cycle_search_time_millis){
                previous_cycle_search_millis = sw.getTime("getAllSimpleCycles");
            }

            const int max_cycles_size = (iter == 1 ? (int) cycles.size() : MAX_CYCLES); // original

            if (cycles.size() > max_cycles_size) {
                StandardUtils::shuffle(cycles);
                sort(ALL(cycles), [](auto &v, auto &w) { return v.size() < w.size(); });
                cycles.resize(max_cycles_size);
            }
            if(cnf.write_logs)  DEBUG(cycles.size());
            if(cycles.empty()){
                max_cycle_length++;
                if(cnf.write_logs) clog << "Increasing cycle length" << endl;
                continue;
            }
            else if( cycles.size() <= max_cycles_size / 2 ) max_cycle_length++; // increasing the value quicker

            all_cycles += cycles;
            break;
        }

        bool find_exact_hs = use_exact_algorithm;
        if(!find_exact_hs){
            if(cnf.write_logs) clog << "*********** !!!!!! ********** "
                    "CAUTION! Using a heuristic HS Solver in DFVSSolverE::solveForStronglyConnectedIterativeHittingSet2" << endl;
        }

        bool found_hs_by_ls = false;

        const bool use_hs_ls_search = ( all_cycles.size() > 100 );
        if(use_hs_ls_search){ // original
            if(cnf.write_logs) clog << "Trying to find 'equivalent' hs using local search" << endl;

            int iters = ( hs.size() >= upper_bound - 4 ? ITR1 : ITR2 );
            if( hs.size() == upper_bound-1 ) iters = ITR3;

            VI hs_ls;
            if(find_exact_hs) hs_ls = Utils::hsImprovementLS2( all_cycles, hs, cnf, iters, hs.size() );
            else hs_ls = Utils::hsImprovementLS2( all_cycles, hs, cnf, iters, (int)hs.size()-3 ); // enabling decreasing size by 3

            if( !hs_ls.empty() ){
                if(find_exact_hs && previous_hs_proved_optimal ){
                    if( hs_ls.size() <= hs.size() ) assert(hs_ls.size() == hs.size()); // this will hold only if we find exact HS
                }
                found_hs_by_ls = true;
                hs = hs_ls;
                if(cnf.write_logs) clog << "\tFOUND HS BY LOCAL SEARCH" << endl;
            }
        }

        auto makeHsAValidHittingSet = [&](){
            if(cnf.write_logs) clog << "Making hs a valid HS" << endl;
            VVI unhit_cycles;
            VB in_hs = StandardUtils::toVB(N, hs);

            for( auto & C : all_cycles ){
                bool hit = false;
                for(int d : C) if( in_hs[d] ){ hit = true; break; }
                if(!hit) unhit_cycles.push_back(C);
            }

            while( !unhit_cycles.empty() ) {
                map<int, int> freqs;
                for (auto &C : unhit_cycles) {
                    for (int d : C) freqs[d]++;
                }

                VPII freqs_sorted(ALL(freqs));
                sort(ALL(freqs_sorted), [](auto a, auto b) { return a.second > b.second; });
                int best_el = freqs_sorted[0].first;
                hs.push_back(best_el);
                in_hs[best_el] = true;

                VVI temp = unhit_cycles;
                unhit_cycles.clear();

                for (auto &C : temp) {
                    bool hit = false;
                    for (int d : C)
                        if (in_hs[d]) {
                            hit = true;
                            break;
                        }
                    if (!hit) unhit_cycles.push_back(C);
                }

                if(!use_exact_algorithm && iter > 1)
                {
                    VI xyz = Utils::hsImprovementLS2(all_cycles, hs, cnf, ITR4);
                    if (!xyz.empty()) {
                        hs = xyz;
                        break;
                    }
                }
            }

            {
                VI xyz = Utils::hsImprovementLS2( all_cycles, hs, cnf, ITR5 );
                if( !xyz.empty() ) hs = xyz;
            }
        };

        // take greedily nodes that hits most unhit cycles to hs - not optimal solution! #CAUTION
        if(!found_hs_by_ls && !find_exact_hs ){
            makeHsAValidHittingSet();
        }

        if(!find_exact_hs){
            if(switch_to_exact_after_hs_is_dfvs && (Utils::isFVS(V, hs) || hs.size() >= upper_bound )){
                switch_to_exact_after_hs_is_dfvs = false;
                use_exact_algorithm = find_exact_hs = true;
                previous_hs_proved_optimal = false;
                previous_optimal_hs_size = -1;
                find_exact_hs = true;
            }
        }

        if(find_exact_hs){
            if(!found_hs_by_ls || !previous_hs_proved_optimal) {

                {
                    makeHsAValidHittingSet();
                    if(previous_hs_proved_optimal && hs.size() == previous_optimal_hs_size){
                        // do nothing - hs is optimal
                    }
                    else{
                        VI ub_hs = hs;
                        if(heur_dfvs.size() < ub_hs.size()){
                            ub_hs = heur_dfvs;
                        }

                        if( cnf.write_logs ){
                            clog << "Searching for MIN HS with LB: " << previous_optimal_hs_size << " and UB: "
                                 << ub_hs.size() << endl;
                        }
                        hs = Utils::findMinHittingSet(all_cycles, ub_hs, previous_optimal_hs_size );
                    }
                }
                previous_hs_proved_optimal = true;
                previous_optimal_hs_size = hs.size();
            }
        }

        if(cnf.write_logs) {
            ENDL(1); DEBUG(find_only_size_of_optimal_solution);DEBUG(find_exact_hs);
            DEBUG(known_upper_bound); DEBUG(upper_bound);
        }

        if( use_exact_algorithm && find_only_size_of_optimal_solution &&
            ((hs.size() == known_upper_bound) || ( hs.size() >= upper_bound )) ){
            optimal_solution_size = hs.size();

            if(Utils::isFVS(V, hs)) return hs;
            else return {};
        }

        if( hs.size() >= upper_bound ){
            if(cnf.write_logs && find_exact_hs) clog << "Heuristic solver found optimum, hs.size() >= upper_bound" << endl;
            hs = heur_dfvs;
            break;
        }

        if(cnf.write_logs) {  DEBUG(hs.size()); ENDL(5); ENDLS(50, "*"); }

    }while( !Utils::isFVS(V,hs) );

    if(cnf.write_logs) clog << "Found " << ( use_exact_algorithm ? "optimal" : "heuristically a " ) << " DFVS of size: " << hs.size() << endl;

    if(!use_exact_algorithm && return_heuristic_solution) return hs;

    if(hs.size() < upper_bound){
        upper_bound = hs.size();
        heur_dfvs = hs;
    }
    if(use_exact_algorithm) optimal_solution_size = hs.size();

    if(!use_exact_algorithm)
    {
        if(!used_ihs) improveHeurDFVSByIHS();
        if( hs.size() > heur_dfvs.size() ){ hs = heur_dfvs; assert(Utils::isFVS(V,hs)); }

        set<PII> cycle_arcs;
        for( VI c : all_cycles ){
            for( int i=0; i<c.size(); i++ ) cycle_arcs.insert({c[i], c[(i + 1) % c.size() ] } );
        }

        {
            VB in_V(N, false);
            VB is_end(N, false);
            VI marker(N, 0);
            VI A = CombinatoricUtils::getRandomPermutation(N);
            VVI H = V;
            VVI revH = GraphUtils::reverseGraph(H);

            VVI nonpiH = Utils::getNonPIGraph(H);
            VVI revnonpiH = GraphUtils::reverseGraph(nonpiH);

            int max_cycle_length = 5;
            VVI cycles = Utils::getAllSimpleCycles2(H, revH, nonpiH, revnonpiH, A, in_V, marker,
                                                    is_end, max_cycle_length);

            if(cnf.write_logs) DEBUG(cycle_arcs.size());
            for( VI c : cycles ){
                for( int i=0; i<c.size(); i++ ) cycle_arcs.insert({c[i], c[(i + 1) % c.size() ] } );
            }
            if(cnf.write_logs) DEBUG(cycle_arcs.size());
        }

        if(cnf.write_logs) {
            ENDL(15);
            clog << "Validating if heuristic HS found optimal DFVS" << endl;
            clog << "Creating graph H induced by arcs in all_cycles" << endl;

            DEBUG(cycle_arcs.size());
            DEBUG(GraphUtils::countEdges(V, true));
        }

        {
            VPII arcs(ALL(cycle_arcs));
            VVI H = GraphUtils::getGraphForEdges(arcs, true);
            if(cnf.write_logs) DEBUG(GraphUtils::countEdges(H,true));

            Reducer red(H,cnf);
            red.cnf.enableAllReductions();
            red.cnf.reducer_use_domination_6inserter = false;
            red.cnf.disableAllRecursiveReductions();
            red.cnf.disableAllConditionalReductions();

            auto reductions = red.reduce();
            VI red_dfvs = Reducer::convertKernelizedReductions(reductions);

            if(cnf.write_logs) red.writeTotals();
            H = red.V;
            if(cnf.write_logs) DEBUG(red_dfvs.size());
            if(cnf.write_logs) DEBUG(GraphUtils::countEdges(H,true));

            InducedGraph hh = GraphInducer::induceByNonisolatedNodes(H);
            H = hh.V;
            if(cnf.write_logs) DEBUG(GraphUtils::countEdges(H,true));

            DFVSSolverE solver(&H, cnf);
            solver.write_logs = write_logs;
            solver.cnf.vc_improver_milliseconds = 500;
            solver.cnf.disableAllRecursiveReductions();
            solver.cnf.disableAllConditionalReductions();
            solver.setKnownUpperBound( min(hs.size(), heur_dfvs.size()) - red_dfvs.size() );
            solver.setFindOnlySizeOfOptimalSolution(true); // we will find only the size of the solution

            VI optH = solver.solveForInputGraph(H);

            for(int& d : optH) d = hh.nodes[d];

            int lb_size;
            if( optH.empty() ) lb_size = solver.getOptimalSolutionSize(); // we do not know the size of optimal solution
            else lb_size = optH.size();

            optH += red_dfvs;
            lb_size += red_dfvs.size();

            if( lb_size == upper_bound ){
                if(cnf.write_logs) clog << "OPTIMAL RESULT WAS FOUND BY HEURISTIC HS" << endl;
            }else{
                if(cnf.write_logs) clog << "NOT PROVED THAT OPTIMAL result was found by heuristic, lower bound: " << lb_size << endl;

                if( optH.size() != red_dfvs.size() && Utils::isFVS(V, optH) ){
                    if(cnf.write_logs) clog << "... but optH is a DFVS of V as well, so we found optimum anyway" << endl;
                    hs = optH;
                }else{

                    if (cnf.write_logs){
                        clog << "Going back to main HS loop, starting with optH as hs that may not be optimal for"
                                << " current cycles, but is a LB for V, because it is optimal DFVS of a subgraph "
                                   "H of V" << endl;
                    }
                    if( optH.size() != red_dfvs.size()){
                        hs = optH; // optH may be empty, hence this condition
                        assert(!Utils::isFVS(V,hs));
                    }else{
                        while( hs.size() > lb_size ) hs.pop_back();
                        if(Utils::isFVS(V,hs)) return hs;
                    }

                    use_exact_algorithm = true;
                    previous_hs_proved_optimal = true; // previous hs may not be optimal, but it is a LB for DFVS of V
                    previous_optimal_hs_size = lb_size;
                    goto label_loop_beginning;
                }
            }
        }
    }

    assert(Utils::isFVS(V,hs));

    return hs;
}
