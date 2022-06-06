//
// Created by sylwester on 12/20/21.
//

#include <CONTESTS/PACE22/heur/DFVSSolverH.h>
#include <CONTESTS/PACE22/heur/AgentFlow.h>
#include <CONTESTS/PACE22/Reducer.h>
#include <graphs/GraphUtils.h>
#include <graphs/VertexCover/VCUtils.h>
#include <utils/TimeMeasurer.h>
#include <combinatorics/CombinatoricUtils.h>
#include <graphs/GraphInducer.h>
#include <CONTESTS/PACE22/heur/VCImprover.h>
#include <CONTESTS/PACE22/heur/DFVSImprover.h>
#include <graphs/scc/StronglyConnectedComponents.h>
#include <graphs/VertexCover/approximation/LibMVC/fastvc.h>
#include <graphs/VertexCover/approximation/NuMVC/NuMVC.h>
#include <CONTESTS/PACE22/heur/SALS2.h>
#include <graphs/VertexCover/kernelization/KernelizerVC.h>
#include <CONTESTS/PACE22/heur/SALS3.h>
#include <utils/Stopwatch.h>
#include <CONTESTS/PACE22/heur/SALS4.h>
#include <CONTESTS/PACE22/heur/SALS5.h>
#include <graphs/cycles/CycleCounter.h>
#include "StandardUtils.h"


VI DFVSSolverH::solveForGraph(VVI V) {
    return solveForGraph2(V);

    VI dfvs_final;

    {
        if(cnf.write_logs) clog << "Setting in DFVSSolverH::solveForGraph all initial reductions to true" << endl;
        Reducer red(V, cnf);

        if(!cnf.solverh_use_reductions_initial){
            red.cnf.disableAllNonbasicReductions();
        }

        red.cnf.write_logs = false;
        red.cnf.disableAllRecursiveReductions();
        red.cnf.disableAllConditionalReductions();

        // #CAUTION! No conditional reductions can be used here
        red.cnf.disableAllConditionalReductions();
        auto reductions = red.reduce();
        VI red_dfvs = Reducer::convertKernelizedReductions(reductions);

        V = red.V;
        if(cnf.write_logs) DEBUG(red_dfvs.size());

        dfvs_final = red_dfvs;
    }

    cnf.reducer_use_twins_merge = false;
    cnf.reducer_use_nonsimple_cycle_arcs = false;

    if(cnf.write_logs){  DEBUG(V.size()); DEBUG(GraphUtils::countEdges(V, true)); }

    StronglyConnectedComponents scc(V);
    scc.createStronglyConnectedComponents();
    VVI comps = scc.getComponents();
    if(cnf.write_logs) DEBUG(comps.size());
    sort(ALL(comps), [](auto & v, auto & w){
        return v.size() < w.size();
    });

    origV_edges = GraphUtils::countEdges(V,true);
    origV_N = V.size();

    VI dfvs_comps;

    for (auto & cmp : comps) {
        if (cmp.size() == 1) continue;
        if (cmp.size() == 2) { // strongly connected component of size two
            dfvs_comps.push_back(cmp[0]);
            continue;
        }

        if(cnf.write_logs) {
            ENDL(3);
            clog << "********************************************************************************" << endl;
            ENDL(3);
            DEBUG(dfvs_comps.size());
        }

        InducedGraph g = GraphInducer::induce(V, cmp);
        if(cnf.write_logs){
            DEBUG(g.V.size()); DEBUG(GraphUtils::countEdges(g.V, true));
            VVI grevV = GraphUtils::reverseGraph(g.V);
            VB helper(grevV.size(),false);
            DEBUG(Utils::getAllPIEdges(g.V, grevV, helper).size());
        }

        // no need to use reductions if graph was reduced at the beginning of solveForGraph()
        VI dfvs = solveForBiconnectedGraph(g.V);
        for( int & d : dfvs ) d = g.nodes[d];

        dfvs_comps += dfvs;
    }

    assert(Utils::isFVS(V,dfvs_comps)); // this should hold, because V was taken from Reducer
    dfvs_final += dfvs_comps;
    assert(Utils::isFVS(V,dfvs_final));

    if(cnf.write_logs){
        DEBUG(dfvs_final.size());
        LL solution_hash = Utils::getSetHash( V.size(), dfvs_final );
        DEBUG(solution_hash);
    }

    return dfvs_final;
}

VI DFVSSolverH::solveForGraph2(VVI V) {
    VI dfvs_final;

    VVI origV = V;

    vector<DFVSReduction*> reductions;

    {
        if(cnf.write_logs) clog << "Setting in DFVSSolverH::solveForGraph all initial reductions to true" << endl;
        Reducer red(V, cnf);
        if(!cnf.solverh_use_reductions_initial) red.cnf.disableAllNonbasicReductions();
        red.cnf.write_logs = false;
        reductions = red.reduce();
        V = red.V;
        if(cnf.write_logs) DEBUG( Reducer::getReductionsSizeDiff(reductions));
    }

    if(cnf.write_logs){  DEBUG(V.size()); DEBUG(GraphUtils::countEdges(V, true)); }

    StronglyConnectedComponents scc(V);
    scc.createStronglyConnectedComponents();
    VVI comps = scc.getComponents();
    if(cnf.write_logs) DEBUG(comps.size());
    sort(ALL(comps), [](auto & v, auto & w){
        return v.size() < w.size();
    });

    origV_edges = GraphUtils::countEdges(V,true);
    origV_N = V.size();


    vector<InducedGraph> comp_graphs;
    vector<VVI> comp_orig_graphs;
    vector<vector<DFVSReduction*>> comp_reductions;
    VVI comp_dfvs;

    const bool use_vc_improvement_straightaway = (origV_edges > 1e5);

    for (auto & cmp : comps) {
        if (cmp.size() == 1) continue;
        if (cmp.size() == 2) { // strongly connected component of size two
            dfvs_final.push_back(cmp[0]);
            continue;
        }

        if(cnf.write_logs) {
            ENDL(3);
            clog << "********************************************************************************" << endl;
            ENDL(3);
        }

        InducedGraph g = GraphInducer::induce(V, cmp);
        VVI orig_gV = g.V;
        if(cnf.write_logs){
            DEBUG(g.V.size()); DEBUG(GraphUtils::countEdges(g.V, true));
            VVI grevV = GraphUtils::reverseGraph(g.V);
            VB helper(grevV.size(),false);
            DEBUG(Utils::getAllPIEdges(g.V, grevV, helper).size());
        }

        // no need to use reductions if graph was reduced at the beginning of solveForGraph()
        bool old_use_reducer_scc = cnf.solverh_use_reductions_for_each_scc;

        vector<DFVSReduction*> cmp_red;
        if(cnf.solverh_use_reductions_for_each_scc){
            Reducer red(g.V, cnf);
            auto E = GraphUtils::countEdges(V, true);
            red.cnf.reducer_max_time_millis = 5 + 1ll * cnf.reducer_max_time_millis * E / origV_edges;
            red.cnf.write_logs = false;
            cmp_red = red.reduce();
            g.V = red.V; // changing structure of g.V to that after reductions

            if(cnf.write_logs) DEBUG(Reducer::getReductionsSizeDiff(reductions));
        }

        cnf.solverh_use_reductions_for_each_scc = false;
        VI dfvs = solveForBiconnectedGraph(g.V, false, use_vc_improvement_straightaway); // use vc improvements for large graphs
//        if( Utils::getPieEdgesPercentage(g.V) < 0.05 && dfvs.size() > 300 && dfvs.size() < 2'000 ){ // original
//            int orig_dst = cnf.agent_flow_max_distance_from_best;
//            int orig_min_dist = cnf.agent_flow_min_distance;
//
//            if(dfvs.size() > 1'000) {
//                VPII configs;
//                if (dfvs.size() < 1'000)
//                    configs = {{4, 2},
//                               {4, 3}};
//                else configs = {{4, 2}};
//
//                for (auto[div, min_dst] : configs) {
//                    cnf.agent_flow_max_distance_from_best = orig_dst / div;
//                    cnf.agent_flow_min_distance = min_dst;
//                    VI dfvs2 = solveForBiconnectedGraph(g.V, false, false); // do not use improvements here
//                    if (dfvs2.size() < dfvs.size()) dfvs = dfvs2;
//                }
//            }
//            else if( dfvs.size() <= 1'000 ) { // #TEST - check also these two alternatives
//                {
//                    cnf.agent_flow_node_selection_type = cnf.agent_flow_remove_largest_flow_node;
//                    cnf.agent_flow_max_distance_from_best = 1;
//                    cnf.agent_flow_min_distance = 2;
//                    VI dfvs2 = solveForBiconnectedGraph(g.V, false, false); // do not use improvements here
//                    if (dfvs2.size() < dfvs.size()) dfvs = dfvs2;
//                }
//
//                {
//                    cnf.agent_flow_node_selection_type = cnf.agent_flow_merge_smallest_flow_node;
//                    cnf.agent_flow_max_distance_from_best = 1;
//                    cnf.agent_flow_min_distance = 2;
//                    VI dfvs2 = solveForBiconnectedGraph(g.V, false, false); // do not use improvements here
//                    if (dfvs2.size() < dfvs.size()) dfvs = dfvs2;
//                }
//            }
//
//            cnf.agent_flow_max_distance_from_best = orig_dst;
//            cnf.agent_flow_min_distance = orig_min_dist;
//        }
        cnf.solverh_use_reductions_for_each_scc = old_use_reducer_scc;

        assert(Utils::isFVS(g.V, dfvs));

        if(!dfvs.empty()) {
            comp_orig_graphs.push_back(orig_gV);
            comp_graphs.push_back(g);
            comp_reductions.push_back(cmp_red);
            comp_dfvs.push_back(dfvs);
        }else{
            // this is the case that g.V was fully reduced and dfvs is empty
            Reducer::liftSolution(g.V.size(), dfvs, cmp_red);
            for(int & d : dfvs) d = g.nodes[d];
            dfvs_final += dfvs;
        }
    }

    int current_dfvs_size = accumulate( ALL(comp_dfvs), 0, [](int s, VI & v){ return v.size()+s; } );

    if(cnf.write_logs){
        ENDL(10);
        ENDLS(50,"*");
        ENDL(5);
        clog << "Initial DFVS created for SCCs, current_dfvs_size: " << current_dfvs_size << endl;
        ENDL(10);
        ENDLS(50,"*");
    }


    { // now use improvements
        for( int i=0; i<comp_dfvs.size(); i++ ){
            bool use_improvements = (  comp_graphs[i].V.size() >= cnf.solverh_min_graph_size_for_improvements );

            if(use_improvements && !cnf.tle()) { // that is the safety check - we may not manage in time to improve all graphs
                if(cnf.write_logs){
                    int all_arcs = GraphUtils::countEdges(comp_graphs[i].V,true);
                    int pi_arcs = Utils::countPiEdges(comp_graphs[i].V);
                    int nonpi_arcs = all_arcs - pi_arcs;
                    ENDL(1);
                    ENDLS(20,"*");
                    ENDL(1);
                    clog << "Improving graph with " << comp_graphs[i].V.size() << " nodes, "
                         << all_arcs << " arcs, from which pi_arcs: " << pi_arcs << ", nonpi_arcs: " << nonpi_arcs
                         << " and reductions_size_diff: " << Reducer::getReductionsSizeDiff(comp_reductions[i]) << endl;
                }

                int improvement_diff = 0;
                if (!comp_dfvs[i].empty()) {
                    DFVSImprover impr(comp_graphs[i].V, cnf);
//                    VI new_dfvs = impr.fullImprovement(comp_dfvs[i], origV_N, origV_edges); // original
                    VI new_dfvs = impr.fullImprovement(comp_dfvs[i], origV_N, origV_edges, !use_vc_improvement_straightaway);
                    improvement_diff = comp_dfvs[i].size() - new_dfvs.size();
                    comp_dfvs[i] = new_dfvs;
                }

                assert(Utils::isFVS(comp_graphs[i].V, comp_dfvs[i]));
                int old_dfvs_size = current_dfvs_size;
                current_dfvs_size -= improvement_diff;
                if (cnf.write_logs) {
                    DEBUG(Reducer::getReductionsSizeDiff(comp_reductions[i]));
                    DEBUG(old_dfvs_size);
                    DEBUG(improvement_diff);
                    DEBUG(current_dfvs_size);
                }
            }

            Reducer::liftSolution( comp_graphs[i].V.size(), comp_dfvs[i], comp_reductions[i] );
            assert(Utils::isFVS(comp_orig_graphs[i], comp_dfvs[i]));
            for( int & d : comp_dfvs[i] ) d = comp_graphs[i].nodes[d]; // remapping nodes to original ids
        }
    }

    for( VI & cmp_dfvs : comp_dfvs ) dfvs_final += cmp_dfvs;
    Reducer::liftSolution(V.size(), dfvs_final, reductions);
    assert(Utils::isFVS(origV,dfvs_final));

    if(cnf.write_logs){
        DEBUG(dfvs_final.size());
        LL solution_hash = Utils::getSetHash( V.size(), dfvs_final );
        DEBUG(solution_hash);
    }

    return dfvs_final;
}

VI DFVSSolverH::solveForBiconnectedGraph(VVI V, bool allow_improvements_here, const bool use_vc_improvement) {
    if( !Utils::hasCycle(V) ) return {};

    VI dfvs;
    VVI sccV = V;

    const bool use_reductions = true;

    vector<DFVSReduction*> reductions;
    if(use_reductions){
        Reducer red(V, cnf);
        auto E = GraphUtils::countEdges(V, true);
        red.cnf.reducer_max_time_millis = 5 + 1ll * cnf.reducer_max_time_millis * E / origV_edges;

        red.cnf.write_logs = false;
        if(!cnf.solverh_use_reductions_for_each_scc) red.cnf.disableAllNonbasicReductions();

        // conditional reductions possible!
        reductions = red.reduce();

        V = red.V;

        if(cnf.write_logs) DEBUG(Reducer::getReductionsSizeDiff(reductions));

        if( !Utils::hasCycle(V) ){
            VI dfvs = {};
            Reducer::liftSolution(V.size(), dfvs, reductions);
            assert( Utils::isFVS(sccV, dfvs) );
            return dfvs;
        }
    }

    const bool is_pi_graph = Utils::isPIGraph(V);
    bool is_almost_pi_graph = false;
    {
        VVI revV = GraphUtils::reverseGraph(V);
        VB helper(V.size(),false);
        is_almost_pi_graph = ( 1.0 * Utils::getAllPIEdges(V,revV,helper).size()
                / GraphUtils::countEdges(V,true) ) > 0.9;
    }
    if(cnf.solverh_use_superpi_vc_ub && is_almost_pi_graph){
        LL FACT = 3;
        if(is_pi_graph) FACT = 10;

        int millis = 5ll + FACT*ceil(1.0 * cnf.vc_improver_milliseconds
                                * GraphUtils::countEdges(V, true) / origV_edges);
        if(Utils::isPIGraph(V)) millis = min(millis, 120'000);
        else millis = min(millis, 30'000);

        VI vc = Utils::getUpperBoundByVCOnSuperPIGraph(V, millis);

        assert(Utils::isFVS(V, vc));

        dfvs = vc;

        if(is_pi_graph){
            Reducer::liftSolution(V.size(), dfvs, reductions);
            assert(Utils::isFVS(sccV, dfvs));
            if( cnf.write_logs ) clog << "Graph is a PI graph, returning found vc of size: " << vc.size() << endl;
//            return red_dfvs + dfvs;
            return dfvs; // with solution lifting
        }
    }


    VI dfvsH;

    LL N = V.size();
    int E  = GraphUtils::countEdges(V,true);
    double density = 2.0*E / ( N * (N-1));

    DFVSSolverH solver(cnf);

    if( V.size() < 10'000 && E > 1e6 && density > 0.05 ){ // #TEST - move it outside DFVSSolverH
        solver.cnf.agent_flow_node_selection_type = Config::agent_flow_merge_smallest_flow_node;
    }

//    dfvsH = solver.solveByAgentFlowAllWithinDistance(V, 2); // original

    { // #TEST #CAUTION- do not use AF at all for large graphs with high percentage of pi-arcs
        if( E > 1'000'000 && Utils::getPieEdgesPercentage(V) > 0.5 ){
            dfvsH = Utils::getUpperBoundByVCOnSuperPIGraph(V, 3'000);
        }else{
            dfvsH = solver.solveByAgentFlowAllWithinDistance(V, 2); // original
        }
    }



    if(cnf.write_logs) clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;

    assert(Utils::isFVS(V, dfvsH));
    if(cnf.write_logs) DEBUG(dfvsH.size());

    if( dfvs.empty() || dfvsH.size() < dfvs.size()) dfvs = dfvsH;

    bool use_improvements = ( V.size() >= cnf.solverh_min_graph_size_for_improvements );
    if(!allow_improvements_here) use_improvements = false;

    if(use_vc_improvement && !use_improvements){ // using VC improver straightaway
        DFVSImprover impr(V, cnf);

        impr.cnf.vc_improver_milliseconds = 5 + ceil(cnf.vc_improver_milliseconds
                                                     * GraphUtils::countEdges(V, true) / origV_edges);

        const int ITERS = cnf.solverh_improvement_iterations;
        if(cnf.write_logs){  DEBUG(ITERS); DEBUG(1ll * ITERS * impr.cnf.vc_improver_milliseconds); }

        dfvs = impr.improve(dfvs, ITERS);
        if(cnf.write_logs) DEBUG(dfvs.size());
    }

    if(use_improvements){
        DFVSImprover impr(V,cnf);
        dfvs = impr.fullImprovement(dfvs, origV_N, origV_edges);
    }


    assert(Utils::isFVS(V,dfvs));

    Reducer::liftSolution(V.size(), dfvs, reductions);
    assert(Utils::isFVS(sccV,dfvs));
    return dfvs;

}

VI DFVSSolverH::solveByAgentFlow(VVI V) {
    VVI _V0 = V;

    int N = V.size();
    VVI revV = GraphUtils::reverseGraph(V);
    VB helper(N,false);

    VI res;
    int E = GraphUtils::countEdges(V);

    const bool use_reductions = true;
    VD tokens;
    int iteration_cnt = 0;


    while( E > 0 ){
        if(use_reductions){ // using reductions
            Reducer red(V, cnf);
            if(!cnf.solverh_use_reductions_AF) red.cnf.disableAllNonbasicReductions();

            // #CAUTION! No conditional reductions can be used here
            red.cnf.disableAllConditionalReductions();
            auto reductions = red.reduce();
            VI dfvs_red = Reducer::convertKernelizedReductions(reductions);

            V = red.V;
            revV = GraphUtils::reverseGraph(V);
            res += dfvs_red;

            E = GraphUtils::countEdges(V,true);
            DEBUG(res.size());
        }

        if(E == 0) break;


        if(iteration_cnt % cnf.agent_flow_node_update_frequency == 0) {
            AgentFlow af(V);
            int iterations = 2*sqrt(N);
            int void_steps = sqrt(N);

            if (cnf.agent_flow_method == Config::agent_flow_continuous)tokens = af.flowContinuous(iterations, void_steps);
            else if (cnf.agent_flow_method == Config::agent_flow_tokens)tokens = af.flowTokens(iterations, void_steps, 3);
            else if (cnf.agent_flow_method == Config::agent_flow_sinkhorn) {
                tokens = sinkhorn(V);
                for (auto &d : tokens) d = 1.0 - d;
            }
        }

        if(cnf.agent_flow_node_selection_type == Config::agent_flow_remove_largest_flow_node) {
            int a = max_element(ALL(tokens)) - tokens.begin();
            tokens[a] = 0;
            Utils::removeNode(V, revV, a, helper);
            res.push_back(a);
        }else if(cnf.agent_flow_node_selection_type == Config::agent_flow_merge_smallest_flow_node){
            int a = -1;
            double least_flow = 1e9;
            VI loop_nodes;
            for(int i=0; i<N; i++){
                if( V[i].empty() ) continue;
                if(!use_reductions && Utils::hasLoop(V,i)){ // if we do no remove loops, then we need to check that
                    loop_nodes.push_back(i);
                    continue;
                }

                if( tokens[i] < least_flow ){
                    a = i;
                    least_flow = tokens[i];
                }
            }

            if(a != -1){
                Utils::merge(V, revV, a, helper);
            }
            else if(!use_reductions){ // if we do not remove loops, the we can merge node only if it has no loop
                res += loop_nodes;
                break;
            }
        }

        if(cnf.agent_flow_alternate_selection_type) cnf.agent_flow_node_selection_type = 1 - cnf.agent_flow_node_selection_type;
        if(cnf.agent_flow_alternate_flow_method){
            cnf.agent_flow_method = (1+cnf.agent_flow_method) % Config::agent_flow_methods_cnt;
        }

        E = GraphUtils::countEdges(V);
        iteration_cnt++;
    }

    reverse(ALL(res));
    VI res2 = Utils::findAndRemoveRedundantNodes(_V0, res);
    res = res2;

    return res;
}


VI DFVSSolverH::solveByAgentFlowFast(VVI V) {
    VVI _V0 = V;

    int N = V.size();
    VVI revV = GraphUtils::reverseGraph(V);
    VB helper(N,false);

    VI res;


    VD tokens;
    VB was(N,false);
    int next_node_smallest = 0;
    int next_node_largest = N-1;

    VI sorted_nodes(N);

    auto initialize_tokens = [&](){
        AgentFlow af(V);
        int iterations = 2*sqrt(N);
        int void_steps = sqrt(N);

        if(cnf.agent_flow_method == Config::agent_flow_continuous) tokens = af.flowContinuous(iterations, void_steps);
        else if(cnf.agent_flow_method == Config::agent_flow_tokens) tokens = af.flowTokens(iterations, void_steps, 3);
        else if(cnf.agent_flow_method == Config::agent_flow_sinkhorn) {
            tokens = sinkhorn(V);
            for( auto& d : tokens ) d = 1.0 - d;
        }

        sorted_nodes = VI(N,0);
        iota(ALL(sorted_nodes),0);
        sort(ALL(sorted_nodes), [&](int a, int b){
            return tokens[a] < tokens[b];
        });

        VI valid_nodes;
        for( int i=0; i<N; i++ ) if( !was[ sorted_nodes[i] ] ) valid_nodes.push_back(sorted_nodes[i]);
        swap(valid_nodes, sorted_nodes);

        next_node_smallest = 0;
        next_node_largest = (int)sorted_nodes.size()-1;
    };

    initialize_tokens();

    int iteration_cnt = 0;

    const bool use_reductions = true;

    const int MAX_EDGES = 5e5;
    int E = 0;

    while(next_node_smallest < sorted_nodes.size() && next_node_largest >= 0 ){

        if(use_reductions && (iteration_cnt % cnf.agent_flow_node_update_frequency == 0) ){ // using reductions
            Reducer red(V, cnf);
            if(!cnf.solverh_use_reductions_AF) red.cnf.disableAllNonbasicReductions();

            // #CAUTION! No conditional reductions can be used here
            red.cnf.disableAllConditionalReductions();
            auto reductions = red.reduce();
            VI dfvs_red = Reducer::convertKernelizedReductions(reductions);

            V = red.V;
            revV = GraphUtils::reverseGraph(V);
            res += dfvs_red;

            E = GraphUtils::countEdges(V);
            if(E == 0) break;

            initialize_tokens();
        }

        if( (iteration_cnt % cnf.agent_flow_node_update_frequency == 0) && cnf.agent_flow_alternate_selection_type ){
            cnf.agent_flow_node_selection_type = Config::agent_flow_remove_largest_flow_node;
        }else{
            cnf.agent_flow_node_selection_type = Config::agent_flow_merge_smallest_flow_node;
        }

        if(cnf.agent_flow_node_selection_type == Config::agent_flow_remove_largest_flow_node) {
            while(next_node_largest >= 0){
                int a = sorted_nodes[next_node_largest];
                if(was[a]){
                    next_node_largest--;
                    continue;
                }
                was[a] = true;

                if( V[a].empty() ){
                    next_node_largest--;
                    continue;
                }
                if( Utils::hasLoop(V,a) ){
                    next_node_largest--;
                    continue;
                }
                break;
            }

            if(next_node_largest >= 0){
                Utils::removeNode(V, revV, sorted_nodes[next_node_largest], helper);
                res.push_back(sorted_nodes[next_node_largest]);
            }

            next_node_largest--;
        }else if(cnf.agent_flow_node_selection_type == Config::agent_flow_merge_smallest_flow_node){
            while(next_node_smallest < sorted_nodes.size()){
                int a = sorted_nodes[next_node_smallest];
                if(was[a]){
                    next_node_smallest++;
                    continue;
                }
                was[a] = true;

                if( V[a].empty() ){
                    next_node_smallest++;
                    continue;
                }
                if( Utils::hasLoop(V,a) ){
                    next_node_smallest++;
                    continue;
                }
                break;
            }

            if(next_node_smallest < sorted_nodes.size()) Utils::merge(V, revV, sorted_nodes[next_node_smallest], helper);

            next_node_smallest++;
        }

        iteration_cnt++;

        if(cnf.solver_improve_alternate_selection_type){ // changing node selection method
            cnf.agent_flow_node_selection_type = 1 - cnf.agent_flow_node_selection_type;
        }
    }

    if(cnf.agent_flow_node_selection_type == Config::agent_flow_merge_smallest_flow_node){
        for(int i=0; i<N; i++) if(Utils::hasLoop(V,i)) res.push_back(i);
    }

    reverse(ALL(res));
    VI res2 = Utils::findAndRemoveRedundantNodes(_V0, res);
    res = res2;

    return res;
}

VI DFVSSolverH::solveByAgentFlowAllWithinDistance(VVI V, const int minimize) {
    VVI _V0 = V;

    int N = V.size();
    VVI revV = GraphUtils::reverseGraph(V);
    VB helper(N,false);

    VI res;

    VD tokens;
    VB was(N,false);
    int next_node_largest = N-1;

    VI sorted_nodes;

    auto initialize_tokens = [&](){
        AgentFlow af(V);
        int iterations = 2*sqrt(N);
        int void_steps = sqrt(N);

        if(cnf.agent_flow_method == Config::agent_flow_continuous) tokens = af.flowContinuous(iterations, void_steps);
        else if(cnf.agent_flow_method == Config::agent_flow_tokens) tokens = af.flowTokens(iterations, void_steps, 3);
        else if(cnf.agent_flow_method == Config::agent_flow_sinkhorn) {
            tokens = sinkhorn(V);
            for( auto& d : tokens ) d = 1.0 - d;
        }

        sorted_nodes = VI(N,0);
        iota(ALL(sorted_nodes),0);
        sort(ALL(sorted_nodes), [&](int a, int b){
            return tokens[a] < tokens[b];
        });

        if(cnf.agent_flow_node_selection_type == Config::agent_flow_merge_smallest_flow_node){
            reverse(ALL(sorted_nodes));
        }

        VI valid_nodes;
        for( int i=0; i<N; i++ ) if( !was[ sorted_nodes[i] ] ) valid_nodes.push_back(sorted_nodes[i]);
        swap(valid_nodes, sorted_nodes);

        next_node_largest = (int)sorted_nodes.size()-1;
    };

    auto findNodesToRemove = [&](){
        VB is_free(N,false);
        for(int d : sorted_nodes) is_free[d] = true;

        VVI undirV = V;
        for( int i=0; i<N; i++ ) undirV[i] += revV[i]; // creating undirected V

        VI to_remove;

        int max_dist = cnf.agent_flow_max_distance_from_best;
        max_dist = min( 1 + (int)sqrt(sorted_nodes.size()), max_dist );

        if(cnf.tle()) max_dist = 1 + 5*(int)sqrt(sorted_nodes.size()); // finish quickly

        auto cond1 = [&]() {
            return next_node_largest >= 0 && sorted_nodes.size() - next_node_largest <= max_dist;
        };

        while( cond1() ){
            while( cond1() && !is_free[sorted_nodes[next_node_largest]] ) next_node_largest--;
            if(!cond1()) break;

            int v = sorted_nodes[next_node_largest];
            to_remove.push_back(v);

            // now run BFS from v in undirV, within distance cnf.agent_flow_min_distance
            int dst = 1;
            VI sources = {v};
            helper[v] = true;
            is_free[v] = false;

            VI visited_nodes = {v};
            VI temp;
            while( dst < cnf.agent_flow_min_distance ){
                temp.clear();
                for( int u : sources ){
                    for( int w : undirV[u] ){
                        if(!helper[w]){
                            helper[w] = true;
                            is_free[w] = false;
                            temp.push_back(w);
                        }
                    }
                }

                visited_nodes += temp;
                swap(sources,temp);
                dst++;
            }

            for(int d : visited_nodes) helper[d] = false; // clearing helper
        }

        return to_remove;
    };

    int iteration_cnt = 0;
    const bool use_reductions = true;

    int E = 0;

    while(next_node_largest >= 0 ){

//        if(cnf.write_logs) clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;

        Reducer red(V, cnf);
        if(!cnf.solverh_use_reductions_AF) red.cnf.disableAllNonbasicReductions();

        // #CAUTION! No conditional reductions can be used here
        red.cnf.disableAllConditionalReductions();
        auto reductions = red.reduce(revV);
        VI dfvs_red = Reducer::convertKernelizedReductions(reductions);

        V = red.V;
        revV = red.revV;

        res += dfvs_red;
        for(int d : dfvs_red) was[d] = true;

        E = GraphUtils::countEdges(V);
        if(E == 0) break;

        initialize_tokens();

//        TimeMeasurer::start("findNodesToRemove");
        VI selected_nodes = findNodesToRemove();
//        TimeMeasurer::stop("findNodesToRemove");

        if(cnf.agent_flow_node_selection_type == Config::agent_flow_merge_smallest_flow_node){
            for( int v : selected_nodes ){
                if( !Utils::hasLoop(V,v) ){
                    was[v] = true;
                    Utils::merge(V, revV, v, helper);
                }
            }
        }else{
            res += selected_nodes;
            for(int v : selected_nodes){
                was[v] = true;
                Utils::removeNode(V, revV, v, helper);
            }
        }

        iteration_cnt++;
    }

    if(cnf.agent_flow_node_selection_type == Config::agent_flow_merge_smallest_flow_node){
        for(int i=0; i<N; i++) if(Utils::hasLoop(V,i)) res.push_back(i);
    }

    if(cnf.tle()) return res; // return without minimalization if time limit exceeded

    if(minimize == 1) { // minimize the result
        reverse(ALL(res));
        VI res2 = Utils::findAndRemoveRedundantNodes(_V0, res);
        res = res2;
    }else if( minimize >= 2 ){
        VI dfvs = res;
        reverse(ALL(dfvs));
        int max_check_length = dfvs.size();
        VI invalid_dfvs = Utils::findAndRemoveRedundantNodes(_V0, dfvs, max_check_length, false);

        VI nodes = CombinatoricUtils::getFullSetDifference(N, invalid_dfvs);
        InducedGraph g = GraphInducer::induce(_V0, nodes);
        VI dfvs_g = solveByAgentFlowAllWithinDistance(g.V, minimize-1); // minimize once

        for (int &d : dfvs_g) d = g.nodes[d];
        dfvs = invalid_dfvs + dfvs_g;
        dfvs = Utils::findAndRemoveRedundantNodes(_V0, dfvs, max_check_length, true);

        res = dfvs;
    }

    return res;
}

VI DFVSSolverH::improve1(VVI V, VI dfvs, int iterations, double alpha) {
    VI best = dfvs;
    using namespace StandardUtils;

    VB helper(V.size(),false);

    clog << "Before improvement dfvs.size(): " << dfvs.size() << endl;

    for( int iter=1; iter<=iterations; iter++ ){

        if(cnf.solver_improve_alternate_selection_type) cnf.agent_flow_node_selection_type = 1 - cnf.agent_flow_node_selection_type;
        shuffle(dfvs);

        int a = ceil(alpha * dfvs.size());
        VI X = slice(dfvs,0,a); // this set of nodes will remain as the DFVS.

        VVI W = V;
        VVI revW = GraphUtils::reverseGraph(W);
        Utils::removeNodes(W,revW,X, helper );

        VI dfvs_red;
        {
            Reducer red(W,cnf);

            // #CAUTION! No conditional reductions can be used here
            red.cnf.disableAllConditionalReductions();
            auto reductions = red.reduce();
            VI dfvs_red = Reducer::convertKernelizedReductions(reductions);

            W = red.V;
        }

        VI dfvsW = solveByAgentFlowAllWithinDistance(W);

        VCImprover impr(cnf);
        impr.cnf.vc_improver_milliseconds = 10 + cnf.vc_improver_milliseconds
                * (1.0 * GraphUtils::countEdges(W,true) / GraphUtils::countEdges(V,true) ) ;

        for( int i=0; i<10; i++ ){
            dfvsW = impr.improveByVertexCoverMinimizeArcs(W, dfvsW, sqrt(W.size()));
            if(cnf.write_logs) DEBUG(dfvsW.size());
            assert( Utils::isFVS(W, dfvsW) );
        }


        dfvsW += dfvs_red;
        VI new_dfvs = X+dfvsW;
        assert( Utils::isFVS(V,new_dfvs ) );
        new_dfvs = Utils::findAndRemoveRedundantNodes(V, new_dfvs);

        assert( Utils::isFVS(V,new_dfvs ) );

        if(cnf.write_logs) DEBUG(new_dfvs.size());

        if( new_dfvs.size() <= best.size() ){
            best = new_dfvs;
            dfvs = best;
        }

        if(cnf.write_logs) clog << "Violating local optimum" << endl;
        return new_dfvs;
    }

    if(cnf.write_logs) clog << "After improvement best.size(): " << best.size() << endl << endl;

    // best is already minimized, we do not need to minimize it now
    return best;
}


VD DFVSSolverH::sinkhorn(VVI &V) {
    int N = V.size();
    int E = GraphUtils::countEdges(V);

    vector<tuple<int,int,double>> entries;
    entries.reserve(2*E);
    for( int i=0; i<N; i++ ){
        entries.emplace_back(i,i,1);
        for( int d : V[i] ) if(d!=i) entries.emplace_back(i,d,1);
    }

    VD rows(N, 0), cols(N, 0);

    auto normalize = [&](){
        fill(ALL(rows), 0);
        fill(ALL(cols), 0);
        for( auto & [r,c,d] : entries ) rows[r] += d;
        for( auto & [r,c,d] : entries ){
            d /= rows[r];
            cols[c] += d;
        }
        for( auto & [r,c,d] : entries ) d /= cols[c];
    };

    int reps = 10 + log(N);

    while(reps--){
        normalize();
    }

    VD res(N,0);
    for( auto & [r,c,d] : entries ) if(r==c) res[r] = d;
    return res;
}

VI DFVSSolverH::solveByLowerBounding(VVI V) {
    const bool debug = false;
    VVI _V0 = V;

    int E = GraphUtils::countEdges(V,true);
    VI dfvs;
    VVI revV = GraphUtils::reverseGraph(V);

    int N = V.size();
    VB helper(N,false);

    while( E > 0 ){
        if(debug){
            ENDL(1);
            clog << "Edges left in the graph: " << E << endl;
            clog << "Current dfvs.size(): " << dfvs.size() << endl;
        }

        {
            Reducer red(V, cnf);

            // #CAUTION! No conditional reductions can be used here
            red.cnf.disableAllConditionalReductions();
            auto reductions = red.reduce(revV);
            VI red_dfvs = Reducer::convertKernelizedReductions(reductions);

            dfvs += red_dfvs;
            V = red.V;
            revV = red.revV;
            E = GraphUtils::countEdges(V,true);
        }


        VVI cycles = Utils::getSmallLowerBoundCycles(V);
        if(cycles.empty()) break;

        VVI cycV = Utils::getCycleAdjacencyGraph(cycles);
        int millis = 500;
        VI lb_cyc = Utils::getLowerBoundCycles(cycV, millis);

        if(debug) DEBUG(lb_cyc.size());

        VD tokens = sinkhorn(V);

        VI to_remove;

        for( int c : lb_cyc ){
            double m = 1e9;
            int best_v = -1;
            for( int v : cycles[c] ){
                if( tokens[v] < m ){
                    m = tokens[v];
                    best_v = v;
                }
            }

            to_remove.push_back(best_v);
        }

        dfvs += to_remove;

        Utils::removeNodes(V, revV, to_remove, helper);

    }

    if(E > 0){
        if(debug) clog << "Graph processing unfinished, left edges: " << E << endl;
        if(debug) clog << "Current dfvs.size(): " << dfvs.size() << endl;

        VI nodes;
        for(int i=0; i<N; i++) if(!V[i].empty()) nodes.push_back(i);
        InducedGraph g = GraphInducer::induce(V,nodes);

        if(debug) clog << "Processing the remaining graph (without cycles of length <=3) using standard AF" << endl;
        if(debug) clog << "g.V.size(): " << g.V.size() << ", g.V.edges: " << E << endl;

        DFVSSolverH solver(cnf);
        VI dfvs_rem = solver.solveByAgentFlowAllWithinDistance(g.V);
        assert(Utils::isFVS(g.V,dfvs_rem));

        for( auto & d : dfvs_rem ) d = g.nodes[d];
        if(debug) clog << "Found DFVS of g.V of size " << dfvs_rem.size() << endl;


        dfvs += dfvs_rem;
    }

    assert(Utils::isFVS(_V0,dfvs));

    if(debug) clog << "Before minimization, dfvs.size(): " << dfvs.size() << endl;
    dfvs = Utils::findAndRemoveRedundantNodes(_V0, dfvs);
    if(debug) clog << "After minimization, dfvs.size(): " << dfvs.size() << endl;

    assert(Utils::isFVS(_V0,dfvs));

    return dfvs;
}

VI DFVSSolverH::solveByVCForAllPIGraph(VVI V, int milliseconds) {
    int N = V.size();
    bool use_fastvc = ( N > 20'000 );

    VI init_vc; // empty
    VI mis;
    VPII edges = GraphUtils::getGraphEdges(V);

    UniformIntGenerator rnd(0,1e6);
    if(use_fastvc){
        int seed = rnd.rand();

        libmvc::FastVC fastvc(edges, N, 0, std::chrono::milliseconds(milliseconds), false, seed, init_vc);
        fastvc.cover_LS();
        mis = fastvc.get_independent_set(false);

        assert(VCUtils::isIndependentSet(V, mis));
    }
    else{ // using NuMVC
        NuMVC numvc;
        mis = numvc.solve(V, 1.0 * milliseconds / 1000, init_vc); // CAUTION! NuMVC indices are from 1 to N
        for (int &d : mis) d--;

        assert(VCUtils::isIndependentSet(V, mis));
    }

    VI vc = CombinatoricUtils::getFullSetDifference(N, mis);
    return vc;
}

VI DFVSSolverH::solveForStronglyConnectedIterativeHittingSet2(VVI V, VI heur_dfvs, VVI & all_cycles, double F) {
    if( cnf.tle() ) return heur_dfvs;

    {
        double millis_done = cnf.sw.getTime("main");
        int millis_left = cnf.sw.getLimit("main") - millis_done;
        LL short_cycles = Utils::countInducedCycles( V, cnf.ihs_init_cycle_length, millis_left );
        LL threshold = 5'000'000;
        if( short_cycles > threshold ) return heur_dfvs;
    }


    int upper_bound = heur_dfvs.size();

    int N = V.size();
    VB helper(N,false);

    VI hs;
    int iter = 0;
    int max_cycle_length = cnf.ihs_init_cycle_length;
    all_cycles.clear();

    VLL hashes(N);
    constexpr int seed = 910238;
    UniformIntGenerator rnd(0, 1'000'000'000ll * 1'000'000'000, seed);
    for(int i=0; i<N; i++) hashes[i] = rnd.rand();
    auto getHSHash = [&](){
        LL hash = 0;
        for(int d : hs) hash ^= hashes[d];
        return hash;
    };


    int ITR_FIRST = 5e4; // checks done trying to make existing hs a valid HS, before adding nodes to hit unhit cycles
    int ITR_INDUCE = 2e5; // original
    int ITR_SINGLE_NODE = 3e4; // checks done after each addition of a single node
    int ITR_SECOND = 3e5; // checks done after making current HS a valid HS
    int ITR_FINAL = 1e6; // checks done if a valid HS-DFVS is found
    int MAX_CYCLES = 10; // original


    auto rescaleITRS = [&](double F) {
        ITR_FINAL *= F;
        ITR_INDUCE *= F;
        ITR_FIRST *= F;
        ITR_SINGLE_NODE *= F;
        ITR_SECOND *= F;
        MAX_CYCLES = min( (int)V.size() / 10, 100);
        if( V.size() > 10'000 ) MAX_CYCLES = min( (int)V.size() / 10, 200);

        MAX_CYCLES = min(MAX_CYCLES, cnf.solverh_ihs_secondary_max_cycles);
    };

    int max_rescaling_times = cnf.solverh_ihs_max_rescaling_times;

//    const bool deviate_perm = (F > 0.2);
    const bool deviate_perm = true;

    rescaleITRS(F);

//    MAX_CYCLES = 300; // this values seems to work very well!!
    MAX_CYCLES = cnf.solverh_ihs_init_max_cycles;
    {
        double perc = Utils::getPieEdgesPercentage(V);
        if( perc <= 0.15 && !heur_dfvs.empty() && heur_dfvs.size() > 300 ){
            MAX_CYCLES = 2 * heur_dfvs.size();
        }
    }

    DEBUG(MAX_CYCLES);

//    int previous_cycle_search_millis = 0;
    int previous_cycle_search_millis = 1e9;
    const int max_cycle_search_time_millis = 1'000; // 300 seems to work ok

    VVI unhit_cycles;


    do{
//        if( cnf.sw.tle("main") ) break;
        if( cnf.tle() ) break;

        iter++;

        if(cnf.write_logs) {
            DEBUG(iter);
            clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;
            DEBUG(hs.size());
            DEBUG(all_cycles.size());
            DEBUG(max_cycle_length);
            DEBUG(getHSHash());
        }

        VB in_V(N,false);
        VB is_end(N,false);
        VI marker(N,0);
        VI A = CombinatoricUtils::getRandomPermutation(N);

        unhit_cycles.clear();

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
                                                    is_end, max_cycle_length, 1500);
            }
            else cycles = Utils::getAllSimpleCycles3(H, max_cycle_length, max_cycle_search_time_millis,false); // original
//            else cycles = Utils::getAllSimpleCycles4(H, max_cycle_length, max_cycle_search_time_millis,true);

            sw.stop("getAllSimpleCycles");
            if(cnf.write_logs) sw.write("getAllSimpleCycles");

            { // this should make it work almost deterministically (minimizes the risk of nondeterministic run)
                int cnt = 0;
                while( cycles.empty() ) {
                    cycles = Utils::getAllSimpleCycles3(H, max_cycle_length, max_cycle_search_time_millis, true); // original
//                    cycles = Utils::getAllSimpleCycles4(H, max_cycle_length, max_cycle_search_time_millis, true);
                    if(cnt++ == 2){
                        max_cycle_length++;
                        cnt = 0;
                    }
                }
            }

            if(previous_cycle_search_millis <= max_cycle_search_time_millis){
                previous_cycle_search_millis = sw.getTime("getAllSimpleCycles");
            }

            const int max_cycles_size = (iter == 1 ? (int) cycles.size() : MAX_CYCLES);

            if (cycles.size() > max_cycles_size) {
                StandardUtils::shuffle(cycles);
                sort(ALL(cycles), [](auto &v, auto &w) { return v.size() < w.size(); });
//                if(iter&1) reverse(ALL(cycles)); // #TEST - take sometimes longest cycles
                cycles.resize(max_cycles_size);
            }
            if(cnf.write_logs) DEBUG(cycles.size());
            if(cycles.empty()){
                max_cycle_length++;
                if(cnf.write_logs) clog << "Increasing cycle length" << endl;
                continue;
            }
            else if( cycles.size() <= ceil(1.0 * max_cycles_size / min(6,max_cycle_length)) ) max_cycle_length++;

            all_cycles += cycles;
            unhit_cycles = cycles;
            break;
        }

        const bool induce_hs_from_heur_dfvs = true;
        if(hs.empty() && induce_hs_from_heur_dfvs){ // induce solution from heur_dfvs
            if(cnf.write_logs) clog << "Inducing hs from heur_dfvs" << endl;
            VB was(N,false);
            for( int d : heur_dfvs ) was[d] = true;
            unordered_set<int> zb;

            for( VI & C : all_cycles ){
                bool hit = false;
                for( int d : C ) if(zb.count(d)) hit = true;
                if(!hit){
                    for(int d : C){
                        if(was[d]){
                            zb.insert(d);
                            break;
                        }
                    }
                }
            }
            hs = VI(ALL(zb));
            if(cnf.write_logs) DEBUG(hs.size());
            VI hsls = Utils::hsImprovementLS2( all_cycles, hs, cnf, ITR_INDUCE );
            if(!hsls.empty()) hs = hsls;
            if(cnf.write_logs) DEBUG(hs.size());

            const bool remove_some_elements = true;
            if(remove_some_elements) {
                hs.resize(max(1, (int) hs.size() - 50));

                { // creating unhit cycles
                    unhit_cycles.clear();
                    VB in_hs = StandardUtils::toVB(N, hs);
                    for (auto &C : all_cycles) {
                        bool hit = false;
                        for (int d : C)
                            if (in_hs[d]) {
                                hit = true;
                                break;
                            }
                        if (!hit) unhit_cycles.push_back(C);
                    }
                }
            }
        }

        bool found_hs_by_ls = false;

        const bool use_hs_ls_search = ( all_cycles.size() > 100 );
        if(use_hs_ls_search){
            if(cnf.write_logs) clog << "Trying to find 'equivalent' hs using local search" << endl;

            VI hs_ls;
            hs_ls = Utils::hsImprovementLS2( all_cycles, hs, cnf, ITR_FIRST, (int)hs.size()-1 ); // enabling decreasing size by 1

            if( !hs_ls.empty() ){
                found_hs_by_ls = true;
                hs = hs_ls;
                if(cnf.write_logs) clog << "\tFOUND HS BY LOCAL SEARCH, hs.size(): " << hs.size() << endl;
            }
        }

        auto makeHsAValidHittingSet = [&](){
            if(cnf.write_logs) clog << "Making hs a valid HS" << endl;
            VB in_hs = StandardUtils::toVB(N, hs);

            int nodes_added = 0;
            if(cnf.write_logs) clog << "Now trying to greedily hit unhit cycles" << endl;

            while( !unhit_cycles.empty() ) {
                VI freqs(N,0);
                for (auto &C : unhit_cycles) {
                    for (int d : C) freqs[d]++;
                }

                int best_el;
                { // taking random best_el from all elements with maximal freqs, not the first one
                    VI best_elts;
                    int val = -1;
                    for (int i = 0; i < N; i++) {
                        if( freqs[i] > val ) best_elts.clear();
                        if(freqs[i] >= val){
                            val = freqs[i];
                            best_elts.push_back(i);
                        }
                    }
                    best_el = best_elts[ rnd.nextInt(best_elts.size()) ];
                }

                hs.push_back(best_el);
                in_hs[best_el] = true;
                nodes_added++;

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


                if(cnf.solverh_use_hsls_after_each_node_addition){ // trying to find a solution after adding just single node
                    VI xyz = Utils::hsImprovementLS2(all_cycles, hs, cnf, ITR_SINGLE_NODE, hs.size()-1);
                    if (!xyz.empty()) {
                        hs = xyz;
                        break;
                    }
                }
            }

            if(cnf.write_logs) clog << "Greedily nodes selected, hs.size(): " << hs.size() << ", improving using LS" << endl;
            if( Utils::isFVS(V,hs) && hs.size() < heur_dfvs.size()){
                if(cnf.write_logs) clog << "Found new best results of size " << hs.size() << " just after greedy node addition" << endl;
                heur_dfvs = hs;
                upper_bound = heur_dfvs.size();
            }

            { // trying to improve solution just after adding nodes, before checking if HS is a DFVS
                VI xyz = Utils::hsImprovementLS2(all_cycles, hs, cnf, ITR_SECOND, (int) hs.size() - 2*nodes_added - 1, deviate_perm );
                if( !xyz.empty() ) hs = xyz;
            }

            if(cnf.write_logs) clog << "Improved, hs.size(): " << hs.size() << endl;
        };

        // take greedily nodes that hits most unhit cycles to hs
        if(!found_hs_by_ls) makeHsAValidHittingSet();

        if(cnf.write_logs) {
            ENDL(1);
            DEBUG(upper_bound);

            DEBUG(hs.size());
            ENDL(5);
            ENDLS(50, "*");
        }

        if(Utils::isFVS(V,hs)){
            if(cnf.write_logs) clog << "Found HS is a valid DFVS, trying to improve it further" << endl;
            if( hs.size() < heur_dfvs.size() ){
                heur_dfvs = hs; // writing this not to lost
                upper_bound = heur_dfvs.size();
            }
            VI xyz = Utils::hsImprovementLS2(all_cycles, hs, cnf, ITR_FINAL );
            if( !xyz.empty() && xyz.size() < hs.size() ) hs = xyz;

            while(Utils::isFVS(V,hs) && max_rescaling_times){
//                if( cnf.sw.tle("main") ) break;
                if( cnf.tle() ) break;

                if(hs.size() < heur_dfvs.size()){
                    heur_dfvs = hs;
                    upper_bound = heur_dfvs.size();
                }
                max_rescaling_times--;
                const double FF = 1.25;
                if(cnf.write_logs){
                    clog << "Rescaling params by factor FF: " << FF << " and trying to improve, rescaling times left: "
                         << max_rescaling_times << endl;
                    clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;
                }
                rescaleITRS(FF);
                VI xyz = Utils::hsImprovementLS2(all_cycles, hs, cnf, ITR_FINAL); // original
//                VI xyz = Utils::hsImprovementLS2(all_cycles, hs, cnf, ITR_FINAL, -1, true,
//                                                 3 ); // using persistent search!
                if( !xyz.empty() && xyz.size() < hs.size() ){
                    hs = xyz;
                }
//                else if(max_rescaling_times && xyz.empty()){
//                    max_rescaling_times = 1;
//                }
                if(cnf.write_logs) clog << "hs.size(): " << hs.size() << endl;
            }
        }

    }while( !Utils::isFVS(V,hs) );

    if( !Utils::isFVS(V,hs) ){
        assert(cnf.tle());
        hs = heur_dfvs;
        assert( Utils::isFVS(V, heur_dfvs) );
    }

    if( Utils::isFVS(V,hs) && hs.size() > heur_dfvs.size() ) hs = heur_dfvs;

    if(cnf.write_logs) clog << "Found heuristically a DFVS of size: " << hs.size() << endl;
    assert(Utils::isFVS(V,hs));

    return hs;
}




