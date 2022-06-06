//
// Created by sylwester on 1/7/22.
//

#include <CONTESTS/PACE22/Reducer.h>
#include <CONTESTS/PACE22/heur/VCImprover.h>
#include <graphs/GraphUtils.h>
#include <combinatorics/CombinatoricUtils.h>
#include <CONTESTS/PACE22/heur/SALS3.h>
#include <CONTESTS/PACE22/heur/SALS4.h>
#include <CONTESTS/PACE22/heur/SALS5.h>
#include <graphs/GraphInducer.h>
#include <CONTESTS/PACE22/heur/AgentFlow.h>
#include "CONTESTS/PACE22/heur/DFVSImprover.h"
#include "CONTESTS/PACE22/heur/DFVSSolverH.h"
#include "StandardUtils.h"

DFVSImprover::DFVSImprover(VVI V, Config c) {
    cnf = c;
    this->V = V;
}

VI DFVSImprover::improve(VI dfvs, const int iterations) {
    int N = V.size();

    best_dfvs = dfvs;

    int iterations_without_improvement = 0;
    int prev_size = dfvs.size();

    VCImprover impr(cnf);

    for(int i=0; i<iterations; i++){
        if(cnf.tle()) break;

        dfvs = impr.improveByVertexCoverMinimizeArcs(V, dfvs, 1 + sqrt(dfvs.size())); // original

        if(dfvs.size() < best_dfvs.size()) best_dfvs = dfvs;

        if(cnf.write_logs){
            clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;
            DEBUG(dfvs.size() );
        }

        if(prev_size <= dfvs.size()) iterations_without_improvement++;
        else iterations_without_improvement = 0;

        if(cnf.write_logs) DEBUG(iterations_without_improvement);

        if( ((i+1) % cnf.dfvsimprover_local_optimum_violation_frequency) == 0
            || iterations_without_improvement == cnf.dfvsimprover_max_iters_without_improvement ){

            iterations_without_improvement = 0;

            if(dfvs.size() < best_dfvs.size()) best_dfvs = dfvs;
            else if( dfvs.size() > best_dfvs.size() + cnf.dfvsimprover_local_opt_max_deviation_from_best_absolute_addition ){
                dfvs = best_dfvs;
                if(cnf.write_logs) clog << "Getting out of local optimum using best_dfvs, not current one" << endl;
            }

            if(cnf.write_logs) clog << "Trying to get out of local optimum" << endl;
            VI temp = getOutOfLocalOptimum(dfvs, cnf.dfvsimprover_alpha);
            int max_size = cnf.dfvsimprover_local_opt_max_deviation_from_best_absolute_addition +
                           ceil( dfvs.size() * (1 + cnf.dfvsimprover_local_opt_max_deviation_from_best_relative) );
            if(cnf.write_logs) DEBUG(max_size);
            if(cnf.write_logs) DEBUG(max_size - dfvs.size());

            if( temp.size() <= max_size ) {
                dfvs = temp;
                if(cnf.write_logs) clog << "Solution size after violating local optimum: " << dfvs.size() << endl;
            }
            else{
                if(cnf.write_logs) clog << "Local optimum of size " << temp.size() << " is to far from current solution, NOT getting out" << endl;
            }

            if(cnf.write_logs) clog << endl;
        }

        prev_size = dfvs.size();
    }

    return best_dfvs;
}

VI DFVSImprover::getOutOfLocalOptimum(VI dfvs, double alpha){
    VI best = dfvs;
    using namespace StandardUtils;

    VB helper(V.size(),false);

    if(cnf.solver_improve_alternate_selection_type) cnf.agent_flow_node_selection_type = 1 - cnf.agent_flow_node_selection_type;
    shuffle(dfvs);

    int a = ceil(alpha * dfvs.size());
    VI X = slice(dfvs,0,a);

    VVI W = V;
    VVI revW = GraphUtils::reverseGraph(W);
    Utils::removeNodes(W,revW,X, helper );

    VI dfvs_red;
    {
        Reducer red(W,cnf);

        red.cnf.disableAllConditionalReductions();
        auto reductions = red.reduce();
        VI red_dfvs = Reducer::convertKernelizedReductions(reductions);

        dfvs_red = red_dfvs;
        W = red.V;
    }

    DFVSSolverH solver(cnf);
    VI dfvsW = solver.solveByAgentFlowAllWithinDistance(W);

    VCImprover impr(cnf);
    impr.cnf.vc_improver_milliseconds = 10 + cnf.vc_improver_milliseconds
                                             * (1.0 * GraphUtils::countEdges(W,true) / GraphUtils::countEdges(V,true) ) ;

    for( int i=0; i<10; i++ ){
        dfvsW = impr.improveByVertexCoverMinimizeArcs( W, dfvsW, sqrt(W.size()) );
        assert( Utils::isFVS(W, dfvsW) );
    }


    dfvsW += dfvs_red;
    VI new_dfvs = X+dfvsW;
    new_dfvs = Utils::findAndRemoveRedundantNodes(V, new_dfvs);

    if( new_dfvs.size() <= best.size() ){
        best = new_dfvs;
        dfvs = best;
    }

    return new_dfvs;

}

VI DFVSImprover::improveBySALSForSparseGraphs( VI dfvs, int iterations, double alpha, double sals_T0,
      double sals_alpha, int sals_maxmvt, int sals_maxfail) {

    VI best_dfvs = dfvs;
    int N = V.size();
    assert( alpha >=0 && alpha <=1 );
    VB helper(N,false);

    auto mergeNodes = [&]( VVI & W, VI & nodes ){
        VVI revW = GraphUtils::reverseGraph(W);
        int merged = 0;
        for( int d : nodes ){
            if( !Utils::hasLoop(W,d) ){
                Utils::merge(W, revW, d, helper);
                merged++;
            }
        }

        {
            VI loops;
            for (int i = 0; i < N; i++) if (Utils::hasLoop(W,i)) loops.push_back(i);
            for(int d : loops) Utils::removeNode(W, revW, d, helper);
            return loops;
        }
    };

    double density = GraphUtils::density(V, true);
    const bool use_sals3 = ( density > cnf.min_density_for_sals3 );

    for(int iter = 1; iter <= iterations; iter++){
        if(cnf.tle()) break;
        if(cnf.write_logs) clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;

        VI order = CombinatoricUtils::getFullSetDifference(N,best_dfvs);
        StandardUtils::shuffle(order);

        order.resize( ceil( alpha * order.size() ) );

        VVI W = this->V;

        VI red_dfvs = mergeNodes(W, order);
        if(cnf.write_logs){
            DEBUG(red_dfvs.size());
            DEBUG(best_dfvs.size());
            DEBUG( Utils::getPieEdgesPercentage(W))
        }

        VI sals_dfvs;

        if(use_sals3) {
            SALS3 sals3( W, best_dfvs, cnf );
            sals3.cnf.write_logs = false;
            double T0 = SALS3::findInitialTemperature(W, best_dfvs, cnf, 10 * W.size());
            if(T0 > 0) sals_dfvs = sals3.localSearch(T0, 0.97, 25*V.size(), 5);
        }else{
            SALS2 sals(W, best_dfvs, cnf, true);
            sals.cnf.write_logs = false;
            sals_dfvs = sals.localSearch(sals_T0, sals_alpha, sals_maxmvt, sals_maxfail);
        }

        DEBUG(sals_dfvs.size());

        VI new_dfvs = sals_dfvs;
        for( int d : new_dfvs ) helper[d] = true;
        for( int d : red_dfvs ) if(!helper[d]) new_dfvs.push_back(d);
        for( int d : new_dfvs ) helper[d] = false;

        DEBUG(new_dfvs.size());
        if( new_dfvs.size() <= best_dfvs.size() ){
            best_dfvs = new_dfvs;
            assert(Utils::isFVS(V, best_dfvs));
        }

        ENDL(10);
        clog << "After iterations #" << iter << ", best_dfvs.size(): " << best_dfvs.size() << endl;
        ENDL(10);
    }

    return best_dfvs;
}


VI DFVSImprover::improveBySALSForSparseGraphsDensifier( VI dfvs, int iterations, double init_density, double step,
                                               double sals_T0, double sals_alpha, int sals_maxmvt_factor, int sals_maxfail) {

    VI best_dfvs = dfvs;
    int N = V.size();
    VB helper(N,false);

    for(int iter = 1; iter <= iterations; iter++){
        if(cnf.tle()) break;
        if(cnf.write_logs) clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;

        VVI W = this->V;

        VI red_dfvs = mergeNodesUntilDensityReached(W, init_density + (iter-1)*step);
        if(cnf.tle()) break;

        InducedGraph g = GraphInducer::induceByNonisolatedNodes(W);
        W = g.V;
        DEBUG(GraphUtils::density(W,true));
        DFVSSolverH solverh(cnf);
        solverh.cnf.write_logs = false;

        solverh.cnf.agent_flow_max_distance_from_best = sqrt(W.size()); // #TEST
        solverh.cnf.agent_flow_node_selection_type = Config::agent_flow_remove_largest_flow_node;
        solverh.cnf.disableAllNonbasicReductions();

        VI mapped_dfvs = solverh.solveByAgentFlowAllWithinDistance(W);
        if(cnf.tle()) break;

        assert(Utils::isFVS(W,mapped_dfvs));

        if(cnf.write_logs){
            DEBUG(red_dfvs.size());
            DEBUG(best_dfvs.size());
            DEBUG( Utils::getPieEdgesPercentage(W))
        }
        clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;



        SALS2 sals(W, mapped_dfvs, cnf, true);
        sals.cnf.write_logs = false;
        int sals_maxmvt = sals_maxmvt_factor * W.size();
        VI sals_dfvs = sals.localSearch(sals_T0, sals_alpha, sals_maxmvt, sals_maxfail);

        DEBUG(sals_dfvs.size());

        for(int & d : sals_dfvs) d = g.nodes[d]; // remapping nodes

        VI new_dfvs = sals_dfvs;
        for( int d : new_dfvs ) helper[d] = true;
        for( int d : red_dfvs ) if(!helper[d]) new_dfvs.push_back(d);
        for( int d : new_dfvs ) helper[d] = false;

        DEBUG(new_dfvs.size());
        if( new_dfvs.size() <= best_dfvs.size() ){
            best_dfvs = new_dfvs;
            assert(Utils::isFVS(V, best_dfvs));
        }

        ENDL(10);
        clog << "After iterations #" << iter << ", best_dfvs.size(): " << best_dfvs.size() << endl;
        clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;
        ENDL(10);
    }

    return best_dfvs;
}

VI DFVSImprover::fullImprovement(VI dfvs, int origV_N, int origV_edges, const bool use_vc_improvement) {
    int N = V.size();

    if(use_vc_improvement){ // improving
        DFVSImprover impr(V, cnf);

        impr.cnf.vc_improver_milliseconds = 5 + ceil(cnf.vc_improver_milliseconds
                                                     * GraphUtils::countEdges(V, true) / origV_edges);

        const int ITERS = cnf.solverh_improvement_iterations;
        if(cnf.write_logs){  DEBUG(ITERS); DEBUG(1ll * ITERS * impr.cnf.vc_improver_milliseconds); }

        dfvs = impr.improve(dfvs, ITERS);
        if(cnf.write_logs) DEBUG(dfvs.size());
    }

    assert(Utils::isFVS(V, dfvs));
    if(cnf.write_logs) clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;

    const bool use_sals = (cnf.solverh_use_sals);
    if(use_sals) {
        int pie_edges;
        const int E = GraphUtils::countEdges(V, true);
        {
            VVI revV = GraphUtils::reverseGraph(V);
            VB helper(revV.size(), false);
            pie_edges = Utils::getAllPIEdges(V, revV, helper).size();
        }

        double threshold = 0.3; // original

        if (pie_edges < E * threshold) {
            VI best = dfvs;

            {
                double density = GraphUtils::density(V, true);
                if(cnf.write_logs) DEBUG(density);

                const bool use_sals3_before_sals2 = ( density > cnf.min_density_for_sals3 );

                if(use_sals3_before_sals2){
                    double alpha = 0.97;
                    int maxMvt = 25 * N;
                    int maxFail = 50;

                    if(cnf.write_logs) clog << "Running SALS3 with maxMvt: " << (maxMvt / N) << " * N" << endl;
                    double T0 = SALS3::findInitialTemperature(V, dfvs, cnf,maxMvt / 5);

                    if( T0 > 0 ) {
                        T0 *= alpha * alpha;
                        if(cnf.write_logs) clog << "Initial temperature for SALS3: " << T0 << endl;

                        SALS3 sals3(V, dfvs, cnf);
                        int max_iters = 7;
                        dfvs = sals3.localSearch(T0, alpha, maxMvt, maxFail, max_iters);
                    }else{
                        if(cnf.write_logs) clog << "Initial temperature not found for SALS3" << endl;
                    }

                    if(dfvs.size() < best.size()) best = dfvs;
                }


                double alpha = 0.99;
                int maxMvt = V.size();
                if (V.size() < 1'000) maxMvt *= 300;
                else if (V.size() < 5'000) maxMvt *= 200;
                else if (V.size() < 10'000) maxMvt *= 100;
                else if (V.size() < 30'000) maxMvt *= 50;
                else maxMvt *= 30;

                if(dfvs.size() > 300 && dfvs.size() < 2'000) maxMvt *= 1.5;


                if( E > 500'000 ) maxMvt = 30 * V.size();


                if(cnf.write_logs) clog << "Running SALS2 with maxMvt: " << (maxMvt / N) << " * N" << endl;

                int maxFail = 35;
                bool cool_down_quickly = false;
                double T0 = SALS2::findInitialTemperature(V, dfvs, cnf,0.3 * maxMvt);
                if (T0 < 0) { // initial temperature not found
                    if(cnf.write_logs) clog << "Initial temperature not found for SALS2" << endl;
                    T0 = 0.33;
                    cool_down_quickly = true;
                }else{
                    if(cnf.write_logs) clog << "Initial temperature for SALS2: " << T0 << endl;
                }
                T0 *= pow(1.0 / alpha, 2);


                SALS2 sals2(V, dfvs, cnf,true);
                sals2.cool_down_quickly_until_improved = cool_down_quickly;

                int max_iters = 1e9;
                constexpr int DENSIFIER_VSIZE_THRESHOLD = 30'000;
                constexpr double MAX_DENSITY_FOR_DENSIFIER = 0.02;

                if( V.size() < DENSIFIER_VSIZE_THRESHOLD
                    && density < MAX_DENSITY_FOR_DENSIFIER
                ){
                    // original 50,30,20
                    if (dfvs.size() > 500) max_iters = 40;
                    if (dfvs.size() > 1'000) max_iters = 20;
                    if (dfvs.size() > 2'000 && dfvs.size() < 10'000 ) max_iters = 10;
                }

                if (cnf.write_logs) {
                    clog << "Running SALS2 with T0: " << T0 << ", maxMvt: " << maxMvt / V.size()
                         << " * N and maxFail: " << maxFail << endl;
                }
                VI sals2_dfvs = sals2.localSearch(T0, alpha, maxMvt, maxFail, max_iters);
                if (sals2_dfvs.size() < dfvs.size()) dfvs = sals2_dfvs;

                const bool use_sals_improver = (cnf.solverh_use_sals_improver);
                const bool use_sals_conditionally = (cnf.solverh_use_conditional_sals_improver
                                                     && dfvs.size() >= cnf.solverh_conditional_sals_improver_min_size);

                constexpr bool use_sals3_after_sals2 = true;
                if(use_sals3_after_sals2){
                    if(cnf.write_logs) clog << "Using SALS3 after SALS2" << endl;
                    SALS3::findInitialTemperature(V, dfvs, cnf, maxMvt / 5);
                    if(dfvs.size() < best.size()) best = dfvs;
                }

                if(use_sals_improver || use_sals_conditionally ){
                    ENDL(10);
                    ENDLS(10,'*');

                    ENDL(10);
                    DFVSImprover impr(V, cnf);
                    int iters = cnf.solverh_sals_improver_iterations;
                    if( pie_edges < E * 0.05 ){
                        iters += 2;
                        if(dfvs.size() > 700) iters += 2;
                    }

                    bool use_densifier = ( V.size() < DENSIFIER_VSIZE_THRESHOLD ) && (density < MAX_DENSITY_FOR_DENSIFIER );

                    if(!use_densifier){
                        if(cnf.write_logs) clog << "improveBySALSForSparseGraphs" << endl;
                        dfvs = impr.improveBySALSForSparseGraphs( dfvs, iters, 0.15, T0, alpha, maxMvt/2, maxFail/2 ); // original
                    }

                    if(use_densifier){
                        if(cnf.write_logs) clog << "improveBySALSForSparseGraphsDensifier" << endl;

                        double init_density = 0.007;
                        if(density > init_density) init_density = 0.11;

                        dfvs = impr.improveBySALSForSparseGraphsDensifier(dfvs, 5, init_density, 0.004,
                                             0.25, 0.99, maxMvt / V.size(),
                                             30);
                    }
                }

                if (dfvs.size() < best.size()) best = dfvs;
            }

            dfvs = best;
            assert(Utils::isFVS(V, dfvs));
        }else {
            constexpr bool use_sals4 = false;
            if(use_sals4) {  // use SALS4
                double alpha = 0.95;
                int maxMvt = 15 * N;
                int maxFail = 50;

                if(cnf.write_logs) clog << "Running SALS4 with maxMvt: " << (maxMvt / N) << " * N" << endl;
                double T0 = SALS4::findInitialTemperature(V, dfvs, cnf, maxMvt / 5);

                if (T0 > 0) {
                    T0 *= alpha * alpha;
                    if (cnf.write_logs) clog << "Initial temperature for SALS4: " << T0 << endl;

                    SALS4 sals4(V, dfvs, cnf);
                    int max_iters = 7;
                    dfvs = sals4.localSearch(T0, alpha, maxMvt, maxFail, max_iters);
                } else {
                    if (cnf.write_logs) clog << "Initial temperature not found for SALS4" << endl;
                }
            }else{ // use SALS5
                double alpha = 0.97;
                int maxMvt = 20 * N;
                if ( origV_edges > 1'500'000 ) maxMvt *= 0.75;
                if(origV_N > 100'000) maxMvt *= 0.85;

                if ( origV_edges < 1'000'000) maxMvt *= 2;
                if ( origV_edges < 500'000){
                    if(pie_edges < E * 0.9) maxMvt *= 1.5;
                    else maxMvt *= 2;
                }
                if ( origV_edges < 100'000) maxMvt *= 2;

                if(cnf.write_logs) clog << "Running SALS5 with maxMvt: " << (maxMvt / N) << " * N" << endl;

                int maxFail = 50;

                if (cnf.write_logs) clog << "Running SALS5" << endl;
                double T0 = SALS5::findInitialTemperature(V, dfvs, cnf, maxMvt / 3);

                if (T0 > 0) {
                    T0 *= alpha * alpha;
                    if (cnf.write_logs) clog << "Initial temperature for SALS5: " << T0 << endl;

                    SALS5 sals5(V, dfvs, cnf);
                    int max_iters = 10;
                    dfvs = sals5.localSearch(T0, alpha, maxMvt, maxFail, max_iters);
                    assert(Utils::isFVS(V, dfvs));
                } else {
                    if (cnf.write_logs) clog << "Initial temperature not found for SALS5" << endl;
                }
            }
        }
    }

    return dfvs;
}


VI DFVSImprover::mergeNodesUntilDensityReached(VVI &V, double density_threshold) {
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
        tokens = DFVSSolverH::sinkhorn(V);
        for( auto& d : tokens ) d = 1.0 - d;

        sorted_nodes = VI(N,0);
        iota(ALL(sorted_nodes),0);
        sort(ALL(sorted_nodes), [&](int a, int b){ return tokens[a] < tokens[b]; });

        reverse(ALL(sorted_nodes));

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

        int max_dist = 1;

        if(V.size() > 10'000 ) max_dist *= 15;
        if(V.size() > 20'000 ) max_dist *= 50;

        max_dist = min( 1 + (int)sqrt(sorted_nodes.size()), max_dist );

        auto cond1 = [&]() {
            return next_node_largest >= 0 && sorted_nodes.size() - next_node_largest <= max_dist;
        };

        while( cond1() ){
            while( cond1() && !is_free[sorted_nodes[next_node_largest]] ) next_node_largest--;
            if(!cond1()) break;

            int v = sorted_nodes[next_node_largest];
            to_remove.push_back(v);

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

            for(int d : visited_nodes) helper[d] = false;
        }

        return to_remove;
    };

    int iteration_cnt = 0;
    const bool use_reductions = true;

    int E = 0;

    int iter_cnt = 0;
    while(next_node_largest >= 0 ){

        iter_cnt++;

        Reducer red(V, cnf);
        red.cnf.disableAllNonbasicReductions();
        auto reductions = red.reduce(revV);
        VI dfvs_red = Reducer::convertKernelizedReductions(reductions);

        V = red.V;
        revV = red.revV;

        res += dfvs_red;
        for(int d : dfvs_red) was[d] = true;

        if(cnf.tle()) break;

        {
            InducedGraph g = GraphInducer::induceByNonisolatedNodes(V);
            double dens = GraphUtils::density(g.V,true);
            if( dens >= density_threshold ) break;
        }

        E = GraphUtils::countEdges(V);
        if(E == 0) break;

        initialize_tokens();

        VI selected_nodes = findNodesToRemove();

        for( int v : selected_nodes ){
            if( !Utils::hasLoop(V,v) ){
                was[v] = true;
                Utils::merge(V, revV, v, helper);
            }
        }
    }

    return res;
}
