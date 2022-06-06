//
// Created by sylwester on 12/26/21.
//

#include <graphs/GraphInducer.h>
#include <graphs/GraphUtils.h>
#include <combinatorics/CombinatoricUtils.h>
#include "CONTESTS/PACE22/heur/VCImprover.h"
#include <graphs/toposort/TopoSort.h>
#include <utils/StandardUtils.h>
#include <graphs/VertexCover/approximation/NuMVC/NuMVC.h>
#include <graphs/VertexCover/approximation/LibMVC/numvc.h>
#include <graphs/VertexCover/approximation/LibMVC/fastvc.h>
#include <CONTESTS/PACE22/Utils.h>
#include <CONTESTS/PACE22/heur/DFVSSolverH.h>
#include <graphs/scc/StronglyConnectedComponents.h>
#include <graphs/VertexCover/VCUtils.h>
#include "CollectionOperators.h"

VI VCImprover::improveByVertexCoverMinimizeArcs(VVI V, VI dfvs, int update_frequency) {
    VI order = createOrderMinimizeArcs(V, dfvs, update_frequency);

    VI new_dfvs;
    UniformIntGenerator rnd(0,10);
    int r = rnd.rand();

    if(r == 0) new_dfvs = improveByVertexCoverForSemiTopologicalOrder(V, order, dfvs);
    else new_dfvs = improveByVertexCoverForSemiTopologicalOrder(V, order);

    if(new_dfvs.size() <= dfvs.size() ) return new_dfvs;
    else return dfvs;
}

VI VCImprover::createOrderMinimizeArcs(VVI V, VI dfvs, int update_frequency) {
    const bool debug = false;

    int round = 0;
    int N = V.size();
    VVI revV = GraphUtils::reverseGraph(V);

    { // removing possible repetitions of elements... ?
        VB was(N, false);
        VI new_dfvs;
        new_dfvs.reserve(dfvs.size());
        for (int d : dfvs) {
            if (!was[d]) {
                new_dfvs.push_back(d);
                was[d] = true;
            }
        }
        swap(dfvs, new_dfvs);
    }

    StandardUtils::shuffle(dfvs);

    VI order = CombinatoricUtils::getFullSetDifference(N, dfvs); // this is just V \setminus dfvs
    {
        // we need the topological order now
        if(debug) DEBUG(order);

        InducedGraph g = GraphInducer::induce(V, order);
        TopoSort ts(g.V);
        order = ts.sortTopologically();
        for(int & d : order) d = g.nodes[d];

        if(debug) DEBUG(order);
    }


    int dummy_index = 1e9;
    VI in_order(N,dummy_index);
    for( int i=0; i<order.size(); i++ ) in_order[order[i]] = i;

    int p = 0;

    while( p < dfvs.size() ){

        {  // mark nodes from [dfvs] to positions
            do{
                // finding place to insert for node dfvs[round]
                int v = dfvs[round];

                if(debug) clog << "Finding a place for node " << v << endl;

                VPII neigh_ind; neigh_ind.reserve(V[v].size() + revV[v].size() );
                for( int i=0; i<V[v].size(); i++ ){
                    int w = V[v][i];
                    if(in_order[w] != dummy_index){
                        neigh_ind.emplace_back( in_order[w], 0 );
                    }
                }
                for( int i=0; i<revV[v].size(); i++ ){
                    if(in_order[revV[v][i]] != dummy_index){
                        neigh_ind.emplace_back( in_order[revV[v][i]], 1 );
                    }
                }
                sort(ALL(neigh_ind));

                if(debug){
                    clog << "neigh_ind: " << neigh_ind << ", corresponding to: ";
                    for(auto [d,mark] : neigh_ind) clog << order[d] << " "; clog << endl;
                }

                int best_ind = -1;
                int cnt = 0, best_cnt = 1e9;
                int q = neigh_ind.size();

                for( int i=0; i<q; i++ ) cnt += neigh_ind[i].second;
                for( int i=0; i<q; i++ ){
                    if(debug){
                        clog << "When placing " << v << " just before " << order[neigh_ind[i].first]
                             << ", cnt: " << cnt << endl;
                    }
                    if( cnt < best_cnt ){
                        best_cnt = cnt;
                        best_ind = neigh_ind[i].first;
                    }

                    if( neigh_ind[i].second == 0 ) cnt++;
                    if( neigh_ind[i].second == 1 ) cnt--;
                }

                if( cnt < best_cnt && !neigh_ind.empty() ){
                    best_cnt = cnt;
                    best_ind = neigh_ind.back().first;
                    best_ind = 1e9;
                }

                if(debug){
                    if(debug){
                        clog << "When placing " << v << " at the end, after " << order[neigh_ind.back().first]
                             << ", cnt: " << cnt << endl;
                    }
                    DEBUG(best_cnt);
                    DEBUG(best_ind);
                }

                in_order[v] = best_ind;

                round++;
            } while( round < dfvs.size() && (round % update_frequency) != 0 );

        }

        if(debug){
            clog << "Proceeding to update order" << endl;
        }

        { // update order
            VPII insertions;
            for( int i=p; i<round; i++ ) insertions.emplace_back( dfvs[i], in_order[ dfvs[i] ] );
            sort( ALL( insertions ), [&](auto a, auto b){
                return a.second < b.second;
            } );

            if( debug ) DEBUG(insertions);

            VI new_order; new_order.reserve( order.size() + round - p + 1 );
            int x = 0, y = 0;

            while( x < order.size() || y < insertions.size() ){ // merging [order] with insertions
                if( x < order.size() ){
                    int a = order[x];
                    if( y == insertions.size() || in_order[a] < in_order[insertions[y].first] ){
                        new_order.push_back(a);
                        x++;
                    }
                }

                if( y < insertions.size() ){
                    if( x == order.size() || in_order[order[x]] >= in_order[insertions[y].first]  ){
                        new_order.push_back(insertions[y].first);
                        y++;
                    }
                }
            }

            order = new_order;
            fill(ALL(in_order),dummy_index);
            for( int i=0; i<order.size(); i++ ) in_order[order[i]] = i;

            if(debug) DEBUG(new_order);
        }

        p = round;
    }

    if(!(order.size() == N)){
        set<int> zb;
        for( int d : order ){
            if(zb.count(d)) clog << "Element " << d << " occurs multiple times in order" << endl;
            else zb.insert(d);
        }
        DEBUG(order.size());
        DEBUG(N);
        assert(order.size() == N);
    }


    return order;
}

VI VCImprover::improveByVertexCoverForSemiTopologicalOrder(VVI V, VI order, VI current_dfvs) {
    const bool debug = false;

    int N = V.size();
    assert(order.size() == N);

    VI in_order(N);
    for( int i=0; i<order.size(); i++ ) in_order[order[i]] = i;

    VVI H(N);
    for( int i=0; i<N; i++ ){
        for(int d : V[i]){
            if( in_order[d] < in_order[i] ){
                GraphUtils::addEdge(H,i,d,false); // add undirected edge
            }
        }
    }

    if(debug) DEBUG(H);

    bool use_numvc = ( N <= 20'000 );
    bool use_fastvc = (!use_numvc);
    bool use_improved_numvc = false; // use NuMVC with improvement by fmoessbauer

    if(!current_dfvs.empty()){
        if(cnf.write_logs) clog << "Initializing with VC" << endl;
    }

    if(use_fastvc){
        int milliseconds = cnf.vc_improver_milliseconds;
        int seed;

        {
            UniformIntGenerator rnd(0,1e6);
            seed = rnd.rand();
        }

        VPII edges = GraphUtils::getGraphEdges(H);

        VI init_vc = current_dfvs;

        libmvc::FastVC fastvc(edges, N, 0, std::chrono::milliseconds(milliseconds), false, seed, init_vc);
        fastvc.cover_LS();

        VI mis = fastvc.get_independent_set(false);
        assert(VCUtils::isIndependentSet(H,mis));

        set<int> zb(ALL(mis));

        VI vc;
        for(int i=0; i<N; i++){
            if( H[i].size() > 0 && zb.count(i) == 0 ){
                vc.push_back(i);
            }
        }
        assert(VCUtils::isVertexCover(H, vc));

        if(minimize_found_vc) vc = Utils::findAndRemoveRedundantNodes(V, vc);
        assert(Utils::isFVS(V, vc));

        return vc;
    }
    else if(use_numvc){ // using NuMVC
        if(!use_improved_numvc) {
            NuMVC numvc;
            int milliseconds = cnf.vc_improver_milliseconds;

            VI init_vc = current_dfvs;

            VI mis = numvc.solve(H, 1.0 * milliseconds / 1000, init_vc); // CAUTION! NuMVC indices are from 1 to N
            for (int &d : mis) d--;
            assert(VCUtils::isIndependentSet(H, mis));

            set<int> zb(ALL(mis));

            VI vc;
            for (int i = 0; i < N; i++) {
                if (H[i].size() > 0 && zb.count(i) == 0) {
                    vc.push_back(i);
                }
            }
            assert(VCUtils::isVertexCover(H, vc));

            if(minimize_found_vc) vc = Utils::findAndRemoveRedundantNodes(V, vc);
            assert(Utils::isFVS(V, vc));

            return vc;
        }
        else{ // using NuMVC improved by fmoessbauer
            int seconds = 1;
            int seed = 1712839;
            VPII edges = GraphUtils::getGraphEdges(H);

            VI init_vc = current_dfvs;

            if(cnf.write_logs) clog << "Initializing improved NuMVC" << endl;
            libmvc::NuMVC numvc(edges, N, 0, std::chrono::seconds(seconds), false, seed, init_vc);

            numvc.cover_LS();
            VI mis = numvc.get_independent_set(false);
            if(cnf.write_logs) clog << "VC found" << endl;
            assert(VCUtils::isIndependentSet(H,mis));

            set<int> zb(ALL(mis));

            VI vc;
            for(int i=0; i<N; i++){
                if( H[i].size() > 0 && zb.count(i) == 0 ){
                    vc.push_back(i);
                }
            }

            assert(VCUtils::isVertexCover(H, vc));
            if(minimize_found_vc) vc = Utils::findAndRemoveRedundantNodes(V, vc);
            assert(Utils::isFVS(V, vc));

            return vc;
        }
    }
    else{
        string msg = "PartitionSVC not available here, use NuMVC or FastVC";
        clog << msg << endl;
        cerr << msg << endl;
        exit(1);
    }

}

VI VCImprover::improveByVertexCoverRandom(VVI V, VI dfvs) {
    const bool debug = false;

    int N = V.size();
    VVI revV = GraphUtils::reverseGraph(V);

    VI order = CombinatoricUtils::getFullSetDifference(N, dfvs); // this is just V \setminus dfvs

    // we need the topological order now
    {
        if(debug) DEBUG(order);

        InducedGraph g = GraphInducer::induce(V, order);
        TopoSort ts(g.V);
        order = ts.sortTopologically();
        for(int & d : order) d = g.nodes[d];

        if(debug) DEBUG(order);
    }

    VI insertion_points = CombinatoricUtils::getRandomSequence( order.size(), dfvs.size());
    VPII insertions = StandardUtils::zip( dfvs,insertion_points );

    insertion_points = VI(order.size(),0);
    iota(ALL(insertion_points),0);
    VPII ins_order = StandardUtils::zip(order,insertion_points);

    insertions += ins_order;
    StandardUtils::shuffle(insertions);
    sort(ALL(insertions), []( auto & a, auto & b ){ return a.second < b.second; });

    VI new_order = StandardUtils::unzip(static_cast<vector<pair<int, int>> &&>(insertions)).first;

    if(debug) DEBUG(new_order);
    VI new_dfvs;
    UniformIntGenerator rnd(0,1);
    int r = rnd.rand();


    if(r == 0) new_dfvs = improveByVertexCoverForSemiTopologicalOrder(V, new_order); // do not initialize NuMVC with dfvs-induced vertex cover
    else new_dfvs = improveByVertexCoverForSemiTopologicalOrder(V, new_order, dfvs); // initialize NuMVC with dfvs-induced vertex cover

    if(new_dfvs.size() <= dfvs.size() ) return new_dfvs;
    else return dfvs;
}

VI VCImprover::createOrderByDFAS(VVI V, VI dfvs) {
    const bool debug = true;
    int N = V.size();
    InducedGraph g = GraphInducer::induce(V, dfvs);
    auto [dfasV, edges] = Utils::transformDFAStoDFVS(g.V);

    if(debug){ DEBUG(edges.size());DEBUG(dfasV.size());DEBUG(GraphUtils::countEdges(dfasV,true)); }

    DFVSSolverH solver(cnf);
    solver.cnf.agent_flow_min_distance = 2;

    StronglyConnectedComponents scc(dfasV);
    scc.createStronglyConnectedComponents();
    VVI comps = scc.getComponents();

    if(debug){
        int cmp_cnt = 0;
        for(auto & v : comps) if(v.size() >= 2) cmp_cnt++;
        clog << "Line graph dfasV has " << cmp_cnt << " strongly connected components of size >= 2" << endl;
    }

    VI dfvs_dfasV; // DFVS of graph dfasV
    for(VI & cmp : comps){
        if( cmp.size() < 2 ) continue;
        if(debug) DEBUG(cmp.size());
        InducedGraph gcmp = GraphInducer::induce(dfasV, cmp);

        solver.cnf.agent_flow_max_distance_from_best = 1 + sqrt( cmp.size() );
        VI gcmp_dfvs = solver.solveByAgentFlowAllWithinDistance(gcmp.V);

        for( int &d : gcmp_dfvs ) d = gcmp.nodes[d];
        // now gcmp_dfvs is a dfvs of component cmp in dfasV
        dfvs_dfasV += gcmp_dfvs;
    }

    if(debug) DEBUG(dfvs_dfasV.size());
    assert(Utils::isFVS(dfasV, dfvs_dfasV));

    VPII arcs; // DFAS of graph V
    for( int d : dfvs_dfasV ) arcs.push_back( edges[d] );

    VVI H = g.V;
    VB helper(V.size(),false);
    Utils::removeEdges( H, arcs, helper );

    assert(!Utils::hasCycle(H));

    VI order1; // permutation of dfvs
    VI order2; // topological order of V \ dfvs
    {
        TopoSort ts(H);
        order1 = ts.sortTopologically();

        { // just to check, #CAUTION, #TEST
            for( auto [a,b] : arcs ) H[a].push_back(b);
            assert(Utils::countBackEdgesForNodeOrder(H, order1) == dfvs_dfasV.size());
            if(debug) DEBUG(Utils::countBackEdgesForNodeOrder(H, order1));
            Utils::removeEdges( H, arcs, helper );
        }

        for( int & d : order1 ) d = g.nodes[d]; // after this order1 is a permutation of dfvs
        assert(order1.size() == dfvs.size());
    }

    {
        VI nodes = CombinatoricUtils::getFullSetDifference(N,dfvs);
        InducedGraph gg = GraphInducer::induce( V, nodes );
        TopoSort ts(gg.V);
        order2 = ts.sortTopologically();
        for(int & d : order2) d = gg.nodes[d];
        assert( order1.size() + order2.size() == N );
    }

    return mergeOrders(V, order1, order2);
}

VI VCImprover::mergeOrders(VVI V, VI ord1, VI ord2) {
    const bool debug = false;

    if(debug) Utils::writeRemainingGraph(V);

    VVI revV = GraphUtils::reverseGraph(V);
    int N = V.size();

    int A = ord1.size(), B = ord2.size();

    /**
     * dp[i][j] = min_k dp[i-1][k] + insertion_cost( ord1[i], j)
     */
    VVI dp(A, VI(B+1,1e9));
    VVI min_ind(A, VI(B+1,1e9)); // min_ind[i][j] is the index k for which dp[i][j] is minimized

    VB in_ord2(N,false);
    for(int d : ord2) in_ord2[d] = true;
    VI ind_in_ord2(N,-1);
    for( int i=0; i<ord2.size(); i++ ) ind_in_ord2[ord2[i]] = i;

    if(debug){ DEBUG(ord1);DEBUG(ord2);DEBUG(in_ord2);DEBUG(ind_in_ord2); }


    VI cost_change(B+1,0);
    for( int i=0; i<A; i++ ){
        int u = ord1[i];
        int ins_cost = 0; // cost of inserting node ord1[i] at position j in ord2. Initially equal to revV[u].size()
        for( int d : V[u] ) if( in_ord2[d] ) cost_change[ind_in_ord2[d]]++;
        for( int d : revV[u] ){
            if( in_ord2[d] ){
                cost_change[ind_in_ord2[d]]--;
                ins_cost++;
            }
        }

        if(debug){ DEBUG(u);DEBUG(cost_change); }

        int prefix_dp_min = 1e9, prefix_dp_min_ind = -1;
        if(i == 0) prefix_dp_min = 0; // this will let take dp[0][j] = insertion_cost of ord1[0] at position j

        for( int j=0; j<=B; j++ ){
            if(debug) clog << "insertion cost of " << u << " at position " << j << ": " << ins_cost << endl;

            if(i > 0) {
                if (dp[i - 1][j] <= prefix_dp_min) {
                    prefix_dp_min = dp[i - 1][j];
                    prefix_dp_min_ind = j;
                }
            }

            dp[i][j] = prefix_dp_min + ins_cost;
            min_ind[i][j] = (i == 0 ? j : prefix_dp_min_ind);

            // now update insertion cost
            if( j < B ) ins_cost += cost_change[j];
        }

        if(debug){
            clog << "dp[" << i << "]: " << dp[i] << endl;
            clog << "min_ind[" << i << "]: " << min_ind[i] << endl;
        }

        for( int d : V[u] ) if( in_ord2[d] ) cost_change[ind_in_ord2[d]] = 0; // clearing costs
        for( int d : revV[u] ) if( in_ord2[d] ) cost_change[ind_in_ord2[d]] = 0; // clearing costs
    }

    VI insertion_points(A);
    int ind = -1, val = 1e9;
    for( int i=0; i<dp.back().size(); i++ ){
        if( dp.back()[i] < val ){
            val = dp.back()[i];
            ind = i;
        }
    }

    if(debug) DEBUG(ind);
    for( int i=A-1; i>=0; i-- ){
        insertion_points[i] = ind;
        ind = min_ind[i][ind];
        if(debug)  DEBUG(ind);
    }

    if(debug) DEBUG(insertion_points);

    for( int i=1; i<ord1.size(); i++ ) assert( insertion_points[i] >= insertion_points[i-1] );

    VI merged_order;
    merged_order.reserve(ord1.size() + ord2.size());

    int p = 0;
    for( int i=0; i<=B; i++ ){
        while( p < ord1.size() && insertion_points[p] == i ){
            merged_order.push_back(ord1[p]);
            p++;
        }

        if(i < B) merged_order.push_back( ord2[i] );
    }

    if(debug) DEBUG(merged_order);

    return merged_order;
}

VI VCImprover::createOrderIterativeMerge( VVI & V, VI init_ord, int iterations) {
    VI ord = init_ord;
    int N = V.size();
    assert(ord.size() == N);

    int backgoing_arcs = Utils::countBackEdgesForNodeOrder(V,ord);

    VI best_ord;
    int best_backgoing_edges = 1e9;

    for( int i=0; i<iterations; i++ ){

        VI subs = CombinatoricUtils::getRandomSubset(N-1,N >> 1);
        sort(ALL(subs));
        for(int & d : subs) d = ord[d];

        VI subs2;
        VB in_subs = StandardUtils::toVB(N,subs);
        for( int d : ord ) if( !in_subs[d] ) subs2.push_back(d);

        ord = mergeOrders(V, subs, subs2);

        int temp = Utils::countBackEdgesForNodeOrder(V, ord);

        if(temp > backgoing_arcs){
            assert(temp <= backgoing_arcs);
        }

        if( temp == backgoing_arcs ){
            rotate(ord.begin(), ord.begin() + (N>>1), ord.end());
            backgoing_arcs = Utils::countBackEdgesForNodeOrder(V,ord);
        }else{
            backgoing_arcs = temp;
        }

        if(backgoing_arcs < best_backgoing_edges){
            best_backgoing_edges = temp;
            best_ord = ord;
        }
    }

    return best_ord;
}


