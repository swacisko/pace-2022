//
// Created by sylwester on 12/20/21.
//

#include <graphs/GraphUtils.h>
#include <utils/RandomNumberGenerators.h>
#include <graphs/GraphInducer.h>
#include <utils/StandardUtils.h>
#include <graphs/VertexCover/kernelization/KernelizerVC.h>
#include "CONTESTS/PACE22/Reducer.h"
#include <ranges>
#include "EDReducer.h"
#include <functional>


using namespace Utils;

Reducer::Reducer(VPII edges, Config c) {
    cnf = c;
    primary_edges = edges;
    primaryN = 0;
    for ( auto [a,b] : edges ) primaryN = max(primaryN, max(a,b)+1);
}


ReducedInstance Reducer::reduce() {
    ReducedInstance reduced_instance;
    reduced_instance.cnf = cnf;
    reduced_instance.primaryN = primaryN;

    // VPII reduced_graph_edges;
    // tie( reduced_graph_edges, reduced_instance.primary_liftables) = primaryReduce(primary_edges);


    sw.setLimit(reducer_str, cnf.reducer_max_time_millis);
    sw.start(reducer_str);

    V = GraphUtils::getGraphForEdges(primary_edges);
    N = V.size();

    auto induceAndRemapToBFSOrder = [&]() {
        VI neigh;
        neigh.reserve(N);
        was = VB(N);

        for ( int i=0; i<N; i++ ) if (!V[i].empty() && !was[i]) {
            neigh.push_back(i);
            was[i] = true;

            for ( int j=(int)neigh.size()-1; j<neigh.size(); j++ ) {
                int v = neigh[j];

                for (int d : V[v]) if (!was[d]) {
                    was[d] = true;
                    neigh.push_back(d);
                }
            }
        }

        assert(neigh.size() <= N);

        int N0 = V.size();

        { // inducing graph
            int indN = neigh.size();
            VVI indV(indN);
            for (int i=0; i<neigh.size(); i++) indV[i].reserve( V[neigh[i]].size() );
            VI mapper(N,-1);
            for (int i=0; i<neigh.size(); i++) mapper[neigh[i]] = i;
            for ( int i=0; i<neigh.size(); i++ ) {
                int v = neigh[i];
                for (int d : V[v]) indV[mapper[v]].push_back(mapper[d]);
            }
            swap(V,indV);
        }

        swap(reduced_instance.primary_indg_nodes,neigh);
        N = V.size();
        was = was2 = helper = helper2 = VB(N);

        clog << "Remapped graph with " << N0 << " nodes to a graph with " << N << " nodes using BFS ordering" << endl;
    };

    induceAndRemapToBFSOrder();


    // { // #TEST - here should be implemented a more efficient graph inducing than the following...
    // auto indg = GraphInducer::induceByNonisolatedNodes(V);
    //     reduced_instance.primary_indg_nodes = indg.nodes;
    //     V = indg.V;
    //     N = V.size();
    //     was = was2 = helper = helper2 = VB(N);
    // }

    reduced_instance.secondaryN = N;
    tie( V, reduced_instance.secondary_liftables) = secondaryReduce();
    {
        auto indg = GraphInducer::induceByNonisolatedNodes(V);
        reduced_instance.secondary_indg_nodes = indg.nodes;
        reduced_instance.coreV = indg.V;
    }


    return reduced_instance;
}


pair<VPII, vector<VCReduction*>> Reducer::primaryReduce(VPII & edges) {
    assert(false && "Implement fast primaryReduce");



    return {};
}

pair<VVI, vector<VCReduction *>> Reducer::primaryReduce(VVI &V) {
    auto edges = GraphUtils::getGraphEdges(V);
    auto [res_edges, liftables] = primaryReduce(edges);
    return {GraphUtils::getGraphForEdges(res_edges), liftables};
}

vector<VCReduction*> Reducer::propagateDeg1RuleSlow(int v) {
    VI removed_nodes;

    VI q;
    q.push_back(v);

    while (!q.empty()) {
        v = q.back();
        q.pop_back();
        if (V[v].empty()) continue;

        for ( int d : V[v] ) {
            GraphUtils::removeEdge(V,d,v);
            if ( V[d].size() == 1 ) q.push_back(V[d][0]);
        }

        V[v].clear();
        removed_nodes.push_back(v);
    }

    return {static_cast<VCReduction *>(new KernelizedNodesReduction(removed_nodes))};
}


pair<VVI, vector<VCReduction*>> Reducer::secondaryReduce() {
    constexpr bool debug = false;
    constexpr bool write_progress_on_the_fly = false;

    bool modified;


    vector<VCReduction*> secondary_reduce_liftables;
    KernelizedNodesReduction * knr = nullptr;

    auto addKNR = [&]( VI nodes ){
        if(nodes.empty()) return;
        if( knr == nullptr ) knr = new KernelizedNodesReduction(nodes);
        else knr->addToKer(nodes);
    };

    auto addLiftables = [&]( vector<VCReduction*> liftables ) {
        if(knr != nullptr){ secondary_reduce_liftables.push_back(knr); knr = nullptr; }
        for(auto *x : liftables) secondary_reduce_liftables.push_back(x);
    };

    auto applyDeg1AndDominationWithMasking = [&]() {
        VI removed_nodes;
        VI affected_nodes; affected_nodes.reserve(N);
        VI deg(N,0);
        VB in_V(N);
        VI index(N,-1);
        VI q; q.reserve(sqrt(N));

        bool changes = true;
        for (int i=0; i<N; i++) {
            deg[i] = V[i].size();
            in_V[i] = (deg[i] > 0);
        }
        while (changes) {

            for (int i=0; i<N; i++) if ( deg[i] == 1 ) {
                int v = -1;
                for (int d : V[i]) if ( in_V[d] ) v=d;
                // assert(v != -1);
                // assert(deg[v] > 0);
                // assert(in_V[v]);

                q.clear();
                q.push_back(v);

                while (!q.empty()) {
                    v = q.back();
                    q.pop_back();
                    if (deg[v] == 0) continue;

                    for ( int d : V[v] ) if (in_V[d]) {
                        deg[d]--;
                        if (deg[d] == 0) in_V[d] = false;

                        if ( deg[d] == 1 ) {
                            for (int dd : V[d]) if (in_V[dd]) { // dd is the only element in V[d] that has set in_V
                                q.push_back(dd);
                                break;
                            }
                        }
                    }

                    deg[v] = 0;
                    in_V[v] = false;
                    removed_nodes.push_back(v);
                    affected_nodes.push_back(v);
                }
            }

            auto checkDomination = [&](){
                // clog << "Checking domination!" << endl;
                // checks domination using 'triangle enumeration approach',
                // this is guaranteed to be efficient for dense graphs, but might be slower than brute force
                // for sparse graphs


                int P = affected_nodes.size();
                for (int d : affected_nodes) was[d] = true;
                for (int i=0; i<P; i++) {// was marks nodes that can be checked for domination
                    int v = affected_nodes[i];
                    for (int d : V[v]) if (!was[d]) { was[d] = true; affected_nodes.push_back(d); }
                }
                int P1 = affected_nodes.size();
                for (int i=P; i<P1; i++) { // repeat the same to get N^2(X) where X was nodes removed
                    int v = affected_nodes[i];
                    for (int d : V[v]) if ( !was[d]) { was[d] = true; affected_nodes.push_back(d); }
                }

                // for ( int i=0; i<N; i++ ) if ( in_V[i] && was[i] ) nodes.push_back(i);
                VI nodes = affected_nodes;
                sort(ALL(nodes), [&](int a, int b){ return deg[a] > deg[b]; });

                // DEBUG(affected_nodes.size()); DEBUG(nodes.size());

                VB is_affected = was;
                for (int d : affected_nodes) was[d] = false;

                affected_nodes.clear();

                if (!nodes.empty()) {

                    for (int i=0; i<nodes.size(); i++) index[nodes[i]] = i;
                    for ( int v : nodes ) if (in_V[v]) {
                        was[v] = true;
                        for ( int d : V[v] ) if (in_V[d]) was[d] = true;
                        // for ( int d : V[v] ) if (in_V[d] && index[d] > index[v]) {
                        for ( int d : V[v] ) if (in_V[d] && index[d] > index[v] && is_affected[d]) {
                            bool is_d_dominated_by_v = true;
                            for ( int dd : V[d] ) if ( in_V[dd] && !was[dd] ){ is_d_dominated_by_v = false; break; }
                            // assert(c <= deg[d]);
                            // if (c == deg[d]) { // node v dominates node d
                            if (is_d_dominated_by_v) { // node v dominates node d
                                total_dominations_done++;
                                // clog << "\t Domination holds!" << endl; exit(4);

                                changes = true;
                                removed_nodes.push_back(v);
                                affected_nodes.push_back(v);
                                for ( int x : V[v] ) if (in_V[x]) {
                                    deg[x]--;
                                    if (deg[x] == 0) in_V[x] = false;
                                }
                                deg[v] = 0;
                                in_V[v] = false;
                                break;
                            }
                        }
                        was[v] = false;
                        for ( int d : V[v] ) was[d] = false;
                    }
                }
            };


            changes = false; // we only mark changes in domination, as otherwise graph is exhaustively "deg-1 reduced"
            checkDomination();

            // assertions to check correctness
            // for ( int i=0; i<N; i++ ) if (in_V[i]) {
            //     int c = 0;
            //     for (int d : V[i]) c += in_V[d];
            //     assert(deg[i] == c);
            // }else assert(deg[i] == 0);
            // assert(ranges::none_of(was,std::identity{}));
        }

        addKNR(removed_nodes);

        GraphUtils::removeNodes(V,removed_nodes,helper);

        return !removed_nodes.empty();
    };


    int basic_red_applied = 0;
    function<bool(bool)> applyBasicReductions = [&](bool allow_crown_and_lp){

        if (basic_red_applied++ & 1) { // hybrid approach - in every second call, used masked checks
            if (!allow_crown_and_lp) {
                Stopwatch s; string opt = ( allow_crown_and_lp ? "basic_kern_lp_crown" : "basic_kern"); s.start(opt);
                auto mod = applyDeg1AndDominationWithMasking();
                s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
                return mod;
            }
        }

        // this is the old verion of basic reducer - just to use crown and LP without specific reimplementation

        VVI Vcp = V;
        KernelizerVC kern;

        Stopwatch s; string opt = ( allow_crown_and_lp ? "basic_kern_lp_crown" : "basic_kern"); s.start(opt);
        kern.use_crown_and_lp_checks = allow_crown_and_lp;
        auto [kern_nodes, edges_removed] = kern.initialKernelization(Vcp);
        s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

        addKNR(kern_nodes);
        GraphUtils::removeNodes(V, kern_nodes,helper);

        return !kern_nodes.empty() || !edges_removed.empty();
    };



    int ed_rules_checked_total = 0;
    int ed_node_removal_rules_checked = 0;
    int ed_with_edge_insertion_rules_checked = 0;
    int ed_with_edge_removal_rules_checked = 0;
    int general_folding_rules_checked = 0;

    int ed_edge_insertion_iterations_without_change = 0;
    int ed_last_edge_insertion_node_cnt = 0;


    do{

        modified = false;
        helper = VB(N,false);

        auto suppr = applyBasicReductions(false); // the graph can be modified, but is it modified exhaustively, so no need to rerun it

        if(modified) continue;
        if (sw.tle(reducer_str)) break;


        if(cnf.reducer_use_folding){
            Stopwatch s; string opt = "folding"; s.start(opt);
            auto folds = folding();
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            if(!folds.empty()) modified = true;
            addLiftables(folds);
            assert( GraphUtils::isSimple(V) );

            if(modified) continue;
        }



        if(cnf.reducer_use_funnel){
            if(write_progress_on_the_fly) clog << "Running funnel" << endl;

            Stopwatch s; string opt = "funnel"; s.start(opt);
            auto liftables = funnel();
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            addLiftables(liftables);

            modified |= !liftables.empty();
            if(modified) continue;
        }

        if(cnf.reducer_use_desk){
            Stopwatch s; string opt = "desk"; s.start(opt);
            vector<VCReduction*> desk_liftables = desk();
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            addLiftables(desk_liftables);

            modified |= !desk_liftables.empty();
            // if(modified) continue;
        }

        if(cnf.reducer_use_twins){
            Stopwatch s; string opt = "twins"; s.start(opt);

            vector<VCReduction*> twin_liftables = twins();
            secondary_reduce_liftables += twin_liftables;

            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            modified |= !twin_liftables.empty();
            // if(modified) continue;
        }

        // running unconfined before funnel, or even folding, can help achieve slightly better reduction ratio,
        // but it can make total reduction time up to 1.5x slower on some instances...
        if(cnf.reducer_use_unconfined){
            Stopwatch s; string opt = "unconfined"; s.start(opt);
            if(write_progress_on_the_fly) clog << "Running unconfined" << endl;
            VI uncon = unconfined();
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            addKNR(uncon);
            total_unconfined_nodes += uncon.size();
            GraphUtils::removeNodes(V, uncon, helper);

            modified |= !uncon.empty();
            if(modified) continue;
        }

        // standard node-removal version
        bool ed_application_cond = ( cnf.ed_application_mode == 0 || (cnf.ed_application_mode == 1 && ed_node_removal_rules_checked == 0) );
        if(cnf.reducer_use_ed && cnf.ed_use_node_removal && ed_application_cond){
            ed_rules_checked_total++;
            ed_node_removal_rules_checked++;
            int edge_cnt = GraphUtils::countEdges(V);
            int node_cnt = ranges::count_if(V,[&](auto & v){ return !v.empty(); });
            clog << "Running ED node removal POINT-1, time: " << sw.getTime(reducer_str) / 1000
                 << ", nodes: " << node_cnt << ", edges: " << edge_cnt << endl;

            Stopwatch s; string opt = "ED node removal"; s.start(opt);
            EDReducer edred(V.size(), cnf);
            edred.resetAllUsedTechniques();
            edred.cnf.ed_use_node_removal = true;
            edred.max_time_millis = sw.getLimit(reducer_str) - sw.getTime(reducer_str);

            auto liftables = edred.reduce(V);
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            // assert(res.size() == edred.last_reduce_nodes_removed);
            ed_nodes_reduced += edred.ed_node_applied_cnt;
            // ed_edges_removed += edred.ed_edge_applied_cnt;
            ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;

            addLiftables(liftables);

            if (edred.madeChangesInLastReduce()) V = edred.getV();
            else assert(liftables.empty());
            modified |= edred.madeChangesInLastReduce();

            if(modified) continue;
        }


        bool lp_and_crown_improved = applyBasicReductions(true); // use crown and LP in addition to degree-1 and domination
        modified |= lp_and_crown_improved;
        if (modified) continue;


        // standard node-removal version
        if(cnf.reducer_use_ed && cnf.ed_use_node_removal && !lp_and_crown_improved){
            ed_rules_checked_total++;
            ed_node_removal_rules_checked++;
            int edge_cnt = GraphUtils::countEdges(V);
            int node_cnt = ranges::count_if(V,[&](auto & v){ return !v.empty(); });
            clog << "Running ED node removal POINT-2, time: " << sw.getTime(reducer_str) / 1000
                 << ", nodes: " << node_cnt << ", edges: " << edge_cnt << endl;

            Stopwatch s; string opt = "ED node removal"; s.start(opt);
            EDReducer edred(V.size(), cnf);
            edred.resetAllUsedTechniques();
            edred.cnf.ed_use_node_removal = true;
            edred.max_time_millis = sw.getLimit(reducer_str) - sw.getTime(reducer_str);

            auto liftables = edred.reduce(V);
            // assert(res.size() == edred.last_reduce_nodes_removed);
            ed_nodes_reduced += edred.last_reduce_nodes_removed;
            // ed_edges_removed += edred.last_reduce_edges_removed;
            ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;

            addLiftables(liftables);

            if (edred.madeChangesInLastReduce()) V = edred.getV();
            else assert(liftables.empty());
            modified |= edred.madeChangesInLastReduce();

            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            if(modified) continue;
        }

        bool first_run_condition = (ed_with_edge_removal_rules_checked == 0);
        bool ed_insert_edges_permanently = ( !cnf.ed_remove_added_t1_constraints_if_kernelized_node_found || !cnf.ed_remove_added_t1_constraints_if_no_kernelized_node_found);
        bool interleaving_cond = (cnf.edge_use_edge_removal_and_insertion_interleaving
            && ed_edge_insertion_iterations_without_change <= cnf.ed_max_edge_removal_and_insertion_iterations_without_change );
        if (cnf.reducer_use_ed && cnf.ed_use_edge_removal
            && (first_run_condition || !ed_insert_edges_permanently || interleaving_cond )
            ) {
            ed_with_edge_removal_rules_checked++;
            ed_rules_checked_total++;
            int edge_cnt = GraphUtils::countEdges(V);
            int node_cnt = ranges::count_if(V,[&](auto & v){ return !v.empty(); });
            clog << "Running ED EDGE REMOVAL, time: " << sw.getTime(reducer_str) / 1000
                 << ", nodes: " << node_cnt << ", edges: " << edge_cnt << endl;

            if (ed_edge_insertion_iterations_without_change > 0) {
                clog << "\tRunning iterative edge-removal and edge-insertion for the "
                     << ed_edge_insertion_iterations_without_change << "-th iteration without success"  << endl << endl;
            }

            Stopwatch s; string opt = "ED edge removal"; s.start(opt);
            EDReducer edred(V.size(), cnf);
            edred.resetAllUsedTechniques();
            edred.cnf.ed_use_node_removal = false;
            edred.cnf.ed_use_edge_removal = true;
            edred.max_time_millis = sw.getLimit(reducer_str) - sw.getTime(reducer_str);

            auto liftables = edred.reduce(V);
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);


            // clog << "last_reduce_edges_removed: " << edred.last_reduce_edges_removed << endl;
            // clog << "last_reduce_nodes_removed (deg1-propagated): " << edred.last_reduce_nodes_removed << endl;
            clog << "#CAUTION! Solution lifting not supported yet for ED-edge-removal rule" << endl;

            // assert(res.size() == edred.last_reduce_nodes_removed);
            ed_nodes_reduced += edred.ed_node_applied_cnt;
            ed_edges_removed += edred.ed_edge_applied_cnt;
            ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;

            addLiftables(liftables);
            // if (edred.last_reduce_edges_removed > 0)
                DEBUG(edred.last_reduce_edges_removed);
            if (edred.madeChangesInLastReduce()) V = edred.getV();
            else assert(liftables.empty());
            modified |= edred.madeChangesInLastReduce();

            if(modified) continue;
        }

        if(cnf.reducer_use_general_folding && ++general_folding_rules_checked <= 2){ // at most two times use general folding
            Stopwatch s; string opt = "general folding"; s.start(opt);
            auto reductions = generalFolding();
            total_general_folds_done += reductions.size();
            if(!modified) modified = (!reductions.empty());

            { // add to resulting kernelization objects
                if(knr != nullptr){ secondary_reduce_liftables.push_back(knr); knr = nullptr; }
                for(auto *x : reductions) secondary_reduce_liftables.push_back(x);
            }

            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            if(modified) continue;
        }

        // standard edge-insertion - those edges that are found using consider(v) for single-node initial sets S
        if (cnf.ed_apply_type1_constraints_on_the_fly)
        if (cnf.reducer_use_ed && cnf.ed_use_edge_insertion) {
            vector<VCReduction*> res_liftables;
            bool made_changes = false;
            bool made_changes_in_iteration = false;
            VPII init_V_edges = GraphUtils::getGraphEdges(V);
            addLiftables({}); // we simply want to flush any remaining knr.
            vector<VCReduction*> liftables;

            {
                int node_cnt = ranges::count_if(V,[&](auto & v){ return !v.empty(); });
                if (node_cnt == ed_last_edge_insertion_node_cnt) {
                    ed_edge_insertion_iterations_without_change++;
                }else ed_edge_insertion_iterations_without_change = 0;
                ed_last_edge_insertion_node_cnt = node_cnt;
            }

            do {
                ed_with_edge_insertion_rules_checked++;
                ed_rules_checked_total++;
                int edge_cnt = GraphUtils::countEdges(V);
                int node_cnt = ranges::count_if(V,[&](auto & v){ return !v.empty(); });
                clog << "Running ED with EDGE INSERTION, time: " << sw.getTime(reducer_str) / 1000
                     << ", nodes: " << node_cnt << ", edges: " << edge_cnt << endl;


                made_changes_in_iteration = false;

                Stopwatch s; string opt = "ED edge insertion"; s.start(opt);
                EDReducer edred(V.size(), cnf);
                edred.resetAllUsedTechniques();
                edred.cnf.ed_use_node_removal = true;
                edred.cnf.ed_apply_type1_constraints_on_the_fly = true;
                edred.cnf.edge_use_edge_removal_and_insertion_interleaving = cnf.edge_use_edge_removal_and_insertion_interleaving;
                edred.max_time_millis = sw.getLimit(reducer_str) - sw.getTime(reducer_str);

                liftables = edred.reduce(V);
                if (!liftables.empty()) res_liftables += liftables;
                // assert(res.size() == edred.last_reduce_nodes_removed);
                ed_nodes_reduced += edred.last_reduce_nodes_removed;
                // ed_edges_removed += edred.last_reduce_edges_removed;
                ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;


                made_changes_in_iteration = edred.madeChangesInLastReduce();
                made_changes |= made_changes_in_iteration;
                s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

                if (made_changes_in_iteration) {
                    VI res;
                    for ( auto r : liftables ) {
                        if (auto* derived = dynamic_cast<KernelizedNodesReduction*>(r)) res += derived->getKer();
                        else assert(false && "at this moment here should be only kernelized nodes, no other reductions");
                    }

                    if (!liftables.empty()) assert(!res.empty());

                    if (!res.empty() && cnf.ed_remove_added_t1_constraints_if_kernelized_node_found) {
                        clog << "\tFound " << res.size() << " kernelized nodes when using edge-insertion mode in ED!";
                        clog << "\tReverting graph state and removing nodes" << endl << endl;
                        V = GraphUtils::getGraphForEdges(N,init_V_edges); // revert changes to the original graph
                        GraphUtils::removeNodes(V,res,helper); // and remove nodes...
                        break;

                        // an alternative to breaking...
                        // we simply remove the nodes and continue with edge insertion mode, but remember new
                        // set of init_edges
                        init_V_edges = GraphUtils::getGraphEdges(V);
                    }else {
                        V = edred.getV();
                    }
                }
                assert( GraphUtils::isSimple(V) );

            }while (liftables.empty() && made_changes_in_iteration);

            if ( cnf.ed_remove_added_t1_constraints_if_no_kernelized_node_found && res_liftables.empty() ) {
                clog << "\tDid not find any kernelized node using ED-edge-insertion, reverting to original state" << endl;
                V = GraphUtils::getGraphForEdges(N,init_V_edges);
                made_changes = false;
            }

            if (made_changes) secondary_reduce_liftables += res_liftables;

            modified |= made_changes;

            if(modified) continue;
        }

    } while(modified);

    if(debug){ DEBUG(V);}

    if(knr != nullptr){ secondary_reduce_liftables.push_back(knr); knr = nullptr;}


    return make_pair(V, secondary_reduce_liftables);
}

vector<VCReduction*> Reducer::twins() {
    vector<VCReduction*> liftables;

    VLL hashes(N);
    IntGenerator rnd;
    for (auto & d : hashes) d = rnd.rand();
    vector<pair<LL,int>> neigh_h;
    neigh_h.reserve(N);

    VB affected(N);

    for ( int i=0; i<N; i++ ) if (!V[i].empty()) {
        LL h = 0;
        for (int d : V[i]) h ^= hashes[d];
        neigh_h.emplace_back(h,i);
    }
    sort(ALL(neigh_h));

    int p = 0, q = p;
    while ( p < neigh_h.size() ) {
        q = p+1;
        while ( q < neigh_h.size() && neigh_h[q].first == neigh_h[p].first ) q++;
        if (q == p+1){ p = q; continue; }

        int deg = V[neigh_h[p].second].size();
        // We have twins in set X and their neighborhood in set T
        // if |T| <= |X|, then we can add T to the solution
        // if |T| = |X|+1 and T is an independent set, we can fold twins
        // additionally, the size of vc in G[T] can be taken into account as well

        int vc_in_T_size = 0;
        if ( q-p >= 2 && V[neigh_h[p].second].size() <= q-p+4 ) {
            VI T = V[neigh_h[p].second];
            auto indg = GraphInducer::induce(V,T);
            vc_in_T_size = Utils::getMinVcCPSAT(indg.V).size();
        }

        if ( q-p >= deg - vc_in_T_size ) { // we have a set X of twins with |N(X)| >= |X|. We denote T = N(X)
            VI T = V[neigh_h[p].second];
            bool is_affected = false;
            for (int t : T) is_affected |= affected[t];
            for ( int i=p; i<q; i++ ) is_affected |= affected[ neigh_h[i].second ]; // check nodes in X for affected
            if (!is_affected) {
                total_twins_done++;
                liftables.push_back(new KernelizedNodesReduction(T));
                for (int t : T) affected[t] = true;
                for ( int i=p; i<q; i++ ) affected[ neigh_h[i].second ] = true;
                GraphUtils::removeNodes(V,T,was);
            }
        }

        if ( q-p+1 == deg - vc_in_T_size ) { // if T is an independent set, we can fold those twins
            VI T = V[neigh_h[p].second];

            bool is_affected = false;
            for (int t : T) is_affected |= affected[t];
            for ( int i=p; i<q; i++ ) is_affected |= affected[ neigh_h[i].second ]; // check nodes in X for affected

            if (!is_affected) {
                for (int t : T) was[t] = true;
                bool is_mis = true;
                for ( int t : T ) for (int d : V[t]) if ( was[d] ){ is_mis = false; break; };
                for (int t : T) was[t] = false;

                if (is_mis) { // we can fold
                    total_twins_done++;
                    VI neigh;
                    for ( int i=p; i<q; i++ ) was[neigh_h[i].second] = true;
                    for (int t : T) for (int d : V[t]) if (!was[d]) {
                        was[d] = true;
                        neigh.push_back(d);
                    }
                    for ( int i=p; i<q; i++ ) was[neigh_h[i].second] = false;
                    for (int d : neigh) was[d] = false;

                    GraphUtils::removeNodes(V,T,helper);

                    VI X; X.reserve(q-p);
                    for (int i=p; i<q; i++) X.push_back(neigh_h[i].second);
                    int merge_to = T.back();

                    for ( int d : neigh ) GraphUtils::addEdge(V,merge_to,d);

                    T.pop_back();
                    // clog << "\tApplying twin folding!! Check the correctness of it, test it properly!" << endl;
                    liftables.push_back(new FoldingTwinReduction(merge_to,T, X));

                    for (int x : X) assert(V[x].empty());
                    for (int x : T) assert(V[x].empty());
                }
            }
        }

        p = q;
    }

    VI rem;
    for (int i=0; i<N; i++) if (V[i].size() == 1) liftables += propagateDeg1RuleSlow(V[i][0]);
    if (!rem.empty()) liftables.push_back(new KernelizedNodesReduction(rem));

    return liftables;
}


vector<VCReduction*> Reducer::folding() {
    vector<VCReduction*> liftables;

    bool changes = true;
    while (changes) {
        changes = false;

        for (int i=0; i<N; i++) if (V[i].size() == 2) {
            int a = V[i][0], b = V[i][1];
            changes = true;

            if( ranges::contains(V[a],b) ) {
                liftables.push_back(new KernelizedNodesReduction({a,b}));
                for (int d : V[a]) GraphUtils::removeEdge(V,d,a);
                for (int d : V[b]) GraphUtils::removeEdge(V,d,b);

                // for (int d : V[a]) if (V[d].size() == 1) liftables += propagateDeg1RuleSlow(d);
                // for (int d : V[b]) if (V[d].size() == 1) liftables += propagateDeg1RuleSlow(d);
                for (int d : V[a]) if (V[d].size() == 1) liftables += propagateDeg1RuleSlow(V[d][0]);
                for (int d : V[b]) if (V[d].size() == 1) liftables += propagateDeg1RuleSlow(V[d][0]);

                V[a].clear();
                V[b].clear();

                continue;
            }
            liftables.push_back(new FoldingReduction(a,b,i));
            GraphUtils::removeNodeFromGraph(V,i);
            total_folds_done++;

            for (int d : V[a] ) was[d] = true;
            for ( int d : V[b] ) if (!was[d]) {
                V[a].push_back(d);
                V[d].push_back(a);
            }
            for (int d : V[a] ) was[d] = false;

            GraphUtils::removeNodeFromGraph(V,b);
        }
    }

    return liftables;
}










void Reducer::writeTotals() {
    DEBUG(total_dominations_done);
    DEBUG(total_folds_done);
    DEBUG(total_funnels_done);
    DEBUG(total_desks_done);
    DEBUG(total_twins_done);
    DEBUG(total_unconfined_nodes);
    DEBUG(total_general_folds_done);

    DEBUG(ed_nodes_reduced);
    DEBUG(ed_edges_removed);
    DEBUG(ed_t1_inference_rules_added);
    DEBUG(ed_total_t2_inference_rules_created);
    DEBUG(ed_t2_inference_rules_added);

    clog << "Reduction times (sec.): " << endl;
    for ( auto [k,v] : reduction_times_millis ) clog << k << " -> " << v / 1000.0 << endl;
    ENDL(1);
}





void Reducer::disableAllNonbasicReductions() {
    cnf.disableAllNonbasicReductions();
}

void Reducer::disableAllConditionalReductions() {
    cnf.disableAllConditionalReductions();
}



vector<VCReduction *> Reducer::applyAlternativeSets(VI A, VI B, bool log) {
    vector<VCReduction*> liftables;


    { // remove nodes from N(A) \cap N(B) and add
        VI to_remove;
        for (int a : A) was[a] = true;
        for (int a : A) for (int d : V[a]) if (!was[d]) helper[d] = true; // N(A) without nodes in A

        for (int b : B) was2[b] = true;
        for (int b : B) for (int d : V[b]) if (!was2[d] && helper[d] && !helper2[d]) {
            to_remove.push_back(d);
            helper2[d] = true;
        }
        for (int d : to_remove) helper2[d] = false;
        for (int b : B) was2[b] = false;

        for (int a : A) for (int d : V[a]) helper[d] = false;
        for (int a : A) was[a] = false;

        if (!to_remove.empty()) {
            if (log) clog << "In alternative sets, found nonempty intersection of NA and NB: " << to_remove << endl;
            liftables.push_back( new KernelizedNodesReduction(to_remove) );
            GraphUtils::removeNodes(V,to_remove,helper);
        }
    }


    VI NA, NB;
    { // add lacking connections
        for (int a : A) was[a] = true;
        for (int a : A) for (int d : V[a]) if (!was[d] && !helper[d]) {
            helper[d] = true;
            NA.push_back(d);
        }
        for (int d : NA) helper[d] = false;
        for (int a : A) was[a] = false;

        for (int b : B) was[b] = true;
        for (int b : B) for (int d : V[b]) if (!was[d] && !helper[d]) {
            helper[d] = true;
            NB.push_back(d);
        }
        for (int d : NB) helper[d] = false;
        for (int b : B) was[b] = false;

        for ( int b : NB ) {
            for (int d : V[b]) was[d] = true;
            for (int d : A) was[d] = true;

            for (int a : NA) if (!was[a]) {
                // assert(!ranges::contains( V[a],b ));
                // assert(!ranges::contains( V[b],a ));
                GraphUtils::addEdge(V,a,b);
            }

            for (int d : V[b]) was[d] = false;
            for (int d : A) was[d] = false;
        }
    }

    if (log){ DEBUG(A); DEBUG(B); DEBUG(NA); DEBUG(NB); }
    NA = StandardUtils::setDifference(NA,B,helper);
    if ( !NA.empty() && !NB.empty() ) liftables.push_back( new AlternativeSetsReduction(NA,B,A) );

    GraphUtils::removeNodes(V,A,helper);
    GraphUtils::removeNodes(V,B,helper);

    return liftables;
}

vector<VCReduction*> Reducer::funnel() {
    vector<VCReduction*> liftables;

    constexpr bool run_correctness_assertions = false;

    if constexpr(run_correctness_assertions) assert(ranges::none_of(was, std::identity{})); // #TEST #CAUTION - just an assertion for tests

    bool changes = true;
    constexpr bool run_exhaustively = true;

    // finds a candidate node x for a funnel edge {v,x}
    auto findIsolatedNodeCandidate = [&](int v) -> pair<int,bool> {
        int cand = -1;

        if constexpr(run_correctness_assertions) assert(ranges::none_of(was, std::identity{})); // #TEST #CAUTION - just an assertion for tests
        if constexpr(run_correctness_assertions) assert(ranges::none_of(was2, std::identity{})); // #TEST #CAUTION - just an assertion for tests

        for ( int d : V[v] ) was[d] = true;
        int a = V[v][0]; // any node from V[v]
        int cnt = 0;
        for (int d : V[a]) cnt += was[d];
        if (cnt >= V[v].size()) { DEBUG(cnt); DEBUG(PII(v,a)); DEBUG(V[v]); DEBUG(V[a]); }
        assert(cnt <= (int)V[v].size()-1);
        for ( int d : V[v] ) was[d] = false;

        if (cnt == 0) cand = a; // a is isolated from N(v)
        else if (cnt < (int)V[v].size()-1 ) { // there exists node in N(v) that is not in N[a], we need to find it
            for (int d : V[a]) was2[d] = true;
            for (int x : V[v]) if (!was2[x] && x != a) cand = x;
            for (int d : V[a]) was2[d] = false;
            assert(cand != -1);
        }else {
            // node a dominates node v
            // clog << "Found a dominating node using funnel, a: " << a << ", V[a]: " << V[a] << endl;
            // clog << "v: " << v << ", V[v]: " << V[v] << endl;
            return {a,true};
        }

        return {cand,false};
    };


    // checks whether V[v] \setminus {cand} is a clique
    auto isClq = [&](int v, int cand)-> bool {
        bool is_clq = true;
        int clq_size = (int)V[v].size()-1;

        for ( int d : V[v] ) if (d != cand) was[d] = true;
        for ( int d : V[v] ) if ( d != cand ) {
            int c = 0;
            for (int dd : V[d]) c += was[dd];
            assert(c+1 <= clq_size);
            if (c+1 < clq_size) { is_clq = false; break; }
        }
        for ( int d : V[v] ) was[d] = false;

        return is_clq;
    };

    constexpr bool allow_separate_domination = true;

    while (changes) {
        changes = false;

        // VI nodes(N); iota(ALL(nodes),0);
        // sort(ALL(nodes), [&](int a, int b){ return V[a].size() < V[b].size(); });
        // for (int v : nodes) if (V[v].size() >= 2) if ( (int)V[v].size()-1 <= cnf.reducer_max_funnel_clique_size ) {
        for (int v=0; v<N; v++) if (V[v].size() >= 2) if ( (int)V[v].size()-1 <= cnf.reducer_max_funnel_clique_size ) {
            auto [cand,dominates] = findIsolatedNodeCandidate(v);

            if constexpr(allow_separate_domination) if (dominates) {
                liftables += propagateDeg1RuleSlow(cand);
                if (run_exhaustively) changes = true;
                continue;
            }

            bool is_clq = isClq(v,cand);
            if (!is_clq) continue;

            if constexpr(allow_separate_domination) if (!dominates) {
                auto neigh = V[v];
                auto cand_neigh = V[cand];
                for (int d : V[v]) if (d != cand) was[d] = true;
                for ( int d : cand_neigh ) if (was[d] && !V[d].empty() && !V[v].empty()) {
                    dominates = true;
                    // clog << "Found a dominating node in funnel!" << endl;
                    liftables += propagateDeg1RuleSlow(d);
                    if (run_exhaustively) changes = true;
                }
                for (int d : neigh) was[d] = false;
                if (dominates) continue;
            }





            // we found a funnel {v,cand}
            // clog << endl << "Found a funnel!" << endl;
            // DEBUG(PII(v,cand));
            // clog << "V[" << v << "]: " << V[v] << endl;
            // for (int d : V[v]) clog << "V[" << d << "]: " << V[d] << endl;
            // clog << "Considered clique: "; for (int d : V[v]) if (d != cand) clog << d << " "; clog << endl;


            auto lft = applyAlternativeSets({v}, {cand});

            if constexpr(run_correctness_assertions) assert(ranges::none_of(was, std::identity{})); // #TEST #CAUTION - just an assertion for tests
            if constexpr(run_correctness_assertions) assert(GraphUtils::isSimple(V)); // #TEST just for debugging

            // DEBUG(lft.size());
            // assert(lft.size() == 1);
            total_funnels_done += lft.size();
            if (run_exhaustively) changes = true;

            for (auto l : lft) {
                if (auto* derived = dynamic_cast<AlternativeSetsReduction*>(l)) derived->red_name = "funnel";
                else {
                    // else this is a kernelized nodes reduction - in case of a funnel it should not exist if
                    // domination was applied before
                    if constexpr(allow_separate_domination)
                    assert(false && "this should not happen, as we distinguish domination case separately,"
                                    " so no kernelized nodes should be created for found alternative sets");
                }
            }
            liftables += lft;
        }
    }

    return liftables;
}







void Reducer::liftSolution(int N, VI &dfvs, vector<VCReduction *> &reductions, bool clear_reductions) {
    VB in_dfvs = StandardUtils::toVB(N, dfvs);

    for( int i = (int)reductions.size()-1; i>=0; i-- ){
        reductions[i]->lift(dfvs, in_dfvs);
    }

    if (clear_reductions) clearReductionObjects(reductions);
}

void Reducer::clearReductionObjects(vector<VCReduction *> &reductions) {
    for(int i=0; i<reductions.size(); i++ ){
        delete reductions[i];
        reductions[i] = nullptr;
    }
}

VI Reducer::convertKernelizedReductions(vector<VCReduction *> &reductions) {
    assert(reductions.size() <= 1);
    VI red_dfvs;
    if(!reductions.empty()){
        KernelizedNodesReduction * knr = (KernelizedNodesReduction*) reductions[0];
        red_dfvs = knr->getKer();
        Reducer::clearReductionObjects(reductions);
    }
    return red_dfvs;
}

int Reducer::getReductionsOffset(vector<VCReduction *> &reductions) {
    int res = 0;
    for(auto * x : reductions) res += x->offset();
    return res;
}

void Reducer::writeReductions(vector<VCReduction *> &reductions) {
    clog << "Reductions: " << endl;
    for(auto * x : reductions) clog << x->toString() << endl;
}

void Reducer::createVFromEdges(VPII & edges) {
    int N = 0;
    for (auto [a,b] : edges) N = max(N, max(a,b)+1);
    V.resize(N);
    for (auto [a,b] : edges) GraphUtils::addEdge(V,a,b);
}


vector<VCReduction*> Reducer::desk(){
    vector<VCReduction*> liftables;

    auto edges = GraphUtils::getGraphEdges(V);

    VI desk, neigh;

    for (auto [a,b] : edges) if (V[a].size() == 3 && V[b].size() == 3) {
        desk.clear();
        bool adj = false;
        for (int x : V[a]) adj |= (x == b);
        if (!adj) continue;

        for (int d : V[a]) was[d] = true;
        for ( int c : V[b] ) if ( c != a && V[c].size() == 3 && !was[c] ) {
            for ( int d : V[c] ) if ( d != a && d != b && was[d] && V[d].size() == 3 ) {
                desk = {a,b,c,d};
                break;
            }
            if (!desk.empty()) break;
        }
        for (int d : V[a]) was[d] = false;

        if (desk.empty()) continue;
        int c = desk[2], d = desk[3];
        for (int x : V[a]) if (x == c) desk.clear();
        for (int x : V[b]) if (x == d) desk.clear();
        if (desk.empty()) continue;

        neigh.clear();
        for (int x : desk) was[x] = true;
        for (int x : desk) for (int y : V[x]) if (!was[y]){ was[y] = true; neigh.push_back(y); }
        for (int x : desk) was[x] = false;
        for (int x : neigh) was[x] = false;


        constexpr bool debug = false;
        assert(neigh.size() <= 4);

        if ( neigh.size() < 4 ) {

            for (int i=(int)neigh.size()-1; i>=0; i--) {
                int x = neigh[i];
                bool remove_x = (ranges::contains(V[x],a) && ranges::contains(V[x],b));
                remove_x |= (ranges::contains(V[x],b) && ranges::contains(V[x],c));
                remove_x |= (ranges::contains(V[x],c) && ranges::contains(V[x],d));
                remove_x |= (ranges::contains(V[x],d) && ranges::contains(V[x],a));
                if (!remove_x) REM(neigh,i);
            }

            if (neigh.empty()) continue;

            if (debug) clog << "Found DESK DOMINATION, need to handle it properly..." << endl;
            total_desks_done++;

            GraphUtils::removeNodes(V,neigh,was);
            liftables.push_back(new KernelizedNodesReduction(neigh));
            continue;
        }

        if (debug) { DEBUG(desk); DEBUG(neigh); for (int x : desk) clog << "V[" << x << "]: " << V[x] << endl; }

        auto lft = applyAlternativeSets({a,c}, {b,d}, debug);
        for (auto l : lft) if (auto* derived = dynamic_cast<AlternativeSetsReduction*>(l)) derived->red_name = "desk";
        liftables += lft;

        if (debug) clog << "\t\tFound a DESK FOLDING!" << endl;
        total_desks_done++;
    }

    VI rem;
    for (int i=0; i<N; i++) if (V[i].size() == 1) liftables += propagateDeg1RuleSlow(V[i][0]);
    if (!rem.empty()) liftables.push_back(new KernelizedNodesReduction(rem));

    return liftables;
}



vector<GeneralFoldingReduction *> Reducer::generalFolding() {
    constexpr bool debug = false;

    int N = V.size();
    VB affected(N,false);
    VB helper(N,false);
    VB was(N,false);

    vector<GeneralFoldingReduction *> res;

    VI order(N); iota(ALL(order),0);
    sort(ALL(order), [&](int a, int b){ return V[a].size() < V[b].size(); });
    reverse(ALL(order)); // starting node selection from largest degree in general_folding

    for( int w : order ){
        if(V[w].size() <= 1) continue;
        if(affected[w]) continue;
        if( V[w].size() > cnf.reducer_max_general_folding_neighborhood_size ) continue;

        VI W = V[w];
        bool aff = false;
        for(int d : W){
            if(affected[d]) aff = true;
            for( int u : V[d] ) if(affected[u]) aff = true;
        }
        if(aff) continue;

        int vc_size;
        InducedGraph g = GraphInducer::induce(V, W);
        vc_size = Utils::getMinVcCPSAT(g.V).size();

        if( vc_size + 2 < W.size() ) continue;

        if(debug){ ENDL(5); clog << "Found W with DFVS(G[W]) >= W.size()-2" << endl; DEBUG(w); DEBUG(W); }


        VPII antiedges;
        {
            VVI comppigv = GraphUtils::getComplimentaryGraph(g.V);
            antiedges = GraphUtils::getGraphEdges(comppigv, false); // we want undirected here
            for( auto & [a,b] : antiedges ){ // remapping antiedges to original ids
                a = g.nodes[a];
                b = g.nodes[b];
            }
        }

        if(debug) DEBUG(antiedges);

        if( antiedges.size() > min((int)W.size(), cnf.reducer_max_general_folding_antiedges) ) continue;

        affected[w] = true;
        for(int d : W){
            affected[d] = true;
            for( int u : V[d] ) affected[u] = true;
        }

        VPII arcs_to_remove;
        {
            for( int d : W ) for( int u : V[d] ) arcs_to_remove.emplace_back(d,u);
            StandardUtils::makeUnique(arcs_to_remove);
        }

        VPII arcs_to_add;
        vector<tuple<int,int,int>> antiedges_tuples;

        { // finding arcs to add
            VI free_ids = (VI(1,w) + W);

            for( int i=0; i<antiedges.size(); i++){
                int id = free_ids[i];
                int a = antiedges[i].first;
                int b = antiedges[i].second;

                VI Npi = StandardUtils::setUnion(V[a], V[b], helper);
                StandardUtils::removeFromArrayPreserveOrderInplace( Npi, free_ids, helper );
                arcs_to_add += StandardUtils::product( VI({id}), Npi );
                arcs_to_add += StandardUtils::product( Npi, VI({id}) );

                antiedges_tuples.emplace_back( id,a,b );
            }

            for( int i=0; i<antiedges.size(); i++ ){
                for( int j=i+1; j<antiedges.size(); j++ ){
                    int id1 = free_ids[i];
                    int id2 = free_ids[j];
                    arcs_to_add.emplace_back(id1,id2);
                    arcs_to_add.emplace_back(id2,id1);
                }
            }
        }

        // GraphUtils::removeEdges(V, arcs_to_remove, helper); // implement this efficiently...
        GraphUtils::removeEdges(V, arcs_to_remove);

        // Utils::addEdges(V, revV, arcs_to_add, helper);
        { // make arcs_to_add unique edges - remove duplicates of same arcs with opposite direction
            int P = arcs_to_add.size();
            for (int i=0; i<P; i++) {
                auto [a,b] = arcs_to_add[i];
                arcs_to_add.emplace_back(b,a);
            }
            StandardUtils::makeUnique(arcs_to_add);
            for (int i=(int)arcs_to_add.size()-1; i>=0; i--) {
                auto [a,b] = arcs_to_add[i];
                if ( a > b ) { // there will also be edge [b,a] in arcs_to_add, we do not wany copies...
                    swap(arcs_to_add[i], arcs_to_add.back());
                    arcs_to_add.pop_back();
                }
            }
        }
        for ( auto [a,b] : arcs_to_add ) GraphUtils::addEdge(V,a,b);
        res.push_back( new GeneralFoldingReduction( w, W, antiedges_tuples ) );

    }
    return res;
}


VI Reducer::unconfined() {
    constexpr bool debug = false;
    VI removed_nodes;

    VB inS(N), inNS(N);
    VI S;
    VI deg_in_S(N,0);
    VI deg_out_NS(N,0);
    int ns_size = 0;
    VB calculated_NS_outdeg(N);

    VI cand_NS;

    constexpr bool check_basic_assertions = true;
    constexpr bool check_expensive_assertions = false;

    auto checkAssertions = [&]() {
        if constexpr (!check_expensive_assertions) return;

        if constexpr (debug) clog << "\t#CAUTION! Checking slow assertions in unconfined" << endl;

        for (int i=0; i<N; i++) if (!inNS[i]) {
            assert(deg_in_S[i] == 0);
            assert(deg_out_NS[i] == 0);
            assert(calculated_NS_outdeg[i] == false);
        }

        for (int s : S) {
            int d_in_s = 0;
            for (int d : V[s]) d_in_s += inS[d];
            assert(d_in_s == deg_in_S[s]);

            int d_out_ns = 0;
            for (int d : V[s]) d_out_ns += !inNS[d];
            if (d_out_ns != 0) DEBUG(PII(s,d_out_ns));
            assert(d_out_ns == 0);

            for (int d : V[s]) if (!inS[d]) {
                assert(inNS[d]);

                d_in_s = 0;
                for (int dd : V[d]) d_in_s += inS[dd];
                if (d_in_s != deg_in_S[d]) DEBUG(PII(d_in_s, deg_in_S[d]));
                assert(d_in_s == deg_in_S[d] );

                if (calculated_NS_outdeg[d]) {
                    d_out_ns = 0;
                    for (int dd : V[d]) d_out_ns += !inNS[dd];
                    if (d_out_ns != deg_out_NS[d]) {
                        DEBUG(d_out_ns);
                        DEBUG(deg_out_NS[d]);
                    }
                    assert(d_out_ns == deg_out_NS[d]);
                }
            }
        }
    };

    auto clearForS = [&]() {
        int ns_size = 0, ns_sumdeg = 0;
        for ( int s : S ) for (int d : V[s]) if (!was[d]) {
            was[d] = true;
            ns_size++;
            ns_sumdeg += V[d].size();
        }
        for ( int s : S ) for (int d : V[s]) was[d] = false;
        // clog << "S.size(): " << S.size() << ", ns_size: " << ns_size << ", ns_sumdeg: " << ns_sumdeg << endl;


        for (int a : S) {
            calculated_NS_outdeg[a] = inS[a] = inNS[a] = deg_in_S[a] = deg_out_NS[a] = 0;
            for (int d : V[a]) calculated_NS_outdeg[d] = inS[d] = inNS[d] = deg_in_S[d] = deg_out_NS[d] = 0;
        }
        S.clear();
        cand_NS.clear();
    };

    auto getLbNSOutdegForUnexpandedNode = [&](int d) {
        // // with 'return 0' uncommented the unconfined rule is stronger, but could potentially be much slower.
        // It should be the same, but there probably must be some hidden bug somewhere...
        return 0;

        if (calculated_NS_outdeg[d]) return deg_out_NS[d];
        return max(0, (int)V[d].size() - deg_in_S[d] - ns_size);
    };

    auto calculateNSOutdeg = [&](int v) {
        calculated_NS_outdeg[v] = true;
        deg_out_NS[v] = 0;
        for (int d : V[v]) deg_out_NS[v] += !inNS[d];
    };

    auto check = [&](int v) {
        if constexpr (debug) clog << "checking node v: " << v << endl;

        S.clear(); S.push_back(v);
        inS[v] = inNS[v] = true;
        ns_size = V[v].size();
        for (int d : V[v]) {
            inNS[d] = true;
            deg_in_S[d] = 1;
        }
        cand_NS.clear();
        for (int d : V[v]) if (getLbNSOutdegForUnexpandedNode(d) <= 1) {
            calculated_NS_outdeg[d] = true;
            for ( int dd : V[d] ) deg_out_NS[d] += !inNS[dd];
            if (deg_out_NS[d] == 1) cand_NS.push_back(d);
        }

        while ( !cand_NS.empty() ) {
            if constexpr (debug) clog << "\tcand_NS: " << cand_NS << endl;
            checkAssertions();

            int u = cand_NS.back(); // we will move node u to S - but here it is the node that has |N(u) \ N(S)| = 1
            cand_NS.pop_back();
            if constexpr (check_basic_assertions) assert(!inS[u]);

            if ( deg_in_S[u] != 1 ) continue;
            if ( deg_out_NS[u] == 0 ) return true;
            if constexpr (check_basic_assertions) assert(deg_out_NS[u] == 1);

            int cnt = 0;
            for (int d : V[u]) cnt += !inNS[d];
            if constexpr (check_basic_assertions) assert(cnt == 1);

            int only_neighbor = -1;
            for ( int d : V[u] ) if ( !inNS[d] ){ only_neighbor = d; break; }
            if constexpr (check_basic_assertions) assert(only_neighbor != -1);

            if constexpr (debug) clog << "\tcand u: " << u << ", only_neighbor: " << only_neighbor << endl;

            // now u should be the only neighbor, so we can move u to S
            u = only_neighbor;
            if constexpr (debug) clog << "\tmoving node " << u << " to S, V[" << u << "]: " << V[u] << endl;

            if constexpr (check_basic_assertions) assert(deg_in_S[u] == 0);
            inS[u] = inNS[u] = true;
            S.push_back(u);
            deg_out_NS[u] = 0;


            bool can_return_true = false;

            for ( int d : V[u] ) if (inNS[d]) {
                // u may have several neighbors in S, but is was the only neighbor outside N[S] for some node in N(S)
                deg_in_S[d]++;
                deg_out_NS[d]--;
            }

            for (int d : V[u]) if (!inNS[d]) { // moving conceptually node d to N(S)
                if constexpr (debug) clog << "\t\tmoving node d: " << d << " to N(S), V[" << d << "]: " << V[d] << endl;
                inNS[d] = true;
                ns_size++;
                if constexpr (check_basic_assertions) assert(deg_out_NS[d] == 0);

                for ( int dd : V[d] ) {
                    if (dd == u) { // this is the only neighbor of d that can be in S
                        // this was not taken into account when iterating over V[u] earlier, as d was not set in inNS
                        deg_out_NS[d]--;
                        deg_in_S[d]++;
                        continue;
                    }

                    if constexpr (check_basic_assertions) assert(!inS[dd]);
                    if ( inNS[dd] ) { // we need to update degrees of all neighbors dd of node d which is moved to N(S)
                        if ( getLbNSOutdegForUnexpandedNode(dd) <= 1 ) {
                            if ( !calculated_NS_outdeg[dd] ) calculateNSOutdeg(dd);
                            else deg_out_NS[dd]--;

                            // if (deg_out_NS[dd] == 0) return true;
                            if (deg_out_NS[dd] == 0 && deg_in_S[dd] <= 1) {
                                if constexpr (check_basic_assertions) assert(deg_in_S[dd] == 1);
                                can_return_true = true;
                            }
                            if (deg_out_NS[dd] == 1) cand_NS.push_back(dd);
                        }
                    }
                    else { // dd is not in NS
                        deg_out_NS[d]++;
                    }
                }

                // if (!calculated_NS_outdeg[d])
                    calculateNSOutdeg(d);
                // checkAssertions();

                if (calculated_NS_outdeg[d] && deg_out_NS[d] == 1) cand_NS.push_back(d);

                // checkAssertions();

                if (can_return_true) return true;


                // checkAssertions();
            }

            checkAssertions();
        }

        checkAssertions();

        return false;
    };

    for ( int v=0; v<N; v++ ) {
        if constexpr (check_expensive_assertions) {
            assert(ranges::none_of(deg_in_S, std::identity{}));
            assert(ranges::none_of(deg_out_NS, std::identity{}));
            assert(ranges::none_of(was, std::identity{}));
            assert(ranges::none_of(inS, std::identity{}));
            assert(ranges::none_of(inNS, std::identity{}));
            assert(ranges::none_of(calculated_NS_outdeg, std::identity{}));
        }

        if (check(v)) {
            clearForS();
            removed_nodes.push_back(v);
            GraphUtils::removeNodeFromGraph(V,v);
        }
        else clearForS();
    }

    return removed_nodes;
}

