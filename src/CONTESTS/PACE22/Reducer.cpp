//
// Created by sylwester on 12/20/21.
//

#include <graphs/GraphUtils.h>
#include <utils/RandomNumberGenerators.h>
#include <graphs/GraphInducer.h>
#include <utils/StandardUtils.h>
#include <combinatorics/CombinatoricUtils.h>
#include <graphs/VertexCover/kernelization/KernelizerVC.h>
#include "CONTESTS/PACE22/Reducer.h"

#include "EDReducer.h"

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


    V = GraphUtils::getGraphForEdges(primary_edges);
    N = V.size();



    { // #TEST - here should be implemented a more efficient graph inducing than the following...
        auto indg = GraphInducer::induceByNonisolatedNodes(V);
        reduced_instance.primary_indg_nodes = indg.nodes;
        V = indg.V;
        N = V.size();
        was = helper = VB(N);
    }

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

    Stopwatch sw;
    string reducer_str = "reducer";
    sw.setLimit(reducer_str, cnf.reducer_max_time_millis);
    sw.start(reducer_str);

    vector<VCReduction*> secondary_reduce_liftables;
    KernelizedNodesReduction * knr = nullptr;

    auto addKNR = [&]( VI nodes ){
        if(nodes.empty()) return;
        if( knr == nullptr ) knr = new KernelizedNodesReduction(nodes);
        else knr->addToKer(nodes);
    };

    auto addLiftables = [&]( auto liftables ) {
        if(knr != nullptr){ secondary_reduce_liftables.push_back(knr); knr = nullptr; }
        for(auto *x : liftables) secondary_reduce_liftables.push_back(x);
    };

    function<void()> applyBasicReductions = [&](){
        VVI Vcp = V;
        KernelizerVC kern;

        Stopwatch s; string opt = "basic_kern"; s.start(opt);
        // kern.use_crown_and_lp_checks = false;
        auto [kern_nodes, edges_removed] = kern.initialKernelization(Vcp);
        s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

        addKNR(kern_nodes);
        GraphUtils::removeNodes(V, kern_nodes,helper);
    };


    int ed_rules_checked = 0;

    do{

        modified = false;
        helper = VB(N,false);

        applyBasicReductions();

        if(modified) continue;
        if (sw.tle(reducer_str)) break;


        if(cnf.reducer_use_folding){
            Stopwatch s; string opt = "folding"; s.start(opt);
            auto folds = folding();
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            if(write_progress_on_the_fly) DEBUG(total_folds_done);

            if(!folds.empty()) modified = true;
            addLiftables(folds);
            assert( GraphUtils::isSimple(V) );

            if(modified) continue;
        }


        if(cnf.reducer_use_unconfined){
            Stopwatch s; string opt = "unconfined"; s.start(opt);
            clog << "Running unconfined" << endl;
            VI uncon = unconfined();
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);


            addKNR(uncon);
            total_unconfined_nodes += uncon.size();
            GraphUtils::removeNodes(V, uncon, helper);

            modified |= !uncon.empty();
            if(modified) continue;
        }

        if (false)
        if(cnf.reducer_use_desk){
            Stopwatch s; string opt = "desk"; s.start(opt);
            vector<VCReduction*> desk_liftables = desk();
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            for (auto l : desk_liftables) {
                if ( l->name() == "desk" ) total_desk_folds++;
                if ( l->name() == "knr" ) total_desk_dominations += l->offset();
            }
            addLiftables(desk_liftables);

            modified |= !desk_liftables.empty();
            if(modified) continue;
        }

        if (false)
        if(cnf.reducer_use_funnel){
            if(!cnf.reducer_use_domination){
                clog << "CAUTION! Calling funnel reduction without domination rule before!" << endl;
            }
            if(write_progress_on_the_fly) DEBUG(total_funnels_done);

            Stopwatch s; string opt = "funnel"; s.start(opt);
            auto funnel_liftables = funnel();
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            total_funnels_done += funnel_liftables.size();
            if(write_progress_on_the_fly) DEBUG(total_funnels_done);

            addLiftables(funnel_liftables);

            modified |= !funnel_liftables.empty();
            if(modified) continue;
        }

        // standard node-removal version
        bool ed_application_cond = ( cnf.ed_application_mode == 0 || (cnf.ed_application_mode == 1 && ed_rules_checked == 0) );
        if( cnf.reducer_use_ed && cnf.ed_use_node_removal && ed_application_cond){
            ed_rules_checked++;
            clog << "Running ED node removal rules in POINT-1" << endl;

            EDReducer edred(V.size(), cnf);
            edred.resetAllUsedTechniques();
            edred.cnf.ed_use_node_removal = true;

            Stopwatch s; string opt = "ED node removal"; s.start(opt);
            VI res = edred.reduce(V);
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            assert(res.size() == edred.last_reduce_nodes_removed);
            ed_nodes_reduced += edred.last_reduce_nodes_removed;
            ed_edges_removed += edred.last_reduce_edges_removed;
            ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;

            addKNR(res);
            if (edred.last_reduce_edges_removed > 0) DEBUG(edred.last_reduce_edges_removed);
            if (edred.madeChangesInLastReduce())  V = edred.getV();
            modified |= edred.madeChangesInLastReduce();

            if(modified) continue;
        }


        // if(cnf.reducer_use_folding_twins) {
        //     Stopwatch s; string opt = "folding twins"; s.start(opt);
        //     auto [twin_folds, to_remove] = foldingTwins();
        //     if (write_progress_on_the_fly) DEBUG(total_twin_folds_done);
        //     total_twin_folds_done += twin_folds.size() + to_remove.size();
        //     if (write_progress_on_the_fly) DEBUG(total_twin_folds_done);
        //     addKNR(to_remove);
        //     Utils::removeNodes(V, revV, to_remove, helper);
        //     if(!modified) modified = (!twin_folds.empty() || !to_remove.empty());
        //     s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
        //
        //     { // add to resulting kernelization objects
        //         if (knr != nullptr) { res.push_back(knr); knr = nullptr; }
        //         for (auto *x : twin_folds) res.push_back(x);
        //     }
        //
        //     if(modified) continue;
        // }


        {
            VVI Vcp = V;
            KernelizerVC kern;
            auto [kern_nodes, edges_removed] = kern.lpDecomposition(Vcp);
            addKNR(kern_nodes);
            GraphUtils::removeNodes(V, kern_nodes,helper);
            if(!modified) modified = (!kern_nodes.empty());
        }


        if(cnf.reducer_use_general_folding){
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

        if (false)
        if(cnf.reducer_use_twins){
            Stopwatch s; string opt = "twins merge"; s.start(opt);
            if(write_progress_on_the_fly) DEBUG(total_twins_merged);

            vector<VCReduction*> liftables = twins();
            secondary_reduce_liftables += liftables;

            if(write_progress_on_the_fly) DEBUG(total_twins_merged);
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            modified |= !liftables.empty();
            if(modified) continue;
        }


        // standard node-removal version
        if( cnf.reducer_use_ed && cnf.ed_use_node_removal){
            ed_rules_checked++;
            Stopwatch s; string opt = "ED node removal"; s.start(opt);
            clog << "Running ED node removal rules in POINT-2" << endl;

            EDReducer edred(V.size(), cnf);
            edred.resetAllUsedTechniques();
            edred.cnf.ed_use_node_removal = true;

            VI res = edred.reduce(V);
            assert(res.size() == edred.last_reduce_nodes_removed);
            ed_nodes_reduced += edred.last_reduce_nodes_removed;
            ed_edges_removed += edred.last_reduce_edges_removed;
            ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;
            addKNR(res);

            if (edred.madeChangesInLastReduce())  V = edred.getV();
            modified |= edred.madeChangesInLastReduce();

            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            if(modified) continue;
        }


        // standard edge-insertion - those edges that are found using consider(v) for single-node initial sets S
        if (cnf.ed_apply_type1_constraints_on_the_fly)
        if ( cnf.reducer_use_ed && cnf.ed_use_edge_insertion) {
            VI res;
            bool made_changes = false;
            do {
                ed_rules_checked++;
                Stopwatch s; string opt = "ED edge insertion"; s.start(opt);
                clog << "Running ED with edge insertion" << endl;

                EDReducer edred(V.size(), cnf);
                edred.resetAllUsedTechniques();
                edred.cnf.ed_use_node_removal = true;
                edred.cnf.ed_apply_type1_constraints_on_the_fly = true;

                res = edred.reduce(V);
                assert(res.size() == edred.last_reduce_nodes_removed);
                ed_nodes_reduced += edred.last_reduce_nodes_removed;
                ed_edges_removed += edred.last_reduce_edges_removed;
                ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;

                addKNR(res);
                if (edred.madeChangesInLastReduce())  V = edred.getV();
                assert( GraphUtils::isSimple(V) );
                made_changes = edred.madeChangesInLastReduce();

                ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;
                modified |= edred.madeChangesInLastReduce();

                s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            }while (res.empty() && made_changes);

            if(modified) continue;
        }

    } while(modified);

    if(debug){ DEBUG(V);}

    if(knr != nullptr){ secondary_reduce_liftables.push_back(knr); knr = nullptr;}


    return make_pair(V, secondary_reduce_liftables);
}

vector<VCReduction*> Reducer::twins() {
    assert(false && "Implement folding twins");
    return {};
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

                for (int d : V[a]) if (V[d].size() == 1) liftables += propagateDeg1RuleSlow(d);
                for (int d : V[b]) if (V[d].size() == 1) liftables += propagateDeg1RuleSlow(d);

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
    DEBUG(total_folds_done);
    DEBUG(total_general_folds_done);
    DEBUG(total_twin_folds_done);
    DEBUG(total_desk_folds);
    DEBUG(total_desk_dominations);
    DEBUG(total_unconfined_nodes);
    DEBUG(total_funnels_done);
    DEBUG(total_twins_merged);

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




vector<FunnelReduction *> Reducer::funnel() {
    assert(false && "funnel not implemented yet, should be incorporated to the main \'basic workflow\'");
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
   assert(false && "Implement desk efficiently for VC - both the desk domination and desk folding");
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
    reverse(ALL(order)); // #TEST - starting node selection from largest degree in general_folding

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

        int vs_size;
        InducedGraph g = GraphInducer::induce(V, W);
        vs_size = Utils::getMinVcCPSAT(g.V).size();

        if( vs_size + 2 < W.size() ) continue;

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
    VI removed_nodes;

    VB inS(N), inNS(N);
    VI S;
    VI deg_in_S(N,0);
    VI deg_out_NS(N,0);
    int ns_size = 0;
    VB calculated_NS_outdeg(N);

    VI cand_NS;

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
        // return 0;

        if (calculated_NS_outdeg[d]) return deg_out_NS[d];
        return max(0, (int)V[d].size() - deg_in_S[d] - ns_size);
    };

    auto calculateNSOutdeg = [&](int v) {
        calculated_NS_outdeg[v] = true;
        deg_out_NS[v] = 0;
        for (int d : V[v]) deg_out_NS[v] += !inNS[d];
    };

    auto check = [&](int v) {
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
            int u = cand_NS.back();
            cand_NS.pop_back();
            if ( deg_in_S[u] != 1 ) continue;
            if ( deg_out_NS[u] == 0 ) return true;
            assert(deg_out_NS[u] == 1);

            inS[u] = inNS[u] = true;
            S.push_back(u);
            deg_out_NS[u] = 0;

            for ( int d : V[u] ) deg_in_S[d]++;
            for (int d : V[u]) if (!inS[d]) { // moving conceptually node d to N(S)
                if ( inNS[d] ) continue; // we can consider only those nodes that are not in NS, as those in NS have deg_in_S > 1
                assert(!inNS[d]);
                inNS[d] = true;
                ns_size++;

                for ( int dd : V[d] ) {
                    assert((dd == u) || !inS[dd]);
                    if ( inNS[dd] ) {
                        if ( getLbNSOutdegForUnexpandedNode(dd) <= 1 ) {
                            if ( !calculated_NS_outdeg[dd] ) calculateNSOutdeg(dd);
                            else deg_out_NS[dd]--;

                            if (deg_out_NS[dd] == 0) return true;
                            if (deg_out_NS[dd] == 1) cand_NS.push_back(dd);
                        }
                    }else { // dd is not in NS
                        deg_out_NS[d]++;
                    }
                }

                if (!calculated_NS_outdeg[d]) calculateNSOutdeg(d);
                if (calculated_NS_outdeg[d] && deg_out_NS[d] == 1) cand_NS.push_back(d);
            }

        }

        return false;
    };

    for ( int v=0; v<N; v++ ) {
        if (check(v)) {
            clearForS();
            removed_nodes.push_back(v);
            GraphUtils::removeNodeFromGraph(V,v);
        }
        else clearForS();
    }

    return removed_nodes;
}

