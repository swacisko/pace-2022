//
// Created by sylwester on 12/20/21.
//

#include <graphs/GraphUtils.h>
#include <graphs/scc/StronglyConnectedComponents.h>
// #include <utils/TimeMeasurer.h>
#include <utils/RandomNumberGenerators.h>
#include <graphs/GraphInducer.h>
#include <CONTESTS/PACE22/exact/DFVSSolverE.h>
#include <utils/StandardUtils.h>
#include <combinatorics/CombinatoricUtils.h>
#include <graphs/cliques/CliqueUtils.h>
#include <graphs/components/ConnectedComponents.h>
#include <graphs/VertexCover/VCUtils.h>
#include <graphs/VertexCover/kernelization/KernelizerVC.h>
#include <graphs/cliques/CliqueExtension.h>
#include <graphs/graphtraversals/BFS.h>
#include "CONTESTS/PACE22/Reducer.h"

#include "EDReducer.h"

using namespace Utils;

Reducer::Reducer(VVI &V, Config c) : origN(V.size()) {
    cnf = c;
    this->V = V;
    N = V.size();

    hashes = VLL(N);
    UniformIntGenerator rnd(0,1'000'000'000ll * 1'000'000'000);
    for(int i=0; i<N; i++) hashes[i] = rnd.rand();
}


vector<DFVSReduction*> Reducer::reduce(VVI _revV) {
    constexpr bool debug = false;
    constexpr bool write_progress_on_the_fly = false;

    VVI prevV = V;
    bool modified;


    auto reducer_start_time = chrono::steady_clock::now();
    VB helper(N,false);

    vector<DFVSReduction*> res;
    KernelizedNodesReduction * knr = nullptr;

    auto addKNR = [&]( VI nodes ){
        if(nodes.empty()) return;
        if( knr == nullptr ) knr = new KernelizedNodesReduction(nodes);
        else knr->addToKer(nodes);
    };

    function<void()> applyBasicReductions = [&](){

    };


    int ed_rules_checked = 0;

    do{

        {
            Stopwatch s; string opt = "basic_kern"; s.start(opt);
            VVI Vcp = V;
            KernelizerVC kern;
            auto [kern_nodes, edges_removed] = kern.initialKernelization(Vcp);
            addKNR(kern_nodes);
            // Utils::removeNodes(V, revV, kern_nodes,helper);
            GraphUtils::removeNodes(V, kern_nodes,helper);
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
        }


        modified = false;
        helper = VB(N,false);

        applyBasicReductions();

        if(modified) continue;

        auto time_total = chrono::duration<double, std::milli >
                (chrono::steady_clock::now() - reducer_start_time ).count();
        if(time_total > cnf.reducer_max_time_millis) break;




        if(cnf.reducer_use_unconfined){
            Stopwatch s; string opt = "unconfined"; s.start(opt);
            VI uncon = unconfined();
            addKNR(uncon);
            total_unconfined_nodes += uncon.size();
            // Utils::removeNodes(V, revV, uncon, helper);
            GraphUtils::removeNodes(V, uncon, helper);
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            if(!modified) modified = (!uncon.empty());
            if(modified) continue;
        }


        // standard node-removal version
        bool ed_application_cond =  ( cnf.ed_application_mode == 0 || (cnf.ed_application_mode == 1 && ed_rules_checked == 0) );
        if( cnf.reducer_use_ed && cnf.ed_use_node_removal && ed_application_cond){
            ed_rules_checked++;
            Stopwatch s; string opt = "ED node removal"; s.start(opt);
            // clog << "Running ED node removal rules in POINT-1" << endl;

            EDReducer edred(V.size(), cnf);
            edred.resetAllUsedTechniques();
            edred.cnf.ed_use_node_removal = true;

            VI res = edred.reduce(V);
            assert(res.size() == edred.last_reduce_nodes_removed);
            ed_nodes_reduced += edred.last_reduce_nodes_removed;
            ed_edges_removed += edred.last_reduce_edges_removed;
            ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;

            addKNR(res);
            if (edred.last_reduce_edges_removed > 0) DEBUG(edred.last_reduce_edges_removed);
            if (edred.madeChangesInLastReduce())  V = edred.getV();
            modified |= edred.madeChangesInLastReduce();

            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            if(modified) continue;
        }


        if(cnf.reducer_use_folding){
            Stopwatch s; string opt = "folding"; s.start(opt);
            auto folds = folding();
            if(write_progress_on_the_fly) DEBUG(total_folds_done);
            total_folds_done += folds.size();
            if(write_progress_on_the_fly) DEBUG(total_folds_done);
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);

            if(!folds.empty()) modified = true;

            { // add to resulting kernelization objects
                if(knr != nullptr){ res.push_back(knr); knr = nullptr; }
                for(auto *x : folds) res.push_back(x);
            }

            assert( GraphUtils::isSimple(V) );

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

        if(cnf.reducer_use_desk){
            Stopwatch s; string opt = "desk"; s.start(opt);
            auto [desk_folds, desk_dominations, arc_diff] = desk();
            total_desk_folds += desk_folds.size();
            total_desk_dominations += desk_dominations.size();
            total_desk_arcs_added += arc_diff;
            addKNR(desk_dominations);
            { // add to resulting kernelization objects
                if(knr != nullptr){ res.push_back(knr); knr = nullptr; }
                for(auto *x : desk_folds) res.push_back(x);
            }
            if(!modified) modified = ( !desk_folds.empty() || !desk_dominations.empty() || arc_diff );

            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            if(modified) continue;
        }

        if(cnf.reducer_use_funnel){
            Stopwatch s; string opt = "funnel"; s.start(opt);
            if(!cnf.reducer_use_domination){
                clog << "CAUTION! Calling funnel reduction without domination rule before!" << endl;
            }
            //TimeMeasurer::start("Reducer::funnel");
            if(write_progress_on_the_fly) DEBUG(total_funnels_done);
            auto funnels = funnel();
            total_funnels_done += funnels.size();
            if(write_progress_on_the_fly) DEBUG(total_funnels_done);
            if(!funnels.empty()) modified = true;

            { // add to resulting kernelization objects
                if(knr != nullptr){ res.push_back(knr); knr = nullptr; }
                for(auto *x : funnels) res.push_back(x);
            }

            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            if(modified) continue;
        }


        {
            VVI Vcp = V;
            KernelizerVC kern;
            auto [kern_nodes, edges_removed] = kern.lpDecomposition(Vcp);
            addKNR(kern_nodes);
            // Utils::removeNodes(V, revV, kern_nodes,helper);
            GraphUtils::removeNodes(V, kern_nodes,helper);
            if(!modified) modified = (!kern_nodes.empty());
        }


        if(cnf.reducer_use_general_folding){
            Stopwatch s; string opt = "general folding"; s.start(opt);
            auto reductions = generalFolding();
            total_general_folds_done += reductions.size();
            if(!modified) modified = (!reductions.empty());

            { // add to resulting kernelization objects
                if(knr != nullptr){ res.push_back(knr); knr = nullptr; }
                for(auto *x : reductions) res.push_back(x);
            }

            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            if(modified) continue;
        }


        if(cnf.reducer_use_twins_merge){
            Stopwatch s; string opt = "twins merge"; s.start(opt);
            if(write_progress_on_the_fly) DEBUG(total_twins_merged);

            bool mod = mergeTwins();

            if(write_progress_on_the_fly) DEBUG(total_twins_merged);
            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            modified |= mod;
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
        // if (false)
        if ( cnf.reducer_use_ed && cnf.ed_use_edge_insertion) {
            ed_rules_checked++;
            Stopwatch s; string opt = "ED edge insertion"; s.start(opt);
            clog << "Running ED with edge insertion" << endl;

            EDReducer edred(V.size(), cnf);
            edred.resetAllUsedTechniques();
            edred.cnf.ed_use_node_removal = true;
            edred.cnf.ed_apply_type1_constraints_on_the_fly = true;

            VI res = edred.reduce(V);
            assert(res.size() == edred.last_reduce_nodes_removed);
            ed_nodes_reduced += edred.last_reduce_nodes_removed;
            ed_edges_removed += edred.last_reduce_edges_removed;
            ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;

            addKNR(res);
            if (edred.madeChangesInLastReduce())  V = edred.getV();
            assert( GraphUtils::isSimple(V) );

            ed_t1_inference_rules_added += edred.last_reduce_inf_rules_1_added;
            modified |= edred.madeChangesInLastReduce();

            s.stop(opt); reduction_times_millis[opt] += s.getTime(opt);
            if(modified) continue;
        }


    }while(modified);

    if(debug){ DEBUG(V);}

    if(knr != nullptr){ res.push_back(knr); knr = nullptr;}
    return res;
}









pair<vector<FoldingTwinReduction*>, VI> Reducer::foldingTwins() {

}

bool Reducer::mergeTwins() {
    assert(false && "Implement twin merging");

}


vector<FoldingReduction*> Reducer::folding() {

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







void Reducer::liftSolution(int N, VI &dfvs, vector<DFVSReduction *> &reductions, bool clear_reductions) {
    VB in_dfvs = StandardUtils::toVB(N, dfvs);

    for( int i = (int)reductions.size()-1; i>=0; i-- ){
        reductions[i]->lift(dfvs, in_dfvs);
    }

    if (clear_reductions) clearReductionObjects(reductions);
}

void Reducer::clearReductionObjects(vector<DFVSReduction *> &reductions) {
    for(int i=0; i<reductions.size(); i++ ){
        delete reductions[i];
        reductions[i] = nullptr;
    }
}

VI Reducer::convertKernelizedReductions(vector<DFVSReduction *> &reductions) {
    assert(reductions.size() <= 1);
    VI red_dfvs;
    if(!reductions.empty()){
        KernelizedNodesReduction * knr = (KernelizedNodesReduction*) reductions[0];
        red_dfvs = knr->getKer();
        Reducer::clearReductionObjects(reductions);
    }
    return red_dfvs;
}

int Reducer::getReductionsOffset(vector<DFVSReduction *> &reductions) {
    int res = 0;
    for(auto * x : reductions) res += x->sizeDiffUB();
    return res;
}

void Reducer::writeReductions(vector<DFVSReduction *> &reductions) {
    clog << "Reductions: " << endl;
    for(auto * x : reductions) clog << x->toString() << endl;
}




tuple<vector<DeskReduction*>,VI, int> Reducer::desk(){
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

        int dfvs_size;
        InducedGraph g = GraphInducer::induce(V, W);
        dfvs_size = Utils::getMinVcCPSAT(g.V).size();

        if( dfvs_size + 2 < W.size() ) continue;

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

    if(debug) clog << "Starting unconfined" << endl;

    int N = V.size();
    VB affected(N,false);
    VB helper(N,false);
    VB was(N,false);
    VB in_S(N,false);
    VB in_NS(N,false);

    VI res;

    VVI & G = V;

    VI order = CombinatoricUtils::getRandomPermutation(N);
    sort(ALL(order), [&](int a, int b){ return G[a].size() > G[b].size(); } );


    auto unconfined = [&](int v){
        VI S = {v};
        VI ws;  // N(u) \ N[S]
        VI NS = G[v];
        in_S[v] = true;
        for(int d : NS) in_NS[d] = true;

        bool can = true;

        while(can) {

            if(debug){DEBUG(S); DEBUG(NS);}

            int best_u = -1, best_val = 1e9;
            VI best_ws;
            for (int u : NS) {
                int cnt = 0;
                ws.clear();
                for( int d : G[u] ){
                    if( in_S[d] ) cnt++;
                    if( !in_S[d] && !in_NS[d] ) ws.push_back(d);
                }
                if(cnt != 1) continue;
                if( ws.size() < best_val ){
                    best_val = ws.size();
                    best_u = u;
                    best_ws = ws;
                }
            }

            if(debug){ DEBUG(best_u); DEBUG(best_ws); }

            if( best_u == -1 ){ can = false;break; }
            if( best_ws.empty() ){
                /* node v is unconfined*/
                if(debug) clog << "Node v is unconfined!" << endl;
                can = true;break;
            }
            if( best_ws.size() == 1 ){
                int w = best_ws[0];
                if(debug) clog << "Pushing node w: " << w << " to S" << endl;
                S.push_back(w);
                in_S[w] = true;
                for( int u : G[w] ){
                    if(!in_S[u] && !in_NS[u]){
                        NS.push_back(u);
                        in_NS[u] = true;
                    }
                }
            }else{ can = false; break; }
        }

        VI X = S + NS;
        for( int d : X ){
            in_S[d] = in_NS[d] = false;
            if( affected[d] ) can = false;
        }

        return can;
    };

    for( int v : order ){
        if( affected[v] ) continue;

        bool aff = false;
        for(int d : G[v]) if(affected[d]) aff = true;
        if(aff) continue;

        bool unconf = unconfined(v);
        if(unconf){
            res.push_back(v);
            affected[v] = true;
            for( int d : G[v] ) affected[d] = true;
        }
    }

    return res;
}

