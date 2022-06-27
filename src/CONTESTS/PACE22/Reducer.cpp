//
// Created by sylwester on 12/20/21.
//

#include <graphs/GraphUtils.h>
#include <graphs/scc/StronglyConnectedComponents.h>
#include <utils/TimeMeasurer.h>
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

using namespace Utils;

Reducer::Reducer(VVI &V, Config c) : origN(V.size()) {
    cnf = c;
    this->V = V;
    N = V.size();

    hashes = VLL(N);
    UniformIntGenerator rnd(0,1'000'000'000ll * 1'000'000'000);
    for(int i=0; i<N; i++) hashes[i] = rnd.rand();
}

VI Reducer::loop() {
    VI res;
    for(int i=0; i<N; i++) if( hasLoop(V,i) ) res.push_back(i);
    return res;
}

VI Reducer::inOut1() {
    VI res;
    VI inDeg(N,0);
    for( int i=0; i<N; i++ ) for(int d : V[i]) inDeg[d]++;
    for(int i=0; i<N; i++){
        if( V[i].size() != 0 || inDeg[i] != 0 ){
            if(V[i].size() <= 1 || inDeg[i] <= 1) res.push_back(i);
        }
    }
    return res;
}

VI Reducer::inOutClique() {
    TimeMeasurer::start("Reducer::inOutClique");
    const bool debug = false;

    int N = V.size();
    VB helper(N,false);

    VPII pie_edges = getAllPIEdges(V, revV, helper);

    if(debug) DEBUG(pie_edges);

    VVI H(N);
    for( auto [a,b] : pie_edges ) H[a].push_back(b);

    VVB neigh_marker(N);
    int E = GraphUtils::countEdges(H);
    for( int i=0; i<N; i++ ){
        if( 1ll * H[i].size() * H[i].size() >= 1ll * E ){
            neigh_marker[i] = VB(N,false);
            for( int d : H[i] ) neigh_marker[i][d] = true;
        }
    }

    auto isDClique = [&](VI & A){
        bool can = true;
        for( int d : A ) helper[d] = true;

        for (int d : A) {
            int cnt = 1;

            if( neigh_marker[d].empty() ){
                for (int v : H[d]) if(helper[v]) cnt++;
            }else{
                for( int dd : A ){
                    if( dd == d ) continue;
                    if( neigh_marker[d][dd] ) cnt++;
                    else{
                        can = false;
                        break;
                    }
                }
            }

            if(cnt != A.size() ) can = false;
            if(!can) break;
        }

        for( int d : A ) helper[d] = false;

        return can;
    };

    VI to_merge;

    VB affected(N,false);

    for( int i=0; i<N; i++){
        if( V[i].empty() && revV[i].empty() ) continue;

        {
            if( affected[i] ) continue;
            bool aff = false;
            VI temp = StandardUtils::setUnion(V[i], revV[i], helper);
            for(int d : temp) if(affected[d]) aff = true;
            if(aff) continue;
        }

        bool isInClique = false;
        if( !revV[i].empty() ) isInClique = isDClique(revV[i]);

        bool isOutClique = false;
        if(!isInClique && !V[i].empty() ) isOutClique = isDClique(V[i]);

        if( isInClique || isOutClique ){
            if(debug) clog << "inOutClique node " << i << " ---> merging node " << i << endl;
            to_merge.push_back(i);

            affected[i] = true;
            VI temp = StandardUtils::setUnion(V[i], revV[i], helper);
            for(int d : temp) affected[d] = true;
        }
    }

    TimeMeasurer::stop("Reducer::inOutClique");

    return to_merge;
}

VI Reducer::core() {
    TimeMeasurer::start("Reducer::core");
    const bool debug = false;

    int N = V.size();
    VB helper(N,false);

    VPII pie_edges = getAllPIEdges(V, revV, helper);

    if(debug) DEBUG(pie_edges);

    VVI H(N);
    for( auto [a,b] : pie_edges ) H[a].push_back(b);

    VI order;
    for(int i=0; i<N; i++){
        if( !H[i].empty() && H[i].size() == V[i].size() ){
            if(cnf.reducer_use_pie) assert( H[i].size() == revV[i].size() ); // this will hold if pie() was called before core()
            order.push_back(i);
        }
    }

    if(debug){
        DEBUG(order.size());
        DEBUG(order);
        Utils::writeRemainingGraph(H);
    }

    sort(ALL(order),[&](int a, int b){
       return V[a].size() < V[b].size();
    });

    VVB neigh_marker(N);
    long long E = GraphUtils::countEdges(H);
    for( int i=0; i<N; i++ ){
        if( H[i].empty() ) continue;
        if( 1ll * H[i].size() * H[i].size() >= E ){
            neigh_marker[i] = VB(N,false);
            for( int d : H[i] ) neigh_marker[i][d] = true;
        }
    }

    auto inducesClique = [&](int u){
        for( int d : H[u] ) if( H[d].size() < H[u].size() ) return false;

        bool can = true;
        for( int d : H[u] ) helper[d] = true;

        sort(ALL(H[u]), [&](int a, int b){
           return H[a].size() < H[b].size();
        });

        for (int d : H[u]) {
            int cnt = 1;

            if( neigh_marker[d].empty() ){
                for (int v : H[d]) if(helper[v]) cnt++;
            }else{
                for( int dd : H[u] ){
                    if( dd == d ) continue;
                    if( neigh_marker[d][dd] ) cnt++;
                    else{
                        can = false;
                        break;
                    }
                }
            }

            if( cnt != H[u].size() ) can = false;
            if(!can) break;
        }

        for( int d : H[u] ) helper[d] = false;
        return can;
    };

    VB was(N,false);
    set<int> zb;
    for( int u : order ){
        if(!was[u]){
            if( inducesClique(u) ){
                if(debug) clog << "Core node " << u << " --> adding H[u] = {" << H[u] << "} to remove" << endl;
                zb.insert(ALL(H[u]));
            }
        }
        for(int d : H[u]) was[d] = true;
        was[u] = true;
    }

    TimeMeasurer::stop("Reducer::core");

    return VI(ALL(zb));
}


VI Reducer::strongly_connected() {
    StronglyConnectedComponents scc(V);
    scc.createStronglyConnectedComponents();

    VI indeg(N,0);
    for(int i=0; i<N; i++) for(int d : V[i]) indeg[d]++;

    VI res;
    VVI comps = scc.getComponents();
    for( VI & v : comps ){
        if( v.size() == 1 ){
            assert(!hasLoop(V,v[0]));
            for(int d : v) if( V[d].size() > 0 || indeg[d] > 0 ) res.push_back(d);
        }
    }

    return res;
}

VVI Reducer::pathCompression(VVI & revV) {
    VB was(N,false);

    VVI res;
    for( int i=0; i<N; i++ ){
        if(was[i] || V[i].size() != 1 || revV[i].size() != 1) continue;
        VI pth_prev, pth_succ;

        int p = i;
        while( !was[p] && V[p].size() == 1 && revV[p].size() == 1 ){
            pth_succ.push_back(p);
            was[p] = true;
            p = V[p][0];
        }
        if(p != i) pth_succ.push_back(p);

        p = i;
        was[i] = false;
        while( !was[p] && V[p].size() == 1 && revV[p].size() == 1 ){
            pth_prev.push_back(p);
            was[p] = true;
            p = revV[p][0];
        }
        pth_prev.push_back(p);
        was[p] = true;

        reverse(ALL(pth_prev));
        pth_prev.pop_back();
        VI pth = pth_prev + pth_succ;

        if( pth.size() >= 2 && pth[0] == pth.back() ) pth.pop_back();

        res.push_back( pth );
    }

    return res;
}

vector<DFVSReduction*> Reducer::reduce(VVI _revV) {
    constexpr bool debug = false;
    constexpr bool write_progress_on_the_fly = false;
    constexpr bool check_correspondings = false; // set to true to run a lot of assertion trying to find hidden bugs

    VVI prevV = V;
    bool modified;

    int nonsimple_arc_full_cnt = 0;
    int domination4_cnt = 0;
    int mixed_domination_full_cnt = 0;
    int domination5_cnt = 0;

    int core_reductions_done = 0;

    if( !_revV.empty() ) revV = _revV;
    else revV = GraphUtils::reverseGraph(V);

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
        bool local_improvement = false;

        {
            VI loops = loop();
            addKNR(loops);

            if(!modified) modified = !loops.empty();
            if(!local_improvement) local_improvement = !loops.empty();

            for (int d : loops) {
                removeNode(V,revV,d,helper);
                if(debug)clog << "Removing node with loop: " << d << endl;
            }
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(debug){ DEBUG(V);DEBUG(revV); }


        if(cnf.reducer_use_strongly_connected){ // strongly_connected
            VI to_remove = strongly_connected();
            if(debug) clog << "Removing nodes: " << to_remove << " - not in strongly_connected components" << endl;
            removeNodes(V,revV, to_remove, helper);
            if(!modified) modified = !to_remove.empty();
            if(!local_improvement) local_improvement = !to_remove.empty();
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(debug){ DEBUG(V);DEBUG(revV); }

        { // path compression
            // path compression could be done using merging nodes with indeg1/outdeg1, but this way is much faster
            VVI paths = pathCompression(revV);
            bool creates_loop = false, path_compression_modified = false;
            for( VI & v : paths ){
                if(debug) clog << "Compressing path: " << v << endl;
                if(v.size() <= 2) continue;
                for( int i=1; i<v.size(); i++ ){
                    int a = v[i-1], b = v[i];
                    removeEdge(V,a,b);
                    removeEdge(revV,b,a);
                    path_compression_modified = true;
                }
                if( v.size() > 1 ){
                    addEdge(V,v[0],v.back());
                    addEdge(revV,v.back(), v[0]);
                }
                if( hasLoop(V,v[0]) ) creates_loop = true;
            }

            if(!modified) modified = !paths.empty();
            if(!local_improvement) local_improvement = !paths.empty();

            if(creates_loop){
                assert(modified);
                assert(local_improvement);
                applyBasicReductions();
                if(check_correspondings) assert(Utils::isCorresponding(V, revV));
                return;
            }
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(debug){ DEBUG(V);DEBUG(revV); }

        {
            VI in1out1 = inOut1();

            for( int d : in1out1 ){
                if( V[d].size() > 1 && revV[d].size() > 1 ) continue;

                if(!hasLoop(V,d)){
                    merge( V,revV, d, helper );
                    modified = local_improvement = true;
                    if(debug) clog << "Merging in1/out1 node: " << d << endl;
                }
            }
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(local_improvement) applyBasicReductions();
    };

    if(check_correspondings) assert(Utils::isCorresponding(V, revV));

    if(Utils::isPIGraph(V, revV, helper)){
        VVI Vcp = V;
        KernelizerVC kern;
        auto [kern_nodes, edges_removed] = kern.initialKernelization(Vcp);
        addKNR(kern_nodes);
        Utils::removeNodes(V, revV, kern_nodes,helper);
    }

    do{

        if(write_progress_on_the_fly){
            StronglyConnectedComponents scc(V, revV);
            scc.createStronglyConnectedComponents();
            auto comps = scc.getComponents();
            sort(ALL(comps),[]( auto & v1, auto & v2 ){
                return v1.size() < v2.size();
            });
            int sccs = comps.size();
            clog << "\tThere are " << sccs << " strongly connected components in the graph, largest with "
                 << comps.back().size() << " nodes" << endl;
        }

        if(debug){ clog << endl << "REDUCER: Starting new reducing iteration" << endl;DEBUG(V); }

        modified = false;
        helper = VB(N,false);
        if(debug) DEBUG(revV);


        if(check_correspondings && !Utils::isCorresponding(V, revV)){
            Utils::writeRemainingGraph(V);
            Utils::writeRemainingGraph(revV);
            assert(Utils::isCorresponding(V, revV));
        }

        applyBasicReductions();

        if(check_correspondings){
            assert(GraphUtils::isSimple(V));
            assert(GraphUtils::isSimple(revV));
            if(!Utils::isCorresponding(V, revV)){
                DEBUG(V);
                DEBUG(revV);
                assert(Utils::isCorresponding(V, revV));
            }
        }

        if(modified) continue;

        auto time_total = chrono::duration<double, std::milli >
                (chrono::steady_clock::now() - reducer_start_time ).count();
        if(time_total > cnf.reducer_max_time_millis) break;


        if(cnf.reducer_use_pie){
            VPII edges_to_remove = pie();
            if(write_progress_on_the_fly) DEBUG(total_pie_edges_removed);
            total_pie_edges_removed += edges_to_remove.size();
            if(write_progress_on_the_fly) DEBUG(total_pie_edges_removed);

            if( !edges_to_remove.empty() ){
                if(debug) clog << "There are " << edges_to_remove.size() << " edges to remove using PIE, out of "
                     << GraphUtils::countEdges(V, true) << " total edges in V" << endl;

                Utils::removeEdges(V, edges_to_remove, helper);
                for( auto & [a,b] : edges_to_remove ) swap(a,b);
                Utils::removeEdges(revV, edges_to_remove, helper);
                modified = true;
            }

            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_core && core_reductions_done <= 1){
            core_reductions_done++;
            VI core_nodes_to_remove = core();

            if(write_progress_on_the_fly) DEBUG(total_core_nodes_removed);
            total_core_nodes_removed += core_nodes_to_remove.size();
            if(write_progress_on_the_fly) DEBUG(total_core_nodes_removed);

            if(!core_nodes_to_remove.empty()){
                if(debug) clog << "There are " << core_nodes_to_remove.size() << " nodes to remove, using CORE" << endl;
            }
            if(!core_nodes_to_remove.empty()){
                removeNodes(V, revV, core_nodes_to_remove, helper);
                addKNR(core_nodes_to_remove);
                modified = true;
            }

            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_inoutclique){
            if(write_progress_on_the_fly) DEBUG(total_inoutclique_nodes_merged);
            VI to_merge = inOutClique();

            for( int d : to_merge ){
                if(!hasLoop(V,d)) {
                    merge(V, revV, d, helper);
                    total_inoutclique_nodes_merged++;
                    modified = true;
                }
            }
            if(write_progress_on_the_fly) DEBUG(total_inoutclique_nodes_merged);

            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_dome){ // dome
            VPII dominated_arcs = dome();
            if(write_progress_on_the_fly) DEBUG(total_dome_edges_removed);
            total_dome_edges_removed += dominated_arcs.size();
            if(write_progress_on_the_fly) DEBUG(total_dome_edges_removed);

            Utils::removeEdges(V, dominated_arcs, helper);
            for( auto & [a,b] : dominated_arcs ) swap(a,b);
            Utils::removeEdges(revV, dominated_arcs, helper);
            if( !dominated_arcs.empty() ){
                if(debug) clog << "There are " << dominated_arcs.size() << " dominated arcs to remove using DOME" << endl;
                modified = true;
            }

            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_domination){
            TimeMeasurer::start("Reducer::domination1");
            VI dominated = domination1();
            addKNR(dominated);

            if(write_progress_on_the_fly) DEBUG(total_dominated_nodes1);
            total_dominated_nodes1 += dominated.size();
            if(write_progress_on_the_fly) DEBUG(total_dominated_nodes1);
            Utils::removeNodes(V, revV, dominated, helper);
            TimeMeasurer::stop("Reducer::domination1");
            if(!dominated.empty()) modified = true;
            if(modified) continue;



            TimeMeasurer::start("Reducer::domination2");
            if(write_progress_on_the_fly) DEBUG(total_dominated_nodes2);
            dominated = domination2();
            if(write_progress_on_the_fly) DEBUG(total_dominated_nodes2);
            addKNR(dominated);
            total_dominated_nodes2 += dominated.size();
            Utils::removeNodes(V, revV, dominated, helper);
            TimeMeasurer::stop("Reducer::domination2");
            if(!dominated.empty()) modified = true;
            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_unconfined && Utils::isPIGraph(V, revV, helper)){
            TimeMeasurer::start("Reducer::unconfined");
            VI uncon = unconfined();
            addKNR(uncon);
            total_unconfined_nodes += uncon.size();
            Utils::removeNodes(V, revV, uncon, helper);
            TimeMeasurer::stop("Reducer::unconfined");
            if(!modified) modified = (!uncon.empty());
            if(modified) continue;
        }

        if(cnf.reducer_use_folding){
            TimeMeasurer::start("Reducer::folding");
            auto folds = folding();
            if(write_progress_on_the_fly) DEBUG(total_folds_done);
            total_folds_done += folds.size();
            if(write_progress_on_the_fly) DEBUG(total_folds_done);
            if(!folds.empty()) modified = true;

            { // add to resulting kernelization objects
                if(knr != nullptr){
                    res.push_back(knr);
                    knr = nullptr;
                }
                for(auto *x : folds) res.push_back(x);
            }

            TimeMeasurer::stop("Reducer::folding");

            assert( GraphUtils::isSimple(V) );
            assert( GraphUtils::isSimple(revV) );

            if(modified) continue;
        }

        if(cnf.reducer_use_folding_twins) {
            TimeMeasurer::start("Reducer::folding_twins");
            auto [twin_folds, to_remove] = foldingTwins();
            if (write_progress_on_the_fly) DEBUG(total_twin_folds_done);
            total_twin_folds_done += twin_folds.size() + to_remove.size();
            if (write_progress_on_the_fly) DEBUG(total_twin_folds_done);
            addKNR(to_remove);
            Utils::removeNodes(V, revV, to_remove, helper);
            if(!modified) modified = (!twin_folds.empty() || !to_remove.empty());

            { // add to resulting kernelization objects
                if (knr != nullptr) { res.push_back(knr); knr = nullptr; }
                for (auto *x : twin_folds) res.push_back(x);
            }

            TimeMeasurer::stop("Reducer::folding_twins");
            if(modified) continue;
        }

        if(cnf.reducer_use_desk){
            TimeMeasurer::start("Reducer::desk");
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

            TimeMeasurer::stop("Reducer::desk");
            if(modified) continue;
        }

        if(cnf.reducer_use_funnel){
            if(!cnf.reducer_use_domination){
                clog << "CAUTION! Calling funnel reduction without domination rule before!" << endl;
            }
            TimeMeasurer::start("Reducer::funnel");
            if(write_progress_on_the_fly) DEBUG(total_funnels_done);
            auto funnels = funnel();
            total_funnels_done += funnels.size();
            if(write_progress_on_the_fly) DEBUG(total_funnels_done);
            if(!funnels.empty()) modified = true;

            { // add to resulting kernelization objects
                if(knr != nullptr){
                    res.push_back(knr);
                    knr = nullptr;
                }
                for(auto *x : funnels) res.push_back(x);
            }

            TimeMeasurer::stop("Reducer::funnel");
            if(modified) continue;
        }

        if(cnf.reducer_use_domination_6){
            TimeMeasurer::start("Reducer::domination6");
            VI to_remove = domination6();
            addKNR(to_remove);
            total_dominated_nodes6 += to_remove.size();
            Utils::removeNodes(V, revV, to_remove, helper);
            TimeMeasurer::stop("Reducer::domination6");
            if(!modified) modified = (!to_remove.empty());
            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if( cnf.reducer_use_pie && cnf.reducer_use_nonsimple_cycle_arcs ){
            TimeMeasurer::start("Reducer::nonsimple_cycle_arcs");
            if(write_progress_on_the_fly) DEBUG(total_nonsimple_cycle_arcs_removed);
            if(!modified) modified = nonSimpleCycleArc2();
            if(write_progress_on_the_fly) DEBUG(total_nonsimple_cycle_arcs_removed);
            TimeMeasurer::stop("Reducer::nonsimple_cycle_arcs");

            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if( cnf.reducer_use_pie && cnf.reducer_use_mixed_domination ){
            TimeMeasurer::start("Reducer::mixed_domination");

            if (write_progress_on_the_fly) clog << "Running mixed domination" << endl;
            bool changes = mixedDomination();
            if(changes) modified = true;
            if (write_progress_on_the_fly) clog << "Mixed domination applied: " << changes << endl;

            TimeMeasurer::stop("Reducer::mixed_domination");

            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if( cnf.reducer_use_pie && cnf.reducer_use_mixed_domination_full ){
            int old_time_millis = cnf.reducer_mixed_domination_full_max_time_millis_total;
            cnf.reducer_mixed_domination_full_max_time_millis_total /= (1 << (mixed_domination_full_cnt));
            mixed_domination_full_cnt++;

            if (write_progress_on_the_fly) clog << "Running mixed domination full" << endl;
            TimeMeasurer::start("Reducer::mixed_domination_full");
            bool changes = mixedDominationFull();
            if(changes) modified = true;
            if (write_progress_on_the_fly) clog << "Mixed domination full applied: " << changes << endl;
            TimeMeasurer::stop("Reducer::mixed_domination_full");

            cnf.reducer_mixed_domination_full_max_time_millis_total = old_time_millis;

            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_domination_3) {
            TimeMeasurer::start("Reducer::domination3");
            VI dominators = domination3();
            addKNR(dominators);

            if (write_progress_on_the_fly) DEBUG(total_dominated_nodes3);
            total_dominated_nodes3 += dominators.size();
            if (write_progress_on_the_fly) DEBUG(total_dominated_nodes3);
            Utils::removeNodes(V, revV, dominators, helper);
            TimeMeasurer::stop("Reducer::domination3");
            if(!dominators.empty()) modified = true;
            if(modified) continue;
        }

        if(cnf.reducer_use_edge_neighborhood_blocker){
            TimeMeasurer::start("Reducer::edge_neighborhood_blocker");
            VPII arcs_to_add = edgeNeighborhoodBlocker();

            int arcs_before = GraphUtils::countEdges(V,true);
            Utils::addEdges(V, arcs_to_add, helper);
            for( auto& [a,b] : arcs_to_add ) swap(a,b);
            Utils::addEdges(revV, arcs_to_add, helper);
            int arcs_after = GraphUtils::countEdges(V,true);
            int arc_diff = arcs_after - arcs_before;

            total_edge_neighborhood_blocker_edges_added += arc_diff;
            TimeMeasurer::stop("Reducer::edge_neighborhood_blocker");

            if(!modified) modified = (arc_diff > 0);
            if(modified) continue;
        }



        if( Utils::isPIGraph(V, revV, helper) ){
            TimeMeasurer::start("Reducer::LP_relaxation");
            VVI Vcp = V;
            KernelizerVC kern;
            auto [kern_nodes, edges_removed] = kern.lpDecomposition(Vcp);
            addKNR(kern_nodes);
            Utils::removeNodes(V, revV, kern_nodes,helper);
            TimeMeasurer::stop("Reducer::LP_relaxation");
            if(!modified) modified = (!kern_nodes.empty());
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_domination_4){
            int old_time_millis = cnf.reducer_domination_3_4_max_time_millis_total;
            cnf.reducer_domination_3_4_max_time_millis_total /= (1 << (domination4_cnt));
            domination4_cnt++;

            TimeMeasurer::start("Reducer::domination4");
            VI dominators = domination4(cnf.reducer_domination4_max_time_millis_per_node);
            addKNR(dominators);

            if(write_progress_on_the_fly) DEBUG(total_dominated_nodes4);
            total_dominated_nodes4 += dominators.size();
            if(write_progress_on_the_fly) DEBUG(total_dominated_nodes4);
            Utils::removeNodes(V, revV, dominators, helper);
            TimeMeasurer::stop("Reducer::domination4");
            cnf.reducer_domination_3_4_max_time_millis_total = old_time_millis;

            if(!dominators.empty()){
                modified = true;
                continue;
            }
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));


        if(cnf.reducer_use_full_bipartite_blocker){
            TimeMeasurer::start("Reducer::full_bipartite_blocker");
            auto reductions = fullBipartiteBlocker();
            total_full_bipartite_blockers += reductions.size();

            { // add to resulting kernelization objects
                if(knr != nullptr){ res.push_back(knr); knr = nullptr; }
                for(auto *x : reductions) res.push_back(x);
            }

            if(!reductions.empty()) modified = true;
            TimeMeasurer::stop("Reducer::full_bipartite_blocker");

            if(modified) continue;
        }

        if(cnf.reducer_use_reverse_triangle_gadgets){
            TimeMeasurer::start("Reducer::reverse_triangle_gadgets");
            auto [nodes_to_remove, arcs_to_add, rev_tr_gadgets, kern_red_dom6] = reverseTriangleGadget();
            total_reverse_triangle_gadgets_applied += rev_tr_gadgets.size();
            total_reverse_triangle_gadget_dom6_cases += kern_red_dom6.size();

            addKNR(kern_red_dom6);
            { // add to resulting kernelization objects
                if(knr != nullptr){ res.push_back(knr); knr = nullptr; }
                for(auto *x : rev_tr_gadgets) res.push_back(x);
            }

            {
                Utils::removeNodes(V, revV, kern_red_dom6, helper);
                Utils::removeNodes(V, revV, nodes_to_remove, helper);

                Utils::addEdges(V, arcs_to_add, helper);
                for (auto &[a, b] : arcs_to_add) swap(a, b);
                Utils::addEdges(revV, arcs_to_add, helper);
            }

            TimeMeasurer::stop("Reducer::reverse_triangle_gadgets");
            if(!modified) modified = (!rev_tr_gadgets.empty() || !kern_red_dom6.empty());
            if(modified) continue;
        }


        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_cycle_folding){
            TimeMeasurer::start("Reducer::cycle-folding");
            auto cycle_folds = cycleFolding();
            total_cycle_folds_done += cycle_folds.size();
            if(!cycle_folds.empty()) modified = true;

            { // add to resulting kernelization objects
                if(knr != nullptr){ res.push_back(knr); knr = nullptr; }
                for(auto *x : cycle_folds) res.push_back(x);
            }

            TimeMeasurer::stop("Reducer::cycle-folding");
            if(modified) continue;
        }

        if(cnf.reducer_use_general_folding){
            TimeMeasurer::start("Reducer::general_folding");
            auto reductions = generalFolding();
            total_general_folds_done += reductions.size();
            if(!modified) modified = (!reductions.empty());

            { // add to resulting kernelization objects
                if(knr != nullptr){ res.push_back(knr); knr = nullptr; }
                for(auto *x : reductions) res.push_back(x);
            }

            TimeMeasurer::stop("Reducer::general_folding");
            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if( cnf.reducer_use_nonsimple_cycle_arcs_full ){
            int old_time_millis = cnf.reducer_nonsimple_cycle_arcs_full_max_time_millis_total;

            cnf.reducer_nonsimple_cycle_arcs_full_max_time_millis_total /= (1<<(nonsimple_arc_full_cnt));
            nonsimple_arc_full_cnt++;

            TimeMeasurer::start("Reducer::nonsimple_cycle_arcs_full");
            VPII arcs_to_remove = nonSimpleCycleArcFull();

            if(write_progress_on_the_fly) DEBUG(total_nonsimple_cycle_arcs_full_removed);
            total_nonsimple_cycle_arcs_full_removed += arcs_to_remove.size();
            if(write_progress_on_the_fly) DEBUG(total_nonsimple_cycle_arcs_full_removed);

            Utils::removeEdges( V, arcs_to_remove, helper );
            for(auto& [a,b] : arcs_to_remove) swap(a,b);
            Utils::removeEdges( revV, arcs_to_remove, helper );

            TimeMeasurer::stop("Reducer::nonsimple_cycle_arcs_full");
            modified = !arcs_to_remove.empty();

            cnf.reducer_nonsimple_cycle_arcs_full_max_time_millis_total = old_time_millis;
            if(modified) continue;
        }

        if(cnf.reducer_use_domination_5){
            int old_time_millis = cnf.reducer_domination_5_max_time_millis_total;
            cnf.reducer_domination_5_max_time_millis_total /= (1 << (domination5_cnt));
            domination5_cnt++;

            if(write_progress_on_the_fly){
                DEBUG(total_dominated_nodes5);
                DEBUG(total_domination5_pi_arcs_added);
            }

            TimeMeasurer::start("Reducer::domination5");
            auto [dominated, arcs_to_add] = domination5(cnf.reducer_domination_5_max_time_millis_total,
                                                        cnf.reducer_domination5_max_time_millis_per_node);
            if( !dominated.empty() ) {
                addKNR(dominated);

                if (write_progress_on_the_fly) DEBUG(total_dominated_nodes5);
                total_dominated_nodes5 += dominated.size();
                if (write_progress_on_the_fly) DEBUG(total_dominated_nodes5);
                Utils::removeNodes(V, revV, dominated, helper);
            }
            else if( !arcs_to_add.empty() ){
                total_domination5_pi_arcs_added += arcs_to_add.size();
                Utils::addEdges(V, arcs_to_add, helper );
                for(auto & [a,b] : arcs_to_add) swap(a, b);
                Utils::addEdges(revV, arcs_to_add, helper );
            }

            if(write_progress_on_the_fly){
                DEBUG(total_dominated_nodes5);
                DEBUG(total_domination5_pi_arcs_added);
            }

            TimeMeasurer::stop("Reducer::domination5");
            cnf.reducer_domination_5_max_time_millis_total = old_time_millis;

            if(check_correspondings) assert(Utils::isCorresponding(V, revV));

            if( !dominated.empty() || !arcs_to_add.empty() ){
                modified = true;
                continue;
            }
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_twins_merge){
            TimeMeasurer::start("Reducer::twins");
            if(write_progress_on_the_fly) DEBUG(total_twins_merged);

            double millis_done = chrono::duration<double, std::milli >
                    (chrono::steady_clock::now() - reducer_start_time ).count();
            int millis_left = cnf.reducer_max_time_millis - millis_done;

            bool mod = mergeTwins(millis_left);

            if(write_progress_on_the_fly) DEBUG(total_twins_merged);
            TimeMeasurer::stop("Reducer::twins");
            if(mod) modified = true;
            if(mod) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_domination_6inserter){
            TimeMeasurer::start("Reducer::domination6inserter");
            auto [nodes_removed, pi_edges_added] = domination6Inserter();
            addKNR(nodes_removed);
            total_domination6inserter_nodes_removed += nodes_removed.size();
            total_domination6inserter_pi_edges_inserted += pi_edges_added.size();
            TimeMeasurer::stop("Reducer::domination6inserter");
            if(!modified) modified = ( pi_edges_added.size() > 0 || nodes_removed.size() > 0 );
            if(modified) continue;
        }

        if(cnf.reducer_use_bottleneck2){
            TimeMeasurer::start("Reducer::bottleneck2");

            VI to_remove = bottleneck2();
            addKNR(to_remove);
            total_bottleneck2_nodes_removed += to_remove.size();
            Utils::removeNodes(V, revV, to_remove, helper);

            TimeMeasurer::stop("Reducer::bottleneck2");
            if(!modified) modified = !to_remove.empty();
            if(modified) continue;
        }

        if(cnf.reducer_use_bottleneck){
            TimeMeasurer::start("Reducer::bottleneck");

            if(write_progress_on_the_fly){
                DEBUG(total_bottlenecks_applied);
                DEBUG(total_bottleneck_nodes);
            }

            double millis_done = chrono::duration<double, std::milli >
                    (chrono::steady_clock::now() - reducer_start_time ).count();
            int millis_left = cnf.reducer_max_time_millis - millis_done;

            auto btl_nodes = bottleneck(millis_left);

            if( !btl_nodes.empty() ){
                modified = true;
                total_bottleneck_nodes += btl_nodes.size();

                Utils::removeNodes(V, revV, btl_nodes, helper);
                addKNR(btl_nodes);
            }
            if(write_progress_on_the_fly){
                DEBUG(total_bottlenecks_applied);
                DEBUG(total_bottleneck_nodes);
            }

            TimeMeasurer::stop("Reducer::bottleneck");
            if(modified) continue;
        }

        if(cnf.reducer_use_recursive_reducer){
            TimeMeasurer::start("Reducer::recursive_reducer");
            VI removed = recursiveReducer();
            addKNR(removed);
            total_recursive_reducer_nodes_removed += removed.size();
            TimeMeasurer::stop("Reducer::recursive_reducer");

            if(!modified) modified = (!removed.empty());
            if(modified) continue;
        }

        if(check_correspondings) assert(Utils::isCorresponding(V, revV));

        if(cnf.reducer_use_spiderweb_gadgets){
            TimeMeasurer::start("Reducer::spiderweb_gadgets");

            if(V.size() > 50){
                ENDL(1);
                DEBUG(V.size());
            }

            auto nonpi_edges_before = GraphUtils::countEdges(V,true) -
                    2 * Utils::getAllPIEdges(V,revV,helper).size();

            if(write_progress_on_the_fly){
                DEBUG(total_spiderweb_gadgets_applied);
                DEBUG(total_spiderweb_nodes_added);
            }

            auto cycle_gadgets = spiderwebGadgets();

            { // add to resulting kernelization objects
                if(knr != nullptr){
                    res.push_back(knr);
                    knr = nullptr;
                }
                for(auto *x : cycle_gadgets) res.push_back(x);
            }

            if(!cycle_gadgets.empty()){
                if(V.size() > 50){
                    DEBUG(cycle_gadgets.size());
                    DEBUG(V.size());
                    ENDLS(10,"*");
                }

                UniformIntGenerator rnd(0,1'000'000'000ll * 1'000'000'000);
                while( hashes.size() < N ){
                    hashes.push_back(rnd.rand());
                    helper.push_back(false);
                }

                modified = true;
            }

            auto nonpi_edges_after = GraphUtils::countEdges(V,true) -
                    2 * Utils::getAllPIEdges(V,revV,helper).size();
            int nonpi_edges_removed = nonpi_edges_before - nonpi_edges_after;

            total_spiderweb_arcs_removed += nonpi_edges_removed;
            total_spiderweb_gadgets_applied += cycle_gadgets.size();

            if(write_progress_on_the_fly){
                DEBUG(total_spiderweb_gadgets_applied);
                DEBUG(total_spiderweb_nodes_added);
            }

            TimeMeasurer::stop("Reducer::spiderweb_gadgets");
            if(modified) continue;
        }

    }while(modified);

    if(debug){ DEBUG(V);}

    if(knr != nullptr){ res.push_back(knr); knr = nullptr;}
    return res;
}


VPII Reducer::pie() {
    TimeMeasurer::start("Reducer::pie");

    int N = V.size();
    VB helper(N,false);

    VVI H = Utils::getNonPIGraph(V, revV);

    VI in_comp(N,-1);

    {
        StronglyConnectedComponents scc(H);
        scc.createStronglyConnectedComponents();
        swap(in_comp, scc.compParent);
    }

    VPII to_remove;
    for(int i=0; i<N; i++){
        for( int d : H[i] ){
            if( in_comp[i] != in_comp[d] ) to_remove.emplace_back(i,d);
        }
    }

    TimeMeasurer::stop("Reducer::pie");

    return to_remove;
}

VPII Reducer::dome() {
    TimeMeasurer::start("Reducer::dome");
    constexpr bool debug = false;

    const int E = GraphUtils::countEdges(V, true);
    const int SQRT_E = sqrt(E);

    if(debug){ DEBUG(E); DEBUG(SQRT_E); }

    VI large_nodes;
    for( int i=0; i<N; i++ ) if( 1ll * V[i].size() + revV[i].size() >= SQRT_E ) large_nodes.push_back(i);
    VB is_large(N,false);
    for(int d : large_nodes) is_large[d] = true;

    if(debug) clog << "Nodes with degree >= " << SQRT_E << ": " << large_nodes << endl;

    VPII dominated;
    VB helper(N,false);

    VVI nonpiV = Utils::getNonPIGraph(V, revV);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    {
        VB is_nonpi_pred_u(N, false);
        VB is_succ_u(N, false);

        for (int u : large_nodes) {
            for (int d : revnonpiV[u]) is_nonpi_pred_u[d] = true;
            for (int d : V[u]) is_succ_u[d] = true;

            if(debug){  ENDL(1); DEBUG(u); DEBUG(revnonpiV[u]); }

            for (int v : V[u]) { // for each arc (u,v)
                if(debug){
                    DEBUG(v);
                }

                bool check = true;
                for( int d : V[v] ) if(d == u) check = false; // we do not check PI-edges
                if(!check){
                    if(debug) clog << "Do not checking PI-edge (" << u << "," << v << ")" << endl;
                    continue;
                }

                {
                    int cnt_pred = 0;
                    for (int d : revV[v]) if (is_nonpi_pred_u[d]) cnt_pred++;
                    if (cnt_pred == revnonpiV[u].size()) {
                        if(debug) clog << "Predecessor condition holds" << endl;
                        dominated.emplace_back(u, v);
                        continue;
                    }
                }


                {
                    if(debug) DEBUG(nonpiV[v]);

                    int cnt_succ = 0;
                    for (int d : nonpiV[v]) if (is_succ_u[d]) cnt_succ++;
                    if (cnt_succ == nonpiV[v].size()){
                        if(debug) clog << "Successor condition holds" << endl;
                        dominated.emplace_back(u, v);
                    }
                }
            }

            for (int d : revnonpiV[u]) is_nonpi_pred_u[d] = false;
            for (int d : V[u]) is_succ_u[d] = false;
        }

        if(debug) clog << "Dominated arcs after large_nodes: " << dominated << endl << endl << endl;
    }

    {
        VB is_nonpi_succ_v(N, false);
        VB is_pred_v(N, false);

        for (int v = 0; v < N; v++) {
            for (int d :  nonpiV[v]) is_nonpi_succ_v[d] = true;
            for (int d : revV[v]) is_pred_v[d] = true;

            if(debug){ ENDL(1); DEBUG(v); DEBUG( nonpiV[v]); }

            for (int u : revV[v]) {
                if (is_large[u]){
                    if(debug) clog << "u: " << u << " is large predecessor!" << endl;
                    continue;
                }

                if(debug){
                    DEBUG(u);
                }

                bool check = true;
                for( int d : revV[u] ) if( d == v ) check = false;
                if(!check){
                    if(debug) clog << "Do not checking PI-edge (" << u << "," << v << ")" << endl;
                    continue;
                }

                {
                    int cnt_succ = 0;
                    for (int d : V[u]) if (is_nonpi_succ_v[d]) cnt_succ++;
                    if (cnt_succ ==  nonpiV[v].size()) {
                        dominated.emplace_back(u, v);
                        if(debug) clog << "Successor condition holds" << endl;
                        continue;
                    }
                }

                {
                    int cnt_pred = 0;
                    if(debug) DEBUG( revnonpiV[u]);

                    for (int d :  revnonpiV[u]) if( is_pred_v[d] ) cnt_pred++;
                    if(cnt_pred ==  revnonpiV[u].size()){
                        if(debug) clog << "Predecessor condition holds" << endl;
                        dominated.emplace_back(u,v);
                    }
                }
            }

            for (int d : revV[v]) is_pred_v[d] = false;
            for (int d :  nonpiV[v]) is_nonpi_succ_v[d] = false;
        }

        if(debug) clog << "Dominated arcs after large_nodes + small_nodes: " << dominated << endl;
    }

    {
        sort(ALL(dominated));
        dominated.resize(unique(ALL(dominated)) - dominated.begin());
    }

    TimeMeasurer::stop("Reducer::dome");
    return dominated;
}


tuple<VI, VI, VI>
Reducer::enhanceTwins(VI &v, VB &in_L, VB &in_Q, VB &in_R, VB &in_Pm, VB &in_Pp, VB &in_Np, VB &in_Nm, VB &in_Npi,
                      VB &in_N2, VI &neigh_marker) {
    const bool debug = false;

    VI L,Q = v,R,Pm,Pp,P;
    for( int d : Q ) in_Q[d] = true;

    Pp = V[Q[0]];
    Pm = revV[Q[0]];
    P = Pp;
    for( int d : Pp ) in_Np[d] = in_Pp[d] = true;
    for( int d : Pm ){
        in_Nm[d] = in_Pm[d] = true;
        if(in_Pp[d]) in_Npi[d] = true;
        else P.push_back(d);
    }

    {
        VI N2;
        for( int p : P ){
            for( int d : V[p] ){
                if( in_Pm[d] || in_Pp[d] || in_Q[d] ) continue;
                if(!in_N2[d]){
                    N2.push_back(d);
                    in_N2[d] = true;
                }
                neigh_marker[d]++;
            }

            for( int d : revV[p] ){
                if( in_Pm[d] || in_Pp[d] || in_Q[d] ) continue;
                if(!in_N2[d]){
                    N2.push_back(d);
                    in_N2[d] = true;
                }
                neigh_marker[d]++;
            }
        }

        if(debug){
            ENDLS(20,"*");
            DEBUG(Q); DEBUG(P); DEBUG(L); DEBUG(R); DEBUG(Pp); DEBUG(Pm); DEBUG(N2);
            ENDL(1);
            for( int d : Q+P ){
                DEBUG(d); DEBUG(V[d]); DEBUG(revV[d]); ENDL(1);
            }

            if(Q.size() > 1){ // just an assertion
                for (int d : Q) {
                    sort(ALL(V[d]));
                    sort(ALL(revV[d]));
                }
                for(int i=1; i<Q.size(); i++){
                    sort(ALL(V[Q[i]]));
                    sort(ALL(revV[Q[i]]));
                    assert(equal(ALL(V[Q[0]]), ALL(V[Q[i]])));
                    assert(equal(ALL(revV[Q[0]]), ALL(revV[Q[i]])));
                }
            }
        }


        VI candidates;
        for( int d : N2 ){
            int dsize = V[d].size() + revV[d].size();
            assert(dsize >= neigh_marker[d]);

            if( neigh_marker[d] <= dsize + 1 ) candidates.push_back(d);
        }

        if(debug){  DEBUG(candidates); ENDL(1); }

        VPII candidates2;
        VI valid_candidates;

        for( int d : candidates ){
            int dsize = V[d].size() + revV[d].size();
            int max_errors;
            if( neigh_marker[d] == dsize ) max_errors = 0;
            else max_errors = 1; // node d has one neighbor outside P

            int paired_node = -1;

            int errors = 0;
            for( int v : V[d] ){
                if( !in_Pp[v] ){
                    if( in_N2[v] ){
                        errors++;
                        paired_node = v;
                    }
                    else errors = 1e9; // invalid, its one neighbor outside P is not in N2
                }
            }
            for( int v : revV[d] ){
                if( !in_Pm[v] ){
                    if( in_N2[v] ){
                        errors++;
                        paired_node = v;
                    }
                    else errors = 1e9;
                }
            }

            if( errors == 0 ) valid_candidates.push_back(d);
            else if( errors <= max_errors ) candidates2.emplace_back(d,paired_node);
        }

        if(debug){
            DEBUG(candidates2);
            DEBUG(valid_candidates);
            ENDL(1);
        }

        for( auto [v, paired_node] : candidates2 ){
            bool in_candidates2 = false;
            for( auto x : candidates2 ){
                if( x.first == paired_node ){
                    assert( x.second == v );
                    valid_candidates.push_back(v);
                }
            }
        }

        if(debug){  DEBUG(valid_candidates); ENDL(1); }

        Q += valid_candidates;
        for( int d : valid_candidates ) in_Q[d] = true;

        for(int d : N2) in_N2[d] = false; // clearing helper array
        for( int p : P ){
            for( int d : V[p] ) neigh_marker[d] = 0;
            for( int d : revV[p] ) neigh_marker[d] = 0;
        }

        if(debug){
            DEBUG(Q);
        }
    }


    {
        auto getNodesToR = [&](){

            auto cond1 = [&](int u){
                for( int v : V[u] ) if( !in_Np[v] && !in_Q[v] && !in_R[v] ) return false;
                for( int v : revV[u] ) if( !in_Pm[v] ) return false;
                return true;
            };

            auto cond2 = [&](int u){
                for( int v : V[u] ) if( !in_Np[v] && !in_Q[v] ) return false;
                for( int v : revV[u] ) if( !in_Pm[v] && !in_R[v] ) return false;
                return true;
            };

            VI res;
            for( int d : Pm ){
                if( in_Npi[d] ) continue;
                if( cond1(d) || cond2(d) ) res.push_back(d);
            }

            return res;
        };

        auto getNodesToL = [&](){
            auto cond1 = [&](int u){
                for( int v : V[u] ) if( !in_Pp[v] ) return false;
                for( int v : revV[u] ) if( !in_Nm[v] && !in_Q[v] && !in_L[v] ) return false;
                return true;
            };

            auto cond2 = [&](int u){
                for( int v : V[u] ) if( !in_Pp[v] && !in_L[v] ) return false;
                for( int v : revV[u] ) if( !in_Nm[v] || !in_Q[v] ) return false;
                return true;
            };

            VI res;
            for( int d : Pp ){
                if( in_Npi[d] ) continue;
                if( cond1(d) || cond2(d) ) res.push_back(d);
            }

            return res;
        };

        bool moved = true;
        while(moved){
            auto nodes_to_R = getNodesToR();
            auto nodes_to_L = getNodesToL();
            moved = !nodes_to_L.empty() || !nodes_to_R.empty();

            VB& helper = in_N2;

            if(debug){
                DEBUG(P); DEBUG(Pp) DEBUG(L); DEBUG(Pm)DEBUG(R);
                DEBUG(nodes_to_L); DEBUG(nodes_to_R);
            }

            for( int u : nodes_to_R ){
                in_Pm[u] = false;
                R.push_back(u);
                in_R[u] = true;
            }
            StandardUtils::removeFromArrayPreserveOrderInplace(Pm, nodes_to_R, helper);
            StandardUtils::removeFromArrayPreserveOrderInplace(P, nodes_to_R, helper);

            for( int u : nodes_to_L ){
                in_Pp[u] = false;
                L.push_back(u);
                in_L[u] = true;
            }
            StandardUtils::removeFromArrayPreserveOrderInplace(Pp, nodes_to_L, helper);
            StandardUtils::removeFromArrayPreserveOrderInplace(P, nodes_to_L, helper);

            if(debug){  DEBUG(P); DEBUG(Pp); DEBUG(L); DEBUG(Pm) DEBUG(R); ENDL(1); }
        }
    }

    VI A,B;
    VI XY = P + L + R + Q;

    if(debug) DEBUG(XY);

    {
        VB& helper = in_N2;

        for( int d : Pp ) helper[d] = true;
        VI cand;
        VI LRQ = L + R + Q;
        if(debug) DEBUG(LRQ);

        for( int a : LRQ ){
            if(helper[a]) continue; // nodes from P cannot be in A nor B
            int cnt = 0;
            for( int d : V[a] ) if( helper[d] ) cnt++;

            if( cnt == Pp.size() ) cand.push_back(a);
            else B.push_back(a);
        }
        for( int d : Pp ) helper[d] = false;

        if(debug){  DEBUG(cand); DEBUG(A); DEBUG(B); }

        for( int d : Pm ) helper[d] = true;
        for( int a : cand ){
            if(helper[a]) continue; // nodes from P cannot be in A nor B
            int cnt = 0;
            for( int d : revV[a] ) if( helper[d] ) cnt++;
            if( cnt == Pm.size() ) A.push_back(a);
            else B.push_back(a);
        }
        for( int d : Pm ) helper[d] = false;
    }

    if(debug){  DEBUG(A); DEBUG(B); }

    for(int d : XY){ // clearing all helper arrays
        in_L[d] = in_R[d] = in_Q[d] = in_Pp[d] = in_Pm[d] = in_Npi[d] = false;
        in_Nm[d] = in_Np[d] = false;
    }

    {
        VB& in_P = in_Pp;
        for( int d : A+B+P ) in_P[d] = true;

        auto X = A+B;
        assert( !Utils::hasCycle(V, X, in_L, in_Q, in_R) ); // checking if A+B is a DAG
        for( int d : X ) for( int u : V[d] ) if( !in_P[u] ){ DEBUG(d); DEBUG(u); assert(in_P[u]); }
        for( int d : X ) for( int u : revV[d] ) assert( in_P[u] );
        for( int d : A+B+P ) in_P[d] = false;
    }

    return {A,B,P};
}

pair<vector<FoldingTwinReduction*>, VI> Reducer::foldingTwins() {
    const bool debug = false;

    bool modified = false;
    int N = V.size();
    VB helper(N,false);
    VB was(N,false), on_path(N,false);

    VB affected(N,false);

    map<LL, VI> twinsInOut;
    UniformIntGenerator rnd(0, 1'000'000'000ll * 1'000'000'000ll);
    LL rand_mod = rnd.rand();

    VVI piV = Utils::getUnderlyingPIGraph(V);

    for (int i = 0; i < N; i++) {
        if (V[i].empty() || revV[i].empty()) continue;
        if(!Utils::isPiNode(V, revV, piV, i)) continue;

        LL hash = 0;
        for (int d : V[i]) hash ^= hashes[d];
        for (int d : revV[i]) hash ^= (hashes[d] + rand_mod );
        twinsInOut[hash].push_back(i);
    }

    VVI to_apply;

    VI to_remove;
    for( auto & [h,v] : twinsInOut ){
        if( v.size() == 1 ) continue;
        constexpr bool test = false;
        if(test){ // just an assertion
            int s1 = V[v[0]].size();
            int s2 = V[v[1]].size();
            for( int d : v ){
                assert( V[d].size() == s1 );
                assert( revV[d].size() == s2 );
            }
        }

        VI A = v;
        VI NA = piV[v[0]];
        if( NA.size() != A.size()+1 ) continue;

        if( !Utils::hasCycle(V, NA, was, on_path, helper) ) to_apply.push_back(v);
        else{
            bool aff = false;
            for( int d : A ) if(affected[d]) aff = true;
            for( int d : NA ) if(affected[d]) aff = true;
            if(aff) continue;

            assert(VCUtils::isIndependentSet(piV, A));
            if(debug) clog << "Found valid twins to MERGE, A: " << A << ", adding NA: " << NA << endl;
            to_remove += NA;
            for( int d : A ) affected[d] = true;
            for( int d : NA ){
                affected[d] = true;
                for( int u : V[d] ) affected[u] = true;
                for( int u : revV[d] ) affected[u] = true;
            }
        }
    }

    assert( to_remove.size() == set<int>(ALL(to_remove)).size() );

    vector<FoldingTwinReduction*> res;

    constexpr bool use_conditional_folding = true; // set to false to disable conditional twin folding
    if(!use_conditional_folding) to_apply.clear();

    for(auto &A : to_apply){
        bool apply = true;
        for( int d : A ) if( affected[d] ) apply = false;
        VI NA = piV[A[0]];
        for( int d : NA ) if(affected[d]) apply = false;
        if(!apply) continue;

        for( int d : A ) affected[d] = true;
        for( int d : NA ) affected[d] = true;

        int if_node = NA.back();
        VI else_nodes = NA; else_nodes.pop_back();
        VI fold_nodes = A;

        if(debug){
            clog << "Before identifying nodes" << endl;

            DEBUG(if_node); DEBUG(V[if_node]); DEBUG(revV[if_node]);
            ENDL(1);

            DEBUG(fold_nodes);
            for( int a : fold_nodes ){  DEBUG(a); DEBUG(V[a]); DEBUG(revV[a]);ENDL(1); }

            DEBUG(else_nodes);
            for( int a : else_nodes ){ DEBUG(a); DEBUG(V[a]); DEBUG(revV[a]); ENDL(1); }

            ENDL(1);
        }


        res.push_back( new FoldingTwinReduction(if_node, else_nodes, fold_nodes) );
        Utils::removeNodes(V, revV, fold_nodes, helper);

        if(debug){
            clog << "After removing fold nodes" << endl;

            DEBUG(if_node); DEBUG(V[if_node]); DEBUG(revV[if_node]);
            ENDL(1);

            DEBUG(fold_nodes);
            for( int a : fold_nodes ){  DEBUG(a); DEBUG(V[a]); DEBUG(revV[a]); ENDL(1); }

            DEBUG(else_nodes);
            for( int a : else_nodes ){  DEBUG(a); DEBUG(V[a]); DEBUG(revV[a]); ENDL(1); }

            ENDL(1);
        }


        for(int u : else_nodes){
            Utils::contractNodeToNode( V, revV, u, if_node, helper);

            if(debug) {
                clog << "After contracting node " << u << " to " << if_node << ":" << endl;
                DEBUG(if_node); DEBUG(V[if_node]); DEBUG(revV[if_node]);
                ENDL(1);
                DEBUG(u); DEBUG(V[u]); DEBUG(revV[u]);
                ENDL(1); ENDLS(10, "*");
            }
        }

        if(debug){
            clog << "After identifying nodes" << endl;

            DEBUG(if_node); DEBUG(V[if_node]); DEBUG(revV[if_node]);
            ENDL(1);

            DEBUG(fold_nodes);
            for( int a : fold_nodes ){  DEBUG(a); DEBUG(V[a]); DEBUG(revV[a]); ENDL(1); }

            DEBUG(else_nodes);
            for( int a : else_nodes ){  DEBUG(a); DEBUG(V[a]); DEBUG(revV[a]); ENDL(1); }

            ENDL(2); ENDLS(20,"*"); ENDL(2);
        }
    }

    return {res,to_remove};
}

bool Reducer::mergeTwins(int max_milliseconds) {
    Stopwatch sw;
    sw.setLimit("twins_merge", max_milliseconds);
    sw.start("twins_merge");

    bool modified = false;

    int N = V.size();
    VB helper(N,false);

    VB was(N,false);
    VI temp_zb;

    auto findLBSizeOfInducedGraph = [&](VI & A ){
        if(A.empty()) return (unsigned long)0;

        InducedGraph g = GraphInducer::induce(V, A);
        auto gvcp = g.V;

        VVI revgV = GraphUtils::reverseGraph(g.V);

        const bool find_only_vc = (A.size()  > 15);
        if(find_only_vc){
            VI vc = Utils::getLowerBoundByVCOnPIGraph(g.V, revgV);
            return vc.size();
        }else {

            assert(g.V == gvcp);

            DFVSSolverE solver(&g.V, cnf);
            solver.cnf.write_logs = false;
            solver.cnf.disableAllRecursiveReductions();
            solver.cnf.disableAllConditionalReductions();

            solver.cnf.solverh_use_superpi_vc_ub = false;
            solver.cnf.solverh_improvement_iterations = false;
            solver.cnf.vc_improver_milliseconds = 5;
            solver.cnf.solverh_improvement_iterations = 0;

            VI dfvs_e = solver.solveForInputGraph(g.V);
            assert(g.V == gvcp);

            if( dfvs_e.size() >= g.V.size() ){
                DEBUG(g.V); DEBUG(gvcp); DEBUG(dfvs_e);
                assert(dfvs_e.size() < g.V.size());
            }

            return dfvs_e.size();
        }
    };


    {
        map<LL, VI> twinsInOut;
        UniformIntGenerator rnd(0, 1'000'000'000ll * 1'000'000'000ll);
        LL rand_mod = rnd.rand();


        for (int i = 0; i < N; i++) {
            if (V[i].empty() || revV[i].empty()) continue;

            LL hash = 0;
            for (int d : V[i]) hash ^= hashes[d];
            for (int d : revV[i]) hash ^= (hashes[d] + rand_mod );
            twinsInOut[hash].push_back(i);
        }

        VB in_L(N,false), in_R(N,false), in_Q(N,false), in_N2(N,false), in_Pp(N,false),
                in_Pm(N,false), in_Npi(N,false), in_Np(N,false), in_Nm(N,false);
        VI neigh_marker(N,0);

        vector<pair<LL,VI>> twins_pairs(ALL(twinsInOut));
        sort(ALL(twins_pairs), [&](auto & a, auto & b){
           int va = a.second[0];
           int vb = b.second[0];
           return V[va].size() + revV[va].size() < V[vb].size() + revV[vb].size();
        });


        int cnt_progress = 0;
        for( auto & [h,v] : twins_pairs){
            if( V[v[0]].empty() || revV[v[0]].empty() ){
                assert( V[v[0]].size() + revV[v[0]].size() == 0 );
                continue;
            }

            if(sw.tle("twins_merge")) break;

            const int max_size = cnf.reducer_max_twin_merge_neighborhood_size;
            if( (v.size() >= 2 && V[v[0]].size() + revV[v[0]].size() <= max_size + 2*v.size() )
                || V[v[0]].size() + revV[v[0]].size() <= max_size ){

                bool loops = false; // if a node from
                bool is_affected = false;
                for( int d : v ){

                    if( Utils::hasLoop(V,d) ){
                        loops = true;
                        break;
                    };

                    if(d == v[0]) { // we need to check this only once, because v contains twins
                        for (int x : V[d])if (Utils::hasLoop(V, x)) { loops = true;break; }
                        if (loops) break;
                        for (int x : revV[d])if (Utils::hasLoop(V, x)) { loops = true;break; }
                        if (loops) break;
                    }
                }

                if(loops) continue;

                constexpr bool test = false;
                if(test && v.size() > 1){ // just an assertion
                    VI t1 = V[v[0]]; sort(ALL(t1));
                    VI t2 = revV[v[0]]; sort(ALL(t2));
                    for( int d : v ){
                        if(d == v[0]) continue;
                        VI tt1 = V[d]; sort(ALL(tt1));
                        VI tt2 = revV[d]; sort(ALL(tt2));
                        assert(equal(ALL(t1), ALL(tt1)));
                        assert(equal(ALL(t2), ALL(tt2)));
                    }
                }


                auto merge_v = [&](VI v){
                    assert(!v.empty());
                    modified = true;
                    for (int d : v){
                        if(!Utils::hasLoop(V, d)){
                            total_twins_merged++;
                            Utils::merge(V, revV, d, helper);
                        }
                    }
                };

                if(v.size() >= 2){ // checking condition for N+(a)
                    int a = v[0];
                    VI A = V[a];
                    int t = findLBSizeOfInducedGraph(A);
                    if(A.size() <= v.size() + t){
                        merge_v(v);
                        continue;
                    }
                }

                if(v.size() >= 2){
                    int a = v[0];
                    VI A = revV[a];
                    int t = findLBSizeOfInducedGraph(A);
                    if(A.size() <= v.size() + t){
                        merge_v(v);
                        continue;
                    }
                }

                constexpr bool use_twin_enhancement = true;

                if(!use_twin_enhancement){ // original twin merge rule based on joined neighborhood
                    int a = v[0];
                    VI A = revV[a] + V[a];
                    set<int> zb(ALL(A));
                    A = VI(ALL(zb));
                    int t = findLBSizeOfInducedGraph(A);

                    if(A.size() <= v.size() + t){
                        merge_v(v);
                        continue;
                    }

                    VI Api; // intersection of N+(a) and N-(a)
                    for( int d : V[a] ) helper[d] = true;
                    for( int d : revV[a] ) if(helper[d]) Api.push_back(d);
                    for( int d : V[a] ) helper[d] = false;

                    int t2 = findLBSizeOfInducedGraph(Api);

                    if( t + 2*v.size() + 2 + t2 > A.size() + Api.size() ) merge_v(v);
                }
                else{ // upgraded twin merge rule based on enhanced twins
                    int a = v[0];
                    temp_zb.clear();

                    for( int u : V[a] ){
                        for( int v : V[u] ){
                            if(was[v]) continue;
                            was[v] = true;
                            temp_zb.push_back(v);

                            if(Utils::hasLoop(V,v) ) loops = true;
                        }
                        for( int v : revV[u] ){
                            if(was[v]) continue;
                            was[v] = true;
                            temp_zb.push_back(v);

                            if(Utils::hasLoop(V,v) ) loops = true;
                        }
                    }
                    for( int u : revV[a] ){
                        for( int v : V[u] ){
                            if(was[v]) continue;
                            was[v] = true;
                            temp_zb.push_back(v);

                            if(Utils::hasLoop(V,v) ) loops = true;
                        }
                        for( int v : revV[u] ){
                            if(was[v]) continue;
                            was[v] = true;
                            temp_zb.push_back(v);

                            if(Utils::hasLoop(V,v) ) loops = true;
                        }
                    }

                    for(int xx : temp_zb) was[xx] = false; // clearing was

                    if(loops) continue;

                    auto[A, B, P] = enhanceTwins(v, in_L, in_Q, in_R, in_Pm, in_Pp, in_Np, in_Nm, in_Npi,
                                                 in_N2, neigh_marker);

                    if( (V[a].empty() || revV[a].empty()) && ( !B.empty() || !P.empty() ) ){
                        DEBUG(v); DEBUG(a); DEBUG(V[a]); DEBUG(revV[a]);
                        DEBUG(A); DEBUG(B); DEBUG(P);
                        assert( !V[a].empty() && !revV[a].empty() );
                    }


                    a = A[0];
                    VI PP = P+B;
                    int t = findLBSizeOfInducedGraph(PP);

                    if(P.size() <= A.size() + t){
                        merge_v(A+B); // this should be he same as merge_v(A) and then apply pie()
                        continue;
                    }

                    VI Api; // intersection of N+(a) and N-(a). This must be a subset of P, so we can get it this way
                    for( int d : V[a] ) helper[d] = true;
                    for( int d : revV[a] ) if(helper[d]) Api.push_back(d);
                    for( int d : V[a] ) helper[d] = false;

                    int t2 = findLBSizeOfInducedGraph(Api);

                    if( t + 2*A.size() + 2 + t2 > P.size() + Api.size() ){
                        merge_v(A+B); // this should be he same as merge_v(A) and then apply pie()
                    }
                }

            }
        }
    }

    return modified;
}


vector<FoldingReduction*> Reducer::folding() {
    VVI underPI = getUnderlyingPIGraph(V, revV);

    VB affected(N,false);
    VB helper(N,false);

    vector<FoldingReduction*> res;

    VVI& piV = underPI;
    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VVI revUnderPI = underPI;

    for( int i=0; i<V.size(); i++ ){
        if(underPI[i].size() == 2 && V[i].size() == 2 && revV[i].size() == 2 ){

            int a = i;
            int b = V[a][0];
            int c = V[a][1];
            assert(b != c);

            if( affected[a] || affected[b] || affected[c] ) continue;

            {
                bool can_be_applied = true;
                bool has_bc_pi_node = (Utils::isPiNode(V, revV, piV, b)
                                       || Utils::isPiNode(V, revV, piV, c) );

                if( !has_bc_pi_node ){

                    {
                        bool partial_domination =
                                Utils::isPartiallyDominated(V, revV, nonpiV, revnonpiV, b, c, helper) ||
                                Utils::isPartiallyDominated(V, revV, nonpiV, revnonpiV, c, b, helper);
                        if (partial_domination) can_be_applied = true;
                        else can_be_applied = false;
                    }
                }

                if(!can_be_applied) continue;
            }

            affected[a] = affected[b] = affected[c] = true;

            if( cnf.reducer_use_domination ) for( int d : piV[b] ) assert( d != c );

            res.push_back( new FoldingReduction( c,b,a ) );
            Utils::contractNodeToNode(V, revV, b,c, helper);
            Utils::removeNode(V, revV, a, helper);

            {
                Utils::contractNodeToNode(underPI, revUnderPI, b,c, helper);
                Utils::removeNode(underPI, revUnderPI, a, helper);

                Utils::contractNodeToNode(nonpiV, revnonpiV, b,c, helper);
                Utils::removeNode(nonpiV, revnonpiV, a, helper);
            }

            constexpr bool test = false;
            if(test) {
                for (int d : underPI[b]) { // just an assertion
                    if (d == c) {
                        DEBUG(a);DEBUG(b);DEBUG(c);
                        DEBUG(V[a]);DEBUG(revV[a]);DEBUG(V[b]);DEBUG(revV[b]);DEBUG(V[c]);DEBUG(revV[c]);
                        DEBUG(underPI[a]);DEBUG(underPI[b]);DEBUG(underPI[c]);
                        assert(d != c);
                    }
                }
            }
        }
    }

    return res;
}


VI Reducer::domination1() {
    int N = V.size();
    VB was1(N,false);
    VB was2(N,false);
    VB affected(N,false);

    VI perm(N);
    iota(ALL(perm),0);
    sort(ALL(perm),[&](int a, int b){
        return V[a].size() + revV[a].size() > V[b].size() + revV[b].size();
    });

    VVI piV = Utils::getUnderlyingPIGraph(V, revV);

    VI res;

    for( int a : perm ){
        if(affected[a]) continue;
        if( V[a].empty() || revV[a].empty() ) continue;

        was1[a] = was2[a] = true;
        for( int d : V[a] ) was1[d] = true;
        for( int d : revV[a] ) was2[d] = true;

        bool check = true;
        for( int d : V[a] ) if( affected[d] ) check = false;
        for( int d : revV[a] ) if( affected[d] ) check = false;

        if(check) {
            for (int b : piV[a]) {
                if(affected[b]) continue;

                if( V[b].size() > V[a].size() || revV[b].size() > revV[a].size() ) continue;

                bool dominated = true;

                for (int d : V[b]) {
                    if (!was1[d]) {
                        dominated = false;
                        break;
                    }
                }

                if (!dominated) continue;

                for (int d : revV[b]) {
                    if (!was2[d]) {
                        dominated = false;
                        break;
                    }
                }

                if (dominated) {
                    res.push_back(a);
                    for (int d : V[a]) affected[d] = true;
                    for (int d : revV[a]) affected[d] = true;
                    break;
                }
            }
        }

        was1[a] = was2[a] = false;
        for( int d : V[a] ) was1[d] = false;
        for( int d : revV[a] ) was2[d] = false;
    }

    return res;
}

VI Reducer::domination2() {
    int N = V.size();
    VB was1(N,false);
    VB was2(N,false);
    VB affected(N,false);

    VI perm(N);
    iota(ALL(perm),0);
    sort(ALL(perm),[&](int a, int b){
        return V[a].size() + revV[a].size() > V[b].size() + revV[b].size();
    });

    VVI piV = Utils::getUnderlyingPIGraph(V, revV);

    VI res;

    for( int a : perm ) {
        if (affected[a]) continue;
        if (V[a].empty() || revV[a].empty()) continue;

        was1[a] = true;
        for( int d : piV[a] ) was1[d] = true;

        bool check = true;
        for( int d : V[a] ) if( affected[d] ) check = false;
        for( int d : revV[a] ) if( affected[d] ) check = false;

        if(check) {
            for (int b : piV[a]) {
                if( V[b].size() > piV[a].size() && revV[b].size() > piV[a].size() ) continue;

                bool dominated = true;

                if( V[b].size() <= piV[a].size()) {
                    for (int d : V[b]) {
                        if (!was1[d]) {
                            dominated = false;
                            break;
                        }
                    }
                }else dominated = false;

                if( !dominated && revV[b].size() <= piV[a].size() ){
                    dominated = true;

                    for( int d : revV[b] ){
                        if( !was1[d] ){
                            dominated = false;
                            break;
                        }
                    }
                }

                if (dominated) {
                    res.push_back(a);
                    for (int d : V[a]) affected[d] = true;
                    for (int d : revV[a]) affected[d] = true;
                    break;
                }
            }
        }

        was1[a] = false;
        for( int d : piV[a] ) was1[d] = false;
    }

    return res;
}

VI Reducer::domination3(bool search_for_simple_cycle, int max_time_millis_per_node) {
    int N = V.size();

    VVI piV = Utils::getUnderlyingPIGraph(V, revV);
    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB is_source(N,false), is_end(N,false);
    VB in_V(N,true);
    VI marker(N,0);

    VI res;
    VB affected(N,false);

    auto start_domination4_total = chrono::steady_clock::now();

    VI perm = CombinatoricUtils::getRandomPermutation(N);

    if(cnf.write_logs && N >= 50){
        clog << "\rdomination " << (search_for_simple_cycle ? "4" : "3") << "                    " << flush;
    }

    for( int a : perm ){
        if(affected[a]) continue;

        bool can = true;
        for(int d : V[a]) if(affected[d]) can = false;
        for(int d : revV[a]) if(affected[d]) can = false;
        if(!can) continue;

        in_V[a] = false;
        for( int d : piV[a] ) in_V[d] = false;

        for( int b : piV[a] ){
            if( V[b].size() == piV[b].size() || revV[b].size() == piV[b].size() ) continue;

            bool npi_contained = true;
            for( int d : piV[b] ) if( in_V[d] ) npi_contained = false;
            if(!npi_contained) continue;

            VI Np; // N+(b) \ Npi(a)
            for( int d : V[b] ) if(in_V[d]) Np.push_back(d);
            VI Nm; // N-(b) \ Npi(a)
            for(int d : revV[b]) if(in_V[d]) Nm.push_back(d);

            bool has_path;
            if(!search_for_simple_cycle){
                has_path = Utils::hasPathFromTo(nonpiV, revnonpiV, in_V, Np, Nm, is_source, is_end );
            }else {
                assert(in_V[b] == false);
                in_V[b] = true;
                bool has_simple_cycle = Utils::containsSimpleCycleWithNode(V, revV, nonpiV, revnonpiV,
                                     in_V, marker, is_end, b, max_time_millis_per_node,
                                     cnf.reducer_simple_cycle_max_branch_depth);
                in_V[b] = false;

                constexpr bool test = false;
                if(test){
                    assert(marker == VI(N,0)); assert(is_end == VB(N,false));

                    bool has_just_path = Utils::hasPathFromTo(nonpiV, revnonpiV, in_V, Np, Nm, is_source, is_end);
                    if (V.size() > 30 && has_just_path != has_simple_cycle) {
                        ENDL(1);DEBUG(V.size());DEBUG(has_just_path);DEBUG(has_simple_cycle);ENDL(2);
                    }
                }

                has_path = has_simple_cycle;
            }

            if(!has_path){
                res.push_back(a);
                affected[a] = true;
                for(int d : V[a]) affected[d] = true;
                for(int d : revV[a]) affected[d] = true;
                break;
            }

            if(chrono::duration<double, std::milli >(chrono::steady_clock::now() - start_domination4_total ).count()
                > cnf.reducer_domination_3_4_max_time_millis_total) break;
        }

        in_V[a] = true;
        for( int d : piV[a] ) in_V[d] = true;

        if(chrono::duration<double, std::milli >(chrono::steady_clock::now() - start_domination4_total ).count()
           > cnf.reducer_domination_3_4_max_time_millis_total) break;
    }

    assert( set<int>(ALL(res)).size() == res.size() );

    return res;
}

VI Reducer::domination4(int max_time_millis_per_node) {
    VI res = domination3( true, max_time_millis_per_node );
    return res;
}

void Reducer::writeTotals() {
    DEBUG(total_pie_edges_removed);
    DEBUG(total_core_nodes_removed);
    DEBUG(total_dome_edges_removed);
    DEBUG(total_dominated_nodes1);
    DEBUG(total_dominated_nodes2);
    DEBUG(total_dominated_nodes3);
    DEBUG(total_dominated_nodes4);
    DEBUG(total_dominated_nodes5);
    DEBUG(total_domination5_pi_arcs_added);
    DEBUG(total_dominated_nodes6);
    DEBUG(total_domination6inserter_nodes_removed);
    DEBUG(total_domination6inserter_pi_edges_inserted);
    DEBUG(total_reverse_triangle_gadgets_applied);
    DEBUG(total_reverse_triangle_gadget_dom6_cases);
    DEBUG(mixed_domination_nodes_excluded);
    DEBUG(mixed_domination_nodes_full_excluded);
    DEBUG(total_inoutclique_nodes_merged);
    DEBUG(total_folds_done);
    DEBUG(total_general_folds_done);
    DEBUG(total_twin_folds_done);
    DEBUG(total_full_bipartite_blockers);
    DEBUG(total_edge_neighborhood_blocker_edges_added);
    DEBUG(total_desk_folds);
    DEBUG(total_desk_dominations);
    DEBUG(total_unconfined_nodes);
    DEBUG(total_desk_arcs_added);
    DEBUG(total_funnels_done);
    DEBUG(total_cycle_folds_done);
    DEBUG(total_twins_merged);
    DEBUG(total_nonsimple_cycle_arcs_removed);
    DEBUG(total_nonsimple_cycle_arcs_full_removed);
    DEBUG(total_spiderweb_gadgets_applied);
    DEBUG(total_spiderweb_nodes_added);
    DEBUG(total_spiderweb_edges_added);
    DEBUG(total_spiderweb_fill_edges_added);
    DEBUG(total_spiderweb_arcs_removed);
    DEBUG(total_bottleneck_nodes);
    DEBUG(total_bottlenecks_applied);
    DEBUG(total_bottleneck2_nodes_removed);
    DEBUG(total_recursive_reducer_nodes_removed);
}

VPII Reducer::nonSimpleCycleArcFull() {
    const bool debug = false;

    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VPII arcs_to_remove;
    int N = V.size();

    VI marker(N,0);
    VB is_end(N,false);
    VB helper(N,false);

    int beg_node = -1;

    const int MAX_TIME_MILLIS = cnf.reducer_nonsimple_cycle_arcs_full_max_time_millis_per_arc;
    int depth = 0;

    vector< set<int> > arcs_in_simple_cycles(N);
    {
        VVI piV = Utils::getUnderlyingPIGraph(V);
        for (int i = 0; i < N; i++) arcs_in_simple_cycles[i].insert(ALL(piV[i]));

        const bool add_triangle_arcs = true;
        if(add_triangle_arcs){

            TimeMeasurer::start("Reducer::arc_full_add_triangles");
            VI order(N);
            iota(ALL(order),0);
            sort(ALL(order), [&]( int a, int b ){ return nonpiV[a].size() > nonpiV[b].size(); } );
            VI inOrder(N);
            for(int i=0; i<N; i++) inOrder[order[i]] = i;

            VB was(N,false);

            for(int a : order){
                for( int d : revnonpiV[a] ) was[d] = true;
                bool found = false;
                for( int b : nonpiV[a] ){
                    if( inOrder[b] < inOrder[a] ) continue;
                    for(int c : nonpiV[b]) {
                        if (was[c] && inOrder[c] > inOrder[a]) {
                            arcs_in_simple_cycles[a].insert(b);
                            arcs_in_simple_cycles[b].insert(c);
                            arcs_in_simple_cycles[c].insert(a);
                            found = true;
                            break;
                        }
                    }
                }
                for( int d : revnonpiV[a] ) was[d] = false;
            }
            TimeMeasurer::stop("Reducer::arc_full_add_triangles");
        }
    }

    auto total_start_time = chrono::steady_clock::now();
    auto arc_start_time = chrono::steady_clock::now();

    auto getMillisFromStart = [&](){
        return chrono::duration<double, std::milli >(chrono::steady_clock::now() - arc_start_time ).count();
    };

    bool branch_depth_exceeded = false;

    function<bool(int,int)> is_in_simple_cycle = [&](int num, int par){
        if(depth == 0) assert(!branch_depth_exceeded);
        if(branch_depth_exceeded) return true;

        for( int d : V[num] ) marker[d]++;
        for( int d : revV[num] ) marker[d]++;
        if( nonpiV[num].size() >= 2 ) depth++;

        bool found_cycle = false;
        for( int d : nonpiV[num] ){
            if( marker[d] == 1 && is_end[d] ){
                found_cycle = true;
                arcs_in_simple_cycles[num].insert(d);
                arcs_in_simple_cycles[d].insert(beg_node);
                break;
            }
        }

        if(depth > cnf.reducer_simple_cycle_max_branch_depth){
            found_cycle = true;
            branch_depth_exceeded = true;
        }

        if( !found_cycle ){
            for( int d : nonpiV[num] ){
                if( marker[d] == 1 ){
                    if((depth % 4 == 0) && getMillisFromStart() > MAX_TIME_MILLIS ){
                        found_cycle = true;
                        break;
                    }

                    if( is_in_simple_cycle( d, num ) ){
                        found_cycle = true; // cycle might not have been found, but search was terminated
                        if(getMillisFromStart() < MAX_TIME_MILLIS && !branch_depth_exceeded){
                            arcs_in_simple_cycles[num].insert(d);
                        }
                        break;
                    }
                }
            }
        }

        for( int d : V[num] ) marker[d]--;
        for( int d : revV[num] ) marker[d]--;
        if( nonpiV[num].size() >= 2 ) depth--;

        if(found_cycle) return true;
        return false;
    };

    int progress_cnt = 0;

    if(cnf.write_logs && N > 50){
        clog << "\rNonsimple cycle arcs full                                     " << flush;
    }

    VI perm = CombinatoricUtils::getRandomPermutation(N);
    for( int a : perm ){

        if( nonpiV[a].size() < 2 ) continue;

        for( int d : revnonpiV[a] ) is_end[d] = true;
        for( int d : V[a] ) marker[d]++; // for the start node i we mark only out-neighbors

        VPII temp;
        beg_node = a;
        for( int b : nonpiV[a] ){
            if(arcs_in_simple_cycles[a].count(b)) continue;

            arc_start_time = chrono::steady_clock::now();
            assert(depth == 0);
            branch_depth_exceeded = false;

            if( !is_in_simple_cycle(b, a ) ){
                temp.emplace_back(a,b);
            }

            auto time_total = chrono::duration<double, std::milli >(chrono::steady_clock::now() - total_start_time ).count();
            if(time_total > cnf.reducer_nonsimple_cycle_arcs_full_max_time_millis_total) break;
        }

        if(debug){  DEBUG(a); DEBUG(temp); ENDL(1); }

        arcs_to_remove += temp;

        for( int d : revnonpiV[a] ) is_end[d] = false;
        for( int d : V[a] ) marker[d]--;


        Utils::removeEdges( nonpiV, temp, helper );
        for(auto & [x,y] : temp) swap(x,y);
        Utils::removeEdges( revnonpiV, temp, helper );

        auto time_total = chrono::duration<double, std::milli >(chrono::steady_clock::now() - total_start_time ).count();
        if(time_total > cnf.reducer_nonsimple_cycle_arcs_full_max_time_millis_total) break;
    }

    return arcs_to_remove;
}

bool Reducer::nonSimpleCycleArc2() {
    const bool debug = false;

    VVI G = Utils::getNonPIGraph(V);
    VVI piV = Utils::getUnderlyingPIGraph(V);

    VVI revG = GraphUtils::reverseGraph(G);
    int N = G.size();

    VB was(N,false);
    VB was2(N,false);
    VB helper(N,false);
    VPII res;

    if(debug){  DEBUG(G); DEBUG(piV); }

    assert(GraphUtils::isSimple(V)); assert(GraphUtils::isSimple(G));

    bool modified = false;

    auto checkForPaths = [&]( VI & pthB2A2, VI & pthA2A1, VI & pthA1B1 ){

        bool is_in_simple_cycle = true;

        if( !pthA2A1.empty() && !pthA1B1.empty() ) {
            assert( pthA2A1.size() > 0 && pthA1B1.size() > 0 );

            for( int d : pthA1B1 ) was2[d] = true;
            int a2 = pthB2A2.back();
            int a1 = pthA1B1[0];

            for (int u : pthB2A2) {
                for (int v : G[u]) {
                    if (was2[v]) {
                        if (u != a2 || v != a1) {
                            is_in_simple_cycle = false;
                            break;
                        }
                    }
                }

                if (!is_in_simple_cycle) break;
            }
            for (int d : pthA1B1) was2[d] = false;
        }

        {
            VI pth_union = pthB2A2 + pthA2A1 + pthA1B1;
            for( int d : pth_union ) was2[d] = true;

            for( int u : pth_union ){
                for( int v : piV[u] ){
                    if( was2[v] ){
                        is_in_simple_cycle = false;
                        break;
                    }
                }
                if(!is_in_simple_cycle) break;
            }

            for( int d : pth_union ) was2[d] = false;
        }

        return is_in_simple_cycle;
    };

    for( int a2 = 0; a2<N; a2++ ){
       if( G[a2].size() < 2 ) continue;

       for( int d : G[a2] ){
           int a1 = d;
           VI pthA2A1 = {a2,a1};
           while( G[a1].size() == 1 && revG[a1].size() == 1 ){
               a1 = G[a1][0];
               pthA2A1.push_back(a1);
           }

           int b1 = a1;
           VI pthA1B1 = {b1};
           while( G[b1].size() == 1 ){
               b1 = G[b1][0];
               pthA1B1.push_back(b1);
           }

           int b2 = a2;
           VI pthB2A2 = {b2};
           while( b2 != b1 && revG[b2].size() == 1 ){
               b2 = revG[b2][0];
               pthB2A2.push_back(b2);
           }
           reverse(ALL(pthB2A2));

           if(debug){  DEBUG(pthB2A2); DEBUG(pthA2A1); DEBUG(pthA1B1); ENDL(2); }

           if( checkForPaths( pthB2A2, pthA2A1, pthA1B1 ) == false ){
               for( int i=1; i<pthA2A1.size(); i++ ){
                   res.emplace_back( pthA2A1[i-1], pthA2A1[i] );
               }
           }

           was[a2] = true;
           for( int i=0; i+1<pthA2A1.size(); i++ ) was[pthA2A1[i]] = true;
           if( G[a1].size() == 1 ){
               for( int i=0; i+1<pthA1B1.size(); i++ ) was[ pthA1B1[i] ] = true;
           }else assert( pthA1B1.size() == 1 );
       }
    }

    for( int i=0; i<N; i++ ){
        if( was[i] || G[i].empty() ) continue;

        if(G[i].size() != 1){
            DEBUG(i); DEBUG(G[i]); DEBUG(revG[i]);
            assert( G[i].size() == 1 );
        }

        if(revG[i].size() != 1){
            DEBUG(i); DEBUG(G[i]); DEBUG(revG[i]);
            assert( revG[i].size() == 1 );
        }

        int a = G[i][0];
        VI pthB2A2 = {a};
        while( G[a].size() == 1 && revG[a].size() == 1 && a != i ){
            a = G[a][0];
            pthB2A2.push_back(a);
        }

        assert(a == i);

        VI pthA2A1, pthA1B1;

        if( checkForPaths( pthB2A2, pthA2A1, pthA1B1 ) == false ){
            for( int i=1; i<pthB2A2.size(); i++ ){
                res.emplace_back( pthB2A2[i-1], pthB2A2[i] );
            }
        }

        for(int d : pthB2A2) was[d] = true;
    }

    if(debug)
        if( !res.empty() )
            clog << "There were " << res.size() << " arcs removed that belonged to nonsimple cycles (v.2)" << endl;

    total_nonsimple_cycle_arcs_removed += res.size();

    assert( set<PII>(ALL(res)).size() == res.size() );

    Utils::removeEdges(V, res, helper);
    for(auto & [a,b] : res) swap(a,b);
    Utils::removeEdges(revV, res, helper);

    return !res.empty();
}


void Reducer::disableAllNonbasicReductions() {
    cnf.disableAllNonbasicReductions();
}

void Reducer::disableAllConditionalReductions() {
    cnf.disableAllConditionalReductions();
}

pair<VI, VPII> Reducer::domination5(int max_time_millis_total, int max_time_millis_per_node) {
    const bool debug = false;
    if(debug) clog << "REDUCER DOMINATION 5!" << endl;

    int N = V.size();

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB is_source(N,false), is_end(N,false);
    VB in_V(N,true);
    VI marker(N,0);
    VB in_Np(N,false);

    VB affected(N,false);

    int cnt = 0;

    auto start_domination5_total = chrono::steady_clock::now();
    VI perm = CombinatoricUtils::getRandomPermutation(N);

    VVI noninclusionGraph(N);

    VI bs; bs.reserve(N);

    if(cnf.write_logs && V.size() >= 50){
        clog << "\rdomination5                                    " << flush;
    }

    for( int a : perm ){
        if(V[a].empty()) continue;
        if(affected[a]) continue;


        bool can = true;
        for(int d : V[a]) if(affected[d]) can = false;
        for(int d : revV[a]) if(affected[d]) can = false;
        if(!can) continue;


        in_V[a] = false;
        for( int d : piV[a] ) in_V[d] = false;

        bs.clear();
        bs = nonpiV[a] + revnonpiV[a] + piV[a]; // checking only neighbors

        auto millis_left = max_time_millis_total - chrono::duration<double, std::milli >
                (chrono::steady_clock::now() - start_domination5_total ).count();

        VI noninclusion_nodes = Utils::getAllNoninclusionNodesForNode(V, revV, nonpiV, revnonpiV,
                         in_V, marker, is_end, in_Np, a, bs, millis_left, max_time_millis_per_node,
                                           cnf.reducer_simple_cycle_max_branch_depth);

        if( !noninclusion_nodes.empty() && debug ){
            DEBUG(a); DEBUG(V[a]); DEBUG(revV[a]);
            ENDL(1);

            DEBUG(noninclusion_nodes);
            for(int b : noninclusion_nodes){  DEBUG(b); DEBUG(V[b]); DEBUG(revV[b]);ENDL(1); }

            ENDLS(15,"*"); ENDL(5);
        }

        noninclusionGraph[a] = noninclusion_nodes;

        in_V[a] = true;
        for( int d : piV[a] ) in_V[d] = true;

        if(chrono::duration<double, std::milli >(chrono::steady_clock::now() - start_domination5_total ).count()
           > max_time_millis_total) break;
    }

    if(debug) Utils::writeRemainingGraph(noninclusionGraph);

    VVI &G = noninclusionGraph;
    VVI revG = GraphUtils::reverseGraph(G);

    VI nodes_to_remove;
    VPII arcs_to_add;

    VB helper1(N,false), helper2(N,false), helper3(N,false);

    {
        VB was(N,false);
        VB was2(N,false);

        for( int x=0; x<N; x++ ){
            VI R = {x};
            for( int i=0; i<R.size(); i++ ){
                int u = R[i];
                was2[u] = true;
                for( int v : G[u] ){
                    if(!was2[v] && !was[v]){
                        was2[v] = true;
                        R.push_back(v);
                    }
                }
            }
            for(int d : R) was2[d] = false;

            bool has_cycle = Utils::hasCycle( V, R, helper1, helper2, helper3 );
            if(has_cycle){
                nodes_to_remove.push_back(x);
                was[x] = true;

                if(debug){
                    DEBUG(x); DEBUG(R);

                    for( int r : R ){
                        DEBUG(r); DEBUG(G[r]); DEBUG(V[r]); DEBUG(revV[r]);
                        ENDL(1);
                    }


                    ENDLS(5, "DOM5 SUCCESS"); ENDL(5);
                }
            }
        }
    }

    constexpr bool add_filling = true;
    if(add_filling) {

        if (nodes_to_remove.empty()) {

            VVI revG = GraphUtils::reverseGraph(G);

            VB was(N, false);
            for (int i = 0; i < N; i++) for (int d : G[i]) was[d] = true;

            VPII pi_edges = Utils::getAllPIEdges(V, revV, helper1);

            for (auto[u, v] : pi_edges) {
                if (was[u] && was[v]) {
                    VPII to_add = StandardUtils::product(revG[u], revG[v]);
                    for (auto[t1, t2] : to_add)
                        assert(t1 != t2);
                    arcs_to_add += to_add;

                    if (debug) DEBUG(to_add);
                }
            }
        }

        {
            int P = arcs_to_add.size();
            for (int i = 0; i < P; i++) arcs_to_add.emplace_back(arcs_to_add[i].second, arcs_to_add[i].first);
            sort(ALL(arcs_to_add));
            arcs_to_add.resize(unique(ALL(arcs_to_add)) - arcs_to_add.begin());

            if (debug) {  DEBUG(arcs_to_add); DEBUG(arcs_to_add.size()); }

            VPII nonpresent_arcs;
            vector<set<int>> zb(N);
            for (int i = 0; i < N; i++) zb[i].insert(ALL(V[i]));
            for (auto[a, b] : arcs_to_add) if (!zb[a].count(b)) nonpresent_arcs.emplace_back(a, b);
            swap(arcs_to_add, nonpresent_arcs);

            if (debug) {  DEBUG(arcs_to_add); DEBUG(arcs_to_add.size()); }
        }
    }

    if(debug) DEBUG(nodes_to_remove);

    return {nodes_to_remove, arcs_to_add};
}

bool Reducer::mixedDomination() {
    const bool debug = false;

    VPII to_remove, to_add;

    VVI piV = Utils::getUnderlyingPIGraph(V, revV);
    VVI nonpiV = Utils::getNonPIGraph(V, revV);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    int N = V.size();
    VB helper(N,false), was(N,false), was2(N,false);

    bool changes_done = false;

    if(cnf.write_logs && V.size() >= 50){
        clog << "\rmixedDomination                                                " << flush;
    }

    for( int u=0; u<N; u++ ){

        bool is_dominated = false;
        int dominator = -1;

        for( int d : piV[u] ) was[d] = true;

        constexpr bool use_full_domination = true;

        if(!use_full_domination) {
            if (nonpiV[u].size() == 1) { // we can check all v's accessible from u
                int v = nonpiV[u][0];
                do {
                    if (v == u) break;
                    if (piV[v].size() >= piV[u].size()) {
                        int cnt = 0;
                        for (int d : piV[v]) if (was[d]) cnt++;
                        if (cnt == piV[u].size()) {
                            is_dominated = true;
                            dominator = v;
                            break;
                        }
                    }

                    if (nonpiV[v].size() == 1) {
                        assert(v != nonpiV[v][0]);
                        v = nonpiV[v][0];
                    }
                } while (v != u && nonpiV[v].size() == 1);
            }

            if (!is_dominated && revnonpiV[u].size() == 1) {
                int v = revnonpiV[u][0];

                do {
                    if (v == u) break;
                    if (piV[v].size() >= piV[u].size()) {
                        int cnt = 0;
                        for (int d : piV[v]) if (was[d]) cnt++;
                        if (cnt == piV[u].size()) {
                            is_dominated = true;
                            dominator = v;
                            break;
                        }
                    }

                    if (revnonpiV[v].size() == 1) {
                        assert(v != revnonpiV[v][0]); // no loops should be possible
                        v = revnonpiV[v][0];
                    }
                } while (v != u && revnonpiV[v].size() == 1);
            }
        }else{
            unordered_set<int> zb;
            VI temp;

            int a = u;
            while( revnonpiV[a].size() == 1){
                a = revnonpiV[a][0];

                if(!helper[a]){
                    temp.push_back(a);
                    helper[a] = true;
                }else break;

                if(a == u) break;
                zb.insert(a);
            }

            int b = u;
            while( nonpiV[b].size() == 1 && b != a ){
                b = nonpiV[b][0];

                if(!helper[b]){
                    temp.push_back(b);
                    helper[b] = true;
                }else break;

                if(b == u) break;
                zb.insert(b);
            }

            for(int d : temp) helper[d] = false;

            zb.erase(u);
            VI D(ALL(zb));

            temp.clear();
            for( int u : D ){
                for( int d : piV[u] ){
                    if( was[d] && !helper[d] ){
                        helper[d] = true;
                        temp.push_back(d);
                    }
                }
            }
            for(int d : temp) helper[d] = false;

            if( !D.empty() && temp.size() == piV[u].size() ) is_dominated = true;
        }

        for( int d : piV[u] ) was[d] = false;

        if(is_dominated){

            if( debug ){
                clog << "partially-merging node in mixeD_domination" << endl;
                DEBUG(u); DEBUG(V[u]); DEBUG(revV[u]);
                ENDL(1);
                DEBUG(dominator); DEBUG(V[dominator]); DEBUG(revV[dominator]);
                ENDLS(10,"*");
            }

            mixed_domination_nodes_excluded++;

            for( int x : revnonpiV[u] ) to_remove.emplace_back( x,u );
            for( int y : nonpiV[u] ) to_remove.emplace_back( u,y );

            VPII temp = StandardUtils::product( revnonpiV[u], nonpiV[u] );
            to_add += temp;

            Utils::partialMerge(V, revV, nonpiV, revnonpiV, piV, u, helper);
            changes_done = true;
        }
    }

    return changes_done;
}

vector<FunnelReduction *> Reducer::funnel() {
    const bool debug = false;

    int N = V.size();
    VB affected(N,false);
    VB helper(N,false);
    VB was(N,false);

    vector<FunnelReduction*> res;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VPII pi_edges = Utils::getAllPIEdges(V, revV, helper);

    sort(ALL(pi_edges), [&](auto a, auto b){
        return 1ll * piV[a.first].size() * piV[a.second].size() > 1ll * piV[b.first].size() * piV[b.second].size(); // most edges
    });

    VI A,B;

    auto isClq = [&]( VI & A ){
        for( int d : A ) was[d] = true;

        bool can_be = true;
        for( int a : A ) if( piV[a].size()+1 < A.size() ) can_be = false;

        if(can_be){
            sort(ALL(A), [&]( auto & x, auto & y ){
               return piV[x].size() < piV[y].size();
            });

            for( int a : A ){
                int cnt = 1;
                for( int d : piV[a] ){
                    if(was[d]) cnt++;
                }

                if(cnt != A.size()){
                    assert(cnt < A.size());
                    can_be = false;
                    break;
                }
            }
        }

        for( int d : A ) was[d] = false;

        return can_be;
    };

    auto check = [&](){
        for( auto & [v,u] : pi_edges ){
            if(affected[u] || affected[v]) continue;
            if( !Utils::isPiNode(V, revV, piV, v) ) continue;
            if( !Utils::isPiNode(V, revV, piV, u) ) continue;

            bool is_affected = false;
            for( int d : V[u] ) if( affected[d] ) is_affected = true;
            for( int d : V[v] ) if( affected[d] ) is_affected = true;

            if(is_affected) continue;

            A.clear();
            for( int d : V[v] ) if(d != u) A.push_back(d);

            bool common_intersection = false;
            for( int a : A ) was[a] = true;
            for( int d : piV[u] ) if(was[d]) common_intersection = true;
            for( int a : A ) was[a] = false;

            if(cnf.reducer_use_domination && common_intersection) continue;

            if( isClq(A) ){ // apply funnel reduction
                if(cnf.reducer_use_domination) assert(!common_intersection);

                if(debug){
                    ENDL(1);
                    clog << "Applying funnel reduction to edge (" << v << "," << u << ")" << endl;
                    DEBUG(v); DEBUG(u); DEBUG(piV[v]); DEBUG(piV[u]);
                    ENDL(1);
                }

                VI if_nodes = A;
                int else_node = u;
                int funnel_node = v;
                res.push_back( new FunnelReduction( A, else_node, funnel_node ) );

                affected[u] = affected[v] = true;
                for( int a : A ) affected[a] = true;

                B.clear(); // N(u) \ v
                for(int d : piV[u]) if(d != v) B.push_back(d);
                for( int b : B ) affected[b] = true;

                VPII arcs_to_add = StandardUtils::product(A,B);
                int P = arcs_to_add.size();
                for( int i=0; i<P; i++ ){
                    int a = arcs_to_add[i].first;
                    int b = arcs_to_add[i].second;
                    assert(a != b && a != u && a != v && b != u && b != v);
                    arcs_to_add.emplace_back(b,a);
                }

                assert(arcs_to_add.size() % 2 == 0);

                if(debug){
                    DEBUG(A); DEBUG(B); DEBUG(arcs_to_add);
                    ENDL(3); ENDLS(10,"*");
                }

                Utils::addEdges( V, arcs_to_add, helper );
                Utils::addEdges( revV, arcs_to_add, helper );

                Utils::addEdges( piV, arcs_to_add, helper );

                Utils::removeNode(V, revV, v, helper);
                Utils::removeNode(V, revV, u, helper);

                GraphUtils::removeNodeFromGraph(piV, v);
                GraphUtils::removeNodeFromGraph(piV, u);
            }
        }
    };

    check();

    return res;
}

vector<CycleFoldingReduction*> Reducer::cycleFolding() {
    const bool debug = false;

    int N = V.size();
    vector<CycleFoldingReduction*> res;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB helper(N,false), was(N,false), affected(N,false);

    StronglyConnectedComponents scc(nonpiV);
    scc.createStronglyConnectedComponents();
    auto comps = scc.getComponents();

    auto isCycle = [&]( VI& cmp ){
        for( int d : cmp ) if( nonpiV[d].size() != 1 || revnonpiV[d].size() != 1 ) return false;
        return true;
    };

    for( auto& cmp : comps ){
        if( cmp.size() <= 2 ) continue;
        if( !isCycle(cmp) ) continue;

        // now we have a cycle
        VI C(1, cmp[0]);
        int a = nonpiV[C[0]][0];
        while(a != C[0]){
            C.push_back(a);
            a = nonpiV[a][0];
        }

        assert( C.size() == cmp.size() );


        bool is_affected = false;
        for( int d : C ) if( affected[d] ) is_affected = true;
        if(is_affected) continue;

        {
            bool common_intersection = false;

            for (int a : C) {
                for (int d : piV[a]) {
                    if (was[d]) {
                        common_intersection = true;
                        break;
                    }
                    was[d] = true;
                }
                if (common_intersection) break;
            }

            for (int a : C) for (int d : piV[a]) was[d] = false; // clearing was array

            if (common_intersection) continue;
        }

        {
            if(debug) clog << "In cycle-folding, checking all cliques" << endl;
            bool all_cliques = true;
            for( int a : C ){
                if( !CliqueUtils::isClique( piV, piV[a], helper ) ){
                    all_cliques = false;
                    break;
                }
            }

            if( !all_cliques ) continue;
        }

        if(debug){
            clog << "Applying cycle-folding reduction!!" << endl;
            DEBUG(C);
            for( int a : C ){
                ENDL(1);
                DEBUG(a);
                DEBUG(V[a]);
                DEBUG(revV[a]);
                DEBUG(piV[a]);
                ENDL(3);
                ENDLS(10,"*");
            }
        }

        VI if_not_nodes, else_nodes;

        for( int a : C ){
            int s = piV[a].size();
            if_not_nodes += piV[a];
            for( int i=0; i<s; i++ ) else_nodes.push_back(a);
        }

        res.push_back( new CycleFoldingReduction( if_not_nodes, else_nodes ) );

        VPII pi_edges_to_add;
        VI NC;
        for( int a : C ){
            affected[a] = true;
            for( int d : piV[a] ) affected[d] = true;

            NC += piV[a];
        }

        for( int i=0; i<NC.size(); i++ ){
            for( int j=i+1; j<NC.size(); j++ ){
                pi_edges_to_add.emplace_back( NC[i], NC[j] );
                pi_edges_to_add.emplace_back( NC[j], NC[i] );
            }
        }

        Utils::addEdges(piV, pi_edges_to_add, helper);

        Utils::addEdges(V, pi_edges_to_add, helper);
        Utils::addEdges(revV, pi_edges_to_add, helper);

        Utils::removeNodes(V, revV, C, helper);
        Utils::removeNodes(nonpiV, revnonpiV, C, helper);
        Utils::removeNodes(piV, piV, C, helper);
    }

    assert( GraphUtils::isSimple(V) );
    assert(Utils::isCorresponding(V, revV));

    return res;
}


vector<CycleGadgetReduction *> Reducer::spiderwebGadgets() {
    const bool debug = false;

    vector<CycleGadgetReduction *> res;

    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB helper(N, false);
    VB in_V(N,false);
    VI marker(N,0);
    VB is_end(N,false);

    StronglyConnectedComponents scc(nonpiV);
    scc.createStronglyConnectedComponents();
    auto comps = scc.getComponents();

    sort(ALL(comps), []( auto& v, auto& w ){return v.size() < w.size();});

    auto isCycle = [&](VI &cmp) {
        for (int d : cmp) if (nonpiV[d].size() != 1 || revnonpiV[d].size() != 1) return false;
        return true;
    };

    int first_free_node = N;

    int gadget_arcs_added = 0;
    int cycle_arcs_removed = 0;

    auto addGadgetForCycle = [&](VI C) {
        for (int i = 0; i < C.size(); i++) {
            V.push_back({});
            revV.push_back({});
            helper.push_back(false);
            N++;
        }

        VPII arcs_to_add;
        for (int i = 0; i < C.size(); i++) { // adding connection of gadget to cycle C
            arcs_to_add.emplace_back(C[i], N - 1 - i);
            arcs_to_add.emplace_back(N - 1 - i, C[i]);
        }

        for (int i = 0; i < C.size(); i++) {
            int a = N - 1 - i;
            for (int j = i + 1; j < C.size(); j++) {
                int b = N - 1 - j;
                arcs_to_add.emplace_back(a, b);
                arcs_to_add.emplace_back(b, a);
            }
        }

        VPII arcs_to_remove;
        arcs_to_remove.reserve(C.size());
        for (int i = 0; i < C.size(); i++) arcs_to_remove.emplace_back(C[i], C[(i + 1) % C.size()]);

        Utils::removeEdges(V, arcs_to_remove, helper);
        for (auto&[a, b] : arcs_to_remove) swap(a, b);
        Utils::removeEdges(revV, arcs_to_remove, helper);
        cycle_arcs_removed += arcs_to_remove.size();

        Utils::addEdges(V, arcs_to_add, helper);
        Utils::addEdges(revV, arcs_to_add, helper);
        gadget_arcs_added += arcs_to_add.size();
    };

    for (auto &cmp : comps) {
        assert( V.size() == N && revV.size() == N );
        if (cmp.size() <= 2) continue;
        if(cmp.size() > cnf.reducer_max_component_size_for_spiderweb_gadgets) continue;

        const bool consider_only_cycle_components = false;
        if(consider_only_cycle_components) {
            if (!isCycle(cmp)) continue;

            VI C(1, cmp[0]);
            int a = nonpiV[C[0]][0];
            while (a != C[0]) {
                C.push_back(a);
                a = nonpiV[a][0];
            }

            assert(C.size() == cmp.size());

            total_spiderweb_nodes_added += C.size();
            total_spiderweb_edges_added += 2*C.size() + C.size() * (C.size()-1);

            VI if_all_nodes(C.size());
            iota(ALL(if_all_nodes), first_free_node);
            first_free_node += C.size();

            int else_node = C[0];
            res.push_back(new CycleGadgetReduction(if_all_nodes, else_node));

            addGadgetForCycle(C);
        }
        else{
            const int max_cycle_length =  2*cnf.reducer_max_component_size_for_spiderweb_gadgets;
            VVI cycles = Utils::getAllSimpleCycles2(V, revV, nonpiV, revnonpiV, cmp, in_V, marker,
                                                   is_end, max_cycle_length);

            unordered_set<int> zb;

            for( auto &C : cycles ) {
                if( C.size() == 2 ) continue;

                total_spiderweb_nodes_added += C.size();
                total_spiderweb_edges_added += 2*C.size() + C.size() * (C.size()-1);

                VI if_all_nodes(C.size());
                iota(ALL(if_all_nodes), first_free_node);
                first_free_node += C.size();

                int else_node = C[0]; // any node from cycle C
                res.push_back(new CycleGadgetReduction(if_all_nodes, else_node));

                addGadgetForCycle(C);
                zb.insert(ALL(C));
            }
        }
    }

    if(debug) {  DEBUG(cycle_arcs_removed); DEBUG(gadget_arcs_added); }

    return res;
}

VI Reducer::bottleneck(int max_millis) {
    const bool debug = false;

    Stopwatch sw;
    sw.setLimit("bottleneck", max_millis);
    sw.start("bottleneck");

    VI res;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB helper(N, false), was(N,false);

    StronglyConnectedComponents scc(nonpiV);
    scc.createStronglyConnectedComponents();
    auto comps = scc.getComponents();

    const int MAX_SIZE = 30;


    auto findLBSizeOfInducedGraph = [&](VI & A ){
        if(A.empty()) return (unsigned long)0;

        InducedGraph g = GraphInducer::induce(V, A);
        auto gvcp = g.V;

        DFVSSolverE solver(&g.V, cnf);
        solver.cnf.write_logs = false;
        solver.cnf.disableAllRecursiveReductions();
        solver.cnf.disableAllConditionalReductions();

        solver.cnf.solverh_use_superpi_vc_ub = false;
        solver.cnf.solverh_improvement_iterations = false;
        solver.cnf.vc_improver_milliseconds = 5;
        solver.cnf.solverh_improvement_iterations = 0;

        VI dfvs_e = solver.solveForInputGraph(g.V);
        assert(g.V == gvcp);

        if( dfvs_e.size() >= g.V.size() ){
            DEBUG(g.V); DEBUG(gvcp); DEBUG(dfvs_e);
            assert(dfvs_e.size() < g.V.size());
        }

        return dfvs_e.size();
    };

    VI marker(N,0);
    VB affected(N,false);

    auto solveForX = [&](VI cmp ) {
        if (cmp.size() > MAX_SIZE) return;


        VI X = cmp;
        VI NX = GraphUtils::getNeighborhoodExclusive(V, X, helper)
                + GraphUtils::getNeighborhoodExclusive(revV, X, helper);
        StandardUtils::makeUnique(NX);


        if( NX.empty() ) return;
        if (X.size() + NX.size() > MAX_SIZE) return;

        VI XNX = X + NX;

        const bool expand_xnx = true;
        if (expand_xnx) {
            for (int d : XNX) helper[d] = true;
            for (int v : NX) {
                for (int d : V[v]) marker[d]++;
                for (int d : revV[v]) marker[d]++;
            }

            for (int v : NX) {
                VI temp = V[v] + revV[v];
                for (int d : temp) {
                    if (!helper[d] && (V[d].size() + revV[d].size() == marker[d]) ) {
                        helper[d] = true;
                        XNX.push_back(d);
                    }
                }
            }

            for (int d : XNX) helper[d] = false;
            for (int v : NX) {
                for (int d : V[v]) marker[d]--;
                for (int d : revV[v]) marker[d]--;
            }
        }

        const int MAX_COMP_SIZE = 20;
        const bool expand_xnx_vastly = true;

        if (expand_xnx_vastly) {
            TimeMeasurer::start("Reducer::bottleneck-expand-vastly");
            const bool use_fast_expansion = true;

            if(!use_fast_expansion){ // slow expansion
                for (int d : XNX) helper[d] = true;
                VI nodes;
                for (int i = 0; i < N; i++) if (!helper[i]) nodes.push_back(i);

                VVI superPiV = Utils::getSuperPIGraph(V);
                InducedGraph g = GraphInducer::induce(superPiV, nodes);
                auto g_cmps = ConnectedComponents::getConnectedComponents(g.V);
                for (auto &vv : g_cmps) for (int &d : vv) d = g.nodes[d];

                for (auto &c : g_cmps) {
                    if (c.size() > MAX_COMP_SIZE) continue;
                    bool has_neighbor_in_nx = false;
                    for (int v : c) {
                        VI temp = V[v] + revV[v];
                        for (int d : temp) if (helper[d]) { has_neighbor_in_nx = true; break; }
                        if (has_neighbor_in_nx) break;
                    }

                    if (has_neighbor_in_nx) {
                        XNX += c;
                        X += c;
                        if (expand_xnx) assert(c.size() >= 2);
                    }
                }
                for (int d : XNX) helper[d] = false;
            }
            else {
                for( int d : XNX ) was[d] = helper[d] = true;

                VI visited = XNX;
                VI basic_XNX = XNX;
                for( int d : basic_XNX ){
                    VI neigh = V[d] + revV[d];
                    for( int u : neigh ){
                        if( !was[u] ){
                            VI temp(1,u);

                            was[u] = true;
                            visited.push_back(u);

                            for( int i=0; i<temp.size(); i++ ){
                                if( temp.size() > MAX_COMP_SIZE ) break;

                                int v = temp[i];
                                VI neigh_v = V[v] + revV[v];
                                for( int w : neigh_v ){
                                    if( !was[w] ){
                                        was[w] = true;
                                        visited.push_back(w);
                                        temp.push_back(w);
                                    }
                                }
                            }

                            if( temp.size() <= MAX_COMP_SIZE ){
                                XNX += temp;
                                X += temp;
                            }
                        }
                    }
                }

                for( int d : visited ) was[d] = helper[d] = false;
            }

            TimeMeasurer::stop("Reducer::bottleneck-expand-vastly");
        }

        if (XNX.size() > MAX_SIZE) return;
        for( int d : XNX ) if(affected[d]) return;

        {
            VI X2, NX2;
            for(int d : XNX) helper[d] = true;
            for( int d : XNX ){
                bool in_boundary = false;
                for( int u : V[d] ) if( !helper[u] ) in_boundary = true;
                for( int u : revV[d] ) if( !helper[u] ) in_boundary = true;
                if( in_boundary ) NX2.push_back(d);
                else X2.push_back(d);
            }
            for(int d : XNX) helper[d] = false;

            X = X2;
            NX = NX2;
        }

        if(NX.empty()) return;

        constexpr bool test = false;
        if(test){ // just an assertion
            for (int d : XNX) helper[d] = true;
            for (int v : X) {
                for (int d : V[v]){
                    if(!helper[d]){
                        ENDL(1);
                        DEBUG(v); DEBUG(V[v]);DEBUG(revV[v]);
                        ENDL(1);
                        DEBUG(d); DEBUG(V[d]); DEBUG(revV[d]);
                        DEBUG(X); DEBUG(NX); DEBUG(XNX);
                        assert(helper[d]);
                    }
                }
                for (int d : revV[v]) assert(helper[d]);
            }
            for (int d : XNX) helper[d] = false;
        }

        auto s1 = findLBSizeOfInducedGraph(XNX);
        auto s2 = findLBSizeOfInducedGraph(X);

        if (s1 >= NX.size() + s2) {
            total_bottlenecks_applied++;
            assert( s1 == NX.size() + s2 );

            if (debug) {
                ENDL(3);  clog << "Found bottleneck!" << endl;
                DEBUG(X);  DEBUG(NX);  DEBUG(XNX);  DEBUG(s1);  DEBUG(NX.size());  DEBUG(s2);
                ENDL(3);  ENDLS(20, "*");
            }

            res += NX;
            for(int d : XNX) affected[d] = true;
        }
    };

    for(auto &cmp : comps) if( cmp.size() >= 3 ) solveForX(cmp);

    const bool check_bottleneck_for_each_node = true;
    if(check_bottleneck_for_each_node) {
        if (res.empty()) {
            int cnt = 0;
            while(cnt < N) {
                if(sw.tle("bottleneck")) break;
                solveForX(VI(1, cnt));
                cnt++;
            }
        }
    }

    sort(ALL(res));
    res.resize(unique(ALL(res)) - res.begin());

    return res;
}

bool Reducer::mixedDominationFull() {
    const bool debug = false;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    int N = V.size();
    VB helper(N,false), is_end(N,false), was(N,false);
    VI marker(N,0), cnt_marker(N,0);
    const int max_cycle_length = cnf.reducer_simple_cycle_max_branch_depth;
    const int millis = cnf.reducer_mixed_domination_full_max_time_millis_per_node;

    VI temp;

    bool changes_done = false;

    auto start_time_node = chrono::steady_clock::now();
    auto getMillisFromStart = [&](){
        return chrono::duration<double, std::milli >(chrono::steady_clock::now() - start_time_node ).count();
    };

    for( int v=0; v<N; v++ ){
        if( nonpiV[v].empty() || revnonpiV[v].empty() ) continue;

        if( getMillisFromStart() > cnf.reducer_mixed_domination_full_max_time_millis_total ) break;

        for( int d : piV[v] ) was[d] = true;

        VI D = Utils::getIntersectionOfAllInducedNonpiCyclesWithNode( V, revV, nonpiV, revnonpiV, v, marker, is_end,
                                                                      cnt_marker, max_cycle_length, millis  );

        bool is_dominated = false;

        temp.clear();
        for( int u : D ){
            if( u == v) continue;
            for( int d : piV[u] ){
                if( was[d] && !helper[d] ){
                    helper[d] = true;
                    temp.push_back(d);
                }
            }
        }

        if( !D.empty() && temp.size() == piV[v].size() ){
            is_dominated = true;
            if(debug){
                clog << "In mixedDominationFull, nodes from D: " << D << " dominate node " << v << endl;
                DEBUG(piV[v]);
                for( int d : D ){
                    if(d == v) continue;
                    ENDL(1);
                    DEBUG(d);
                    DEBUG(piV[d]);
                }
                ENDL(5);
            }
        }

        for(int d : temp) helper[d] = false;
        for( int d : piV[v] ) was[d] = false;

        if(is_dominated){

            if( debug ){
                clog << "partially-merging node in mixedDominationFull for node " << v << endl;
                DEBUG(v);DEBUG(V[v]);DEBUG(revV[v]);ENDL(1);ENDLS(10,"*");
            }

            mixed_domination_nodes_full_excluded++;

            Utils::partialMerge(V, revV, nonpiV, revnonpiV, piV, v, helper );
            changes_done = true;
        }
    }

    return changes_done;
}

void Reducer::liftSolution(int N, VI &dfvs, vector<DFVSReduction *> &reductions) {
    VB in_dfvs = StandardUtils::toVB(N, dfvs);

    for( int i = (int)reductions.size()-1; i>=0; i-- ){
        reductions[i]->lift(dfvs, in_dfvs);
    }

    clearReductionObjects(reductions);
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

int Reducer::getReductionsSizeDiff(vector<DFVSReduction *> &reductions) {
    int res = 0;
    for(auto * x : reductions) res += x->sizeDiffUB();
    return res;
}

void Reducer::writeReductions(vector<DFVSReduction *> &reductions) {
    clog << "Reductions: " << endl;
    for(auto * x : reductions) clog << x->toString() << endl;
}

VI Reducer::domination6() {
    constexpr bool debug = false;

    if(debug) clog << "Starting domination6" << endl;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB helper(N,false), affected(N,false), was(N,false);
    VI res;

    VVI Np(N);
    VVI Nm(N);

    VVI in_Np(N);
    VVI in_Nm(N);

    VI in_Np_cnt(N,0);
    VI in_Nm_cnt(N,0);

    VI X;

    for( int w=0; w<N; w++ ){
        if( affected[w] || V[w].empty() || revV[w].empty() ) continue; // admit one pi-edge (self-domination possible)

        bool can_be_added = true;
        for( int d : piV[w] ) if(affected[d]) can_be_added = false;

        VI W = piV[w];
        in_Np[w].clear();
        for( int u : W ){
            in_Np[u].clear();
            for( int d : V[u] ){
                if(affected[d]) can_be_added = false;
                in_Np[d].clear();
            }

            for( int d : revV[u] ){
                if(affected[d]) can_be_added = false;
                in_Nm[d].clear();
            }
        }

        if(!can_be_added) continue;

        for( int i=0; i<W.size(); i++ ){
            int u =  W[i];

            Np[u] = V[u];
            W.push_back(w);
            StandardUtils::removeFromArrayPreserveOrderInplace( Np[u], W, helper );
            W.pop_back();
            for( int d : Np[u] ) in_Np[d].push_back(u);

            Nm[u] = revV[u];
            W.push_back(w);
            StandardUtils::removeFromArrayPreserveOrderInplace( Nm[u], W, helper );
            W.pop_back();
            for( int d : Nm[u] ) in_Nm[d].push_back(u);
        }

        auto writeMe = [&]() {
            DEBUG(W);
            DEBUG(w); writeNeighborhood(V, revV, w);
            ENDL(1);
            DEBUG(piV[w]); DEBUG(W); ENDL(1);
            for(int u : W){
                ENDL(1); DEBUG(u); DEBUG(piV[u]);
                writeNeighborhood(V, revV, u);

                set<int> zb(ALL(V[u]));
                zb.insert(ALL(revV[u]));
                for(int d : zb){  DEBUG(d); writeNeighborhood(V, revV, d); }
            }
            ENDL(5);

            for (int i = 0; i < W.size(); i++) {  int u = W[i]; DEBUG(u); DEBUG(Np[u]); DEBUG(Nm[u]); }
        };

        if(debug) writeMe();

        int valid_uj = -1;
        X.clear();
        for( int j=0; j<W.size(); j++ ){
            int uj = W[j];
            if(debug){ DEBUG(w);DEBUG(uj);}
            X.clear();

            VI zb;
            {
                was[w] = true; for (int d : W) was[d] = true;
                VI vrevv = V[uj] + revV[uj];
                for( int d : vrevv ){
                    if(!was[d]){
                        was[d] = true;
                        zb.push_back(d);
                    }
                }
                was[w] = false; for (int d : W) was[d] = false;
                for(int d : zb) was[d] = false;
            }


            if(debug) DEBUG(zb);

            if( zb.size() > cnf.reducer_domination6_max_neigh_size ) continue;

            VI temp;
            for( int a : zb ){
                temp.clear();
                bool condition_met_for_a = false;

                for( int d : piV[a] ){
                    for( int y : in_Np[d] ){
                        temp.push_back(y);
                        in_Np_cnt[y]++;
                        if( in_Np_cnt[y] == Np[y].size() ){
                            condition_met_for_a = true;
                            if(debug)
                                clog << "\tNp[" << y << "] \\subset piV[" << a << "]" << endl;
                        }
                    }
                    if(condition_met_for_a) break;

                    for( int y : in_Nm[d] ){
                        temp.push_back(y);
                        in_Nm_cnt[y]++;
                        if( in_Nm_cnt[y] == Nm[y].size() ){
                            condition_met_for_a = true;
                            if(debug)
                                clog << "\tNm[" << y << "] \\subset piV[" << a << "]" << endl;
                        }
                    }
                    if(condition_met_for_a) break;
                }

                if(debug){ DEBUG(a); DEBUG(condition_met_for_a); }

                for( int d : temp ) in_Np_cnt[d] = in_Nm_cnt[d] = 0;
                temp.clear();

                if(condition_met_for_a) X.push_back(a);
            }


            if(debug) DEBUG(X);

            if(!X.empty()) {
                for (int d : X) was[d] = true;

                int cnt = 0;
                for (int d : Np[uj]) if (was[d]) cnt++;
                if (cnt == Np[uj].size()){
                    valid_uj = uj;
                    if(debug){
                        DEBUG(valid_uj); DEBUG(Np[uj]);
                    }
                }

                cnt = 0;
                for (int d : Nm[uj]) if (was[d]) cnt++;
                if (cnt == Nm[uj].size()){
                    valid_uj = uj;
                    if(debug){
                        DEBUG(valid_uj); DEBUG(Nm[uj]);
                    }
                }

                for (int d : X) was[d] = false;
            }


            if(valid_uj != -1) break;
        }

        if(debug){ DEBUG(valid_uj); }

        if(valid_uj != -1){
            res.push_back(w);
            if(debug)
                clog << "domination6 applied to node w: " << w << "! valid_uj: " << valid_uj << endl;

            affected[w] = true;
            for(int u : W){
                affected[u] = true;
                for(int d : V[u]) affected[d] = true;
                for(int d : revV[u]) affected[d] = true;
            }
        }

        in_Np[w].clear();
        in_Nm[w].clear();
        for( int u : W ){
            in_Np[u].clear();
            in_Nm[u].clear();
            for( int d : V[u] ){
                in_Np[d].clear(); in_Nm[d].clear();
            }
            for( int d : revV[u] ){
                in_Np[d].clear(); in_Nm[d].clear();
            }
        }
        if(debug){ ENDL(5); ENDLS(30, "*"); ENDL(5); }
    }

    if(debug) clog << "domination6() finished, res: " << res << endl;

    return res;
}


pair<VI,VPII> Reducer::domination6Inserter() {
    constexpr bool debug = false;

    if(debug)
        clog << "Starting domination6inserter" << endl;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI piVcp = piV;

    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB helper(N,false), affected(N,false), was(N,false);
    VPII pi_edges_added;
    VPII nonpi_arcs_transformed_to_pi_edges;

    VVI V0 = V; // this is the initial state of graph V
    VVI revV0 = revV; // this is the initial state of graph V
    VVI piV0 = piV;

    VVI Np(N);
    VVI Nm(N);

    VVI in_Np(N);
    VVI in_Nm(N);

    VI in_Np_cnt(N,0);
    VI in_Nm_cnt(N,0);

    VI X;

    auto findValidUj = [&](VI & W, VI & W0){

        int valid_uj = -1;
        X.clear();
        for (int j = 0; j < W0.size(); j++) {
            int uj = W0[j];
            X.clear();

            VI zb = StandardUtils::setUnion(V[uj], revV[uj], helper);
            zb = StandardUtils::setDifference(zb, W, helper);

            if (debug) {
                ENDL(1);
                DEBUG(uj);
                DEBUG(zb);
            }

            if (zb.size() > cnf.reducer_domination6inserter_max_neigh_size) continue;

            VI temp;
            for (int a : zb) {
                temp.clear();

                bool condition_met_for_a = false;

                for (int d : piV[a]) {
                    for (int y : in_Np[d]) {
                        temp.push_back(y);
                        in_Np_cnt[y]++;
                        if (in_Np_cnt[y] == Np[y].size()) {
                            condition_met_for_a = true;
                            if (debug) clog << "\tNp[" << y << "] \\subset piV[" << a << "]" << endl;
                        }
                    }
                    if (condition_met_for_a) break;

                    for (int y : in_Nm[d]) {
                        temp.push_back(y);
                        in_Nm_cnt[y]++;
                        if (in_Nm_cnt[y] == Nm[y].size()) {
                            condition_met_for_a = true;
                            if (debug) clog << "\tNm[" << y << "] \\subset piV[" << a << "]" << endl;
                        }
                    }
                    if (condition_met_for_a) break;
                }

                if (debug) {
                    DEBUG(a);
                    clog << "condition met for a: " << (condition_met_for_a ? "true" : "false") << endl;
                }

                for (int d : temp) in_Np_cnt[d] = in_Nm_cnt[d] = 0;
                temp.clear();

                if (condition_met_for_a) X.push_back(a);
            }


            if (debug) DEBUG(X);
            if (debug) DEBUG(valid_uj);

            {
                for (int d : X) was[d] = true;

                int cnt = 0;
                for (int d : Np[uj]) if (was[d]) cnt++;
                if (cnt == Np[uj].size()) {
                    valid_uj = uj;
                    if (debug) {
                        DEBUG(valid_uj);
                        DEBUG(Np[uj]);
                    }
                }

                cnt = 0;
                for (int d : Nm[uj]) if (was[d]) cnt++;
                if (cnt == Nm[uj].size()) {
                    valid_uj = uj;
                    if (debug) {
                        DEBUG(valid_uj);
                        DEBUG(Nm[uj]);
                    }
                }

                for (int d : X) was[d] = false;
            }

            if(valid_uj != -1) return valid_uj;
        }

        assert(valid_uj == -1);
        return valid_uj;
    };

    int reps = 1;
    while(reps--) {
        for (int w1 = 0; w1 < N; w1++) {
            if (piV[w1].empty()) continue;

            VVI layers = BFS::getBfsLayers(V, w1);

            VI w2s_to_check;
            for (int i = 1; i < min(cnf.reducer_domination6inserter_distance, (int) layers.size()); i++) w2s_to_check += layers[i];
            unordered_set<int> pi_neighbors(ALL(piV[w1]));

            for (int w2 : w2s_to_check) {
                if (piV[w2].empty()) continue;
                if (pi_neighbors.count(w2)) continue;

                VI W1 = piV[w1];
                VI W2 = piV[w2];

                VI W = W1 + W2;
                StandardUtils::makeUnique(W);

                for (int u : W) {
                    in_Np[u].clear();
                    in_Nm[u].clear();
                    for (int d : V[u]) in_Np[d].clear();
                    for (int d : revV[u]) in_Nm[d].clear();
                }

                VI W0;
                {
                    VI inters = StandardUtils::setIntersection(W1, W2, helper);
                    W0 = StandardUtils::setDifference(W, inters, helper);
                }

                bool has_nonpi_connection = false;
                for (int u : W0) {
                    for (int d : nonpiV[u]) if (d == w1 || d == w2) has_nonpi_connection = true;
                    for (int d : revnonpiV[u]) if (d == w1 || d == w2) has_nonpi_connection = true;
                }
                // if there is a nonpi arc between W0 and node w1 or w2, than we cannot execute the reduction
                if (has_nonpi_connection) continue;

                W += VI({w1, w2}); // #CAUTION from now on w1 nad w2 are in W
                for (int u : W0) {
                    Np[u] = V[u];
                    StandardUtils::removeFromArrayPreserveOrderInplace(Np[u], W, helper);
                    for (int d : Np[u]) in_Np[d].push_back(u);

                    Nm[u] = revV[u];
                    StandardUtils::removeFromArrayPreserveOrderInplace(Nm[u], W, helper);
                    for (int d : Nm[u]) in_Nm[d].push_back(u);
                }

                if (debug) {
                    ENDL(3);
                    ENDLS(50, "*");
                    ENDL(3);
                    DEBUG(w1);
                    DEBUG(w2);
                    DEBUG(W);
                    DEBUG(W1);
                    DEBUG(W2);
                    DEBUG(W0);
                }

                int valid_uj = findValidUj( W, W0 );

                if (debug) DEBUG(valid_uj);


                if (valid_uj != -1) {
                    if (StandardUtils::find(nonpiV[w1], w2)) {
                        nonpi_arcs_transformed_to_pi_edges.emplace_back(w1, w2);
                    } else if (StandardUtils::find(revnonpiV[w1], w2)) {
                        nonpi_arcs_transformed_to_pi_edges.emplace_back(w2, w1);
                    }

                    VPII arcs = {{w1, w2},
                                 {w2, w1}};
                    Utils::addEdges(V, revV, arcs, helper);
                    Utils::addEdges(piV, piVcp, arcs, helper);

                    pi_edges_added.emplace_back(w1, w2);
                    pi_edges_added.emplace_back(w2, w1);

                    if (debug) {
                        clog << "---------------------------> pi-edge (" << w1 << "," << w2 << ") can be added!"
                             << endl;
                    }
                }

                for (int u : W) {
                    in_Np[u].clear();
                    in_Nm[u].clear();
                    for (int d : V[u]) {
                        in_Np[d].clear();
                        in_Nm[d].clear();
                    }
                    for (int d : revV[u]) {
                        in_Np[d].clear();
                        in_Nm[d].clear();
                    }
                }

                if (debug) {
                    ENDL(5);
                    ENDLS(30, "*");
                    ENDL(5);
                }

            }
        }
    }

    VI nodes_removed;

    /**
     * #CAUTION! IT SEEMS THAT THERE IS A BUG SOMEWHERE WITH TRANS_DOMINATIONS!!!
     */
    constexpr bool use_trans_dominations = true; // set to true to enable trans-dominations

    constexpr bool use_trans_domination1 = true;
    constexpr bool use_trans_domination2 = true;

    constexpr bool use_trans_domination6 = false;
    constexpr bool safety_mode_return_after_first_removed_node = true;

    constexpr bool add_pi_edges_if_no_nodes_removed = true;

    int pi_edges_added_permanently = 0;

    if(use_trans_dominations){

        VVI trans_V(N);
        Utils::addEdges(trans_V, pi_edges_added, helper);

        VVI tempV = V0;
        VVI temprevV = revV0;


        /**
         * Removes node a on from V, revV, tempV, temprevV, piV
         */
        auto removeNodeFromAll = [&](int a){
            Utils::removeNode(V, revV, a, helper);
            Utils::removeNode(tempV, temprevV, a, helper);
            Utils::removeNode(piV, piVcp, a, helper);
            GraphUtils::removeNodeFromGraph(trans_V,a);

            Utils::removeNode(nonpiV, revnonpiV, a, helper);
        };

        //****************************************************************** BEG OF TRANS DOMINATION 1

        if(use_trans_domination1 || use_trans_domination2) {
            if (debug) {
                ENDL(2);
                clog << "Starting trans-domination1" << endl;
            }
            fill(ALL(affected), false);

            auto checkTransDomination1 = [&](int a, int b) {
                VI prevVa = tempV[a];
                VI prevrevVa = temprevV[a];

                tempV[a] = V[a];
                temprevV[a] = revV[a];

                bool fully_dominated = Utils::isFullyDominated(tempV, temprevV, b, a, helper, was);

                tempV[a] = prevVa;
                temprevV[a] = prevrevVa;

                return fully_dominated;
            };

            auto checkTransDomination2 = [&](int a, int b) {
                VI prevVa = tempV[a];
                VI prevrevVa = temprevV[a];

                tempV[a] = V[a];
                temprevV[a] = revV[a];

                bool v2_dominated = false;
                {
                    was[a] = true;
                    for (int d : piV[a]) was[d] = true;

                    int cnt = 0;
                    for (int d : tempV[b]) if(was[d]) cnt++;
                    if( cnt == tempV[b].size() ) v2_dominated = true;
                    if(!v2_dominated){
                        cnt = 0;
                        for (int d : temprevV[b]) if(was[d]) cnt++;
                        if( cnt == temprevV[b].size() ) v2_dominated = true;
                    }

                    was[a] = false;
                    for (int d : piV[a]) was[d] = false;
                }

                tempV[a] = prevVa;
                temprevV[a] = prevrevVa;

                return v2_dominated;
            };


            for (int a = 0; a < N; a++) {
                if (affected[a]) continue;
                if (piV[a].empty()) continue;

                int bb = -1;

                bool can_remove_a = false;
                for (int b : piV[a]) { // we consider here all pi-neighbors, including those that are in trans_V
                    if (affected[b]) continue;
                    bool a_trans_dominates_b = false;
                    if(use_trans_domination1) a_trans_dominates_b = checkTransDomination1(a, b);
                    if(!a_trans_dominates_b && use_trans_domination2){
                        a_trans_dominates_b = checkTransDomination2(a, b);
                    }


                    if (a_trans_dominates_b) {
                        can_remove_a = true;
                        bb = b;
                        break;
                    }
                }

                if (can_remove_a) {
                    affected[a] = affected[bb] = true;
                    for (int d : V[a] + revV[a] + V[bb] + revV[bb]) affected[d] = true;

                    nodes_removed.push_back(a);
                    removeNodeFromAll(a);
                }
            }
        }

        //****************************************************************** END OF TRANS DOMINATION 1


        //****************************************************************** BEG OF TRANS DOMINATION 6

        if(use_trans_domination6) {
            if (debug) {
                ENDL(2);
                clog << "Starting trans-domination6" << endl;
            }
            fill(ALL(affected), false);

            /**
            * Checks if there can be applied domination6 rule to some node
            */
            auto checkTransDomination6 = [&](int a, VPII &required_edge_fill) {
                VI W = piV[a];
                W.push_back(a);
                VI W0 = piV[a];

                VI W2;
                for (int d : W) was[d] = true;
                for (int u : W) {
                    for (int d : tempV[u] + temprevV[u]) {
                        if (!was[d]) {
                            was[d] = true;
                            W2.push_back(d);
                        }
                    }
                }
                for (int d : W + W2) was[d] = false;

                if (debug) {
                    DEBUG(W);
                    DEBUG(W0);
                    DEBUG(W2);
                }


                for (int d : W2) was[d] = true;
                VPII edge_mods;
                for (int d : trans_V[a]) { // adding pi-edges incident to node a
                    edge_mods.emplace_back(a, d);
                    edge_mods.emplace_back(d, a);
                }
                for (int u : W2) {
                    for (int d : trans_V[u]) {
                        if (was[d]) { // adding pi-edges (x,y) with both x and y in W2
                            edge_mods.emplace_back(u, d);
                            edge_mods.emplace_back(d, u);
                        }
                    }
                }
                for (int d : W2) was[d] = false;
                StandardUtils::makeUnique(edge_mods);
                assert(edge_mods.size() % 2 == 0); // there should be only pi-edges incident to W2 or a



                map<int, pair<VI, VI>> prevs; // current state of neighborhoods, necessary to bring back applied changes
                for (int u : W + W2) {
                    prevs[u] = {V[u], revV[u]};
                    V[u] = tempV[u];
                    revV[u] = temprevV[u];
                }

                Utils::addEdges(V, revV, edge_mods, helper);


                {
                    for (int u : W + W2) { // clearing required fields
                        in_Np[u].clear();
                        in_Nm[u].clear();
                        for (int d : V[u]) in_Np[d].clear();
                        for (int d : revV[u]) in_Nm[d].clear();
                    }

                    for (int u : W0) { // creating Np, Nm, in_Np, in_Nm
                        Np[u] = V[u];
                        StandardUtils::removeFromArrayPreserveOrderInplace(Np[u], W, helper);
                        for (int d : Np[u]) in_Np[d].push_back(u);

                        Nm[u] = revV[u];
                        StandardUtils::removeFromArrayPreserveOrderInplace(Nm[u], W, helper);
                        for (int d : Nm[u]) in_Nm[d].push_back(u);
                    }
                }

                int valid_uj = findValidUj(W, W0);

                for (int u : W + W2) {
                    V[u] = prevs[u].first;
                    revV[u] = prevs[u].second;
                }

                if (valid_uj != -1) {
                    VI neigh = StandardUtils::setUnion(tempV[valid_uj], temprevV[valid_uj], helper);
                    neigh = StandardUtils::setDifference(neigh, W, helper);

                    if (debug) DEBUG(edge_mods);

                    for (int d : neigh) was[d] = true;
                    for (auto &[a, b] : edge_mods) {
                        if (was[a] || was[b]) required_edge_fill.emplace_back(a, b);
                    }
                    for (int d : neigh) was[d] = false;

                    assert(required_edge_fill.size() % 2 == 0); // here should be only pi-edges
                }
                return valid_uj != -1;
            };


            for (int a = 0; a < N; a++) {
                if (debug) clog << endl << "Checking trans-dom6 for a: " << a << endl;

                {
                    bool aff = false;
                    VI A = VI({a}) + V[a] + revV[a];
                    int P = A.size();
                    for (int i = 0; i < P; i++) A += V[A[i]] + revV[A[i]];
                    for (int d : A) if (affected[d]) aff = true;
                    if (aff) continue;
                }

                VPII required_edge_fill;
                bool trans_dom6 = checkTransDomination6(a, required_edge_fill);
                if (trans_dom6) {
                    { // marking affected nodes before we update structures
                        VI A = VI({a}) + V[a] + revV[a];
                        int P = A.size();
                        for (int i = 0; i < P; i++) A += V[A[i]] + revV[A[i]];
                        for (int d : A) affected[d] = true;
                    }

                    if (debug) {
                        ENDL(1);
                        clog << "trans-domination6 applies to node: " << a << endl;
                        DEBUG(a);
                        DEBUG(V0[a]);
                        DEBUG(tempV[a]);
                        DEBUG(V[a]);
                        ENDL(1);
                        DEBUG(revV0[a]);
                        DEBUG(temprevV[a]);
                        DEBUG(revV[a]);
                        DEBUG(trans_V[a]);
                        DEBUG(required_edge_fill);
                        ENDL(1);
                    }

                    {
                        VPII incident_edges;
                        for (int u : piV[a]) {
                            for (int d : trans_V[u]) {
                                incident_edges.emplace_back(u, d);
                                incident_edges.emplace_back(d, u);
                            }
                        }
                        StandardUtils::makeUnique(incident_edges);

                        Utils::removeEdges(piV, incident_edges, helper);
                        Utils::removeEdges(trans_V, incident_edges, helper);

                        for (int i = (int) incident_edges.size() - 1; i >= 0; i--) {
                            int x = incident_edges[i].first;
                            int y = incident_edges[i].second;
                            if (StandardUtils::find(V0[x], y)) {
                                swap(incident_edges[i], incident_edges.back());
                                incident_edges.pop_back();
                            }
                        }
                        Utils::removeEdges(V, revV, incident_edges, helper);

                        if (debug) {
                            DEBUG(incident_edges);
                            clog << "trans_V: " << endl;
                            GraphUtils::writeGraphHumanReadable(trans_V);
                            clog << "V: " << endl;
                            GraphUtils::writeGraphHumanReadable(V);
                        }
                    }

                    Utils::addEdges(tempV, temprevV, required_edge_fill, helper);
                    GraphUtils::removeEdges(trans_V, required_edge_fill); // removing undirected edges
                    pi_edges_added_permanently += required_edge_fill.size();


                    nodes_removed.push_back(a);
                    removeNodeFromAll(a);

                    if (safety_mode_return_after_first_removed_node) {
                        V = tempV, revV = temprevV;
                        return {{a}, {}};
                    }

                }
            }
        }
        //****************************************************************** END OF TRANS DOMINATION 6

        if( !nodes_removed.empty() ){
            V = tempV;
            revV = temprevV;
            pi_edges_added.clear();
        }
    }


    if(debug)
        clog << "domination6Inserter() finished, nodes_removed: " << nodes_removed.size() << ", pi_edges_added: "
             << pi_edges_added.size() << ", pi_edges_added_permanently: " << pi_edges_added_permanently << endl;


    assert( nodes_removed.empty() || pi_edges_added.empty() );

    if(nodes_removed.empty() && !add_pi_edges_if_no_nodes_removed) pi_edges_added.clear();

    return {nodes_removed,pi_edges_added};
}

tuple<VI,VPII, vector<ReverseTriangleGadgetReduction*>, VI>
        Reducer::reverseTriangleGadget(bool use_only_when_mixed_domination_applies) {
    constexpr bool debug = false;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB helper(N,false), affected(N,false), was(N,false);

    bool modified = false;

    VB is_pi_node(N,false);
    for( int i=0; i<N; i++ ){
        if( piV[i].size() == V[i].size() && piV[i].size() == revV[i].size() ){
            is_pi_node[i] = true;
        }
    }

    VVI triangles;
    {
        VI pi_nodes_deg3;
        for (int i = 0; i < N; i++) if (is_pi_node[i] && piV[i].size() == 3) pi_nodes_deg3.push_back(i);

        if(debug) DEBUG(pi_nodes_deg3);

        InducedGraph g = GraphInducer::induce(piV, pi_nodes_deg3);

        for (int a = 0; a < g.V.size(); a++) {
            for (int d : g.V[a]) was[d] = true;
            for (int b : g.V[a]) {
                if (b < a) continue;
                for (int c : g.V[b]) if (b < c && was[c]) triangles.push_back(VI({a, b, c}));
            }
            for (int d : g.V[a]) was[d] = false;
        }

        for (VI &tr : triangles) { // remapping triangles
            for (int &d : tr) d = g.nodes[d];
        }
    }

    if(debug) DEBUG(triangles);

    VI nodes_to_remove;
    VPII arcs_to_add;
    vector<ReverseTriangleGadgetReduction*> reductions;
    VI kern_red_dom6;

    VI marker(N,0);
    VB in_triangle(N,false);

    for( VI & tr : triangles ){
        int a = tr[0];
        int b = tr[1];
        int c = tr[2];

        int a2,b2,c2;
        for( int d : piV[a] ) if( d != b && d != c ) a2 = d;
        for( int d : piV[b] ) if( d != a && d != c ) b2 = d;
        for( int d : piV[c] ) if( d != a && d != b ) c2 = d;

        if(cnf.reducer_use_domination) assert( a2 != b2 && a2 != c2 && b2 != c2 );

        VI nds = {a,b,c,a2,b2,c2};
        bool aff = false;
        for(int d : nds) if(affected[d]) aff = true;
        if(aff) continue;


        if(!use_only_when_mixed_domination_applies) {
            if(debug){
                clog << "Valid gadget: " << endl;
                DEBUG(VI({a,b,c,a2,b2,c2}));
                for(int d : nds){
                    DEBUG(d);
                    DEBUG(V[d]);
                    DEBUG(revV[d]);
                }
                ENDL(1);
            }

            if (!is_pi_node[a2] || !is_pi_node[b2] || !is_pi_node[c2]) continue;


            nodes_to_remove += VI({a, b, c});
            arcs_to_add += VPII({{a2, b2}, {b2, c2}, {c2, a2}});
            VI if_nodes = {a2, b2, c2};
            VPII else_nodes;
            else_nodes.emplace_back(b, c);
            else_nodes.emplace_back(c, a);
            else_nodes.emplace_back(a, b);

            reductions.push_back(new ReverseTriangleGadgetReduction(if_nodes, else_nodes));

            for (int d : nds) affected[d] = true;
        }else{

            int dom6_node_to_remove = -1;
            {
                using namespace StandardUtils;
                if( find(piV[a2],b2) ) dom6_node_to_remove = c;
                else if( find(piV[b2],c2) ) dom6_node_to_remove = a;
                else if( find(piV[c2],a2) ) dom6_node_to_remove = b;
            }

            if(dom6_node_to_remove != -1){
                kern_red_dom6.push_back(dom6_node_to_remove);
                for (int d : nds) affected[d] = true;
                if(debug) clog << "dom6 special case node: " << dom6_node_to_remove << endl;
            }else {

                for( int d : piV[a2] ) marker[d]++;
                for( int d : piV[b2] ) marker[d]++;
                for( int d : piV[c2] ) marker[d]++;
                in_triangle[a] = in_triangle[b] = in_triangle[c] = true;

                auto isMixedDominated = [&](int x) {
                    int cnt = 0;
                    for (int d : piV[x]) if (in_triangle[d] || marker[d] > 1) cnt++;
                    return cnt == piV[x].size();
                };

                int x = -1, y = -1, z = -1;
                if( is_pi_node[a2] && isMixedDominated(a2) ){ x = a2, y = b2; z = c2; }
                else if( is_pi_node[b2] && isMixedDominated(b2) ){ x = b2; y = a2; z = c2; }
                else if( is_pi_node[c2] && isMixedDominated(c2) ){ x = c2; y = a2; z = b2; }

                if( x != -1 ){
                    if(debug){
                        clog << "Valid reverse-triangle-gadget with mixed domination: " << endl;
                        DEBUG(VI({a,b,c,a2,b2,c2}));
                        for(int d : nds){
                            DEBUG(d);
                            DEBUG(V[d]);
                            DEBUG(revV[d]);
                        }
                        DEBUG(VI({x,y,z}));
                        ENDL(1);
                    }

                    nodes_to_remove += VI({a, b, c});
                    arcs_to_add += VPII({{y,z}, {z,y}});
                    VI if_nodes = {a2, b2, c2};
                    VPII else_nodes;
                    else_nodes.emplace_back(b, c);
                    else_nodes.emplace_back(c, a);
                    else_nodes.emplace_back(a, b);

                    reductions.push_back(new ReverseTriangleGadgetReduction(if_nodes, else_nodes));

                    for (int d : nds) affected[d] = true; // mark this only if rule was applied
                }

                for( int d : piV[a2] ) marker[d] = 0;
                for( int d : piV[b2] ) marker[d] = 0;
                for( int d : piV[c2] ) marker[d] = 0;
                in_triangle[a] = in_triangle[b] = in_triangle[c] = false;
            }
        }

        if(debug){
            DEBUG(nodes_to_remove);
            DEBUG(arcs_to_add);
            ENDL(2);
        }
    }


    return {nodes_to_remove, arcs_to_add, reductions, kern_red_dom6};
}

vector<FullBipartiteBlockerReduction*> Reducer::fullBipartiteBlocker() {
    constexpr bool debug = false;

    int N = V.size();
    VB affected(N,false);
    VB helper(N,false);
    VB was(N,false);

    int applications = 0;
    vector<FullBipartiteBlockerReduction*> res;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI piV_cp = piV;
    VPII pi_edges = Utils::getAllPIEdges(V, revV, helper);

    auto isDisjoint = [&](VI &X, VI &Y){
        for(int d : Y) was[d] = true;
        bool disjoint = true;
        for(int d : X) if( was[d] ) {disjoint = false; break;}
        for(int d : Y) was[d] = false;
        return disjoint;
    };

    auto hasAllEdges = [&](VI &X, VI & Y){
        for( int x : X ) if( piV[x].size() < Y.size() ) return false;

        for(int d : Y) was[d] = true;

        bool can = true;
        for( int x : X ){
            int cnt = 0;
            for( int d : piV[x] ) if( was[d] ) cnt++;
            if( cnt != Y.size() ){
                can = false;
                break;
            }
        }

        for(int d : Y) was[d] = false;

        return can;
    };

    VB is_pi_node(N);
    for(int i=0; i<N; i++) if( piV[i].size() == V[i].size() && piV[i].size() == revV[i].size() ) is_pi_node[i] = true;

    constexpr int MAX_AB_DEG = 7; // for efficiency reasons

    for( auto [a,b] : pi_edges ){
        if( a > b ) continue;
        if( affected[a] || affected[b] ) continue;
        if( V[a].size() > MAX_AB_DEG || V[b].size() > MAX_AB_DEG ) continue;
        if( revV[a].size() > MAX_AB_DEG || revV[b].size() > MAX_AB_DEG ) continue;

        {
            bool aff = false;
            for(int d : V[a]) if(affected[d]) aff = true;   if(aff) continue;
            for(int d : revV[a]) if(affected[d]) aff = true;   if(aff) continue;
            for(int d : V[b]) if(affected[d]) aff = true;   if(aff) continue;
            for(int d : revV[b]) if(affected[d]) aff = true;   if(aff) continue;
        }

        VI rem = {a,b};

        if(debug){ ENDL(2); DEBUG(VI({a,b})); }

        vector<pair<VI,VI>> to_check = {
                {V[a], V[b]},
                {V[a], revV[b]},
                {revV[a], V[b]},
                {revV[a], revV[b]}
        };

        VI XX, YY;

        for( auto [X,Y] : to_check ){
            StandardUtils::removeFromArrayPreserveOrderInplace(X, rem, helper);
            StandardUtils::removeFromArrayPreserveOrderInplace(Y, rem, helper);
            if( isDisjoint(X,Y) && hasAllEdges(X,Y) ){
                XX = X; YY = Y;
                break;
            }
        }

        if(debug){ DEBUG(XX);DEBUG(YY); }

        if( !XX.empty() ){
            assert(!YY.empty());
            applications++;

            VI temp = V[a] + revV[a] + V[b] + revV[b];

            if(debug){
                ENDL(10);
                clog << "Found valid sets XX: " << XX << "    and YY: " << YY << endl;
                DEBUG(VI({a,b}));
                Utils::writeNeighborhood(V, revV, a);
                Utils::writeNeighborhood(V, revV, b);
                ENDL(1);
            }

            { // remove node a nd b from the graph and add reduction
                VI nodes = {a};
                Utils::removeNodes( V, revV, nodes, helper );
                Utils::removeNodes( piV, piV_cp, nodes, helper );

                nodes = {b}; Utils::removeNodes( V, revV, nodes, helper );
                nodes = {b}; Utils::removeNodes( piV, piV_cp, nodes, helper );

                res.push_back( new FullBipartiteBlockerReduction( XX, b, a ) );
            }

            for(int d : temp) affected[d] = true;
        }
    }


    return res;
}

VPII Reducer::edgeNeighborhoodBlocker(){
    constexpr bool debug = false;

    int N = V.size();
    VB affected(N,false);
    VB helper(N,false);
    VB was(N,false);

    VPII res;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI piV_cp = piV;
    VPII pi_edges = Utils::getAllPIEdges(V, revV, helper);

    constexpr int MAX_AB_DEG = 7;

    auto isDClique = [&](VI & A){
        for( int d : A ) if( piV[d].size() + 1 < A.size() ) return false;

        bool can = true;
        for( int d : A ) helper[d] = true;
        for( int a : A ){
            int cnt = 1;
            for( int d : piV[a] ) if(helper[d]) cnt++;
            if( cnt != A.size() ) can = false;
            if(!can) break;
        }
        for( int d : A ) helper[d] = false;

        return can;
    };


    for( auto [a,b] : pi_edges ){
        if(a > b) continue;
        if( V[a].size() > MAX_AB_DEG || revV[a].size() > MAX_AB_DEG ) continue;
        if( V[b].size() > MAX_AB_DEG || revV[b].size() > MAX_AB_DEG ) continue;
        if(affected[a] || affected[b]) continue;
        if( piV[a].size() == 1 || piV[b].size() == 1)  continue; // there will be no edge to add

        if(debug){ DEBUG(VI({a,b})); }


        VI rem = {a,b};
        if(debug){ ENDL(2); DEBUG(VI({a,b})); }

        vector<pair<VI,VI>> to_check = {
                {V[a], V[b]},
                {V[a], revV[b]},
                {revV[a], V[b]},
                {revV[a], revV[b]}
        };

        bool can = false;
        for( auto [X,Y] : to_check ){
            StandardUtils::removeFromArrayPreserveOrderInplace(X, rem, helper);
            StandardUtils::removeFromArrayPreserveOrderInplace(Y, rem, helper);
            if(isDClique(X) && isDClique(Y)){
                can = true;
                if(debug){
                    clog << "Found valid edge neighborhood blocker:" << endl;
                    DEBUG(X); DEBUG(Y);
                }
                break;
            }
        }

        if(can){
            if(debug){ DEBUG(VI({a,b})); }

            VI A = piV[a];
            StandardUtils::removeFromArrayPreserveOrderInplace(A, rem, helper);
            VI B = piV[b];
            StandardUtils::removeFromArrayPreserveOrderInplace(B, rem, helper);

            if( !A.empty() && !B.empty() ){
                VPII prod = StandardUtils::product(A,B);
                int P = prod.size();
                for(int i=0; i<P; i++) {
                   prod.emplace_back(prod[i].second, prod[i].first);
                }
                res += prod;

                VI tmp = V[a] + revV[a] + V[b] + revV[b];
                for(int d : tmp) affected[d] = true;
            }
        }
    }

    StandardUtils::makeUnique(res);

    return res;
}

tuple<vector<DeskReduction*>,VI, int> Reducer::desk(){
    constexpr bool debug = false;

    int N = V.size();
    VB affected(N,false);
    VB helper(N,false);
    VB was(N,false);

    vector<DeskReduction*> desk_folds;
    VI desk_dominations;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI piV_cp = piV;
    VPII pi_edges = Utils::getAllPIEdges(V, revV, helper);

    auto isDClique = [&](VI & A){
        for( int d : A ) if( piV[d].size() + 1 < A.size() ) return false;

        bool can = true;
        for( int d : A ) helper[d] = true;
        for( int a : A ){
            int cnt = 1;
            for( int d : piV[a] ) if(helper[d]) cnt++;
            if( cnt != A.size() ) can = false;
            if(!can) break;
        }
        for( int d : A ) helper[d] = false;

        return can;
    };

    auto findInducedC4s = [&](){
        constexpr int MAX_DEG = 8;

        VI nodes;
        for( int i=0; i<N; i++ ){
            if(V[i].size() <= MAX_DEG && revV[i].size() <= MAX_DEG){
                if( piV[i].size() >= 2 ){
                    nodes.push_back(i);
                }
            }
        }

        sort(ALL(nodes),[&](int a, int b){
            return piV[a].size() > piV[b].size();
        });
        InducedGraph g = GraphInducer::induce(V, nodes);

        VVI gpiV = Utils::getUnderlyingPIGraph(g.V);
        for( int i=0; i<gpiV.size(); i++ ) sort(ALL(gpiV[i]));
        VVI revgV = GraphUtils::reverseGraph(g.V);

        VVI c4s;

        VB was2(N,false);

        for( int a=0; a < g.V.size(); a++ ){
            for(int u : gpiV[a]) was[u] = true;
            for(int u : g.V[a] + revgV[a]) was2[u] = true;

            for( int b : gpiV[a] ){
                if(b < a) continue;

                for( int u : g.V[b] + revgV[b] ) helper[u] = true;

                for( int c : gpiV[b] ){
                    if( c <= a ) continue;
                    if( was2[c] ) continue;

                    for( int d : gpiV[c] ){
                        if( was[d] && b < d && !helper[d] ) {
                            c4s.push_back({a, b, c, d});
                        }
                    }
                }

                for( int u : g.V[b] + revgV[b] ) helper[u] = false;
            }

            for(int u : gpiV[a]) was[u] = false;
            for(int u : g.V[a] + revgV[a]) was2[u] = false;
        }

        for( auto & dsk : c4s ){
            for( int & d : dsk ) d = g.nodes[d];
        }

        return c4s;
    };

    VVI desks = findInducedC4s();

    VPII arcs_to_add;

    for( auto dsk : desks ){
        int a = dsk[0];
        int b = dsk[1];
        int c = dsk[2];
        int d = dsk[3];

        if(affected[a] || affected[b] || affected[c] || affected[d]) continue;
        bool aff = false;
        VI tmp = V[a] + revV[a] + V[b] + revV[b] + V[c] + revV[c] + V[d] + revV[d];
        for(int x : tmp) if(affected[x]) aff = true;
        if(aff) continue;

        VVI Xs;
        for( int i=0; i<4; i++ ){
            VI X = V[dsk[i]];
            StandardUtils::removeFromArrayPreserveOrderInplace( X, dsk, helper );
            if( isDClique(X) ) Xs.push_back(X);
            else{
                X = revV[dsk[i]];
                StandardUtils::removeFromArrayPreserveOrderInplace( X, dsk, helper );
                if( isDClique(X) ) Xs.push_back(X);
            }

            if(Xs.size() < i+1) break;
        }

        if(Xs.size() == 4){

            if(debug){
                ENDL(5);  clog << "Found C4 with valid Xs:" << endl;  DEBUG(VI({a,b,c,d})); DEBUG(Xs);
                for(int u : dsk){
                    DEBUG(u);
                    DEBUG(V[u]);
                    DEBUG(revV[u]);
                }
            }

            VI pi_inters;

            auto addInters = [&was, &pi_inters]( VI & x, VI & y ){
                for( int d : y ) was[d] = true;
                for(int d : x) if(was[d]) pi_inters.push_back(d);
                for( int d : y ) was[d] = false;
            };

            addInters( piV[a], piV[b] );
            addInters( piV[b], piV[c] );
            addInters( piV[c], piV[d] );
            addInters( piV[d], piV[a] );

            if(debug) DEBUG(pi_inters);
            while(pi_inters.size() > 1) pi_inters.pop_back();
            if(debug) DEBUG(pi_inters);


            if(!pi_inters.empty()){
                desk_dominations += pi_inters;

                if(debug){
                    clog << "Found valid desk domination" << endl;
                    DEBUG(VI({a,b,c,d})); DEBUG(Xs); DEBUG(pi_inters);
                }

                Utils::removeNodes(V, revV, pi_inters, helper);
                Utils::removeNodes(piV, piV_cp, pi_inters, helper);

                for(int x : tmp) affected[x] = true;
            }else{
                bool can_apply_desk_reduction = true;
                if( !Utils::isPiNode(V, revV, piV, a) || piV[a].size() < 3 ) can_apply_desk_reduction = false;
                if( !Utils::isPiNode(V, revV, piV, b) || piV[b].size() < 3) can_apply_desk_reduction = false;
                if( !Utils::isPiNode(V, revV, piV, c) || piV[c].size() < 3) can_apply_desk_reduction = false;
                if( !Utils::isPiNode(V, revV, piV, d) || piV[d].size() < 3) can_apply_desk_reduction = false;

                if(can_apply_desk_reduction) {
                    VI if_nodes = piV[a] + piV[c];
                    StandardUtils::removeFromArrayPreserveOrderInplace(if_nodes, dsk, helper);
                    PII then_nodes = {b, d};
                    PII else_nodes = {a, c};

                    desk_folds.push_back(new DeskReduction(if_nodes, then_nodes, else_nodes));

                    if (debug) clog << "Applying desk folding: " << desk_folds.back()->toString() << endl;

                    VI if_nodes2 = piV[b] + piV[d];
                    StandardUtils::removeFromArrayPreserveOrderInplace(if_nodes2, dsk, helper);
                    VPII arcs_to_add = StandardUtils::product(if_nodes, if_nodes2);
                    int P = arcs_to_add.size();
                    for (int i = 0; i < P; i++) arcs_to_add.emplace_back(arcs_to_add[i].second, arcs_to_add[i].first);

                    { // remove dsk from graph
                        Utils::removeNodes(V, revV, dsk, helper);
                        Utils::removeNodes(piV, piV_cp, dsk, helper);
                    }

                    { // add neccessary pi-connections
                        Utils::addEdges(V, arcs_to_add, helper);
                        // no need to reverse - there are only pi-edges in arcs_to_add
                        Utils::addEdges(revV, arcs_to_add, helper);

                        Utils::addEdges(piV, arcs_to_add, helper);
                    }

                    for (int x : tmp) affected[x] = true;
                }else{
                    VPII to_add;
                    auto addArcs = [&](int a, int b){
                        if( piV[a].size() == 2 || piV[b].size() == 2 ) return;

                        VI Xa = piV[a];
                        StandardUtils::removeFromArrayPreserveOrderInplace(Xa, dsk, helper);
                        VI Xb = piV[b];
                        StandardUtils::removeFromArrayPreserveOrderInplace(Xb, dsk, helper);

                        VPII prod = StandardUtils::product(Xa,Xb);
                        to_add += prod;
                        for(auto & [a,b] : prod) swap(a,b);
                        to_add += prod;
                    };

                    addArcs(a,b);
                    addArcs(b,c);
                    addArcs(c,d);
                    addArcs(d,a);

                    if(!to_add.empty()){
                        arcs_to_add += to_add;
                        if(debug){
                            DEBUG(VI({a,b,c,d}));
                            DEBUG(Xs);
                            clog << "In desk() found pi-arcs to add: " << to_add << endl;
                        }
                    }
                }
            }
        }
    }

    int arc_diff = 0;
    { // adding pi-arcs
        int arcs_before = GraphUtils::countEdges(V,true);

        StandardUtils::makeUnique(arcs_to_add);
        Utils::addEdges(V, arcs_to_add, helper);
        Utils::addEdges(revV, arcs_to_add, helper);

        Utils::addEdges(piV, arcs_to_add, helper);

        int arcs_after = GraphUtils::countEdges(V,true);
        arc_diff = arcs_after - arcs_before;

        if( arc_diff > 0 && V.size() > 50 ){
                clog << "There were " << arc_diff << " arcs added in desk() reduction" << endl;
        }
    }

    return {desk_folds, desk_dominations, arc_diff};
}

VI Reducer::bottleneck2() {
    constexpr bool debug = false;

    int N = V.size();
    VB affected(N,false);
    VB helper(N,false);
    VB was(N,false);

    VI res;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI piV_cp = piV;
    VI temp;

    for(int w=0; w<N; w++){
        if(affected[w]) continue;
        VI W = piV[w];
        if( W.size() > cnf.reducer_max_bottleneck2_neighborhood_size) continue;

        bool aff = false;
        for(int d : W) if(affected[d]) aff=true;
        if(aff) continue;

        temp.clear();
        temp.push_back(w); temp += W;
        VI P;

        for(int d : temp) was[d] = true;
        for( int d : W ){
            VI neigh = V[d] + revV[d];
            for( int dd : neigh ){
                if (!was[dd]) {
                    temp.push_back(dd);
                    P.push_back(dd);
                    was[dd] = true;
                }
            }
        }
        for(int d : temp) was[d] = false;

        for( int d : P ) if(affected[d]) aff=true;
        if(aff) continue;

        if(  P.size() > cnf.reducer_max_bottleneck2_neighborhood_size) continue;

        int y = 0, z = 0;
        {
            InducedGraph g = GraphInducer::induce(V,W);
            auto dfvs = Utils::getMinDFVSDefault(g.V);
            y = dfvs.size();
        }

        {
            InducedGraph g = GraphInducer::induce(V,P);
            auto dfvs = Utils::getMinDFVSDefault(g.V);
            z = dfvs.size();
        }

        if( P.size() - z < W.size() - y ){
            if(debug){
                ENDL(2);
                clog << "Found valid bottleneck2:" << endl;
                DEBUG(w);
                DEBUG(W); DEBUG(y);
                DEBUG(P); DEBUG(z);
                for( int d : W ) Utils::writeNeighborhood(V, revV, d);
            }

            res.push_back(w);

            for( int d : P ) affected[d] = true;
            for( int d : W ) affected[d] = true;
            affected[w] = true;
        }

    }


    return res;
}

vector<GeneralFoldingReduction *> Reducer::generalFolding() {
    constexpr bool debug = false;

    int N = V.size();
    VB affected(N,false);
    VB helper(N,false);
    VB was(N,false);

    vector<GeneralFoldingReduction *> res;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI piVcp = piV;
    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB is_pi_node(N,false);
    for( int i=0; i<N; i++ ) is_pi_node[i] = Utils::isPiNode(V, revV, piV,i);

    constexpr int MAX_SIZE = 10;

    auto partiallyDominates = [&]( int a, int b, VI & W ){
        if( is_pi_node[b] ) return true;

        for(int d : W) helper[d] = true;

        bool can = true;
        {
            for( int d : V[a] ) was[d] = true;
            for( int d : nonpiV[b] ){
                if(!helper[d]){
                    if( !was[d] ){
                        can = false;
                        break;
                    }
                }
            }
            for( int d : V[a] ) was[d] = false;
        }

        if(can){
            for( int d : revV[a] ) was[d] = true;
            for( int d : revnonpiV[b] ){
                if(!helper[d]){
                    if( !was[d] ){
                        can = false;
                        break;
                    }
                }
            }
            for( int d : revV[a] ) was[d] = false;
        }

        for(int d : W) helper[d] = false;

        return can;

    };

    VI order(N); iota(ALL(order),0);
    sort(ALL(order), [&](int a, int b){
        return piV[a].size() < piV[b].size();
    });
    reverse(ALL(order)); // #TEST - starting node selection from largest degree in general_folding

    for( int w : order ){
        if(V[w].size() <= 1 || revV[w].size() <= 1) continue;
        if( !is_pi_node[w] ) continue;
        if(affected[w]) continue;
        if( piV[w].size() > cnf.reducer_max_general_folding_neighborhood_size ) continue;

        VI W = piV[w];
        bool aff = false;
        for(int d : W){
            if(affected[d]) aff = true;
            for( int u : V[d] ) if(affected[u]) aff = true;
            for( int u : revV[d] ) if(affected[u]) aff = true;
        }
        if(aff) continue;

        int dfvs_size;
        InducedGraph g = GraphInducer::induce(V, W);
        dfvs_size = Utils::getMinDFVSDefault(g.V).size();


        if( dfvs_size + 2 < W.size() ) continue;

        if(debug){
            ENDL(5);
            clog << "Found W with DFVS(G[W]) >= W.size()-2" << endl;
            DEBUG(w);
            DEBUG(W);
            Utils::writeNeighborhood(V, revV, w);
        }


        VPII antiedges;
        {
            VVI pigv = Utils::getUnderlyingPIGraph(g.V);
            VVI comppigv = GraphUtils::getComplimentaryGraph(pigv);
            antiedges = GraphUtils::getGraphEdges(comppigv, false); // we want undirected here
            for( auto & [a,b] : antiedges ){ // remapping antiedges to original ids
                a = g.nodes[a];
                b = g.nodes[b];
            }
        }

        if(debug) DEBUG(antiedges);

        if( antiedges.size() > min((int)W.size(), cnf.reducer_max_general_folding_antiedges) ){
            if(debug) clog << "There are too many antiedges: " << antiedges.size() << endl;
            continue;
        }

        bool can_apply = true;

        for( auto [a,b] : antiedges ){
            if( !partiallyDominates(a,b,W) && !partiallyDominates(b,a,W) ){
                if(debug) clog << "Partial domination condition not met for antiedge (" << a << "," << b << ")" << endl;
                can_apply = false;
            }
            if(!can_apply) break;
        }

        if(!can_apply) continue;

        affected[w] = true;
        for(int d : W){
            affected[d] = true;
            for( int u : V[d] ) affected[u] = true;
            for( int u : revV[d] ) affected[u] = true;
        }

        VPII arcs_to_remove;
        {
            for( int d : W ){
                for( int u : V[d] ) arcs_to_remove.emplace_back(d,u);
                for( int u : revV[d] ) arcs_to_remove.emplace_back(u,d);
            }
            StandardUtils::makeUnique(arcs_to_remove);
        }

        VPII arcs_to_add;
        vector<tuple<int,int,int>> antiedges_tuples;

        { // finding arcs to add
            VI free_ids = (VI(1,w) + W);
            if(debug) DEBUG(free_ids);

            for( int i=0; i<antiedges.size(); i++){
                int id = free_ids[i];
                int a = antiedges[i].first;
                int b = antiedges[i].second;

                if(debug){ ENDL(1);DEBUG(VI({id,a,b})); }

                VI Np = StandardUtils::setUnion( nonpiV[a], nonpiV[b], helper );
                StandardUtils::removeFromArrayPreserveOrderInplace( Np, W, helper );
                arcs_to_add += StandardUtils::product( VI({id}), Np );

                VI Nm = StandardUtils::setUnion( revnonpiV[a], revnonpiV[b], helper );
                StandardUtils::removeFromArrayPreserveOrderInplace( Nm, W, helper );
                arcs_to_add += StandardUtils::product( Nm, VI({id}) );

                VI Npi = StandardUtils::setUnion(piV[a], piV[b], helper);
                StandardUtils::removeFromArrayPreserveOrderInplace( Npi, free_ids, helper );
                arcs_to_add += StandardUtils::product( VI({id}), Npi );
                arcs_to_add += StandardUtils::product( Npi, VI({id}) );

                antiedges_tuples.emplace_back( id,a,b );

                if(debug){
                    DEBUG(VI({id,a,b}));
                    DEBUG(Np);DEBUG(Nm);DEBUG(Npi);
                    sort(ALL(arcs_to_add)); DEBUG(arcs_to_add);
                    ENDL(1);
                }
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

        if(debug){ DEBUG(arcs_to_add);DEBUG(arcs_to_remove); }

        if(debug){
            clog << "State before application:" << endl;
            Utils::writeNeighborhood(V, revV, w);
        }

        Utils::removeEdges(V, revV, arcs_to_remove, helper);

        Utils::addEdges(V, revV, arcs_to_add, helper);
        res.push_back( new GeneralFoldingReduction( w, W, antiedges_tuples ) );

        if(debug){
            clog << "State after application:" << endl;
            Utils::writeNeighborhood(V, revV, w);
        }

    }
    return res;
}

VI Reducer::recursiveReducer() {
    const bool debug = (V.size() > 50);

    clog << "Starting recursive reducer" << endl;

    int N = V.size();
    VB affected(N,false);
    VB helper(N,false);
    VB was(N,false);

    VI res;

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI piVcp = piV;
    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB is_pi_node(N,false);
    for( int i=0; i<N; i++ ) is_pi_node[i] = Utils::isPiNode(V, revV, piV,i);

    VVI semi_alternatives;

    VI marker(N,0);

    auto checkABCorrectness = [&](VI &A, VI& B){
        for(int b : B) was[b] = true;
        for( int a : A ){
            int cnt = 0;
            for(int d : piV[a]) if(was[d]) cnt++;
            if( cnt != B.size() ) return false;
        }
        for(int b : B) was[b] = false;
        return true;
    };


    auto maximizeABs = [&]( VI & A, VI & B ){
        assert(!A.empty() && !B.empty());

        bool can = true;
        int swaps = 0;
        auto check = [&piV, &marker,&was,&can]( VI & A, VI & B ) {
            for (int b : B) was[b] = true;
            for (int a : A) for (int d : piV[a]) if (!was[d]) marker[d]++;
            int el = -1, val = 0;
            for (int a : A) {
                StandardUtils::shuffle(piV[a]);
                for (int d : piV[a]) {
                    if (marker[d] == A.size() && piV[d].size() > val) {
                        val = piV[d].size();
                        el = d;
                    }
                }
            }
            for (int a : A) for (int d : piV[a]) marker[d] = 0;
            for (int b : B) was[b] = false;

            if (el != -1) {
                B.push_back(el);
                can = true;
            }
        };


        while(can){
            can = false;
            if( A.size() <= B.size() ) check(A,B);
            if( !can || A.size() >= B.size() ) check(B,A);
        }


        if(swaps & 1) swap(A,B);
        assert(checkABCorrectness(A,B));
    };

    auto constructSetAB0 = [&](int u, int v){
        VI A = {u};
        VI B = {v};
        maximizeABs(A,B);
        semi_alternatives = {A,B};
    };

    auto constructSetAB1 = [&](int u, int v = -1){
        VI A = {u}, B = piV[u];

        int s1 = A.size(), s2 = B.size();
        semi_alternatives.clear();
        if(piV[u].empty()) return;

        maximizeABs( A,B );
        semi_alternatives = VVI( {A,B} );
    };

    auto constructSetAB2 = [&](int u, int v = -1){
        VI nodes = piV[u] + VI({u});
        InducedGraph g = GraphInducer::induce(piV, nodes);
        assert( g.nodes[g.perm[u]] == u );
        VI clq = CliqueExtension::maximizeCliqueGreedy(g.V, VI({ g.perm[u] }));
        for(int & d : clq) d = g.nodes[d];

        assert(CliqueUtils::isClique(piV,clq));

        semi_alternatives.clear();
        if(clq.size() < 4) return;

        StandardUtils::shuffle(clq);
        VI A = StandardUtils::slice( clq, 0, (clq.size()+1)/2 );
        VI B = StandardUtils::slice( clq, (clq.size()+1)/2, clq.size() );

        maximizeABs(A,B);
        semi_alternatives = VVI( {A,B} );
    };

    auto constructSetAB3 = [&](int u, int v = -1){
        VI nodes = piV[u] + VI({u});
        InducedGraph g = GraphInducer::induce(piV, nodes);
        assert( g.nodes[g.perm[u]] == u );
        VI clq = CliqueExtension::maximizeCliqueGreedy(g.V, VI({ g.perm[u] }));
        for(int & d : clq) d = g.nodes[d];

        assert(CliqueUtils::isClique(piV,clq));

        semi_alternatives.clear();
        if(clq.size() < 3) return;

        for( int i=0; i<clq.size(); i++ ){
            VI A;
            for( int j=0; j<clq.size(); j++ ) if(i != j) A.push_back(clq[j]);
            semi_alternatives.push_back(A);
        }
    };

    auto constructSetAB4 = [&](int u, int v = -1){
        constexpr int MIN_DEG = 3, MAX_DEG = 17;
        semi_alternatives.clear();
        if( piV[u].size() < MIN_DEG || piV[u].size() > MAX_DEG ) return;

        int T = piV[u].size();
        StandardUtils::shuffle(piV[u]);
        VI A,B;
        VI best_A, best_B;

        for( int i=1; i<(1<<T); i++ ){
            if( __builtin_popcount(i) != MIN_DEG ) continue;

            A.clear();
            B.clear();
            for( int j=0; j<T; j++ ){
                if( i & (1<<j) ) B.push_back( piV[u][j] );
            }

            for( int b : B ){
                for( int d : piV[b] ){
                    marker[d]++;
                    if(marker[d] == B.size()) A.push_back(d);
                }
            }
            for( int b : B ) for( int d : piV[b] ) marker[d] = 0;

            maximizeABs(A,B);
            if( A.size() * B.size() > best_A.size() * best_B.size() ){
                best_A = A;
                best_B = B;
            }
        }

        semi_alternatives = {best_A, best_B};
    };

    auto getAllNodesInSomeOptDFVSForSet = [&]( VI A ){
        for(int d : A) was[d] = true;
        VI nodes; nodes.reserve(N);
        for(int i=0; i<N; i++) if( !V[i].empty() && !revV[i].empty() && !was[i] ) nodes.push_back(i);
        for(int d : A) was[d] = false;
        InducedGraph g = GraphInducer::induce(V,nodes);

        constexpr bool use_folding = true;

        Reducer red(g.V,cnf);
        red.cnf.write_logs = false;
        red.cnf.enableAllReductions();
        red.cnf.disableAllConditionalReductions();
        red.cnf.disableAllRecursiveReductions();
        if(use_folding) red.cnf.reducer_use_folding = true;

        red.cnf.reducer_use_recursive_reducer = false;
        red.cnf.reducer_use_domination_5 = false;
        red.cnf.reducer_use_domination_6inserter = false; // this would be just too time-consuming, but it would work
        red.cnf.reducer_use_nonsimple_cycle_arcs_full = false;
        red.cnf.reducer_use_mixed_domination_full = false;

        auto reductions = red.reduce();
        VI red_dfvs;
        if(!use_folding) red_dfvs = Reducer::convertKernelizedReductions(reductions); // this is here
        else {
            VB in_red_dfvs(N, false);
            for(int i=(int)reductions.size()-1; i>=0; i--) {
                DFVSReduction* redu = reductions[i];
                KernelizedNodesReduction *kernred = dynamic_cast<KernelizedNodesReduction *>(redu);
                if (kernred != nullptr) {
                    red_dfvs += kernred->getKer();
                    for (int d : kernred->getKer()) in_red_dfvs[d] = true;
                } else {
                    // we have folding
                    FoldingReduction *foldred = dynamic_cast<FoldingReduction *>(redu);
                    assert(foldred != nullptr);
                    int if_node = foldred->getIfNode();
                    int else_node = foldred->getElseNode();
                    int folding_node = foldred->getFoldingNode();

                    if (in_red_dfvs[if_node]) {
                        red_dfvs.push_back(else_node);
                        in_red_dfvs[else_node] = true;
                    }
                }
            }
        }
        for(int & d : red_dfvs) d = g.nodes[d];

        return A + red_dfvs;
    };

    auto proceedForSemiAlternatives = [&](int i){
        TimeMeasurer::start("proceedForSemiAlternatives");
        if(semi_alternatives.empty()) return VI();
        VI inters(N,0);

        for( int j=0; j<semi_alternatives.size(); j++ ){
            VI A = semi_alternatives[j];
            VI red_dfvs = getAllNodesInSomeOptDFVSForSet(A);
            VI zb = StandardUtils::setUnion(A,red_dfvs,helper);
            if(zb.size() != red_dfvs.size()){
                ENDL(2);
                DEBUG(zb);
                DEBUG(A);
                DEBUG(red_dfvs);
            }
            assert(zb.size() == red_dfvs.size());
            for( int d : zb ) inters[d]++;
        }

        VI to_remove;
        int mx = semi_alternatives.size();
        for( int j=0; j<N; j++ ) if( inters[j] == mx ) to_remove.push_back(j);

        if( !to_remove.empty() ){
            if(debug){
                ENDL(5);
                clog << "Found nonempty intersection over all recursive reductions: " << endl;
                DEBUG(i);
                clog << "semi_alternatives:" << endl;
                for(auto & sa : semi_alternatives) clog << sa << endl;
                clog << "semi_alternatives sizes: ";
                for(auto & sa : semi_alternatives) clog << sa.size() << " ";
                clog << endl;
                DEBUG(to_remove);
            }

            Utils::removeNodes( V, revV, to_remove, helper );
            Utils::removeNodes( piV, piVcp, to_remove, helper );
        }

        TimeMeasurer::stop("proceedForSemiAlternatives");
        return to_remove;
    };

    constexpr bool return_on_first_success = false;
    const bool return_early = false;
    constexpr int res_size_to_return = 10;

    VI order = CombinatoricUtils::getRandomPermutation(N);

    int cnt = 0;
    for( int i : order ){
        if(debug){
            clog << "\rRecursive Reducer called for #" << (++cnt) << " / " << N
                 << ", res.size(): " << res.size() << ", total_res: " << total_recursive_reducer_nodes_removed << flush;
        }
        if(V[i].empty() || revV[i].empty()) continue;

        if(return_early && (res.size() >= res_size_to_return)) return res;

        VI removed;
        {
            constructSetAB1(i);
            if (semi_alternatives.empty()) continue;
            removed = proceedForSemiAlternatives(i);
            res += removed;
            if (!removed.empty()){
                if(return_on_first_success) return res;
                else continue;
            }
        }

        {
            constructSetAB2(i);
            if (semi_alternatives.empty()) continue;
            removed = proceedForSemiAlternatives(i);
            res += removed;
            if (!removed.empty()){
                if(return_on_first_success) return res;
                else continue;
            }
        }

        {
            constructSetAB3(i);
            if (semi_alternatives.empty()) continue;
            removed = proceedForSemiAlternatives(i);
            res += removed;
            if (!removed.empty()){
                if(return_on_first_success) return res;
                else continue;
            }
        }

        {
            constructSetAB4(i);
            if (semi_alternatives.empty()) continue;
            removed = proceedForSemiAlternatives(i);
            res += removed;
            if (!removed.empty()){
                if(return_on_first_success) return res;
                else continue;
            }
        }

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

    VVI piV = Utils::getUnderlyingPIGraph(V);
    VVI superPiV = Utils::getSuperPIGraph(V);

    VVI nonpiV = Utils::getNonPIGraph(V);
    VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

    VB is_pi_node(N,false);
    for( int i=0; i<N; i++ ) is_pi_node[i] = Utils::isPiNode(V, revV, piV,i);

    VVI & G = superPiV;

    VI order = CombinatoricUtils::getRandomPermutation(N);
    sort(ALL(order), [&](int a, int b){ return G[a].size() > G[b].size(); } );


    auto unconfined = [&](int v){
        VI S = {v};
        VI ws;  // N(u) \ N[S]
        VI NS = G[v];
        in_S[v] = true;
        for(int d : NS) in_NS[d] = true;

        if(debug){
            ENDL(2);
            clog << "Checking node v: " << v << endl;
            Utils::writeNeighborhood(G, G, v);
        }

        bool can = true;

        while(can) {

            if(debug){DEBUG(S); DEBUG(NS);}

            int best_u = -1, best_val = 1e9;
            VI best_ws;
            for (int u : NS) {
                int cnt = 0;
                ws.clear();
                for( int d : G[u] ){
                    if(in_S[d]) cnt++;
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
            if(best_ws.empty()){
                /* node v is unconfined*/
                if(debug) clog << "Node v is unconfined!" << endl;
                can = true;break;
            }
            if( best_ws.size() == 1 ){
                int w = best_ws[0];
                if(!is_pi_node[w]){can = false; break;}
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
            if( !is_pi_node[d] || affected[d] ) can = false;
            for( int dd : G[d] ) if( !is_pi_node[dd] ) can = false;
        }

        return can;
    };

    for( int v : order ){
        if( !is_pi_node[v] || affected[v] ) continue;

        bool aff = false;
        for(int d : G[v]) if(affected[d] || !is_pi_node[d]) aff = true;
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

