//
// Created by sylwester on 12/20/21.
//

#include <graphs/GraphInducer.h>
#include <graphs/toposort/TopoSort.h>
#include <combinatorics/CombinatoricUtils.h>
#include <utils/StandardUtils.h>
#include <graphs/cycles/CycleCounter.h>
#include <graphs/VertexCover/VCUtils.h>
#include <graphs/VertexCover/approximation/NuMVC/NuMVC.h>
#include <graphs/VertexCover/approximation/LibMVC/fastvc.h>
#include <graphs/GraphWriter.h>
#include <utils/TimeMeasurer.h>
#include <CONTESTS/PACE22/Config.h>
#include <CONTESTS/PACE22/Reducer.h>
#include <graphs/scc/StronglyConnectedComponents.h>
#include <CONTESTS/PACE22/heur/DFVSSolverH.h>
#include <graphs/VertexCover/kernelization/KernelizerVC.h>
#include <CONTESTS/PACE22/heur/HittingSetLS.h>
#include <CONTESTS/PACE22/exact/DFVSSolverE.h>
#include "CONTESTS/PACE22/Utils.h"
#include "graphs/GraphUtils.h"

namespace Utils{


    bool addEdge(VVI & V, int a, int b, bool just_push_back){
        if(just_push_back){
            V[a].push_back(b);
            return true;
        }else{
            bool add = true;
            for( int d : V[a] ) if( d == b ){ add = false;break; }
            if(add) V[a].push_back(b);
            return add;
        }
    }

    void addEdges(VVI & V, VPII & edges, VB & helper){
        sort(ALL(edges));
        int p=0, q=0;
        while( p < edges.size() ){
            q = p+1;
            while( q < edges.size() && edges[q].first == edges[p].first ) q++;

            int a = edges[p].first;
            for( int d : V[a] ) helper[d] = true;
            for( int i=p; i<q; i++ ){
                int b = edges[i].second;
                if( !helper[b] ){
                    helper[b] = true;
                    V[a].push_back(b);
                }
            }
            for( int d : V[a] ) helper[d] = false;

            p = q;
        }
    }

    void addEdges(VVI &V, VVI &revV, VPII arcs, VB &helper) {
        addEdges(V, arcs, helper);
        for(auto & [a,b] : arcs) swap(a,b);
        addEdges(revV, arcs, helper);
    }

    void removeEdge( VVI & V, int a, int b ){
        GraphUtils::removeEdge(V,a,b, true); // removes directed edge a->b
    }

    void removeEdges( VVI & V, VPII & edges, VB & helper ){
        sort(ALL(edges));
        int p=0, q=0;
        while( p < edges.size() ){
            q = p+1;
            while( q < edges.size() && edges[q].first == edges[p].first ) q++;

            int a = edges[p].first;
            for( int i=p; i<q; i++ ) helper[ edges[i].second ] = true;
            for( int i=(int)V[a].size()-1; i>=0; i-- ){
                int d = V[a][i];
                if( helper[d] ){
                    swap( V[a][i], V[a].back() );
                    V[a].pop_back();
                }
            }
            for( int i=p; i<q; i++ ) helper[ edges[i].second ] = false;
            p = q;
        }
    }

    void removeEdges(VVI &V, VVI &revV, VPII arcs, VB &helper) {
        removeEdges(V, arcs, helper);
        for(auto & [a,b] : arcs) swap(a,b);
        removeEdges(revV, arcs, helper);
    }


    void merge( VVI & V, VVI & revV, int a, VB & helper ){
        if(V[a].empty() && revV[a].empty()) return;

        VPII add_edges;
        add_edges.reserve( V[a].size() * revV[a].size() );
        for( int pred : revV[a] ){
            for( int succ : V[a] ){
                add_edges.push_back({pred,succ});
            }
        }

        VPII remove_edges;
        remove_edges.reserve( V[a].size() + revV[a].size() );

        for( int d : V[a] ) remove_edges.push_back({a,d});
        for( int d : revV[a] ) remove_edges.push_back({d,a});

        removeEdges(V, remove_edges, helper);
        addEdges(V, add_edges, helper);

        for( auto & [a,b] : add_edges ) swap(a,b);
        for( auto & [a,b] : remove_edges ) swap(a,b);

        removeEdges(revV, remove_edges, helper);
        addEdges(revV, add_edges, helper);

        assert(V[a].empty());
        assert(revV[a].empty());
    }

    void writeNeighborhood(VVI & V, VVI & revV, int v){
        DEBUG(v);
        DEBUG(V[v]);
        DEBUG(revV[v]);
        ENDL(1);

        set<int> zb(ALL(V[v])); zb.insert(ALL(revV[v]));
        for(int d : zb){
            DEBUG(d);
            DEBUG(V[d]);
            DEBUG(revV[d]);
            ENDL(1);
        }
        ENDL(1);
    }

    void partialMerge(VVI &V, VVI &revV, VVI &nonpiV, VVI &revnonpiV, VVI &piV, int v, VB &helper) {
        VPII to_remove;
        for( int x : revnonpiV[v] ) to_remove.emplace_back( x,v );
        for( int y : nonpiV[v] ) to_remove.emplace_back( v,y );

        VPII to_add = StandardUtils::product( revnonpiV[v], nonpiV[v] );

        VI Np = nonpiV[v];
        VI Nm = revnonpiV[v];

        {
            // making the changes in graph V and revV
            removeEdges(V, to_remove, helper);
            for(auto& [a,b] : to_remove) swap(a,b);
            removeEdges(revV, to_remove, helper);

            addEdges(V, to_add,helper);
            for(auto& [a,b] : to_add) swap(a,b);
            addEdges(revV, to_add,helper);

            { // updating nonpi and pi structures
                for( int d : Nm + Np ){
                    nonpiV[d] = Utils::getNonPISuccessors( V, revV, d, helper );
                    revnonpiV[d] = Utils::getNonPIPredecessors( V, revV, d, helper );
                }

                for(int d : Nm + Np){
                    piV[d].clear();
                    for(int dd : V[d]) helper[dd] = true;
                    for(int dd : revV[d]) if(helper[dd]) piV[d].push_back(dd);
                    for(int dd : V[d]) helper[dd] = false;
                }


                nonpiV[v].clear();
                revnonpiV[v].clear();
            }
        }
    }


    bool hasLoop( VVI & V, int a ){
        for( int d : V[a] ) if(d == a) return true;
        return false;
    }

    VPII getPIEdges( VVI & V, VVI & revV, int a, VB & helper){
        VPII pie_edges;
        for( int d : V[a] ) helper[d] = true;
        for( int d : revV[a] ) if(helper[d]) pie_edges.emplace_back(a,d);
        for( int d : V[a] ) helper[d] = false;
        return pie_edges;
    }

    VPII getAllPIEdges(VVI &V, VVI &revV, VB &helper) {
        VPII pie_edges;
        for( int i=0; i<V.size(); i++ ) pie_edges += getPIEdges(V, revV, i, helper);
        return pie_edges;
    }

    VI getNonPISuccessors( VVI & V, VVI & revV, int a, VB & helper ){
        VI res; res.reserve(V[a].size());
        for( int d : revV[a] ) helper[d] = true;
        for( int d : V[a] ) if(!helper[d]) res.push_back(d);
        for( int d : revV[a] ) helper[d] = false;
        return res;
    }

    VI getNonPIPredecessors( VVI & V, VVI & revV, int a, VB & helper ){
       VI res; res.reserve( revV[a].size() );
       for(int d : V[a]) helper[d] = true;
       for( int d : revV[a] ) if(!helper[d]) res.push_back(d);
       for(int d : V[a]) helper[d] = false;
       return res;
    }

    void removeNode(VVI &V, VVI &revV, int a, VB & helper) {
        for( int d : V[a] ) removeEdge( revV, d,a );
        for( int d : revV[a] ) removeEdge( V, d,a );
        V[a].clear();
        revV[a].clear();
    }

    void removeNodes(VVI &V, VVI &revV, VI &nodes, VB &helper) {
        VPII edges;
        for( int a : nodes ){
            for(int d : V[a]) edges.push_back({a,d});
            for(int d : revV[a]) edges.push_back({d,a});
        }

        sort(ALL(edges));
        edges.resize(unique(ALL(edges)) - edges.begin());

        removeEdges( V, edges, helper );
        for( auto & [a,b] : edges ) swap(a,b);
        removeEdges( revV, edges, helper );
    }

    bool isCorresponding(VVI V, VVI revV) {
        if( V.size() != revV.size() ) return false;
        revV = GraphUtils::reverseGraph(revV);
        for(int i=0; i<V.size(); i++) if( V[i].size() != revV[i].size() ) return false;
        for(int i=0; i<V.size(); i++){
            sort(ALL(V[i]));
            sort(ALL(revV[i]));
            if( !equal(ALL(V[i]), ALL(revV[i])) ) return false;
        }
        return true;
    }

    bool hasCycle(VVI &V) {
        int N = V.size();
        VB was(N,false), on_path(N,false);

        function<bool(int)> dfs = [&](int a){
            if(on_path[a]) return true;
            on_path[a] = true;
            was[a] = true;

            for(int d : V[a]){
                if(on_path[d]) return true;

                if(!was[d]){
                    if( dfs(d) ) return true;
                }
            }

            on_path[a] = false;
            return false;
        };

        for( int i=0; i<N; i++ ) if(!was[i] && dfs(i)) return true;
        return false;
    }

    bool hasCycle( VVI & V, VI & A, VB& was, VB& on_path, VB & helper ){
        for(int a : A) helper[a] = true;

        int N = V.size();
        VI temp;

        function<bool(int)> dfs = [&](int a){
            if(on_path[a]) return true;
            on_path[a] = true;
            was[a] = true;
            temp.push_back(a);

            for(int d : V[a]){
                if(!helper[d]) continue;
                if(on_path[d]) return true;

                if(!was[d]){
                    if( dfs(d) ) return true;
                }
            }

            on_path[a] = false;
            return false;
        };

        bool res = false;
        for( int i : A ){
            if(!was[i] && dfs(i)){
                res = true;
                break;
            }
        }

        for(int a : temp) was[a] = on_path[a] = false;
        for(int a : A) helper[a] = false;

        return res;
    }

    void writeRemainingGraph(VVI &V) {
        clog << "Graph V:" << endl;
        for(int i=0; i<V.size(); i++){
            if(!V[i].empty()){
                clog << i << " --> ";
                for(int d : V[i]) clog << d << " ";
                clog << endl;
            }
        }
    }

    bool isFVS(VVI &V, VI S) {
        int N = V.size();
        VB was(N,false), on_path(N,false);
        for(int d : S) was[d] = true;

        function<bool(int)> dfs = [&](int a){
            on_path[a] = true;
            was[a] = true;

            for(int d : V[a]){
                if(on_path[d]) return true;

                if(!was[d]){
                    if( dfs(d) ) return true;
                }
            }

            on_path[a] = false;
            return false;
        };

        for( int i=0; i<N; i++ ) if(!was[i] && dfs(i)) return false;
        return true;
    }

    bool isFVS2(VVI &V, VI &S) {
        VI nodes = CombinatoricUtils::getFullSetDifference(V.size(), S);
        InducedGraph g = GraphInducer::induce( V, nodes );
        return hasCycle(g.V) == false;
    }


    VI getRedundantNodes(VVI V, VI dfvs, int max_check_length, const bool minimize, const bool return_on_first) {
        int N = V.size();
        VB in_dfvs(N,false);
        for(int d : dfvs) in_dfvs[d] = true;

        VB is_end(N,false);
        VB is_source(N,false);

        VVI revV = GraphUtils::reverseGraph(V);

        //*************** AUXILIARY GRAPH H FOR SPEEDUP - BEG SECTION
        VVI H;
        VVI revH;
        VI order;
        VI in_order(N,-1);
        VI temp; temp.reserve(N);
        //*************** AUXILIARY GRAPH H FOR SPEEDUP - END SECTION

        /**
         * Adds node u to the auxiliary graph H
         */
        auto updateAuxiliaryGraph = [&](int i){
            for( int d : V[i] ){
                if( !in_dfvs[d] ){
                    H[i].push_back(d);
                    revH[d].push_back(i);
                }
            }

            for( int d : revV[i] ){
                if( !in_dfvs[d] ){
                    H[d].push_back(i);
                    revH[i].push_back(d);
                }
            }
        };

        /**
         * Creates a graph that preserves only edges in the topological order - no edges incident to nodes in DFVS
         */
        auto createAuxiliaryGraph = [&](){
            H = VVI(N);
            revH = VVI(N);
            for(int i=0; i<N; i++){
                if( !in_dfvs[i] ){
                    for( int d : V[i] ){
                        if( !in_dfvs[d] ){
                            H[i].push_back(d);
                            revH[d].push_back(i);
                        }
                    }
                }
            }
        };

        /**
        * Checks if [order] corresponds to nodes sorted topologically. Just for debugging
        */
        auto checkOrder = [&](){
            for( int u : order ){
                for( int d : H[u] ) if( in_order[d] <= in_order[u] ) return false;
            }
            return true;
        };

        auto createOrder = [&](){
            VI nodes; nodes.reserve(N);
            for(int i=0; i<N; i++) if(!in_dfvs[i]) nodes.push_back(i);
            InducedGraph g = GraphInducer::induce(V, nodes);
            TopoSort ts(g.V);
            order = ts.sortTopologically();
            for(int & d : order) d = g.nodes[d];
            for( int i=0; i<order.size(); i++ ) in_order[order[i]] = i;
        };

        /**
         * Inserts node u to topological order [order], then makes necessary shifts to preserve topological order.
         */
        function<void(int)> insertIntoOrder = [&]( int u ){
            int m = 1e9, M = -1;
            for( int d : H[u] ) m = min(m, in_order[d]);
            for( int d : revH[u] ) M = max(M, in_order[d]);

            if(m == M){
                DEBUG(hasLoop(V,u));
                DEBUG(PII(m,M));
                assert(m != M);
            }

            if( H[u].empty() ){
                in_order[u] = order.size();
                order.push_back(u);
                return;
            }else if( revH[u].empty() || m > M ){
                for( int i=m; i<order.size(); i++ ) in_order[order[i]]++;
                order.insert(order.begin()+m,u);
                in_order[u] = m;
                return;
            }else{
                VB reachable(N,false);
                for( int d : revH[u] ) reachable[d] = true;

                {
                    temp.clear();
                    temp += revH[u];
                    for( int i=0; i<temp.size(); i++ ){
                        int v = temp[i];
                        for( int d : revH[v] ){
                            if( in_order[d] < m ) continue;
                            if(!reachable[d]){
                                reachable[d] = true;
                                temp.push_back(d);
                            }
                        }
                    }
                    for(int i=0; i<m; i++) reachable[order[i]] = true;
                }


                for(int d : order) in_order[d] = -1; // clearing order

                VI new_order; new_order.reserve(N);
                for( int i=0; i<=M; i++ ){
                    int v = order[i];
                    if( reachable[v] ) new_order.push_back(v);
                }

                for( int i=m; i<=M; i++ ){
                    int v = order[i];
                    if( !reachable[v] ) new_order.push_back(v);
                }

                for( int i=M+1; i<order.size(); i++ ) new_order.push_back(order[i]);
                swap(order,new_order);
                for(int i=0; i<order.size(); i++) in_order[order[i]] = i;

                insertIntoOrder(u);
            }
        };



        /**
         * Runs simultaneously two BFSs from sources and ends (only in the graph induced by nodes
         * order[min_source_order_ind], ..., order[max_end_order_ind])
         */
        auto hasCycleFrom = [&](VI &sources, VI &ends, int min_source_order_ind, int max_end_order_ind ){
            if( sources.empty() || ends.empty() ) return false;

            for(int d : sources) is_source[d] = true;
            for(int d : ends) is_end[d] = true;

            for(int d : sources) if( is_end[d] ) return true;

            int p1 = 0, p2 = sources.size();
            int q1 = 0, q2 = ends.size();
            while( p1 < sources.size() || q1 < ends.size() ){

                for( int i = p1; i<p2; i++ ){ // expand sources
                    int u = sources[i];
                    for( int d : H[u] ){
                        if( is_end[d] ) return true;
                        else if( !is_source[d] && in_order[d] <= max_end_order_ind ){
                            is_source[d] = true;
                            sources.push_back(d);
                        }
                    }
                }
                p1 = p2;
                p2 = sources.size();

                for( int i=q1; i<q2; i++ ){ // expand ends
                    int u = ends[i];
                    for( int d : revH[u] ){
                        if( is_source[d] ) return true;
                        else if( !is_end[d] && in_order[d] >= min_source_order_ind ){
                            is_end[d] = true;
                            ends.push_back(d);
                        }
                    }
                }
                q1 = q2;
                q2 = ends.size();
            }

            return false;
        };

        VI res;
        VI sources, ends;
        sources.reserve(N); ends.reserve(N);
        createAuxiliaryGraph();
        createOrder();

        for( int i=0; i<dfvs.size(); i++ ){
            if( i >= max_check_length && i < (int)dfvs.size() - max_check_length ) continue;

            int u = dfvs[i];

            sources.clear(); ends.clear();
            for( int d : V[u] ) if(!in_dfvs[d]) sources.push_back(d);
            for( int d : revV[u] ) if(!in_dfvs[d]) ends.push_back(d);

            int min_source_order_ind = 1e9, max_end_order_ind = -1;
            for( int d : sources ) min_source_order_ind = min(min_source_order_ind, in_order[d]);
            for( int d : ends ) max_end_order_ind = max(max_end_order_ind, in_order[d]);

            if( !hasCycleFrom(sources, ends, min_source_order_ind, max_end_order_ind ) ){
                res.push_back(u);
                if(return_on_first) return res;

                if(minimize) {
                    in_dfvs[u] = false;
                    updateAuxiliaryGraph(u);
                    insertIntoOrder(u);
                }
            }

            for(int d : sources) is_source[d] = false;
            for(int d : ends) is_end[d] = false;
        }

        return res;
    }

    VI removeRedundantNodes(VI &dfvs, VI &redundant, VB &helper) {
        for(int d : redundant) helper[d] = true;
        VI res;
        res.reserve(dfvs.size());
        for( int d : dfvs ) if(!helper[d]) res.push_back(d);
        for(int d : redundant) helper[d] = false;
        return res;
    }

    VI findAndRemoveRedundantNodes(VVI &V, VI &dfvs, int max_check_length, const bool minimize) {
        VI redundant = getRedundantNodes(V,dfvs, max_check_length, minimize);
        VB helper(V.size());
        return removeRedundantNodes(dfvs,redundant,helper);
    }

    tuple<VVI, VPII> transformDFAStoDFVS(VVI &V) {
        VPII edges = GraphUtils::getDirectedGraphEdges(V);
        int E = edges.size();
        VVI L(E);
        int N = V.size();
        vector<unordered_map<int,int>> rev_mapa(N);

        for(int i=0; i<E; i++) rev_mapa[edges[i].first][edges[i].second] = i;

        for(int i=0; i<E; i++){
            int a = edges[i].second;
            for( int b : V[a] ) L[i].push_back( rev_mapa[a][b] );
        }

        return make_tuple(L, edges);
    }

    int countBackEdgesForNodeOrder(VVI &V, VI &ord) {
        int N = V.size();
        VI ind_in_ord(N,-1);
        for(int i=0; i<ord.size(); i++) ind_in_ord[ord[i]] = i;
        int back_edges = 0;
        for(int u : ord) for(int d : V[u]) if(ind_in_ord[d] < ind_in_ord[u]) back_edges++;
        return back_edges;
    }

    VVI readGraph(istream &str) {
        string s;
        int N = -1, M = -1;

        int id = 1;
        VVI V;

        while(getline(str,s, '\n')){
            if(s.empty()){
                if(N>=0) id++;
                continue;
            }

            if(s[0] == '%')continue;
            auto l = StandardUtils::split(s, " ");

            if(N == -1){
                N = stoi(l[0]);
                M = stoi(l[1]);
                V = VVI(N);
            }else{
                int a = id;

                for( string s : l ){
                    int b = stoi(s);
                    assert(b>=1 && b <= N);

                    V[a-1].push_back(b-1);
                }

                id++;
            }
        }

        return V;
    }

    VVI getCycleAdjacencyGraph(VVI &cycles) {
        int N = cycles.size();
        VVI V(N);

        int M = -1;
        for(VI & c : cycles) M = max(M, *max_element(ALL(c)));
        VVI in_cycles(M+1);
        for(int i=0; i<N; i++) for( int d : cycles[i] ) in_cycles[d].push_back(i);

        VB helper(N,false);

        VI adj_cycles;
        adj_cycles.reserve(N);
        for( int i=0; i<N; i++ ){
            adj_cycles.clear();
            for( int v : cycles[i] ){
                for( int c : in_cycles[v] ){
                    if(!helper[c] && c != i){
                        helper[c] = true;
                        adj_cycles.push_back(c);
                    }
                }
            }

            V[i] = adj_cycles;
            for(int c : adj_cycles) helper[c] = false;
        }

        return V;
    }

    VVI getSmallLowerBoundCycles(VVI V) {
        VVI res;
        int N = V.size();
        VB helper(N,false);
        VVI revV = GraphUtils::reverseGraph(V);
        VPII pi_edges = getAllPIEdges(V, revV, helper);
        for( auto & [a,b] : pi_edges ) if(a<b) res.push_back( VI({a,b}) );

        removeEdges(V, pi_edges, helper);

        VVI triangles = CycleCounter::getAllTrianglesInDirectedGraph(V);
        res += triangles;

        return res;
    }

    VI getLowerBoundCycles(VVI cycV, int milliseconds) {
        VVI & V = cycV;
        int N = cycV.size();

        VI mis;

        const bool use_exact_mis = ( milliseconds == 1e9 );

        VI kern_vc;
        const bool use_kernelization = (!use_exact_mis);
        if(use_kernelization) {
            KernelizerVC kern;
            auto[x, y] = kern.initialKernelization(V);
            kern_vc = x;
            int E = GraphUtils::countEdges(V);
            if(E == 0) return CombinatoricUtils::getFullSetDifference(V.size(), kern_vc);
        }


        if(!use_exact_mis) {
            bool use_fastvc = ( N > 20'000 );

            VI init_vc; // empty
            VPII edges = GraphUtils::getGraphEdges(V);
            UniformIntGenerator rnd(0,1e6);

            if (use_fastvc) {
                int seed = rnd.rand();

                libmvc::FastVC fastvc(edges, N, 0, std::chrono::milliseconds(milliseconds), false, seed, init_vc);
                fastvc.cover_LS();
                mis = fastvc.get_independent_set(false);

                assert(VCUtils::isIndependentSet(V, mis));
            } else { // using NuMVC
                NuMVC numvc;
                mis = numvc.solve(V, 1.0 * milliseconds / 1000, init_vc); // CAUTION! NuMVC indices are from 1 to N
                for (int &d : mis) d--;

                assert(VCUtils::isIndependentSet(V, mis));
            }
        }else{
            VI vc = getLowerBoundByVCOnPIGraph(V,V);
            mis = CombinatoricUtils::getFullSetDifference(V.size(),vc);

            assert(VCUtils::isIndependentSet(V, mis));
        }

        if(use_kernelization) {
            VI vc = CombinatoricUtils::getFullSetDifference(V.size(), mis);
            vc += kern_vc;
            mis = CombinatoricUtils::getFullSetDifference(V.size(), vc);
        }

        return mis;
    }

    VI getLowerBoundByVCOnPIGraph(VVI &V, VVI &revV) {
        int N = V.size();
        VB helper(N,false);

        VVI V2 = V;
        for(int i=0; i<N; i++) V2[i].clear();

        for( int i=0; i<N; i++ ){
            for( int d : V[i] ) helper[d] = true;
            for( int d : revV[i] ) if(helper[d]) V2[i].push_back(d);
            for( int d : V[i] ) helper[d] = false;
        }

        const bool use_wgyc = true;
        VI vc;

        if(use_wgyc) {
            ofstream out("subgraph_wgyc.txt");
            GraphWriter::writeGraphDIMACS(V2, out, false, 1);
            out.close();

            TimeMeasurer::start("Running WGYC");
            auto suppress = system("./vc_solver subgraph_wgyc.txt > temp_wgyc_bin_log_file.txt");
            TimeMeasurer::stop("Running WGYC");

            ifstream wgyc("subgraph_wgyc.txt.vc");
            string s;
            int n, sz;
            wgyc >> s >> s >> n >> sz;

            for (int i = 0; i < sz; i++) {
                wgyc >> n;
                vc.push_back(n - 1);
            }
            wgyc.close();
        }else{
            ofstream out( "subgraph_peaty.txt" );
            GraphWriter::writeGraphDIMACS(V2, out, false, 1);
            out.close();

            TimeMeasurer::start("Running Peaty");
            auto suppress = system("./peaty < subgraph_peaty.txt > subgraph_peaty.vc");
            TimeMeasurer::stop("Running Peaty");

            ifstream in("subgraph_peaty.vc");
            string s;
            while( getline(in, s) ){
                if( s[0] == 'c' || s[0] == 's' ) continue;
                else{
                    vc.push_back( stoi(s) - 1 );
                }
            }
        }

        return vc;
    }

    VI getUpperBoundByVCOnSuperPIGraph(VVI V, int milliseconds, bool use_reductions, bool use_extensive_reductions) {
        int N = V.size();
        VVI H = getSuperPIGraph(V);

        bool use_numvc = ( N <= 20'000 );
        bool use_fastvc = (!use_numvc);
        VI init_vc, mis;

        VI kern_vc;
        vector<DFVSReduction*> reductions;
        const bool use_kernelization = use_reductions;
        const bool use_reducer = use_extensive_reductions;

        if(use_kernelization) {
            if(use_reducer){
                Config cnf;
                cnf.disableAllRecursiveReductions();
                Reducer red(H, cnf);
                reductions = red.reduce();
                H = red.V;
                int E = GraphUtils::countEdges(H);
                if(E == 0){
                    kern_vc.clear();
                    Reducer::liftSolution(V.size(), kern_vc,reductions);
                    assert(VCUtils::isVertexCover(V, kern_vc));
                    return kern_vc;
                }
            }
            else {
                KernelizerVC kern;
                auto[x, y] = kern.initialKernelization(H);
                kern_vc = x;
                int E = GraphUtils::countEdges(H);
                if(E == 0){
                    assert(VCUtils::isVertexCover(V, kern_vc));
                    return kern_vc;
                }
            }
        }

        if(use_fastvc){
            UniformIntGenerator rnd(0,1e6);
            int seed = rnd.rand();
            VPII edges = GraphUtils::getGraphEdges(H);
            libmvc::FastVC fastvc(edges, N, 0, std::chrono::milliseconds(milliseconds), false, seed, init_vc);
            fastvc.cover_LS();

            mis = fastvc.get_independent_set(false);
            assert(VCUtils::isIndependentSet(H,mis));
        }
        else{
            NuMVC numvc;
            mis = numvc.solve(H, 1.0 * milliseconds / 1000, init_vc); // CAUTION! NuMVC indices are from 1 to N
            for (int &d : mis) d--;
            assert(VCUtils::isIndependentSet(H, mis));
        }

        if(use_kernelization) {
            VI vc = CombinatoricUtils::getFullSetDifference(V.size(), mis);
            if(!use_reducer) vc += kern_vc;
            else Reducer::liftSolution(V.size(), vc, reductions);
            mis = CombinatoricUtils::getFullSetDifference(V.size(), vc);
        }

        VI vc = CombinatoricUtils::getFullSetDifference(N,mis);
        assert(VCUtils::isVertexCover(V, vc));
        return vc;
    }

    bool isPIGraph(VVI V) {
        int N = V.size();

        VVI revV = GraphUtils::reverseGraph(V);
        V = GraphUtils::reverseGraph(revV);

        for( int i=0; i<N; i++ ) if(!equal(ALL(V[i]), ALL(revV[i]))) return false;

        return true;
    }

    bool isPIGraph(VVI &V, VVI & revV, VB & helper){
        bool is = true;
        for(int i=0; i<V.size(); i++){
            for(int d : V[i]) helper[d] = true;
            for( int d : revV[i] ) if( !helper[d] ) is = false;
            for(int d : V[i]) helper[d] = false;
            if(!is) return false;
        }
        return true;
    }

    VVI getUnderlyingPIGraph(VVI V) {
        int N = V.size();
        VVI H(N);

        VVI revV = GraphUtils::reverseGraph(V);
        VB helper(N,false);
        VPII pi_arcs = getAllPIEdges(V, revV, helper);
        for(auto & [a,b] : pi_arcs) H[a].push_back(b);

        return H;
    }

    VVI getUnderlyingPIGraph(VVI &V, VVI & revV) {
        int N = V.size();
        VVI H(N);

        VB helper(N,false);
        VPII pi_arcs = getAllPIEdges(V, revV, helper);
        for(auto & [a,b] : pi_arcs) H[a].push_back(b);

        return H;
    }

    VVI getSuperPIGraph(VVI V) {
        int N = V.size();
        VVI revV = GraphUtils::reverseGraph(V);
        VB helper(N,false);
        VPII arcs = GraphUtils::getDirectedGraphEdges(V);

        for( auto & [a,b] : arcs ) if(a>b) swap(a,b);
        sort(ALL(arcs));
        arcs.resize( unique(ALL(arcs)) - arcs.begin() );

        VVI H(N);
        for( int i=0; i<N; i++ ) H[i].reserve(V[i].size());
        for( auto & [a,b] : arcs ){
            H[a].push_back(b);
            H[b].push_back(a);
        }

        return H;
    }

    VVI getNonPIGraph(VVI V) {
        VVI revV = GraphUtils::reverseGraph(V);
        int N = V.size();
        VVI H(N);
        VB helper(N,false);
        for(int i=0; i<N; i++){
            for(int d : revV[i]) helper[d] = true;
            for(int d : V[i]) if( !helper[d] ) H[i].push_back(d);
            for(int d : revV[i]) helper[d] = false;
        }
        return H;
    }

    VVI getNonPIGraph(VVI &V, VVI &revV) {
        int N = V.size();
        VVI H(N);
        VB helper(N,false);
        for(int i=0; i<N; i++){
            for(int d : revV[i]) helper[d] = true;
            for(int d : V[i]) if( !helper[d] ) H[i].push_back(d);
            for(int d : revV[i]) helper[d] = false;
        }
        return H;
    }


    void contractNodeToNode(VVI &V, VVI &revV, int b, int c, VB & helper) {

        VPII arcs_to_add;

        {
            VI Nb = V[b];
            VI Nc = V[c];

            sort(ALL(Nb));
            sort(ALL(Nc));

            VI diff;
            set_difference(ALL(Nb), ALL(Nc), back_inserter(diff));

            for( int d : diff ) if(d != c) arcs_to_add.push_back( {c,d} );
        }

        {
            VI Nb = revV[b];
            VI Nc = revV[c];

            sort(ALL(Nb));
            sort(ALL(Nc));

            VI diff;
            set_difference(ALL(Nb), ALL(Nc), back_inserter(diff));

            for( int d : diff ) {
                if (d != c) {
                    arcs_to_add.push_back({d, c});
                    assert(d != b);
                }
            }
        }

        removeNode(V, revV, b, helper);

        addEdges(V, arcs_to_add, helper);
        for( auto & [a,b] : arcs_to_add ) swap(a,b);
        addEdges(revV, arcs_to_add, helper);
    }

    bool hasPathFromTo(VVI &V, VVI &revV, VB &in_V, VI sources, VI ends, VB &is_source, VB &is_end) {
        const bool debug = false;

        if( sources.empty() || ends.empty() ) return false;

        for( int d : sources ) assert( in_V[d] );
        for( int d : ends ) assert( in_V[d] );

        for(int d : sources) is_source[d] = true;
        for(int d : ends) is_end[d] = true;

        bool has_path = false;

        for(int d : sources) if( is_end[d] ) has_path = true;

        if(debug) DEBUG(has_path);

        if(!has_path) {
            if(debug) clog << "Entering while" << endl;

            int p1 = 0, p2 = sources.size();
            int q1 = 0, q2 = ends.size();

            while (p1 < sources.size() || q1 < ends.size()) {

                for(int i = p1; i < p2; i++) {
                    int u = sources[i];
                    for( int d : V[u]) {
                        if( !in_V[d] ) continue;

                        if(is_end[d]){
                            has_path = true;
                            break;
                        }
                        else if(!is_source[d]) {
                            is_source[d] = true;
                            sources.push_back(d);
                        }
                    }

                    if(has_path) break;
                }
                if(has_path) break;

                p1 = p2;
                p2 = sources.size();

                for(int i = q1; i < q2; i++) { // expand ends
                    int u = ends[i];
                    for(int d : revV[u]) {
                        if( !in_V[d] ) continue;

                        if(is_source[d]){
                            has_path = true;
                            break;
                        }
                        else if(!is_end[d]) {
                            is_end[d] = true;
                            ends.push_back(d);
                        }
                    }
                    if(has_path) break;
                }
                if(has_path) break;
                q1 = q2;
                q2 = ends.size();
            }
        }

        for(int d : sources) is_source[d] = false;
        for(int d : ends) is_end[d] = false;

        if(debug){
            DEBUG(has_path);
            ENDL(2);
        }

        return has_path;
    }

    bool containsSimpleCycleWithNode(VVI &V, VVI &revV, VVI& nonpiV, VVI & revnonpiV,
                                      VB &in_V, VI& marker, VB& is_end, int v, const int max_time_millis,
                                      int max_branch_depth) {
        int depth = 0;

        auto start_time = chrono::steady_clock::now();
        auto getMillisFromStart = [&](){
            return chrono::duration<double, std::milli >(chrono::steady_clock::now() - start_time ).count();
        };

        bool branch_depth_exceeded = false;

        function<bool(int,int)> is_in_simple_cycle = [&](int num, int par){
            if(depth == 0) assert(!branch_depth_exceeded);
            if(branch_depth_exceeded) return true;

            if( nonpiV[num].size() >= 2 ) depth++;

            bool found_cycle = false;
            for( int d : nonpiV[num] ){
                if( !in_V[d] ) continue;
                if( marker[d] == 1 && is_end[d] ){
                    found_cycle = true;
                    break;
                }
            }

            for( int d : V[num] ) marker[d]++;
            for( int d : revV[num] ) marker[d]++;

            if( depth > max_branch_depth ){
                branch_depth_exceeded = true;
                found_cycle = true;
            }

            if( !found_cycle ){
                for( int d : nonpiV[num] ){
                    if( !in_V[d] ) continue;

                    if( marker[d] == 1 ){
                        if((depth % 4 == 0) && getMillisFromStart() > max_time_millis ){
                            found_cycle = true;
                            break;
                        }

                        if( is_in_simple_cycle( d, num ) ){
                            found_cycle = true;
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

        bool has_simple_cycle = false;
        { // here we check if there exists a pi-edge containing node v
            for (int d : revV[v]) is_end[d] = true;
            for (int d : V[v]) if (in_V[d] && is_end[d]) has_simple_cycle = true;
            for (int d : revV[v]) is_end[d] = false;
            if (has_simple_cycle) return true;
        }

        { // if there is no pi-edge containing v, then we can use DFS traversing using nonpiV
            for (int d : revnonpiV[v]) if (in_V[d]) is_end[d] = true;
            branch_depth_exceeded = false;
            assert(depth == 0);
            has_simple_cycle = is_in_simple_cycle(v, v);
            for (int d : revnonpiV[v]) if (in_V[d]) is_end[d] = false;
        }

        return has_simple_cycle;
    }


    VI getAllNoninclusionNodesForNode(VVI &V, VVI &revV, VVI &nonpiV, VVI &revnonpiV, VB &in_V, VI &marker, VB &is_end,
                            VB & in_Np, int a, VI bs, int max_time_millis_total, int max_time_millis_per_node,
                            int max_branch_depth) {

        int depth = 0;

        auto start_time_node = chrono::steady_clock::now();
        auto getMillisFromStart = [&](){
            return chrono::duration<double, std::milli >(chrono::steady_clock::now() - start_time_node ).count();
        };

        bool branch_depth_exceeded = false;

        function<bool(int,int)> is_in_prohibited_simple_cycle = [&](int num, int par){
            if(depth == 0) assert(!branch_depth_exceeded);
            if(branch_depth_exceeded) return true;

            if( in_Np[num] ) for( int d : revV[a] ) marker[d]++;
            if(nonpiV[num].size() >= 2) depth++;

            bool found_cycle = false;
            for( int d : nonpiV[num] ){
                if( !in_V[d] ) continue;
                if( marker[d] == 1 && is_end[d] ){
                    found_cycle = true;
                    break;
                }
            }

            for( int d : V[num] ) marker[d]++;
            for( int d : revV[num] ) marker[d]++;

            if(depth > max_branch_depth){
                branch_depth_exceeded = true;
                found_cycle = true;
            }

            if( !found_cycle ){
                for( int d : nonpiV[num] ){
                    if( !in_V[d] ) continue;

                    if( marker[d] == 1 ){
                        if((depth % 4 == 0) && getMillisFromStart() > max_time_millis_per_node ){
                            found_cycle = true;
                            break;
                        }

                        if( is_in_prohibited_simple_cycle(d, num ) ){
                            found_cycle = true;
                            break;
                        }
                    }
                }
            }

            for( int d : V[num] ) marker[d]--;
            for( int d : revV[num] ) marker[d]--;
            if( in_Np[num] ) for( int d : revV[a] ) marker[d]--;
            if(nonpiV[num].size() >= 2) depth--;

            if(found_cycle) return true;

            return false;
        };

        auto start_time_total = chrono::steady_clock::now();
        VI res;

        for(int d : V[a]) in_Np[d] = true;

        for(int b : bs) {
            if(b == a) continue;

            bool old_inv = in_V[b];
            in_V[b] = true;

            bool old_inNp = in_Np[b];
            in_Np[b] = false;

            bool is_in_bad_cycle = false;
            {
                for (int d : revV[b]) is_end[d] = true;
                for (int d : V[b]) if (in_V[d] && is_end[d]) is_in_bad_cycle = true;
                for (int d : revV[b]) is_end[d] = false;
            }

            if(!is_in_bad_cycle){
                for (int d : revnonpiV[b]) if (in_V[d]) is_end[d] = true;
                start_time_node = chrono::steady_clock::now();
                branch_depth_exceeded = false;
                assert(depth == 0);
                is_in_bad_cycle = is_in_prohibited_simple_cycle(b, b);
                for (int d : revnonpiV[b]) if (in_V[d]) is_end[d] = false;
            }

            if(!is_in_bad_cycle) res.push_back(b);

            in_V[b] = old_inv;
            in_Np[b] = old_inNp;

            if( chrono::duration<double, std::milli >(chrono::steady_clock::now() - start_time_total ).count()
                > max_time_millis_total ) break;
        }

        for(int d : V[a]) in_Np[d] = false;

        return res;
    }

    bool isPiNode(VVI &V, VVI &revV, VVI &piV, int a) {
        if( piV[a].size() != V[a].size() || piV[a].size() != revV[a].size() ) return false;
        else return true;
    }

    VVI getAllSimpleCycles(VVI &V, VVI &revV, VVI &nonpiV, VVI &revnonpiV, VI A, VB &in_V, VI &marker, VB &is_end,
                           int max_cycle_length) {
        for(int d : A) in_V[d] = true;
        int cycle_length = 0;

        int N = V.size();
        VLL hashes = VLL(N);
        UniformIntGenerator rnd(0,1'000'000'000ll * 1'000'000'000);
        for(int i=0; i<N; i++) hashes[i] = rnd.rand();

        /**
         * Stores only one copy of each cycle by keeping the cycles hash.
         */
        unordered_map<LL, VI> found_cycles;

        function<void(int,VI&, LL)> findSimpleCycles = [&](int num, VI & C, LL hash){
            if( cycle_length+1 >= max_cycle_length) return;

            hash ^= hashes[num];
            cycle_length++;
            C.push_back(num);

            for( int d : nonpiV[num] ){
                if( !in_V[d] ) continue;
                if( marker[d] == 1 && is_end[d] ){
                    C.push_back(d);
                    found_cycles[hash^hashes[d]] = C;
                    C.pop_back();
                }
            }

            for( int d : V[num] ) marker[d]++;
            for( int d : revV[num] ) marker[d]++;

            for( int d : nonpiV[num] ){
                if( !in_V[d] ) continue;
                if( marker[d] == 1 && !is_end[d] ) findSimpleCycles(d, C, hash);
            }

            for( int d : V[num] ) marker[d]--;
            for( int d : revV[num] ) marker[d]--;

            cycle_length--;
            C.pop_back();
        };


        for( int v : A) {
            if(!in_V[v]) continue;

            {
                for (int d : revV[v]) is_end[d] = true;
                for (int d : V[v]) if( in_V[d] && is_end[d]) found_cycles[ hashes[v] ^ hashes[d] ] = {v,d};
                for (int d : revV[v]) is_end[d] = false;
            }

            {
                for (int d : revnonpiV[v]) if (in_V[d]) is_end[d] = true;
                VI C;
                assert(cycle_length == 0);
                findSimpleCycles( v, C, 0 );
                for (int d : revnonpiV[v]) if (in_V[d]) is_end[d] = false;
            }
        }

        for(int d : A) in_V[d] = false;
        VVI cycles;
        for( auto & [h,v] : found_cycles ) cycles.push_back(v);
        return cycles;
    }

    VVI getAllSimpleCycles2(VVI &V, VVI &revV, VVI &nonpiV, VVI &revnonpiV, VI A, VB &in_V, VI &marker, VB &is_end,
                           int max_cycle_length, int millis) {
        if(max_cycle_length < 2) return {};
        int cycle_length = 0;
        int N = V.size();
        VVI found_cycles;

        auto start_time = chrono::steady_clock::now();
        auto getMillisFromStart = [&](){
            return chrono::duration<double, std::milli >(chrono::steady_clock::now() - start_time ).count();
        };

        function<void(int,VI&)> findSimpleCycles = [&](int num, VI & C){
            if( cycle_length+1 >= max_cycle_length) return;
            if(getMillisFromStart() > millis) return;

            cycle_length++;
            C.push_back(num);

            for( int d : nonpiV[num] ){
                if( !in_V[d] ) continue;
                if( marker[d] == 1 && is_end[d] ){
                    // found cycle   C+d   of length cycle_length+1
                    C.push_back(d);
                    found_cycles.push_back(C);
                    C.pop_back();
                }
            }

            for( int d : V[num] ) marker[d]++;
            for( int d : revV[num] ) marker[d]++;

            for( int d : nonpiV[num] ){
                if( !in_V[d] ) continue;
                if( marker[d] == 1 && !is_end[d] ) findSimpleCycles(d, C);
            }

            for( int d : V[num] ) marker[d]--;
            for( int d : revV[num] ) marker[d]--;

            cycle_length--;
            C.pop_back();
        };


        for( int v : A) {
            in_V[v] = true;

            { // here we check if there exists a cycle of length 2
                for (int d : revV[v]) is_end[d] = true;
                for (int d : V[v]) if( in_V[d] && is_end[d]) found_cycles.push_back({v,d});
                for (int d : revV[v]) is_end[d] = false;
            }

            { // if there is no pi-edge containing v, then we can use DFS traversing using nonpiV
                for (int d : revnonpiV[v]) if (in_V[d]) is_end[d] = true;
                VI C;
                assert(cycle_length == 0);
                findSimpleCycles( v, C );
                for (int d : revnonpiV[v]) if (in_V[d]) is_end[d] = false;
            }
        }

        for(int d : A) in_V[d] = false;

        return found_cycles;
    }

    VVI getAllSimpleCycles3(VVI V, int max_cycle_length, int millis, bool random_A) {
        StronglyConnectedComponents scc(V);
        scc.createStronglyConnectedComponents();
        auto comps = scc.getComponents();
        VI in_comp = StandardUtils::layersToPartition(comps);
        VPII arcs_to_remove;
        VPII all_arcs = GraphUtils::getGraphEdges(V, true);
        for( auto & [a,b] : all_arcs ) if( in_comp[a] != in_comp[b] ) arcs_to_remove.emplace_back(a,b);

        int N = V.size();
        VB helper(N,false);
        removeEdges( V, arcs_to_remove, helper );


        VB in_V(N,false);
        VB is_end(N,false);
        VI marker(N,0);
        VI A;
        if(random_A) A = CombinatoricUtils::getRandomPermutation(N);
        else{ A = VI(N); iota(ALL(A),0); }

        VVI revV = GraphUtils::reverseGraph(V);
        VVI nonpiV = Utils::getNonPIGraph(V);
        VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);

        return getAllSimpleCycles2( V, revV, nonpiV, revnonpiV, A, in_V, marker, is_end, max_cycle_length, millis );
    }

    VI findMinHittingSet(VVI &cycles, VI upper_bound_solution, int lower_bound) {
        { // checking if we have the vertex-cover instance
            bool vc_instance = true;
            for(auto & c : cycles) if(c.size() != 2) vc_instance = false;
            if(vc_instance) {
                int N = -1;
                for( auto & c : cycles ) for(int d : c) N = max(N,d);
                N++;
                VVI V(N);
                for( auto & c : cycles ){
                    int a = c[0], b = c[1];
                    GraphUtils::addEdge(V,a,b,false);
                }
                return getLowerBoundByVCOnPIGraph(V, V);
            }
        }

        string infile = "hs_in.txt";
        ofstream str(infile.c_str());
        int N = 0;
        for( VI & C : cycles ) for( int d : C ) N = max(d,N);
        str << N+1 << " " << cycles.size() << endl;
        for( VI & C : cycles ){
            str << C.size();
            for( int d : C ) str << " " << d;
            str << endl;
        }
        str.close();


        string settings_file = "hs_settings.txt";
        str.open( settings_file.c_str() );
        ostringstream sett_str;

        sett_str << "{\n" <<
                          "    \"enable_local_search\": true,\n" << // original
//                          "    \"enable_local_search\": false,\n" << // #TEST
                          "    \"enable_max_degree_bound\": true,\n" <<
                          "    \"enable_sum_degree_bound\": false,\n" <<
                          "    \"enable_efficiency_bound\": true,\n" <<
                          "    \"enable_packing_bound\": true,\n" <<
                          "    \"enable_sum_over_packing_bound\": true,\n" <<
//                          "    \"packing_from_scratch_limit\": 5,\n" << // original
                          "    \"packing_from_scratch_limit\": 3,\n" << // #TEST
//                          "    \"greedy_mode\": \"AlwaysBeforeBounds\"\n" <<
//                          "    \"greedy_mode\": \"AlwaysBeforeExpensiveReductions\""; // original
                          "    \"greedy_mode\": \"Once\""; // #TEST
        if(!upper_bound_solution.empty()){
            sett_str << ",\n    \"initial_hitting_set\": [";
            int cnt = 0;
            for(int d : upper_bound_solution){
                if(cnt++ > 0) sett_str << ",";
                sett_str << d;
            }
            sett_str << "]";
        }

        if(lower_bound > 0){
            sett_str << ",\n    \"stop_at\": " << lower_bound;
        }

        sett_str << "\n}";
        string settings = sett_str.str();

        str << settings << endl;
        str.close();

        string report_file = "hs_report.txt";
        string out_file = "hs_out.txt";
        string log_file = "hs_log.txt";
        string command = "./findminhs solve " + infile + " " + settings_file + " > " + out_file;
        auto suppress = system(command.c_str());


        ifstream istr(out_file.c_str());
        string s;
        vector<string> lines;
        while( istr >> s ) lines.push_back(s);
        istr.close();

        VI res;
        for( string s : lines ) res.push_back( stoi(s) );
        return res;
    }

    VI hsImprovementLS2(VVI sets, VI hs, Config cnf, int iters, int lower_bound, bool deviate_perm,
                        int persistent_search_iterations_left) {
        HittingSetLS hsls(sets, hs, cnf);
        hsls.persistent_search_iterations_left = persistent_search_iterations_left;
        return hsls.hsImprovementLS2(iters, lower_bound, cnf.ihs_hsls_perm_deviation_frequency );
    }

    VI getIntersectionOfAllInducedNonpiCyclesWithNode(VVI &V, VVI &revV, VVI &nonpiV, VVI &revnonpiV, int v,
                                   VI &marker, VB &is_end, VI & cnt_marker, int max_cycle_length, int millis) {
        int cycle_length = 0;
        int N = V.size();

        VI res;
        bool found_all_cycles = true;

        auto start_time_node = chrono::steady_clock::now();
        auto getMillisFromStart = [&](){
            return chrono::duration<double, std::milli >(chrono::steady_clock::now() - start_time_node ).count();
        };

        VI visited;
        function<void(int,VI&)> findAllCyclesIntersection = [&](int num, VI & C){
            if(!found_all_cycles) return;
            if( cycle_length+1 >= max_cycle_length){ found_all_cycles = false; return; }
            if( ((cycle_length & 7) == 0) && getMillisFromStart() > millis){ found_all_cycles = false; return; }

            cycle_length++;
            C.push_back(num);

            for( int d : nonpiV[num] ){
                if( marker[d] == 1 && is_end[d] ){
                    // found cycle   C+d   of length cycle_length+1
                    C.push_back(d);

                    // we found cycle C
                    for( int x : C ){
                        cnt_marker[x]++;
                        visited.push_back(x);
                    }
                    int inters_cnt = 0;
                    res.clear();
                    for( int x : C ){
                        if( cnt_marker[x] == cnt_marker[v] ){
                            res.push_back(x);
                            inters_cnt++;
                        }
                    }

                    // if at some point of recursion the intersection is only {v}, then we can terminate further search
                    if(inters_cnt == 1) found_all_cycles = false;

                    C.pop_back();
                }
            }

            for( int d : V[num] ) marker[d]++;
            for( int d : revV[num] ) marker[d]++;

            for( int d : nonpiV[num] ){
                if( marker[d] == 1 && !is_end[d] ) findAllCyclesIntersection(d, C);
            }

            for( int d : V[num] ) marker[d]--;
            for( int d : revV[num] ) marker[d]--;

            cycle_length--;
            C.pop_back();
        };



        { // if there is no pi-edge containing v, then we can use DFS traversing using nonpiV
            for (int d : revnonpiV[v]) is_end[d] = true;
            VI C;
            assert(cycle_length == 0);
            findAllCyclesIntersection(v, C );
            for (int d : revnonpiV[v]) is_end[d] = false;
        }

        for(int d : visited) cnt_marker[d] = 0;

        if(!found_all_cycles) return {}; // we did not find all cycles - return empty intersection
        else return res;
    }

    double getPieEdgesPercentage(VVI &V) {
        VVI revV = GraphUtils::reverseGraph(V);
        return getPieEdgesPercentage(V, revV);
    }

    double getPieEdgesPercentage(VVI &V, VVI &revV) {
        VB helper(V.size(),false);
        double val = getAllPIEdges(V, revV, helper).size();
        val /= GraphUtils::countEdges(V,true);
        return val;
    }

    LL countInducedCycles(VVI &V, int max_cycle_length, int max_millis) {
        // this function is an almost 1-1 copy from getAllSimlpeCycles2

        Stopwatch sw;
        sw.setLimit("CIC", max_millis);

        int cycle_length = 0;
        int N = V.size();
        LL found_cycles = 0;

        VVI revV = GraphUtils::reverseGraph(V);
        VVI nonpiV = getNonPIGraph(V);
        VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);
        VVI piV = getUnderlyingPIGraph(V);

        VB in_V(N,false), is_end(N,false);
        VI marker(N,0);

        function<void(int,VI&)> findSimpleCycles = [&](int num, VI & C){
            if( cycle_length+1 >= max_cycle_length) return;

            cycle_length++;
            C.push_back(num);

            for( int d : nonpiV[num] ){
                if( !in_V[d] ) continue;
                if( marker[d] == 1 && is_end[d] ){
                    // found cycle   C+d   of length cycle_length+1
                    found_cycles++;
                }
            }

            for( int d : V[num] ) marker[d]++;
            for( int d : revV[num] ) marker[d]++;

            for( int d : nonpiV[num] ){
                if( !in_V[d] ) continue;
                if( marker[d] == 1 && !is_end[d] ) findSimpleCycles(d, C);
            }

            for( int d : V[num] ) marker[d]--;
            for( int d : revV[num] ) marker[d]--;

            cycle_length--;
            C.pop_back();
        };

        VI A(N);
        iota(ALL(A),0);

        for( int v : A) {
            if(sw.tle("CIC")) return 1'000'000'000ll * 1'000'000'000;
            in_V[v] = true;

            { // here we check if there exists a cycle of length 2
                for (int d : revV[v]) is_end[d] = true;
                for (int d : V[v]) if( in_V[d] && is_end[d]) found_cycles++;
                for (int d : revV[v]) is_end[d] = false;
            }

            { // if there is no pi-edge containing v, then we can use DFS traversing using nonpiV
                for (int d : revnonpiV[v]) if (in_V[d]) is_end[d] = true;
                VI C;
                assert(cycle_length == 0);
                findSimpleCycles( v, C );
                for (int d : revnonpiV[v]) if (in_V[d]) is_end[d] = false;
            }
        }

        for(int d : A) in_V[d] = false;

        return found_cycles;
    }

    LL getSetHash(int N, VI &s, int seed) {
        VLL hashes(N);
        UniformIntGenerator rnd(0ll, 1'000'000'000ll * 1'000'000'000, seed);
        for( int i=0; i<N; i++ ) hashes[i] = rnd.rand();
        return accumulate(ALL(s), 0ll, [&](LL h, int b){ return h ^ hashes[b]; });
    }

    int countPiEdges(VVI &V) {
        VVI revV = GraphUtils::reverseGraph(V);
        VB helper(V.size(),false);
        return getAllPIEdges(V, revV, helper).size();
    }

    VI getMinDFVSDefault(VVI &V) {
        Config cnf;
        {
            cnf.sw.setLimit("main", 1e9); // this is just to suppress warning logs
            cnf.sw.start("main");
        }
        DFVSSolverE solver(&V, cnf);
        solver.cnf.write_logs = false;
        solver.cnf.disableAllRecursiveReductions();
        solver.cnf.disableAllConditionalReductions();
        solver.cnf.reducer_use_nonsimple_cycle_arcs_full = false;
        solver.cnf.reducer_use_domination_6 = false;
        solver.cnf.reducer_use_domination_6inserter = false;

        solver.cnf.solverh_use_superpi_vc_ub = false;
        solver.cnf.solverh_min_graph_size_for_improvements = 1e9;
        solver.cnf.vc_improver_milliseconds = 5;
        solver.cnf.solverh_improvement_iterations = 0;

        VI dfvs_e = solver.solveForInputGraph(V);

        return dfvs_e;
    }

    void emergencyExit(VVI &V, VI &dfvs, const bool write_communicate) {
        if( Utils::isFVS(V, dfvs) ) return;

        if(write_communicate) clog << "EMERGENCY EXIT! Something went WRONG, fixing solution..." << endl;
        if(write_communicate) clog << "Wrong solution size: " << dfvs.size() << endl;

        int N = V.size();
        VVI G = V;
        VVI revG = GraphUtils::reverseGraph(G);
        VB helper(N,false);
        removeNodes(G, revG, dfvs, helper);

        InducedGraph g = GraphInducer::induceByNonisolatedNodes(G);

        Config cnf;
        cnf.write_logs = false;
        cnf.disableAllNonbasicReductions();
        Reducer red(g.V, cnf);

        auto reductions = red.reduce();

//        VI vc = getUpperBoundByVCOnSuperPIGraph(red.V, 500);
        DFVSSolverH solverh(cnf);
        solverh.cnf.solverh_min_graph_size_for_improvements = 1e9; // disabling improvements
        solverh.cnf.agent_flow_node_selection_type = Config::agent_flow_remove_largest_flow_node;
        solverh.cnf.agent_flow_max_distance_from_best = 1; // #TEST - this may slow down emergencyExit
        VI new_dfvs = solverh.solveForGraph(red.V);

        Reducer::liftSolution( red.V.size(), new_dfvs, reductions );
        for( int & d : new_dfvs ) d = g.nodes[d];

        // now vc should be a DFVS of G
        dfvs += new_dfvs;

        if(write_communicate) clog << "Fixed solution size: " << dfvs.size() << endl;
        assert(Utils::isFVS(V,dfvs));
    }

    vector<Triple<int>> getLengthOfShortestInducedCycleWithArc(VVI V, int max_length, int millis) {
        if(max_length < 2) return {};
        constexpr bool debug = true;

        int N = V.size();
        VPII arcs = GraphUtils::getGraphEdges(V,true);
        int E = arcs.size();

        /**
         * min_length[i] is minimum length of induced cycle that contains arcs arcs[i]
         */
        VI min_length(E,1e9);

        VVI revV = GraphUtils::reverseGraph(V);
        VVI nonpiV = Utils::getNonPIGraph(V);
        VVI revnonpiV = GraphUtils::reverseGraph(nonpiV);
        VI marker(N,0);
        VB helper(N,false), is_end(N,false);
        int current_length = 0;

        int current_arc_id = 0;

        map<PII,int> arc_to_id;
        for( int i=0; i<E; i++ ) arc_to_id[arcs[i]] = i;
        VI current_cycle;
        int current_length_arcs_determined = 0;
        constexpr bool use_mapping = false;

        function<bool(int,int,int)> findMinInducedCycleLength = [&](int num, int par, int max_l){
            if(current_length + 2 >= max_l ) return false;

            // increasing markers here is OK (unlike in functions in Utils, where it needs to be increased after checking
            // if neighbors are in is_end), because here we start from an edge (a,b)
            // so we start dfs from node b, hence we did not increase markers for revnonpiV[a]
            for( int d : V[num] ) marker[d]++;
            for( int d : revV[num] ) marker[d]++;
            current_length++;
            if(use_mapping) current_cycle.push_back(num);

            bool found_cycle = false;
            for( int d : nonpiV[num] ){
                if( marker[d] == 1 && is_end[d] ){
                    // there is a cycle of length current_length + 2
                    found_cycle = true;

                    // this must hold, because we proceed in an iterative way - there was no cycle of length max_l-1
                    // that contained this arc
//                    if(debug) clog << "Found cycle of length " << max_l << " for arcs " << arcs[current_arc_id]
//                         << ", current_length: " << current_length << endl;
                    assert( current_length+2 == max_l );

                    if(!use_mapping) min_length[current_arc_id] = max_l;

                    if(use_mapping) {
                        current_cycle.push_back(d);
                        assert(current_cycle.size() == max_l);
                        for (int i = 0; i < current_cycle.size(); i++) {
                            int a = current_cycle[i];
                            int b = current_cycle[(i + 1) % current_cycle.size()];
                            int id = arc_to_id[PII(a, b)];
                            if (min_length[id] > max_l) {
                                current_length_arcs_determined++;
                                min_length[id] = max_l;
                            }
                        }
                        current_cycle.pop_back();
                    }

                    break;
                }
            }

            if( !found_cycle ){
                for( int d : nonpiV[num] ){
                    if( marker[d] == 1 ){
                        if( findMinInducedCycleLength( d, num, max_l ) ){
                            found_cycle = true;
                            break;
                        }
                    }
                }
            }

            for( int d : V[num] ) marker[d]--;
            for( int d : revV[num] ) marker[d]--;
            current_length--;
            if(use_mapping) current_cycle.pop_back();

            return found_cycle;
        };

        int total_arcs_determined = 0;

        { // case: length = 2
            VPII pi_arcs = Utils::getAllPIEdges(V, revV, helper);
            set<PII> zb(ALL(pi_arcs));
            for( int i=0; i<E; i++ ){
                if(zb.count(arcs[i])){
                    min_length[i] = 2;
                    current_length_arcs_determined++;
                }
            }

            total_arcs_determined += current_length_arcs_determined;
            if(debug) clog << "There are " << current_length_arcs_determined << " in duced cycles of length 2" << endl;
        }


        for( int l = 3; l <= max_length; l++ ){
            if(debug){
                ENDL(1);
                clog << "Checking length l: " << l << endl;
            }

            int arcs_to_check = 0;
            for( int i=0; i<E; i++ ) if( min_length[i] == 1e9 ) arcs_to_check++;
            if(debug) clog << "There are " << arcs_to_check << " arcs to check for given length" << endl;

            current_length_arcs_determined = 0;
            for( int i=0; i<E; i++ ){
                if( min_length[i] < max_length+1 ) continue;

                int a = arcs[i].first;
                int b = arcs[i].second;

                for( int d : revnonpiV[a] ) is_end[d] = true;
                for( int d : V[a] ) marker[d]++; // for the start node i we mark only out-neighbors

                if(use_mapping) current_cycle = {a};
                assert(current_length == 0);
                current_arc_id = i;
                int cnt = findMinInducedCycleLength( b,a,l );
                if(!use_mapping) current_length_arcs_determined += cnt;

                for( int d : revnonpiV[a] ) is_end[d] = false;
                for( int d : V[a] ) marker[d] = 0;
            }

            total_arcs_determined += current_length_arcs_determined;

            if(debug){
                clog << "total_arcs_determined: " << total_arcs_determined << ",\tcurrent_length_arcs_determined: "
                     << current_length_arcs_determined << endl;
            }
        }

        vector<Triple<int>> res; res.reserve(E);
        for( int i=0; i<E; i++ ) res.emplace_back( arcs[i].first, arcs[i].second, min_length[i] );
        return res;
    }

    bool isPartiallyDominated(VVI &V, VVI &revV, VVI &nonpiV, VVI &revnonpiV, int b, int c, VB &helper) {
        bool dominated = true;
        helper[c] = true;

        for( int d : V[c] ) helper[d] = true;
        for(int d : nonpiV[b]) if(!helper[d]) dominated = false;
        for( int d : V[c] ) helper[d] = false;

        for( int d : revV[c] ) helper[d] = true;
        for(int d : revnonpiV[b]) if(!helper[d]) dominated = false;
        for( int d : revV[c] ) helper[d] = false;

        helper[c] = false;

        return dominated;
    }

    bool isFullyDominated(VVI &V, VVI &revV, int b, int c, VB &helper1, VB &helper2) {
        helper1[c] = helper2[c] = true;
        for( int d : V[c] ) helper1[d] = true;
        for( int d : revV[c] ) helper2[d] = true;

        bool fully_dominated = true;

        for( int d : V[b] ) if(!helper1[d]){ fully_dominated = false; break;}
        if(fully_dominated){
            for( int d : revV[b] ) if(!helper2[d]){ fully_dominated = false; break;}
        }

        helper1[c] = helper2[c] = false;
        for( int d : V[c] ) helper1[d] = false;
        for( int d : revV[c] ) helper2[d] = false;

        return fully_dominated;
    }



}