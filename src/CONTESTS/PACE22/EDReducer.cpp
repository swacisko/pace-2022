//
// Created by sylwe on 07/08/2025.
//

#include "EDReducer.h"

#include "GraphUtils.h"
#include "StandardUtils.h"
#include "CONTESTS/PACE22/Utils.h"


VI EDReducer::reduce(VVI V0) {
    V = V0;

    VI nodes(N);
    iota(ALL(nodes),0);
    if ( cnf.ed_node_sorting_mode == 1 ) sort(ALL(nodes), [&]( int a, int b ) { return V[a].size() > V[b].size(); });
    if ( cnf.ed_node_sorting_mode == 2 ) sort(ALL(nodes), [&]( int a, int b ) { return V[a].size() < V[b].size(); });

    VI reducible_nodes;

    // clog << "EDReducer - reducing graph with " << V0.size() << " nodes and " << GraphUtils::countEdges(V0) << " edges" << endl;

    if (write_logs) DEBUG(V);

    // if (V.size() < 10) { clog << "Considering small graph" << endl; }

    { // clear all data for new call to [reduce]
        last_reduce_nodes_removed = last_reduce_edges_removed = 0;
        last_reduce_inf_rules_1_added = last_reduce_inf_rules_2_created = 0;
        inf_rules_1.clear();
        inf_rules_2.clear();
        all_inf_rules_1_found.clear();
        all_inf_rules_2_found.clear();
    }

    bool changes = true;
    constexpr bool use_exhaustively = false;

    while (changes) {
        changes = false;

        if (cnf.ed_use_node_removal) {
            for (int v : nodes) if (!V[v].empty()) {
                // clog << "\rConsidering node " << v << flush;
                if (consider({v})) {
                    if (write_logs)
                        clog << "\t\tNode " << v << " is ED-reducible!   final W.size(): " << W.size() << endl << endl << endl;
                    // reducible_nodes.push_back(v);
                    // GraphUtils::removeNodeFromGraph(V,v);
                    // last_reduce_nodes_removed++;
                    auto temp = propagateDeg1RuleSlow(v);
                    reducible_nodes += temp;
                    last_reduce_nodes_removed += temp.size();
                    ed_node_applied_cnt++;
                    if (use_exhaustively) changes = true;
                }
                else if ( cnf.ed_apply_type1_constraints_on_the_fly && !inf_rules_1.empty() ) {
                    for (int d : V[v]) was[d] = true;
                    StandardUtils::makeUnique(inf_rules_1);
                    for ( auto d : inf_rules_1 ) { // in inf_rules we will have N(v), as we add v to S using moveToS(v).
                        assert(d != v);
                        if (!was[d]) {
                            GraphUtils::addEdge(V,v,d);
                            last_reduce_inf_rules_1_added++;
                            all_inf_rules_1_found.emplace_back(v,d);
                            if (use_exhaustively) changes = true;
                        }
                    }
                    for (int d : V[v]) was[d] = false;
                }
            }
        }

        if (cnf.ed_use_edge_removal) {
            VPII edges = GraphUtils::getGraphEdges(V);
            sort(ALL(edges), [&]( PII e1, PII e2 ) {
                auto [a,b] = e1;
                auto [c,d] = e2;
                return V[a].size() + V[b].size() > V[c].size() + V[d].size();
            });
            if ( cnf.ed_node_sorting_mode == 2 ) reverse(ALL(edges));

            // for (int u : nodes)  if (!V[u].empty()) for ( int v : V[u] ) {
            for ( auto [u,v] : edges ){
                // clog << "Considering edge " << PII(u,v) << " for ED-edge-removal" << endl;
                if ( consider({u,v}) ) {
                    clearAllForConsider();
                    GraphUtils::removeEdge(V,u,v);
                    if (write_logs)
                        clog << "\tRemoving edge " << PII(u,v) << " using ED for edge removal" << endl;
                    last_reduce_edges_removed++;

                    int t = reducible_nodes.size();
                    if ( V[u].size() == 1 ) reducible_nodes += propagateDeg1RuleSlow(V[u][0]);
                    if ( V[v].size() == 1 ) reducible_nodes += propagateDeg1RuleSlow(V[v][0]);
                    last_reduce_nodes_removed += t - reducible_nodes.size();

                    if (use_exhaustively) changes = true;
                }
            }
        }
    }

    return reducible_nodes;
}

bool EDReducer::madeChangesInLastReduce() {
    if ( last_reduce_nodes_removed > 0 ) return true;
    if ( last_reduce_edges_removed > 0 ) return true;
    if ( last_reduce_inf_rules_1_added > 0 ) return true;
    return false;
}

void EDReducer::resetAllUsedTechniques() {
    cnf.ed_use_edge_removal = false;
    cnf.ed_apply_type1_constraints_on_the_fly = false;
    cnf.ed_use_edge_insertion = false;
    cnf.ed_use_extended_edges_insertion = false;
}

int EDReducer::getNonWNeighborhoodSize(int u) {
    return accumulate(ALL(V[u]), 0, [&](int s, auto v) {return s + !inW[v];} );
}

int EDReducer::getSIntersection(int u) {
    return accumulate(ALL(V[u]), 0, [&](int s, auto v) {return s + inS[v];} );
}

bool EDReducer::consider(VI initS) {
    clearAllForConsider();

    if (write_logs) clog << "Considering initS: " << initS << endl;

    if (initS.size() == 1) moveToS(initS[0]);
    else {
        S = W = initS;
        for (int d : S) was[d] = true;
        for (int s : S) for (int d : V[s]) if (!was[d]) temp.push_back(d);
        for (int d : S) was[d] = false;


        for (int u : S) inS[u] = inW[u] = true;
        StandardUtils::makeUnique(temp);
        for (int u : temp) moveToU(u);
        temp.clear();

        updateU1();
    }

    while ( true ) {
        if (write_logs) clog << "\tContinuing ED, next step..." << endl;

        int status = nextStep();
        if (status == 1) return true;
        if (status == -1) return false;
        else{} // nothing to do, wait for the next step
    }

    if (write_logs) clog << "\tStopping ED, nothing more to be done..." << endl;
}

bool EDReducer::existsExtDominator(VB & marked) {
    bool exists_ext_dominator = false;

    for ( int u : U1 ) {
        bool all_marked = true;
        int marked_cnt = 0;

        for ( int w : V[u] ) if (!inW[w]) {
            all_marked &= marked[w];
            marked_cnt++;
        }
        if (all_marked) {
            if (write_logs) {
                clog << "\t\tNode " << u << " is an ext-dominator!!! "
                     << "All of its " << marked_cnt << " non-W neighbors are marked" << endl;
            }
            exists_ext_dominator = true;
            break;
        }
    }

    return exists_ext_dominator;
}



void EDReducer::moveToU(int u) {
    assert(!inW[u]);

    U.push_back(u);
    W.push_back(u);
    U1.push_back(u);
    inU[u] = inW[u] = inU1[u] = true;

    inf_rules_1.push_back(u);

    if (cnf.ed_U_nodes_sorting_mode == 1) sort(ALL(U1), [&](int a, int b){ return V[a].size() > V[b].size(); });
    if (cnf.ed_U_nodes_sorting_mode == 2) sort(ALL(U1), [&](int a, int b){ return V[a].size() > V[b].size(); });
}

void EDReducer::moveToS(int u) {
    assert(!inW[u]);

    S.push_back(u);
    W.push_back(u);
    inS[u] = inW[u] = true;

    for (int d : V[u]) if ( !inW[d] ) {
        auto t = write_logs;
        write_logs = false;

        if (write_logs) clog << "\t\tmoving neighbor " << d << " of node u=" << u << " to U" << endl;
        moveToU(d);

        write_logs = t;
    }


    updateU1();

    // for ( int d : U ) { // marking in helper all nodes that should be removed from U1
    //     int c = 0;
    //     for ( int dd : V[d] ) if (inS[dd]) c++;
    //     if (c > 1) helper[d] = true;
    //     // assert(c >= 1);
    // }
    //
    // // removing now nodes from U1, if necessary
    // for (int i=(int)U1.size()-1; i>=0; i--) if (helper[U1[i]]) {
    //     inU1[U1[i]] = false;
    //     swap(U1[i], U1.back());
    //     U1.pop_back();
    // }
    //
    // // clearing
    // for (int d : V[u]) helper[d] = false;
    // for (int d : U) helper[d] = false;

    inf_rules_2.push_back(u);

    if (cnf.ed_U_nodes_sorting_mode == 1) sort(ALL(U1), [&](int a, int b){ return V[a].size() > V[b].size(); });
    if (cnf.ed_U_nodes_sorting_mode == 2) sort(ALL(U1), [&](int a, int b){ return V[a].size() > V[b].size(); });
}

void EDReducer::updateU1() {
    for ( int d : U ) { // marking in helper all nodes that should be removed from U1
        int c = 0;
        for ( int dd : V[d] ) if (inS[dd]) c++;
        if (c > 1) helper[d] = true;
        // assert(c >= 1);
    }

    // removing now nodes from U1, if necessary
    for (int i=(int)U1.size()-1; i>=0; i--) if (helper[U1[i]]) {
        inU1[U1[i]] = false;
        swap(U1[i], U1.back());
        U1.pop_back();
    }

    // clearing
    for (int d : U) helper[d] = false;
}

void EDReducer::markDominationNodes(VB& marked, bool check_double_ed) {

    // for (int u : W) for (int w : V[u]) for ( int d : V[w] ) {  // clearing marked and cnt arrays - should be clear
    // for (int u : U1) {
    for (int u : U1) if (hasNonWIntersectionAtMost(u,cnf.ed_ext_dom_max_node_neigh)) {
        for (int w : V[u]) for ( int d : V[w] ) {  // clearing marked and cnt arrays - should be clear - perhaps this will be enough clearing...
            marked[u] = marked[w] = marked[d] = false;
            cnt[u] = cnt[w] = cnt[d] = 0;
        }
    }

    checkEmptyArraysAssertions(check_double_ed, !check_double_ed); // we check all arrays, including the marked array, which should be empty here

    if constexpr(Config::use_ed_domination) { // the standard concept, used always
        for (int u : U1) {
            if (V[u].empty()) continue;

            int nonw_neigh_size = getNonWNeighborhoodSize(u);
            if (nonw_neigh_size > cnf.ed_ext_dom_max_node_neigh) continue;

            for (int w : V[u]) if (!inW[w]) for (int d : V[w]) cnt[d]++;

            for ( int w : V[u] ) if (!inW[w]) {
                for ( int d : V[w] ) if ( !inW[d] && cnt[d] == nonw_neigh_size ) {
                    if (write_logs) clog << "\t\tmarking node d = " << d << " for u = " << u << ", _c: " << nonw_neigh_size << endl;
                    marked[d] = true;
                }
                break; // we need only 1 node to mark the intersection, so we can break here
            }

            for (int w : V[u]) for (int d : V[w]) cnt[d] = 0; // clearing cnt array for next node u
        }
    }

    if(cnf.ed_use_same_neigh_domination) {
        // here we use ``same neighborhood domination'' approach.
        // we find all nodes w \in N(u) \setminus W such that N(u) \setminus W \subseteq N[w]
        // we can consider only nodes u \in U_0, that is nodes for which N(u) \cap S = \emptyset
        for ( int u : U1 ) {
            bool empty_S_inters = true;
            for ( int d : V[u] ) empty_S_inters &= !inS[d];
            if (!empty_S_inters) continue;

            int nonw_neigh_size = getNonWNeighborhoodSize(u);
            if (nonw_neigh_size > cnf.ed_ext_dom_max_node_neigh) continue;

            for ( int w : V[u] ) if ( !inW[w] ) helper[w] = true; // marking N(u) \setminus W

            for ( int w : V[u] ) if (!inW[w]) {
                int c = 1;
                for (int d : V[w]) if ( helper[d] ) c++;
                assert(c <= nonw_neigh_size);
                if (c == nonw_neigh_size) {
                    if (write_logs) clog << "\t\t\tmarking node " << w << " to move to U using same-neighborhood-rule for node u = " << u  << endl;
                    marked[w] = true;
                }
            }
            for ( int w : V[u] ) if ( !inW[w] ) helper[w] = false; // clearing
        }
    }

    if(cnf.ed_use_deficit1_domination) {
        // generalization of the ``same neighborhood'' domination. Finds all mirrors for nodes in U0
        VI& T = temp;
        VI& vis = temp2;

        for (int u : U1) {
            bool empty_S_inters = true;
            for ( int d : V[u] ) empty_S_inters &= !inS[d];
            if (!empty_S_inters) continue;

            int nonw_neigh_size = getNonWNeighborhoodSize(u);
            if (nonw_neigh_size > cnf.ed_ext_dom_max_node_neigh) continue;

            T.clear(); vis.clear();
            for (int w : V[u]) if (!inW[w]) T.push_back(w); // now T contains N(u) \ W

            for ( int w : T ) for (int d : V[w]) if (!inW[d]) {
                cnt[d]++;
                if ( cnt[d] >= T.size()-1 && !was[d] ) {
                    was[d] = true;
                    vis.push_back(d);
                }
            }
            for (int d : vis) was[d] = false;

            for ( int w : T ) { // we exclude each node in T in turn and update counters
                for (int d : V[w]) if (!inW[d]) cnt[d]--;

                for ( int i=(int)vis.size()-1; i>=0; i-- ) {
                    int d = vis[i];
                    if ( cnt[d] >= T.size()-1 ) {
                        // if after excluding node w for node d=vis[i] it still covers all but for one node in T
                        // then we can add it to U
                        marked[d] = true;
                        if (write_logs) {
                            clog << "\t\tmarking node d: " << d << " for node w: " << w
                                 << " and u: " << u
                                 << ", where u is in U0, using deficit1 approach"
                                 << ", cnt[" << d << "]: " << cnt[d]
                                 << "\n\t\t\t T.size(): " << T.size()
                                 << "\n\t\t\t T: " << T
                                 << "\n\t\t\t V[" << d << "]: " << V[d]
                                 << "\n\t\t\t V[" << w << "]: " << V[w]
                                 << "\n\t\t\t V[" << u << "]: " << V[u]
                                 << endl;
                        }
                        swap(vis[i], vis.back());
                        vis.pop_back();
                    }
                }

                for (int d : V[w]) if (!inW[d]) cnt[d]++;
            }

            for ( int w : T ) for (int d : V[w]) cnt[d] = 0;
        }

        vis.clear();
        T.clear();
    }

    if(this->check_double_ed && check_double_ed) {
        // ``double-ED domination'' - might be considerably slower than other approches,
        // but addresses some of the cases that the other approaches do not

        checkEmptyArraysAssertions(false,check_double_ed);

        VI candidates;
        {
            // create candidates here...
            for ( int u : U1 ) {
                int nonw_neigh_size = getNonWNeighborhoodSize(u);
                if (nonw_neigh_size > min(cnf.ed_ext_dom_max_node_neigh,7)  ) continue;

                // enable candidates from N(W)
                for ( int w : V[u] ) if (!inW[w]) if (!inW[w] && !was[w]) { was[w] = true; candidates.push_back(w); }

                for ( int w : V[u] ) if (!inW[w]) for ( int d : V[w] ) if ( !inW[d] && !was[d] ) {
                    was[d] = true;
                    candidates.push_back(d);
                }
            }
            for (int d : candidates) was[d] = false;

            // now sort candidates by their |N(d) \cap N(U1)|
            for (int u : U1) for (int w : V[u]) if(!inW[w]) was[w] = true;
            for ( int d : candidates ) for (int w : V[d]) cnt[d] += was[w];
            sort(ALL(candidates), [&](int a, int b){ return cnt[a] > cnt[b]; });
            for ( int d : candidates ) cnt[d] = 0;
            for (int u : U1) for (int w : V[u]) if(!inW[w]) was[w] = false;

            if (candidates.size() > cnf.ed_double_ed_max_candidates) candidates.resize(cnf.ed_double_ed_max_candidates);
        }

        if (write_logs) if ( !candidates.empty() ) {
            clog << "\t\t\t There are " << candidates.size() << " candidates to check for double-ED" << endl;
        }

        auto clearMarked = [&](VB & marked) {
            // for (int u : W)
            for (int u : W) if (hasNonWIntersectionAtMost(u,cnf.ed_ext_dom_max_node_neigh)) {
                for (int w : V[u]) for ( int d : V[w] ) marked[u] = marked[w] = marked[d] = false;
            }
        };

        VI neigh_not_in_W;
        VI nodes_removed_from_U; // now we remove nodes from the set U - they will be added back later
        VI nodes_removed_from_U1; // now we remove nodes from the set U - they will be added back later

        VI prevS, prevU, prevU1, prevW, nodes_to_remove_from_W;
        // VB prevInS, prevInU, prevInU1, prevInW;

        constexpr int option = 2; // options 2 seems to work, while option 1 seems to be buggy...

        for ( int w : candidates ) {
            // we assume that w is not in any optimal solution. Thus we can move it to S
            // if this leads to an existence of an ext-dominator, then we can move w to U.
            // we do not check recursively - we do not consider moving nodes when assuming that w is in S.
            // we only check the current state using the existsExtDominator function.
            // It is time-consuming enough to check that...

            assert(!inW[w]);
            clearMarked(marked2);


            if constexpr (option == 1){ // removing N[w] from the graph by marking all nodes in the inW bitvector and removing N(w) from U and U1
                prevS = S; prevU = U; prevU1 = U1; prevW = W;
                // prevInS = inS; prevInU = inU; prevInW = inW;

                {
                    S.push_back(w); W.push_back(w); inS[w] = inW[w] = true; // marking node w as in S
                    nodes_to_remove_from_W.clear();

                    for (int d : V[w]) if ( !inW[d] ) { // equivalent to moveToU, but not adding inference rules
                        auto t = write_logs; write_logs = false;
                        U.push_back(d); W.push_back(d); U1.push_back(d); inU[d] = inW[d] = inU1[d] = true;
                        nodes_to_remove_from_W.push_back(d); // we need that to unmark inW entries without full copies
                        write_logs = t;
                    }
                }

                // // now removing V[w] from U and U1 - it seems that without it the rule is incorrect...
                for (int d : V[w]) was[d] = true;

                // VB inU0(N); for (int u : U1){ int c = 0; for (int d : V[u]) c += inS[d]; inU0[u] = (c == 0); }
                // for (int i=(int)U.size()-1; i>=0; i--) if ( was[U[i]] && !inU0[U[i]] ) { swap(U[i], U.back()); U.pop_back(); }
                // for (int i=(int)U1.size()-1; i>=0; i--) if ( was[U1[i]] && !inU0[U[i]] ) { swap(U1[i], U1.back()); U1.pop_back(); }

                for (int i=(int)U.size()-1; i>=0; i--) if ( was[U[i]] ) { swap(U[i], U.back()); U.pop_back(); }
                for (int i=(int)U1.size()-1; i>=0; i--) if ( was[U1[i]] ) { swap(U1[i], U1.back()); U1.pop_back(); }

                for (int d : V[w]) was[d] = false;

                updateU1();
            }
            else if constexpr (option == 2){ // removing N[w] from the graph by marking all nodes in the inW bitvector and removing N(w) from U and U1
                neigh_not_in_W.clear();
                nodes_removed_from_U.clear();
                nodes_removed_from_U1.clear();

                neigh_not_in_W.push_back(w);
                inW[w] = true;
                for ( int d : V[w] ) if (!inW[d]) {
                    neigh_not_in_W.push_back(d);
                    inW[d] = true;
                    // marking inW[d] way we effectively remove node d from the graph,
                    // but we still need to remove it from U and U1
                }

                // now removing V[w] from U and U1
                for (int d : V[w]) was[d] = true;
                for (int i=(int)U.size()-1; i>=0; i--) if ( was[U[i]] ) {
                    // int c = 0; for ( int d : V[U1[i]] ) c += inW[d];
                    // if (c == 0) continue; // #TEST - do not remove from U1 nodes that are in U0 - check if it is correct!!

                    nodes_removed_from_U.push_back(U[i]);
                    swap(U[i], U.back()); U.pop_back();
                }
                for (int i=(int)U1.size()-1; i>=0; i--) if ( was[U1[i]] ) {
                    // int c = 0; for ( int d : V[U1[i]] ) c += inW[d];
                    // if (c == 0) continue; // #TEST - do not remove from U1 nodes that are in U0 - check if it is correct!!

                    nodes_removed_from_U1.push_back(U1[i]);
                    swap(U1[i], U1.back()); U1.pop_back();
                }
                for (int d : V[w]) was[d] = false;
            }

            constexpr bool check_double_ed_again = false; // we do not want to end in endless loop
            markDominationNodes(marked2, check_double_ed_again);
            bool exists_ext_dominator = existsExtDominator(marked2);


            if (option == 1){ // apply changes back to bring the original state and get ready for the next node
                clearMarked(marked2);
                inS[w] = inW[w] = false;
                for (int d : nodes_to_remove_from_W) inU[d] = inW[d] = inU1[d] = false;
                S = prevS; U = prevU; U1 = prevU1; W = prevW;
                // inS = prevInS; inU = prevInU; inW = prevInW;
            }
            else if (option == 2){ // apply changes back to bring the original state and get ready for the next node
                U += nodes_removed_from_U;
                U1 += nodes_removed_from_U1;
                for (int d : neigh_not_in_W) inW[d] = false;
            }

            if (exists_ext_dominator) {
                if (write_logs) {
                    clog << "\n\t\t\t using double-ext-domination, candidate w: " << w << ", W[" << w << "]: " << V[w]
                         // << "\n\t\t\t nodes_removed: " << nodes_removed
                         // << "\n\t\t\t edges_removed: " << edges_to_remove
                         // << "\n\t\t\t exists_ext_dominator: " << exists_ext_dominator
                         << "\n\t\t\t marking candidate node w: " << w
                         << endl;
                }
                marked[w] = true;
            }
        }

        clearMarked(marked2);
    }
}

VI EDReducer::findNodesToMoveToU() {
    VI nodes_to_move_to_U;
    nodes_to_move_to_U.reserve(W.size());

    if(cnf.ed_consider_nodes_to_move_outside_NW) {
        // for ( int u : U1 ) {
        for ( int u : U1 ) if (hasNonWIntersectionAtMost(u,cnf.ed_ext_dom_max_node_neigh)) {
            for (int w : V[u]) {
                if (!inW[w] && marked[w]) nodes_to_move_to_U.push_back(w);
                for( int d : V[w] ) if(!inW[d] && marked[d]) nodes_to_move_to_U.push_back(d);
            }
        }
    }else {
        // for ( int u : U1 ) {
        for ( int u : U1 ) if (hasNonWIntersectionAtMost(u,cnf.ed_ext_dom_max_node_neigh)) {
            for (int w : V[u]) if (!inW[w] && marked[w]) nodes_to_move_to_U.push_back(w);
        }
    }

    StandardUtils::makeUnique(nodes_to_move_to_U);
    return nodes_to_move_to_U;
}

int EDReducer::nextStep() {
    if (write_logs) {
        clog << "\t\tS: " << S << endl;
        clog << "\t\tU: " << U << endl;
        clog << "\t\tU1: " << U1 << endl;
    }

    checkEmptyArraysAssertions(true, true);

    markDominationNodes(marked);

    if (write_logs) {
        VI temp;
        for ( int u : U ) for (int w : V[u]) if (!inW[w] && marked[w]) temp.push_back(w);
        StandardUtils::makeUnique(temp);
        clog << "\t\tnodes in N(U) marked: " << temp << endl;
    }

    if (existsExtDominator(marked)) {
        return 1;
    }else if (write_logs) clog << "\t\tdominator does not exist" << endl;

    VI nodes_to_move_to_U = findNodesToMoveToU();
    for (int u : nodes_to_move_to_U) marked[u] = false;

    if ( !nodes_to_move_to_U.empty() ) {
        check_double_ed = false; // if a move is possible, do not check double-ED in next iteration

        if (!cnf.ed_move_nodes_to_U_simultaneously) nodes_to_move_to_U.resize(1);
        if (write_logs) clog << "\t\tMoving nodes " << nodes_to_move_to_U << " to U (and creating type-1 constraints)" << endl;
        for (int u : nodes_to_move_to_U) moveToU(u);
        return 0;
    }

    VI nodes_to_move_to_S;
    for ( int u : U1 ) {
        int c = 0, id = -1;
        for (int w : V[u]) if (!inW[w]){ c++; id=w; }
        if (c == 1) nodes_to_move_to_S.push_back(id);
    }
    if (!nodes_to_move_to_S.empty()) {
        check_double_ed = false; // if a move is possible, do not check double-ED in next iteration

        StandardUtils::makeUnique(nodes_to_move_to_S);

        bool move_all_simultanously = cnf.ed_move_nodes_to_S_simultaneously;

        // move all nodes in the same time - correct, but don't we lose some possible domination situations here?
        // by moving node to S we might remove some nodes from U1, which might contribute to ext-domination otherwise
        if (move_all_simultanously){
            if (write_logs) clog << "\t\tMoving nodes " << nodes_to_move_to_S << " to S (and creating type-2 constraints)" << endl;
            for (int u : nodes_to_move_to_S) moveToS(u);
        }else{
            // move to S only the single node from nodes_to_move_to_S fow which the intersection N(w) \cap U' is smallest
            int id = -1;
            int m = 1e9;
            for (int w : nodes_to_move_to_S) {
                int c = 0;
                for (int d : V[w]) if (inU1[d]) c++;
                assert(c > 0);
                bool cond = (c < m);
                cond |= ( c == m && ( id != -1 && V[c].size() < V[id].size() ) );
                if (cond) {
                    m = c;
                    id = w;
                }
            }

            if (write_logs) clog << "\t\tMoving single node " << id << ", with " << m << " neighbors in U', to S (and creating type-2 constraint)" << endl;
            moveToS(id);
        }

        return 0;
    }

    if ( cnf.ed_use_double_ed_checks && !check_double_ed ) {
        // we allow checking the time-consuming double-ED rule only if no other moves are possible...
        check_double_ed = true;
        return 0;
    }

    return -1;
}


void EDReducer::clearAllForConsider() {
    temp.clear();
    temp2.clear();
    for (int d : W) {
        inS[d] = inU[d] = inW[d] = inU1[d] = was[d] = helper[d] = marked[d] = marked2[d] = false;
    }
    for (int d0 : W) for (int d : V[d0]) {
        inS[d] = inU[d] = inW[d] = inU1[d] = was[d] = helper[d] = marked[d] = marked2[d] = false;
    }

    if (cnf.ed_consider_nodes_to_move_outside_NW) {
        // for (int d0 : W) {
        for (int d0 : W) if (hasNonWIntersectionAtMost(d0,cnf.ed_ext_dom_max_node_neigh)) {
            for (int d1 : V[d0]) for (int d : V[d1]) {
                inS[d] = inU[d] = inW[d] = inU1[d] = was[d] = helper[d] = marked[d] = marked2[d] = false;
            }
        }
    }
    inf_rules_1.clear();
    inf_rules_2.clear();
    S.clear();
    U.clear();
    U1.clear();
    W.clear();

    check_double_ed = false;
}

void EDReducer::checkEmptyArraysAssertions(bool check_marked, bool check_marked2) {
    constexpr bool check_slow_assertions_for_correctness = false;
    if constexpr ( check_slow_assertions_for_correctness ) {
        clog << endl << "CAUTION!! RUNNING VERY SLOW ASSERTIONS TO CHECK IF ARRAYS ARE CORRECTLY CLEARED!!" << endl;
        assert(ranges::all_of(cnt, [&](auto b){return !b;}));
        if (check_marked) {
            for (int i=0; i<N; i++) if (marked[i]) clog << "marked[" << i << "] = " << marked[i] << endl;
            assert(ranges::all_of(marked, [&](auto b){return !b;}));
        }
        if (check_marked2) assert(ranges::all_of(marked2, [&](auto b){return !b;}));
        assert(ranges::all_of(was, [&](auto b){return !b;}));
        assert(ranges::all_of(helper, [&](auto b){return !b;}));
        assert(temp.empty());
        assert(temp2.empty());
    }
}

VI EDReducer::propagateDeg1RuleSlow(int v) {
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

    return removed_nodes;
}
