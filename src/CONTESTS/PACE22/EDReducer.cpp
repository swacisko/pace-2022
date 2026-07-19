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

    Stopwatch sw;
    string ed = "ed";
    sw.setLimit(ed,max_time_millis);
    sw.start(ed);

    while (changes) {
        if (sw.tle(ed)) break;
        changes = false;

        if (cnf.ed_use_node_removal) {
            for (int v : nodes) if (!V[v].empty()) {
                if (sw.tle(ed)) break;

                // clog << "\rConsidering node " << v << flush;
                if (consider({v})) {
                    clearAllForConsider();
                    if (write_logs)
                        clog << "\t\tNode " << v << " is ED-reducible!   final W.size(): " << W.size() << endl << endl << endl;
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

            for ( auto [u,v] : edges ) if ( V[u].size() >= 2 && V[v].size() >= 2 ) {
                if (sw.tle(ed)) break;
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
                    last_reduce_nodes_removed += reducible_nodes.size() - t;

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

int EDReducer::getSIntersectionSize(int u) {
    return accumulate(ALL(V[u]), 0, [&](int s, auto v) {return s + inS[v];} );
}

bool EDReducer::consider(VI initS) {
    clearAllForConsider();

    if (write_logs) clog << "Considering initS: " << initS << endl;

    checkEmptyArraysAssertions(true,true);

    if (initS.size() == 1) {
        int v = initS[0];
        moveToS(v, true);
        for (int d : V[v]) if (!excluded[d] && inU1[d]) to_consider_in_next_step[d] = true;
    }
    else {
        S = W = initS;
        for (int d : S) was[d] = true;
        for (int s : S) for (int d : V[s]) if (!was[d]) temp.push_back(d);
        for (int d : S) was[d] = false;


        for (int u : S) inS[u] = inW[u] = true;
        StandardUtils::makeUnique(temp);
        int W_size_after_all_simultaneous_moves = temp.size();
        for (int u : temp) {
            // bool excl =  (getLowerBoundOnNonWNeighbors(u) >= cnf.ed_min_nonw_deg_to_exclude_node);
            bool excl = (max(0, (int)V[u].size() - W_size_after_all_simultaneous_moves  - 1) >= cnf.ed_min_nonw_deg_to_exclude_node);
            moveToU(u,excl);
        }
        temp.clear();

        updateU1();

        for (int u : S) for (int d : V[u]) if (!excluded[d]) {
            if (inU1[d]) to_consider_in_next_step[d] = true;
            if constexpr (keep_track_of_degrees) deg_in_S[d]++;
        }
    }

    // clog << "Initial sizes: W: " << W.size() << ", U: " << U.size() << ", U1: " << U1.size() << endl;
    checkEmptyArraysAssertions(true,true);

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

    for ( int u : U1 ) if ( to_consider_in_next_step[u] ) {
    // for ( int u : U1 ) {
        bool all_marked = true;
        int marked_cnt = 0;

        for ( int w : V[u] ) if (!inW[w]) {
            all_marked &= marked[w];
            if (!all_marked) break;
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



void EDReducer::moveToU(int u, const bool exclude) {
    assert(!inW[u]);

    if (exclude) {
        excluded[u] = true;
        U.push_back(u);
        W.push_back(u);
        inU[u] = inW[u] = true;
        inU1[u] = inU0[u] = false;
        if constexpr (keep_track_of_degrees) deg_notin_W[u] = deg_in_S[u] = 0;
        if constexpr (keep_track_of_degrees) for (int d : V[u]) deg_notin_W[d] -= inW[d]; // we may  exclude the node, but we need to update neighbors!! Get rid of the ''deg_notin_shit'' completely!
        return;
    }

    U.push_back(u);
    W.push_back(u);
    U1.push_back(u);
    inU[u] = inW[u] = inU1[u] = true;
    inU0[u] = true;

    inf_rules_1.push_back(u);

    if constexpr (keep_track_of_degrees) deg_notin_W[u] = deg_in_S[u] = 0;
    for (int d : V[u]) {
        if constexpr (keep_track_of_degrees) deg_notin_W[d] -= inW[d]; // we decrease value for each neighbor d, because u is moved to W
        if constexpr (keep_track_of_degrees) deg_notin_W[u] += !inW[d]; // we calculate value for the moved node u
        if constexpr (keep_track_of_degrees) deg_in_S[u] += inS[d];
        if (inS[d]) inU0[u] = false;

        if (inU1[d]) to_consider_in_next_step[d] = true;
    }

    to_consider_in_next_step[u] = true;

}

void EDReducer::moveToS(int u, const bool init_exclude) {
    assert(!inW[u]);

    S.push_back(u);
    W.push_back(u);
    inS[u] = inW[u] = true;
    inU[u] = inU1[u] = inU0[u] = false;

    if constexpr (keep_track_of_degrees) deg_in_S[u] = 0;
    for (int d : V[u]) if (!excluded[d]) {
        if constexpr (keep_track_of_degrees) deg_in_S[u] += inS[d];
        inU0[d] = false;
    }

    for (int d : V[u]) if (!excluded[d]) {

        if ( !inW[d] ) {
            auto t = write_logs;
            write_logs = false;

            if (write_logs) clog << "\t\tmoving neighbor " << d << " of node u=" << u << " to U" << endl;
            bool exclude =  (getLowerBoundOnNonWNeighbors(d) >= cnf.ed_min_nonw_deg_to_exclude_node);
            if (init_exclude) {
                int W_size_after_all_simultaneous_moves = V[u].size();
                exclude = (max(0, (int)V[d].size() - W_size_after_all_simultaneous_moves - 1) >= cnf.ed_min_nonw_deg_to_exclude_node);
            }
            moveToU(d, exclude);

            write_logs = t;
        }else { // d is not excluded and is already in W
            if constexpr (keep_track_of_degrees) deg_notin_W[d]--;
            if constexpr (keep_track_of_degrees) deg_in_S[d]++;
        }
    }

    if constexpr (keep_track_of_degrees) deg_notin_W[u] = 0;

    updateU1();

    // for (int d : V[u]) if (!excluded[d] && inU1[d]) to_consider_in_next_step[d] = true;

    inf_rules_2.push_back(u);
}

void EDReducer::updateU1() {

    for (int i=(int)U1.size()-1; i>=0; i--) {
        int u = U1[i];
        int c = 0;
        for (int d : V[u]) c += inS[d];
        if (c > 1) {
            inU1[u] = inU0[u] = false;
            to_consider_in_next_step[u] = false;
            REM(U1,i);
        }
    }
}

void EDReducer::markDominationNodes(VB& marked, bool check_double_ed) {

    for (int u : U1) if (hasNonWIntersectionAtMost(u,cnf.ed_ext_dom_max_node_neigh)) {
        // for (int w : V[u]) for ( int d : V[w] ) {  // clearing marked and cnt arrays - should be clear - perhaps this will be enough clearing...
        for (int w : V[u]) if (!inW[w]) for ( int d : V[w] ) {  // clearing marked and cnt arrays - should be clear - perhaps this will be enough clearing...
            marked[u] = marked[w] = marked[d] = false;
            cnt[u] = cnt[w] = cnt[d] = 0;
        }
    }

    checkEmptyArraysAssertions(check_double_ed, !check_double_ed); // we check all arrays, including the marked array, which should be empty here

    if constexpr(Config::use_ed_domination) { // the standard concept, used always
        // for (int u : U1) if (!V[u].empty()) {
        for (int u : U1) if (!V[u].empty() && to_consider_in_next_step[u]) {
            if (!has_cnf_bounded_neigh[u]) continue;
            if ( cnf.ed_use_deficit1_domination && inU0[u] ) continue; // do not duplicate search

            int nonw_neigh_size = 0;
            for (int w : V[u]) {
                nonw_neigh_size += !inW[w];
                if (!inW[w]) for (int d : V[w]) cnt[d]++;
            }

            if constexpr (keep_track_of_degrees) {
                if(nonw_neigh_size != deg_notin_W[u]) { DEBUG(nonw_neigh_size); DEBUG(deg_notin_W[u]); }
                assert(nonw_neigh_size == deg_notin_W[u]);
            }

            for ( int w : V[u] ) if (!inW[w]) {
                for ( int d : V[w] ) if ( !inW[d] && cnt[d] == nonw_neigh_size ) {
                    if (write_logs) clog << "\t\tmarking node d = " << d << " for u = " << u << ", _c: " << nonw_neigh_size << endl;
                    marked[d] = true;
                }
                break; // we need only 1 node to mark the intersection, so we can break here
            }

            for (int w : V[u]) if (!inW[w]) for (int d : V[w]) cnt[d] = 0; // clearing cnt array for next node u
        }
    }

    if(cnf.ed_use_same_neigh_domination) {
        // here we use ``same neighborhood domination'' approach.
        // we find all nodes w \in N(u) \setminus W such that N(u) \setminus W \subseteq N[w]
        // we can consider only nodes u \in U_0, that is nodes for which N(u) \cap S = \emptyset
        // for ( int u : U1 ) if (!V[u].empty()) {
        for ( int u : U1 ) if (!V[u].empty() && to_consider_in_next_step[u]) {
            if (!has_cnf_bounded_neigh[u]) continue;
            if (!inU0[u]) continue;
            if ( cnf.ed_use_deficit1_domination ) continue; // do not duplicate search

            // int nonw_neigh_size = deg_notin_W[u];
            int nonw_neigh_size = 0;

            for ( int w : V[u] ) if ( !inW[w] ) {
                helper[w] = true; // marking N(u) \setminus W
                nonw_neigh_size++;
            }

            if constexpr (keep_track_of_degrees) {
                if(nonw_neigh_size != deg_notin_W[u]) { DEBUG(nonw_neigh_size); DEBUG(deg_notin_W[u]); }
                assert(nonw_neigh_size == deg_notin_W[u]);
            }

            for ( int w : V[u] ) if (!inW[w]) {
                int c = 1;
                for (int d : V[w]) c += helper[d];
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
        // generalization of the ``same neighborhood'' domination. Finds all mirrors for nodes in U0 that have at
        // most 1 node outside N(u) \setminus W
        VI& T = temp;
        VI& vis = temp2;

        // for (int u : U1) if (!V[u].empty()) {
        for (int u : U1) if (!V[u].empty() && to_consider_in_next_step[u]) {
            if (!has_cnf_bounded_neigh[u]) continue;
            if (!inU0[u]) continue;

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

            if (!vis.empty())
            for ( int w : T ) { // we exclude each node in T in turn and update counters

                for (int d : V[w]) cnt[d] -= !inW[d];

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

                for (int d : V[w]) cnt[d] += !inW[d];
            }

            for ( int w : T ) for (int d : V[w]) cnt[d] = 0;
        }

        vis.clear();
        T.clear();
    }

    if(this->check_double_ed && check_double_ed) {
        // ``double-ED domination'' - might be considerably slower than other approches,
        // but addresses some of the cases that the other approaches do not
        // assert(false && "Modify correctly degrees in deg_in_S and deg_notin_W!");

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
            for (int u : W) if (!excluded[u] && hasNonWIntersectionAtMost(u,cnf.ed_ext_dom_max_node_neigh)) {
                marked[u] = false;
                for (int w : V[u]) if (!inW[w]) for ( int d : V[w] ) marked[w] = marked[d] = false;
            }
        };

        VI neigh_not_in_W;
        VI nodes_removed_from_U; // now we remove nodes from the set U - they will be added back later
        VI nodes_removed_from_U1; // now we remove nodes from the set U - they will be added back later
        VI nodes_removed_from_U0; // now we remove nodes from the set U - they will be added back later
        vector<pair<int,bool>> has_cnf_bndn_to_restore; has_cnf_bndn_to_restore.reserve(U1.size());

        // VI prevS, prevU, prevU1, prevW, nodes_to_remove_from_W;
        // VB prevInS, prevInU, prevInU1, prevInW;

        for ( int w : candidates ) {
            // we assume that w is not in any optimal solution. Thus we can move it to S
            // if this leads to an existence of an ext-dominator, then we can move w to U.
            // we do not check recursively - we do not consider moving nodes when assuming that w is in S.
            // we only check the current state using the existsExtDominator function.
            // It is time-consuming enough to check that...

            assert(!inW[w]);
            clearMarked(marked2);

            { // removing N[w] from the graph by marking all nodes in the inW bitvector and removing N(w) from U and U1
                neigh_not_in_W.clear();
                nodes_removed_from_U.clear();
                nodes_removed_from_U1.clear();
                nodes_removed_from_U0.clear();
                has_cnf_bndn_to_restore.clear();

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
                for (int i=(int)U1.size()-1; i>=0; i--) if ( was[U1[i]] ) {
                    int x = U1[i];

                    has_cnf_bndn_to_restore.emplace_back(x,has_cnf_bounded_neigh[x]);

                    if (inU0[x]) { // if x is in U0, then we do note remove it from U1, we remove it from U0 only
                        nodes_removed_from_U0.push_back(x);
                        inU0[x] = false;
                        continue;
                    }

                    nodes_removed_from_U1.push_back(x);
                    inU1[x] = false;
                    REM(U1,i);
                }
                for (int d : V[w]) was[d] = false;
            }

            constexpr bool check_double_ed_again = false; // we do not want to end in endless loop
            markDominationNodes(marked2, check_double_ed_again);
            bool exists_ext_dominator = existsExtDominator(marked2);


            { // apply changes back to bring the original state and get ready for the next node
                U += nodes_removed_from_U;
                U1 += nodes_removed_from_U1;
                for (int d : nodes_removed_from_U1) inU1[d] = true;
                for (int d : nodes_removed_from_U0) inU0[d] = true;
                for (int d : neigh_not_in_W) inW[d] = false;
                for (auto [d,val] : has_cnf_bndn_to_restore) has_cnf_bounded_neigh[d] = val;
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

int EDReducer::getLowerBoundOnNonWNeighbors(int u) {
    return max(0, (int)V[u].size() - (int)U.size() - 1);
}

VI EDReducer::findNodesToMoveToU() {
    VI nodes_to_move_to_U;
    nodes_to_move_to_U.reserve(U1.size());

    if(cnf.ed_consider_nodes_to_move_outside_NW) {
        for ( int u : U1 ) if (hasNonWIntersectionAtMost(u,cnf.ed_ext_dom_max_node_neigh)) {
            // for (int w : V[u]) {
            //     if (!inW[w] && marked[w]) nodes_to_move_to_U.push_back(w);
            //     for( int d : V[w] ) if(!inW[d] && marked[d]) nodes_to_move_to_U.push_back(d);
            // }
            for (int w : V[u]) if (!inW[w]) {
                if (marked[w]) nodes_to_move_to_U.push_back(w);
                for( int d : V[w] ) if(!inW[d] && marked[d]) nodes_to_move_to_U.push_back(d);
            }
        }
    }else {
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
        for ( int u : U1 ) for (int w : V[u]) if (!inW[w] && marked[w]) temp.push_back(w);
        StandardUtils::makeUnique(temp);
        clog << "\t\tnodes in N(U) marked: " << temp << endl;
    }

    if (existsExtDominator(marked)) {
        return 1;
    }else if (write_logs) clog << "\t\tdominator does not exist" << endl;


    for (int d : U1) to_consider_in_next_step[d] = false;

    VI nodes_to_move_to_U = findNodesToMoveToU();
    for (int u : nodes_to_move_to_U) marked[u] = false;

    if ( !nodes_to_move_to_U.empty() ) {
        check_double_ed = false; // if a move is possible, do not check double-ED in next iteration

        if (!cnf.ed_move_nodes_to_U_simultaneously) nodes_to_move_to_U.resize(1);
        if (write_logs) clog << "\t\tMoving nodes " << nodes_to_move_to_U << " to U (and creating type-1 constraints)" << endl;
        int W_size_after_all_simultaneous_moves = W.size() + nodes_to_move_to_U.size();
        for (int u : nodes_to_move_to_U) {
            bool excl = (max(0, (int)V[u].size() - W_size_after_all_simultaneous_moves  - 1) >= cnf.ed_min_nonw_deg_to_exclude_node);
            moveToU(u, excl);
        }
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

            // for some reason moving nodes simultaneously fails on some tests...
            // for (const int u : nodes_to_move_to_S) moveToS(u);
            // for (const int u : nodes_to_move_to_S) for (int d : V[u]) if (inU1[d]) to_consider_in_next_step[d] = true;

            for (const int u : nodes_to_move_to_S) if (!inW[u]) moveToS(u);
            for (const int u : nodes_to_move_to_S) if (inS[u]) for (int d : V[u]) if (inU1[d]) to_consider_in_next_step[d] = true;
        }else{
            // move to S only the single node from nodes_to_move_to_S fow which the intersection N(w) \cap U' is smallest
            int id = -1;
            int m = 1e9;
            for (int w : nodes_to_move_to_S) {
                assert(!inW[w]);
                int c = 0;
                for (int d : V[w]) if (inU1[d]) c++;
                assert(c > 0 && "this might fail if we design detection of nodes that can be moved to S that are not in N(W)");
                bool cond = (c < m);
                cond |= ( c == m && ( id != -1 && V[c].size() < V[id].size() ) );
                if (cond) {
                    m = c;
                    id = w;
                }
            }

            if (write_logs) clog << "\t\tMoving single node " << id << ", with " << m << " neighbors in U', to S (and creating type-2 constraint)" << endl;
            moveToS(id);
            for (int d : V[id]) if (inU1[d]) to_consider_in_next_step[d] = true;
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

    if (cnf.ed_consider_nodes_to_move_outside_NW) {
        for (int d0 : W) if (has_cnf_bounded_neigh[d0] && !excluded[d0]) {
            for (int d1 : V[d0]) if (!inW[d1] && !clearing_helper[d1]) {
                clearing_helper[d1] = true;
                for (int d : V[d1]) if (!inW[d]) {
                    inS[d] = inU[d] = inW[d] = false;
                    inU1[d] = inU0[d] = was[d] = false;
                    helper[d] = marked[d] = marked2[d] = false;
                    if constexpr (keep_track_of_degrees) deg_in_S[d] = deg_notin_W[d] = 0;
                }
            }
        }
    }

    for (int d0 : W) if (has_cnf_bounded_neigh[d0] && !excluded[d0]) for (int d : V[d0]) if (!inW[d]) {
        inS[d] = inU[d] = inW[d] = false;
        inU1[d] = inU0[d] = was[d] = false;
        helper[d] = marked[d] = marked2[d] = false;
        if constexpr (keep_track_of_degrees) deg_in_S[d] = deg_notin_W[d] = 0;
        clearing_helper[d] = false;
    }

    for (int d : W) {
        inS[d] = inU[d] = inW[d] = false;
        inU1[d] = inU0[d] = was[d] = false;
        helper[d] = marked[d] = marked2[d] = false;
        has_cnf_bounded_neigh[d] = excluded[d] = false;
        if constexpr (keep_track_of_degrees) deg_in_S[d] = deg_notin_W[d] = 0;
        to_consider_in_next_step[d] = false;
    }


    inf_rules_1.clear();
    inf_rules_2.clear();
    S.clear();
    U.clear();
    U1.clear();
    W.clear();

    check_double_ed = false;
}

// void EDReducer::excludeHighDegreeNodesFromU(int max_nonw_deg) {
//     for ( int i=(int)U.size()-1; i>=0; i-- ) {
//         int u = U[i];
//         if (deg_notin_W[u] > max_nonw_deg) {
//             inU[u] = inU0[u] = false;
//             REM(U,i);
//         }
//     }
// }

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
        // assert(ranges::all_of(has_cnf_bounded_neigh, [&](auto b){return !b;}));
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
