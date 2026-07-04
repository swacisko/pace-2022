//
// Created by sylwe on 07/08/2025.
//

#include "EDReducer.h"

#include "GraphUtils.h"
#include "StandardUtils.h"
#include "CONTESTS/PACE22/Utils.h"


VI EDReducer::reduce(VVI& V0) {
    V = V0;

    VI nodes(N);
    iota(ALL(nodes),0);
    if ( cnf.ed_node_sorting_mode == 1 ) sort(ALL(nodes), [&]( int a, int b ) { return V[a].size() > V[b].size(); });
    if ( cnf.ed_node_sorting_mode == 2 ) sort(ALL(nodes), [&]( int a, int b ) { return V[a].size() < V[b].size(); });

    VI reducible_nodes;

    DEBUG(V);

    for (int v : nodes) {
        if (considerNode(v)) {
            if (write_logs) clog << "Node " << v << " is reducible!!!" << endl << endl << endl;
            reducible_nodes.push_back(v);
            GraphUtils::removeNodeFromGraph(V,v);
        }
    }


    return reducible_nodes;
}

bool EDReducer::considerNode(int v) {
    clearAll();

    if (write_logs) clog << "Considering node " << v << endl;

    moveToS(v);
    while ( true ) {
        if (write_logs) clog << "\tContinuing ED, next step..." << endl;

        int status = nextStep();
        if (status == 1) return true;
        if (status == -1) return false;
        else{} // nothing to do, wait for the next step
    }

    if (write_logs) clog << "\tStopping ED, nothing more to be done..." << endl;
}

bool EDReducer::existsExtDominator() {
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
        clog << "\t\tmoving neighbor " << d << " of node u=" << u << " to U" << endl;
        moveToU(d);
    }

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
    for (int d : V[u]) helper[d] = false;
    for (int d : U) helper[d] = false;

    inf_rules_2.push_back(u);

    if (cnf.ed_U_nodes_sorting_mode == 1) sort(ALL(U1), [&](int a, int b){ return V[a].size() > V[b].size(); });
    if (cnf.ed_U_nodes_sorting_mode == 2) sort(ALL(U1), [&](int a, int b){ return V[a].size() > V[b].size(); });
}

void EDReducer::markDominationNodes() {

    for (int u : U1) for (int w : V[u]) for ( int d : V[w] ) marked[d] = cnt[d] = 0; // clearing marked array

    for (int u : U1) {
        if (V[u].empty()) continue;

        int _c = 0;
        for (int w : V[u]) if (!inW[w]) {
            _c++;
            for (int d : V[w]) cnt[d]++;
        }

        for ( int w : V[u] ) if (!inW[w]) {
            for ( int d : V[w] ) if ( !inW[d] && cnt[d] == _c ) {
                if (write_logs) clog << "\t\tmarking node d = " << d << " for u = " << u << ", _c: " << _c << endl;
                marked[d] = true;
            }
            break; // we need only 1 node to mark the intersection, so we can break here
        }

        for (int w : V[u]) for (int d : V[w]) cnt[d] = 0; // clearing cnt array for next node u
    }
}

int EDReducer::nextStep() {
    if (write_logs) {
        clog << "\t\tS: " << S << endl;
        clog << "\t\tU: " << U << endl;
        clog << "\t\tU1: " << U1 << endl;
    }

    markDominationNodes();

    if (write_logs) {
        VI temp;
        for ( int u : U ) for (int w : V[u]) if (!inW[w] && marked[w]) temp.push_back(w);
        StandardUtils::makeUnique(temp);
        clog << "\t\tnodes in N(U) marked: " << temp << endl;
    }

    if (existsExtDominator()) {
        return 1;
    }else if (write_logs) clog << "\t\tdominator does not exist" << endl;

    VI nodes_to_move_to_U;
    for ( int u : U1 ) for (int w : V[u]) if (!inW[w] && marked[w]) nodes_to_move_to_U.push_back(w);
    { // here we find all nodes w \in N(u) \setminus W such that N(u) \setminus W \subseteq N[w]
        for ( int u : U1 ) {
            bool empty_S_inters = true;
            for ( int d : V[u] ) empty_S_inters &= !inS[d];

            if ( empty_S_inters ){
                int neigh_size = 0;
                for ( int w : V[u] ) if ( !inW[w] ) {
                    helper[w] = true; // marking N(u) \setminus W
                    neigh_size++;
                }
                for ( int w : V[u] ) if (!inW[w]) {
                    int c = 1;
                    for (int d : V[w]) if ( helper[d] ) c++;
                    assert(c <= neigh_size);
                    if (c == neigh_size) {
                        if (write_logs) clog << "\t\t\tmarking node " << w << " to move to U using same-neighborhood-rule for node u = " << u  << endl;
                        nodes_to_move_to_U.push_back(w);
                    }
                }
                for ( int w : V[u] ) if ( !inW[w] ) helper[w] = false; // clearing
            }
        }
    }

    if ( !nodes_to_move_to_U.empty() ) {
        StandardUtils::makeUnique(nodes_to_move_to_U);
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
        StandardUtils::makeUnique(nodes_to_move_to_S);

        bool move_all_simultanously = false;

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

    return -1;
}


void EDReducer::clearAll() {
    temp.clear();
    temp2.clear();
    for (int d : W) {
        inS[d] = inU[d] = inW[d] = inU1[d] = was[d] = helper[d] = false;
    }
    for (int d0 : W) for (int d : V[d0]) {
        inS[d] = inU[d] = inW[d] = inU1[d] = was[d] = helper[d] = false;
    }

    inf_rules_1.clear();
    inf_rules_2.clear();
    S.clear();
    U.clear();
    U1.clear();
    W.clear();
}
