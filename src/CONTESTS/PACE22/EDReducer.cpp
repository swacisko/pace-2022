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

    for (int v : nodes) {
        if (considerNode(v)) {
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
            all_marked |= marked[w];
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

    for (int d : V[u]) if ( !inW[d] ) moveToU(d);
    for ( int d : U ) {
        int c = 0;
        for ( int dd : V[d] ) if (inS[dd]) c++;
        if (c > 1) helper[d] = true;
        assert(c >= 1);
    }

    // removing now nodes from U1, if necessary
    for (int i=(int)U1.size()-1; i>=0; i--) if (helper[U1[i]]) {
        inU1[U1[i]] = false;
        swap(U1[i], U1.back());
        U1.pop_back();
    }

    for (int d : V[u]) helper[d] = false;

    inf_rules_2.push_back(u);

    if (cnf.ed_U_nodes_sorting_mode == 1) sort(ALL(U1), [&](int a, int b){ return V[a].size() > V[b].size(); });
    if (cnf.ed_U_nodes_sorting_mode == 2) sort(ALL(U1), [&](int a, int b){ return V[a].size() > V[b].size(); });
}

void EDReducer::markDominationNodes() {

    for (int u : U1) for (int w : V[u]) for ( int d : V[w] ) marked[d] = cnt[d] = 0; // clearing marked array

    for (int u : U1) {
        if (V[u].empty()) continue;

        int c = 0;
        for (int w : V[u]) if (!inW[w]) {
            c++;
            for (int d : V[w]) cnt[d]++;
        }

        for ( int w : V[u] ) if (!inW[w]) for ( int d : V[w] ) if ( cnt[d] == c ) marked[d] = true;

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
        for ( int u : U ) for (int w : V[u]) if (!inW[w]) temp.push_back(w);
        StandardUtils::makeUnique(temp);
        clog << "\t\tnodes marked: " << temp << endl;
    }

    if (existsExtDominator()) {
        return 1;
    }

    VI nodes_to_move_to_U;
    for ( int u : U1 ) for (int w : V[u]) if (!inW[w] && marked[w]) nodes_to_move_to_U.push_back(w);
    if ( !nodes_to_move_to_U.empty() ) {
        StandardUtils::makeUnique(nodes_to_move_to_U);
        clog << "\t\tMoving nodes " << nodes_to_move_to_U << " to U (and creating type-1 constraints)" << endl;
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
        clog << "\t\tMoving nodes " << nodes_to_move_to_S << " to S (and creating type-2 constraints)" << endl;
        for (int u : nodes_to_move_to_S) moveToS(u);
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
