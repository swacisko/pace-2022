//
// Created by sylwester on 1/5/22.
//
#include <CONTESTS/PACE22/heur/SALS.h>
#include <combinatorics/CombinatoricUtils.h>
#include <CONTESTS/PACE22/Utils.h>
#include <graphs/toposort/TopoSort.h>
#include <graphs/GraphInducer.h>
#include <CONTESTS/PACE22/heur/DFVSSolverH.h>
#include "graphs/GraphUtils.h"
#include "StandardUtils.h"
#include "CollectionOperators.h"

SALS::SALS(VVI V, VI dfvs) : rnd(0,1'000'000ll * 1'000'000ll), rnd_d(0,1) {
    this->V = V;
    revV = GraphUtils::reverseGraph(V);
    N = V.size();
    this->dfvs = dfvs;
    in_order = VI(N,-1);
    delta = VPII(N);
    is_valid_delta = VB(N,false);
    best_dfvs = dfvs;

    helper = VB(N,false);

    {
        order = CombinatoricUtils::getFullSetDifference(N, dfvs);
        {
            InducedGraph g = GraphInducer::induce(V, order);
            TopoSort ts(g.V);
            order = ts.sortTopologically();
            for(int & d : order) d = g.nodes[d];
        }

        updateOrder();
        for(int i=0; i<N; i++) if( in_order[i] == -1 ) updateNode(i);
    }


    { // creating probabs
        Config cnf;
        DFVSSolverH solver(cnf);
        probabs = solver.sinkhorn(V);
//        DEBUG(probabs.size());
        {
            for (double &d : probabs) d = 1.0 / d;
            SCALE = 10;
        }
//        {
//            for( double & d : probabs ) d = pow(1-d,2);
//            SCALE = 100;
//        }
    }

    createDfvsPrefSum();
}

double SALS::findInitialTemperature(VVI V, VI &dfvs, int max_mvt) {
    double T0 = 0.33;
    double alpha = 0.99; // anything here, it does not matter

    int N = V.size();

    while( T0 >= 0.1 ){

        SALS sals( V, dfvs );
        auto temp = sals.localSearch( T0, alpha, 10*N, 1, 1 );
        if( temp.size() < dfvs.size() ) return T0;

        T0 -= 0.01;
    }

    return -1;
}

VI SALS::localSearch(double T0, double alpha, int maxMvt, int maxFail, int max_iters) {
    if(write_logs) clog << "Starting SALS local search, T0: " << T0 << endl;
    T = T0;
    nbFail = 0;
    last_improvement_T = 2.0;
    dfvs_check_marker = dfvs.size();

    int iters_done = 0;
    while( nbFail < maxFail ){
        if( iters_done >= max_iters ) break;
        iters_done++;
        if(write_logs) clog << endl << "Starting iteration #" << iters_done << endl;


        nbMvt = 0;
        bool failure = true;


        pair<int,pair<int,bool>> invalid_move = {1e9, {-1,-1}};
        pair<int,pair<int,bool>> best_move = invalid_move;

        int time_since_change = 0;
        int frequency;

//            frequency = dfvs.size();


//        int opt = rnd.nextInt(5);
//        if( opt <= 0 ) frequency = sqrt(dfvs.size());
//        else if( opt <= 1 ) frequency = dfvs.size() / 3;
//        else if( opt <= 2 ) frequency = dfvs.size() / 2;
//        else if( opt <= 3 ) frequency = 2*dfvs.size() / 3;
//        else if( opt <= 4 ) frequency = dfvs.size();


        int opt = rnd.nextInt(4);
        if( opt <= 0 ) frequency = dfvs.size() / 3;
        else if( opt <= 1 ) frequency = dfvs.size() / 2;
        else if( opt <= 2 ) frequency = 2*dfvs.size() / 3;
        else if( opt <= 3 ) frequency = dfvs.size();

        while( nbMvt < maxMvt ){ // original
//        while( nbMvt < ( (last_improvement_T > 1.0) ? maxMvt / 3 : maxMvt) ){ // the condition is to faster find initial good value of T

            int v = getNodeForMove();

            /**
             * If before is true, then v will be inserted jest before its first out-neighbors.
             * Otherwise it will be inserted just after its last in-neighbor
             */
            bool before = rnd.nextInt(2);

            int D_before = evaluateNode(v, before);
            int D_after = evaluateNode(v, !before);

            time_since_change++;
            if (D_before < best_move.first) {
                best_move = {D_before, {v, before}};
            }
            if (D_after < best_move.first) {
                best_move = {D_after, {v, !before}};
            }

            const bool use_metropolis = true;
            if(use_metropolis){
                double r = rnd_d.rand();
                bool cond1 = (D_before <= 0 || exp(-1.0 * D_before / T) >= r);
                bool cond2 = (D_after <= 0 || exp(-1.0 * D_after / T) >= r);
                if (cond1 || cond2 || dfvs_check_marker == 0) { // metropolis
                    if( !cond1 ) before = !before;

                    if(dfvs_check_marker == 0){ // if all nodes from [dfvs] where checked, then we apply the best move
//                        clog << "----------------------------------------- MARKER == 0" << endl;
                        v = best_move.second.first;
                        before = best_move.second.second;
                    }

                    applyMove(v, before);
                    nbMvt++;

                    if (dfvs.size() < best_dfvs.size()) {
                        last_improvement_T = T;
                        failure = false;
                        best_dfvs = dfvs;
                        assert(Utils::isFVS(V, best_dfvs));
                        if(write_logs) clog << "\rFound new DFVS of size " << dfvs.size() << ", T: " << T << flush;
                    }

                    time_since_change = 0;
                    best_move = invalid_move;
                }

            }else {

                if (D_before <= 0 ||
                    (time_since_change % frequency) == 0) { // selecting the best move from last [frequency] moves
                    v = best_move.second.first;
                    before = best_move.second.second;

                    applyMove(v, before);
                    nbMvt++;
                    if (dfvs.size() < best_dfvs.size()) {
                        failure = false;
                        best_dfvs = dfvs;
                        clog << "Found new DFVS of size " << dfvs.size() << ", T: " << T << endl;
                        DEBUG(dfvs.size());
                    }

                    time_since_change = 0;
                    best_move = invalid_move;
                }
            }
        }

        if( failure ) nbFail++;
        else nbFail = 0;

        if( last_improvement_T < 1.0 ) T *= alpha; // #TEST!!! #CAUTION - oscillate near the best temperature
        else{
            if(cool_down_quickly_until_improved) T -= 0.0075; // this is to quickly find a good initial temperature
            else T *= alpha;
        }

        if(T<0.07) break;

        const int max_steps_without_improvement = 10;
        if( (last_improvement_T < 1.0 && T < last_improvement_T * pow(alpha,max_steps_without_improvement))
            || (last_improvement_T > 1.0 && T < T0 * pow(alpha,max_steps_without_improvement) && !cool_down_quickly_until_improved)
                ){
            if( last_improvement_T < 1.0 ) T = last_improvement_T;
            const int steps_back = 7;
            // move [steps_back] temperature change steps back
            for( int i=0; i<steps_back; i++ ) T *= 1.0 / alpha;
            if(write_logs) clog << "************** Reverting temperature just above last improvement T: " << T << endl;
        }

        if(write_logs) clog << "Updating T: " << T << ", dfvs.size(): " << dfvs.size() << ", best_dfvs.size(): "
                            << best_dfvs.size() << endl;


        bool start_from_best = false;

        if(start_from_best) { // #TEST
            dfvs = best_dfvs;
            order = CombinatoricUtils::getFullSetDifference(N, dfvs);
            {
                InducedGraph g = GraphInducer::induce(V, order);
                TopoSort ts(g.V);
                order = ts.sortTopologically();
                for(int & d : order) d = g.nodes[d];
            }
            updateOrder();
            fill(ALL(is_valid_delta),false);
            for (int i = 0; i < N; i++) if (in_order[i] == -1) updateNode(i);

            createDfvsPrefSum();
            if(write_logs)DEBUG(best_dfvs.size());
        }
    }

    return best_dfvs;
}

int SALS::getNodeForMove() {
    if(use_uniform_node_selection){ // uniform selection
        if(dfvs_check_marker == 0) return dfvs[0]; // returning anything
        int r = rnd.nextInt(dfvs_check_marker);
        int v = dfvs[r];
        swap( dfvs[r], dfvs.back() );
        dfvs_check_marker--;
        return v;
    }
    else{ // nonuniform selection
        int r = rnd.nextInt( dfvs_pref_sum.back() );
        int ind = lower_bound(ALL(dfvs_pref_sum), r) - dfvs_pref_sum.begin();
        return dfvs[ind];
    }
}

void SALS::updateOrder() {
    fill(ALL(in_order),-1);
    for(int i=0; i<order.size(); i++) in_order[order[i]] = i;
}

void SALS::updateOrder(int pos, VI &conflicts) {
    // should be the same as updateOrder(), but this should be a bit faster
    int ind = pos;
    for(int d : conflicts){
        if(in_order[d] != -1) ind = min(ind, in_order[d]);
    }

    for(int i=ind; i<order.size(); i++) in_order[order[i]] = i;
    for( int d : conflicts ) in_order[d] = -1;

}


int SALS::evaluateNode(int v, int before) {
    auto old_delta = PII(1e9,1e9);
    bool valid = false;
    if(is_valid_delta[v]){
        valid = true;
        old_delta = delta[v];
    }

    if( !is_valid_delta[v] ) updateNode(v);

    if( valid ) assert(delta[v] == old_delta);

    if(!before) return delta[v].first;
    else return delta[v].second;
}

void SALS::applyMove(int v, int before) {
    int pos;

    if(before){
        pos = 1e9;
        for( int d : V[v] ){
            if( in_order[d] != -1 ) pos = min(pos, in_order[d]);
        }

        if( pos == 1e9 ) pos = rnd.nextInt(order.size()+1);
    }else{
        pos = -1;
        for( int d : revV[v] ){
            if( in_order[d] != -1 ) pos = max(pos, in_order[d]+1);
        }

        if( pos == -1 ) pos = rnd.nextInt(order.size()+1);
    }


    /**
     * Vector of conflicting nodes.
     */
    conflicts.clear();

//    for( int d : V[v] ) if( in_order[d] < pos ) conflicts.push_back(d);
    for( int d : V[v] ) if( in_order[d] != -1 && in_order[d] < pos ) conflicts.push_back(d);
//    for( int d : revV[v] ) if( in_order[d] >= pos ) conflicts.push_back(d);
    for( int d : revV[v] ) if( in_order[d] != -1 && in_order[d] >= pos ) conflicts.push_back(d);


    order.insert(order.begin()+pos, v);
    delta[v] = {1e9,1e9}; // updating delta for node v. It should not be necessary...
    is_valid_delta[v] = false;

    // removes from order all elements in [conflicts]
    StandardUtils::removeFromArrayPreserveOrderInplace(order, conflicts, helper);

    {
        // removing v from dfvs and adding [conflicts] to dfvs.
        // This should be the same as     dfvs = CombinatoricUtils::getFullSetDifference(N, order)
        for( int i=0; i<dfvs.size(); i++ ){
            if( dfvs[i] == v ){
                swap(dfvs[i], dfvs.back());
                dfvs.pop_back();
                break;
            }
        }
        for(int d : dfvs) helper[d] = true;
        for( int d : conflicts ) if(!helper[d]) dfvs.push_back(d);
        for(int d : dfvs) helper[d] = false;
    }

    updateOrder(pos, conflicts); // this updates only indices at positions >= pos

    conflicts.push_back(v);
    updateValidDelta(conflicts);


    dfvs_check_marker = dfvs.size();
    if(!use_uniform_node_selection) createDfvsPrefSum();
}

void SALS::updateNode(int v) {
    int m  = 1e9, M = -1;
    for( int d : V[v] ) if( in_order[d] != -1 ) m = min( m, in_order[d] );
    for( int d : revV[v] ) if( in_order[d] != -1 ) M = max( M, in_order[d] );

    /**
     * cnt_before is the number of conflicting arcs that will be added if we insert v just before its first
     * out-neighbor
     *
     * cnt_after is the number of conflicting arcs that will be added if we insert v just after its last
     * in-neighbor
     */
    int cnt_before = 0, cnt_after = 0;

    for( int d : revV[v] ) if( in_order[d] != -1 && in_order[d] > m ) cnt_before++;
    for( int d : V[v] ) if( in_order[d] != -1 && in_order[d] < M ) cnt_after++;

    int delta_after = 0;
    int delta_before = 0;

    VI conflicts;
    auto calculateDeltaForConflicts = [&](){
        for(int d : conflicts) helper[d] = true;

        int removed_conflict_arcs = 0;
        for( int u : conflicts ){
            for( int d : V[u] ) if( in_order[d] != -1 && in_order[d] < in_order[u] && helper[d] ) removed_conflict_arcs++;
            for( int d : revV[u] ) if( in_order[d] != -1 && in_order[d] > in_order[u] && helper[d] ) removed_conflict_arcs++;
        }

        for(int d : conflicts) helper[d] = false;

        return removed_conflict_arcs;
    };

//    { // calculate delta_before
//        conflicts.clear();
//        for( int d : revV[v] ) if( in_order[d] != -1 && in_order[d] >= m ) conflicts.push_back(d);
//        auto removed_arcs = calculateDeltaForConflicts();
//        delta_before = -1 + cnt_before - removed_arcs;
//    }
//
//    { // calculate delta_after
//        conflicts.clear();
//        for( int d : V[v] ) if( in_order[d] != -1 && in_order[d] < M ) conflicts.push_back(d);
//        auto removed_arcs = calculateDeltaForConflicts();
//        delta_after = -1 + cnt_after - removed_arcs;
//    }

    { // original version, works very well!
        delta_before = -1 + cnt_before;
        delta_after = -1 + cnt_after;
    }

    delta[v] = { delta_after, delta_before }; // correct
    is_valid_delta[v] = true;
}

void SALS::createDfvsPrefSum() {
    dfvs_pref_sum = VLL(dfvs.size(), 0);
    for( int i=0; i<dfvs.size(); i++ ){
        if( i > 0 ) dfvs_pref_sum[i] = dfvs_pref_sum[i-1];
        dfvs_pref_sum[i] += ceil(probabs[ dfvs[i] ] * SCALE);
    }
}

void SALS::updateValidDelta(VI &conflicts) {

    int operations_to_do = 1 + conflicts.size();
    for (int d : conflicts) operations_to_do += V[d].size() + revV[d].size();

    // this if here is avoid peak complexity
    if(operations_to_do < N * sqrt(N) ) {

        for (int d : conflicts) is_valid_delta[d] = false;
        for (int d : conflicts) {
            for (int u : V[d]) {
                if (in_order[u] == -1) {
                    is_valid_delta[u] = false;
                }
            }

            for (int u : revV[d]) {
                if (in_order[u] == -1) {
                    is_valid_delta[u] = false;
                }
            }
        }
    }
    else{
        for( int i=0; i<N; i++ ) if(in_order[i] == -1) is_valid_delta[i] = false;
    }

}




