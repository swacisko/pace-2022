//
// Created by sylwester on 3/25/22.
//

//
// Created by sylwester on 1/15/22.
//
#include <CONTESTS/PACE22/heur/SALS3.h>
#include <combinatorics/CombinatoricUtils.h>
#include <CONTESTS/PACE22/Utils.h>
#include <graphs/toposort/TopoSort.h>
#include <graphs/GraphInducer.h>
#include <CONTESTS/PACE22/heur/DFVSSolverH.h>
#include "graphs/GraphUtils.h"
#include "StandardUtils.h"
#include "CollectionOperators.h"

unsigned long long SALS3::x=123456789, SALS3::y=362436069, SALS3::z=521288629;

SALS3::SALS3(VVI V, VI dfvs, Config c)
        : prs_order(V.size()),
          prs_dfvs( 2*V.size()),
          rnd(0,1'000'000ll * 1'000'000ll), rnd_d(0,1),
          cnf(c)
{
    assert(GraphUtils::isSimple(V));
    this->V = V;
    revV = GraphUtils::reverseGraph(V);
    N = V.size();
    best_dfvs = dfvs;
    helper = VB(N,false);

    initialize(dfvs);

}

void SALS3::initialize(VI dfvs) {

    dfvs_size = dfvs.size();
    in_order = VI(N,-1);
    node_on_pos = VI(N * SPARSITY, -1 );
    delta = VPII(N);
    is_valid_delta = VB(N,false);
    insertion_point = VPII(N,{-1,-1});

    VI order = CombinatoricUtils::getFullSetDifference(N, dfvs);
    {
        InducedGraph g = GraphInducer::induce(V, order);
        TopoSort ts(g.V);
        order = ts.sortTopologically();
        for(int & d : order) d = g.nodes[d];
    }
    for( int i=0; i<order.size(); i++){
        in_order[order[i]] = i;
        node_on_pos[i] = order[i];
    }
    sparsifyOrder(); // make the in_order sparser

    for(int i=0; i<N; i++) if( in_order[i] == -1 ) updateNode(i);

    backgoing_arcs_for_node = VI(N,0);
    total_backgoing_arcs = 0;
//    largest_backgoing_arcs_value = 0;
//    nodes_with_violating_arcs = vector<unordered_set<int>>(N);
    for( int i=0; i<N; i++ ){
        backgoing_arcs_for_node[i] = 0;
//        if(in_order[i] != -1) nodes_with_violating_arcs[0].insert(i);
    }

    prs_order = RandomSelectionSet(N);
    moveRandomNodeToOrder();

    if( !use_metropolis_sa ){
        prs_dfvs = RandomSelectionSet(2*N);
        for( int i=0; i<N; i++ ){
            if(in_order[i] == -1) {
                updateDFVSNodePrsProbabilities(i);
            }
        }
    }
}


VI SALS3::localSearch(double T0, double alpha, int maxMvt, int maxFail, int max_iters) {
    if(cnf.write_logs) clog << "Starting SALS3 local search, T0: " << T0 << endl;
    T = T0;
    nbFail = 0;
    last_improvement_T = 2.0;
    int iters_done = 0;

    int nodes_checked = 0;
    int applications_done = 0;

    while( nbFail < maxFail ){
        if( iters_done >= max_iters ) break;
        iters_done++;
        if(cnf.write_logs){
            clog << endl << "Starting iteration #" << iters_done << endl;
            clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;
        }

        nbMvt = 0;
        bool failure = true;

        sparsifyOrder();

//        DEBUG(total_backgoing_arcs);
//        DEBUG(dfvs_size);

        if(cnf.write_logs) clog << "applications_done / nodes_checked = " << (1.0 * applications_done / nodes_checked) << endl;

        LL safety_counter = 0;
        LL cnt = 0;

        while( nbMvt < maxMvt ){

//            if( (nbMvt & 15) == 0 && cnf.sw.tle("main") ) return best_dfvs;
            if( (cnt++ & 255) == 0 && cnf.tle() ) return best_dfvs;

            int v;
            bool before, cond1, cond2;

            if(use_metropolis_sa) {
                v = getNodeForMove();
                assert(in_order[v] == -1);

                nodes_checked++;

                /**
                 * If before is true, then v will be inserted jest before its first out-neighbors.
                 * Otherwise it will be inserted just after its last in-neighbor
                 */
                before = getRandomInt(2);

                int D_before = evaluateNode(v, before);
                int D_after = evaluateNode(v, !before);

                double r = getRandomDouble();
                cond1 = (D_before <= 0 || exp(-1.0 * D_before / T) >= r);
                cond2 = (D_after <= 0 || exp(-1.0 * D_after / T) >= r);
            }else{
                int ind = prs_dfvs.getRandomElement();
//                DEBUG(ind);
                v = (ind/2);
                before = (ind%2);
                cond1 = true; // this is true to make the swap AND do not change value of [before]
                nodes_checked++;
//                DEBUG(v);
//                DEBUG(total_backgoing_arcs);

                const bool check_assertions = false;
                if(check_assertions) {
                    assert(in_order[v] == -1);

                    if (v == -1) {
                        for (int i = 0; i < N; i++) assert(prs_dfvs.getProbab(i) == 0);
                    }


                    for (int i = 0; i < N; i++) { // #TEST - just an assertion
                        assert(delta[i].first >= -1 && delta[i].second >= -1);

                        for (int j = 0; j < N; j++) {
                            if (in_order[i] == -1 && in_order[j] == -1) {
                                {
                                    int d1 = delta[i].first;
                                    int d2 = delta[j].first;
                                    bool cmp = (d1 < d2);

                                    LL p1 = prs_dfvs.getProbab(2 * i);
                                    LL p2 = prs_dfvs.getProbab(2 * j);

                                    bool cmp2 = (p1 > p2);

                                    if (cmp != cmp2) {
                                        DEBUG(d1);
                                        DEBUG(d2);
                                        DEBUG(p1);
                                        DEBUG(p2);
                                        assert(cmp == cmp2);
                                    }
                                }

                                {
                                    int d1 = delta[i].second;
                                    int d2 = delta[j].second;
                                    bool cmp = (d1 < d2);

                                    LL p1 = prs_dfvs.getProbab(2 * i + 1);
                                    LL p2 = prs_dfvs.getProbab(2 * j + 1);

                                    bool cmp2 = (p1 > p2);
                                    if (cmp != cmp2) {
                                        DEBUG(d1);
                                        DEBUG(d2);
                                        DEBUG(p1);
                                        DEBUG(p2);
                                        assert(cmp == cmp2);
                                    }
                                }
                            }
                        }
                    }


                    if (v == -1) {
                        for (int i = 0; i < 2 * N; i++) assert(prs_dfvs.get(i) == 0); // #TEST - just an assertion
                        v = getNodeForMove();
                    } else {
                        if (in_order[v] != -1) {
                            if (before) assert(delta[v].second == 0);
                            else
                                assert(delta[v].first == 0);
                        }
                        DEBUG(prs_dfvs.get(ind));
                    }
                }

            }

            if(cond1 || cond2 ) { // double check metropolis
                if( !cond1 ) before = !before;
                applyMove(v, before);
                nbMvt++;
                applications_done++;

                if (total_backgoing_arcs == 0) {
                    last_improvement_T = T;
                    failure = false;

                    // we have no access to [dfvs], so we need to create the order anew.
                    best_dfvs.clear();
                    best_dfvs = getCurrentDFVS();
                    assert( Utils::isFVS(V,best_dfvs) );

//                    if(cnf.write_logs) clog << "\rFound new DFVS of size " << dfvs_size << ", T: " << T << flush;
                    if(cnf.write_logs){
                        clog << "\rFound new DFVS of size " << dfvs_size << ", T: " << T
                             << ", main elapsed time: " << cnf.sw.getTime("main") / 1'000 << endl;
                    }

                    moveRandomNodeToOrder();
//                    DEBUG(total_backgoing_arcs);
//                    DEBUG(dfvs_size);
//                    exit(1);
                }
            }else{
                // this ought to disable getting stuck into unbreaking solutions in small graphs
                safety_counter++;
                if(safety_counter > 1ll * N * maxMvt) break;
            }
        }

        if( failure ) nbFail++;
        else nbFail = 0;

        if( last_improvement_T < 1.0 ) T *= alpha; // oscillate near the best temperature
        else{
            if(cool_down_quickly_until_improved) T -= 0.0075; // this is to quickly find a good initial temperature
            else T *= alpha;
        }

        if(T<0.07) break;

        const int max_steps_without_improvement = 20;
        if( (last_improvement_T < 1.0 && T < last_improvement_T * pow(alpha,max_steps_without_improvement))
            || (last_improvement_T > 1.0 && T < T0 * pow(alpha,max_steps_without_improvement) && !cool_down_quickly_until_improved)
                ){
            if( last_improvement_T < 1.0 ) T = last_improvement_T;
            const int steps_back = max_steps_without_improvement / 2;
            // move [steps_back] temperature change steps back
            for( int i=0; i<steps_back; i++ ) T *= 1.0 / alpha;
            if(cnf.write_logs) clog << "************** Reverting temperature just above last improvement T: " << T << endl;
        }

        if(cnf.write_logs) clog << "Updating T: " << T << ", dfvs.size(): " << dfvs_size << ", best_dfvs.size(): "
                            << best_dfvs.size() << endl;

        // initialize() cannot be done here as it was e.g. in SALS2
        if(nbFail % 3 == 2) initialize(best_dfvs); // #TEST #CAUTION! - will create an order anew in each iteration from best_order
    }

    return best_dfvs;
}

int SALS3::getNodeForMove() {
    int v = getRandomInt(N);

    while( in_order[v] != -1 ) v = getRandomInt(N);
    assert(in_order[v] == -1);

    return v;
}

void SALS3::insertIntoOrder(int v, int pos) {
    dfvs_size--;

    assert(in_order[v] == -1);

    if(pos == node_on_pos.size()) node_on_pos.push_back(-1);
    if( node_on_pos[pos] != -1 ) shiftRight(pos); // should make free space for node [v]

    node_on_pos[pos] = v;
    in_order[v] = pos;

    { // updating deltas
        for( int u : revV[v] ){
            if( in_order[u] != -1 || !is_valid_delta[u] ) continue;

            int first_out_neigh = insertion_point[u].first;
            if(first_out_neigh != -1) {
                if (in_order[first_out_neigh] > in_order[v]) {
                    // we will need to scan whole revV[u] to find that proper delta, so we mark it as invalid
                    is_valid_delta[u] = false;
                    continue;
                }
            }else{
                is_valid_delta[u] = false;
                continue;
            }

            int last_in_neigh = insertion_point[u].second;
            if(last_in_neigh != -1) {
                if (in_order[last_in_neigh] > in_order[v]) {
                    delta[u].first++; // arc(u,v) is a back-going arc
                }
            }else{
                is_valid_delta[u] = false;
                continue;
            }
        }

        for( int u : V[v] ){
            if( in_order[u] != -1 || !is_valid_delta[u] ) continue;

            int first_out_neigh = insertion_point[u].first;
            if(first_out_neigh != -1){
                if( in_order[v] > in_order[first_out_neigh] ){
                    delta[u].second++;  // arc (v,u) is a back-going arc
                }
            }else {
                is_valid_delta[u] = false;
                continue;
            }

            int last_in_neigh = insertion_point[u].second;
            if(last_in_neigh != -1){
                if( in_order[last_in_neigh] < in_order[v] ){
                    is_valid_delta[u] = false;
                    continue;
                }
            }else{
                is_valid_delta[u] = false;
                continue;
            }
        }

        is_valid_delta[v] = false;
    }

    auto addBackgoingArc = [&]( int u, int v ){
        // apply changes in backgoing arcs
        backgoing_arcs_for_node[u]++;
        total_backgoing_arcs++;

        backgoing_arcs_for_node[v]++;
        total_backgoing_arcs++;
    };

    { // update backgoing arcs
        assert(backgoing_arcs_for_node[v] == 0);

        for( int u : revV[v] ){
            if(in_order[u] == -1) continue;
            if(in_order[u] < in_order[v]) continue;
            addBackgoingArc(u,v);
            prs_order.set(u, score_fun_order(backgoing_arcs_for_node[u]) );
        }

        for( int u : V[v] ){
            if(in_order[u] == -1) continue;
            if( in_order[u] > in_order[v] ) continue;
            addBackgoingArc(u,v);
            prs_order.set(u, score_fun_order(backgoing_arcs_for_node[u]) );
        }
    }

    // update RandomSelectionSet prs_order
    prs_order.set(v, score_fun_order(backgoing_arcs_for_node[v]) );


    if(!use_metropolis_sa){ // updating all invalidated nodes and inserting new values to prs_dfvs
        VI temp = {v};
        temp += V[v] + revV[v];
        for( int u : temp ){
//            if( in_order[u] == -1 && !is_valid_delta[u] ){
//            if( !is_valid_delta[u] ){
                updateNode(u);
                updateDFVSNodePrsProbabilities(u);
//            }
        }
    }
}

void SALS3::removeFromOrder(int v) {
    dfvs_size++;
    assert(in_order[v] != -1);

    { // removing backgoing arcs
        auto removeBackgoingArc = [&](int u, int v){
            backgoing_arcs_for_node[u]--;
            total_backgoing_arcs--;

            backgoing_arcs_for_node[v]--;
            total_backgoing_arcs--;
        };

        { // update backgoing arcs
            for( int u : revV[v] ){
                if(in_order[u] == -1) continue;
                if(in_order[u] < in_order[v]) continue;
                removeBackgoingArc(u,v);
                prs_order.set(u, score_fun_order(backgoing_arcs_for_node[u]) );
            }

            for( int u : V[v] ){
                if(in_order[u] == -1) continue;
                if( in_order[u] > in_order[v] ) continue;
                removeBackgoingArc(u,v);
                prs_order.set(u, score_fun_order(backgoing_arcs_for_node[u]) );
            }
        }
    }

    // update RandomSelectionSet prs_order
    prs_order.set(v, 0 );

    { // updating deltas

        for( int u : revV[v] ){
            if( in_order[u] != -1 || !is_valid_delta[u] ) continue;

            int first_out_neigh = insertion_point[u].first;
            assert(first_out_neigh != -1); // v is out-neigh of u, so first_out_neigh cannot be -1 here

            if( first_out_neigh == v ){
                is_valid_delta[u] = false;
                continue;
            }

            int last_in_neigh = insertion_point[u].second;
            if( last_in_neigh != -1 ){
                if( last_in_neigh == v ){
                    is_valid_delta[u] = false;
                    continue;
                }
                else if( in_order[last_in_neigh] > in_order[v] ){
                    delta[u].first--; // arc (u,v) is a back-going arc here
                }
            }else{} // if there is no arc (x,u) with x in virtual order, then we do not modify anything
        }

        for( int u : V[v] ){
            if( in_order[u] != -1 || !is_valid_delta[u] ) continue;

            int first_out_neigh = insertion_point[u].first;
            if( first_out_neigh != -1 ){
                if(first_out_neigh == v){
                    is_valid_delta[u] = false;
                    continue;
                }
                else if( in_order[first_out_neigh] < in_order[v] ){
                    delta[u].second--;
                }
            }else{} // nothing to do if there is no arc (u,x) with  in virtual order

            int last_in_neigh = insertion_point[u].second;
            assert(last_in_neigh != -1); // v is in-neigh of u, so last_in_neigh cannot be -1 here
            if( last_in_neigh == v ){
                is_valid_delta[u] = false;
                continue;
            }
        }

        is_valid_delta[v] = false;
    }

    node_on_pos[in_order[v]] = -1;
    in_order[v] = -1;


    if(!use_metropolis_sa){ // updating all invalidated nodes and inserting new values to prs_dfvs
        VI temp = {v};
        temp += V[v] + revV[v];
        for( int u : temp ){
//            if( in_order[u] == -1 && !is_valid_delta[u] ){
//            if( !is_valid_delta[u] ){
                updateNode(u);
                updateDFVSNodePrsProbabilities(u);
//            }
        }
    }
}

void SALS3::shiftRight(int pos) {
    int v = node_on_pos[pos];
    if( v == -1 ) return;
    if( pos+1 == node_on_pos.size() ) node_on_pos.push_back(-1);

    if( node_on_pos[pos+1] != -1 ) shiftRight(pos+1);

    node_on_pos[pos] = -1;
    node_on_pos[pos+1] = v;
    in_order[v] = pos+1;
}

void SALS3::shiftLeft(int pos) {
    int v = node_on_pos[pos];
    if( v == -1 ) return;
    if( pos == 0 ){
        if(cnf.write_logs) clog << "In shiftLeft, pos == 0, should do something here!" << endl;
        int t = in_order[-1]; // causing a crash
    }

    if( node_on_pos[pos-1] != -1 ) shiftLeft(pos-1);

    node_on_pos[pos] = -1;
    node_on_pos[pos-1] = v;
    in_order[v] = pos-1;
}


void SALS3::sparsifyOrder() {
    temp_order.clear();
    for( int i=0; i<N; i++ ) if( in_order[i] != -1 ) temp_order.emplace_back( i, in_order[i] );
    sort(ALL(temp_order), []( auto & a, auto & b ){
        return a.second < b.second;
    });

    for( auto & [v,ind] : temp_order) node_on_pos[ind] = in_order[v] = -1;
    int block_size = (0.9 * node_on_pos.size()) / (1 + temp_order.size() );

    for( int i=1; i<=temp_order.size(); i++ ){
        int v = temp_order[i-1].first;
//        in_order[v] = i * block_size; // original, but shift left is dispensible now
        in_order[v] = (i-1) * block_size;
        node_on_pos[in_order[v]] = v;
    }
}

int SALS3::evaluateNode(int v, int before) {
    // UNCOMMENT
//    auto old_delta = delta[v];
//    auto old_ins_point = insertion_point[v];
//    bool valid = is_valid_delta[v];

    if( !is_valid_delta[v] )
        updateNode(v);

//    if( valid ){
//        assert(delta[v] == old_delta);
//        if(insertion_point[v] != old_ins_point){
//            DEBUG(v);
//            DEBUG(insertion_point[v]);
//            DEBUG(old_ins_point);
//        }
//        assert(insertion_point[v] == old_ins_point);
//    }

    if(!before) return delta[v].first;
    else return delta[v].second;
}

void SALS3::applyMove(int v, int before) {
    int pos = -1;

    if(before){
        int w = insertion_point[v].first;
        if(w == -1) pos = getRandomInt(node_on_pos.size() );
        else pos = in_order[w];
    }else{
        int w = insertion_point[v].second;
        if(w == -1) pos = getRandomInt(node_on_pos.size() );
        else pos = in_order[w]+1;
    }

    // we do not need to remove conflicting nodes from order - we remove only one conflicting node in each
    // move application

    int to_remove = getNodeToRemoveFromOrder();
    removeFromOrder(to_remove);

    insertIntoOrder(v, pos);
}

void SALS3::updateNode(int v) {
    int m  = 1e9, M = -1;
    insertion_point[v] = {-1,-1};

    for( int d : V[v] ){
        if( in_order[d] != -1 ){
            if (in_order[d] < m) {
                m = in_order[d];
                insertion_point[v].first = d;
            }
        }
    }
    for( int d : revV[v] ){
        if( in_order[d] != -1 ){
            if (in_order[d] > M) {
                M = in_order[d];
                insertion_point[v].second = d;
            }

        }
    }

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

    { // original version, works very well!
        delta_before = -1 + cnt_before;
        delta_after = -1 + cnt_after;
    }

    delta[v] = { delta_after, delta_before }; // correct
    is_valid_delta[v] = true;
}


void SALS3::updateValidDelta(VI &conflicts) {
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

unsigned long long SALS3::xorshf96() {          //period 2^96-1
    unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;

    return z;
}

int SALS3::getRandomInt(int mod) {
    return xorshf96() % mod;
}

double SALS3::getRandomDouble() {
    return (double)(xorshf96() & ( (1ll<<31)-1 )) / ((1ll<<31)-1); // this is roughly number of the form   a / 2^31
}


VI SALS3::getCurrentDFVS() {
    VI dfvs;
    for( int i=0; i<N; i++ ) if( in_order[i] == -1 ) dfvs.push_back(i);
    if(dfvs.size() != dfvs_size){
        DEBUG(dfvs.size());
        DEBUG(dfvs_size);
        assert(dfvs.size() == dfvs_size);
    }
    return dfvs;
}

double SALS3::findInitialTemperature(VVI V, VI& dfvs, Config cnf, int max_mvt) {
//    double T0 = 0.65;
    double T0 = 0.7; // #TEST
    double alpha = 0.99; // anything here, it does not matter

    while( T0 >= 0.1 ){
        if(cnf.tle()) break;

        SALS3 sals3( V, dfvs, cnf );
        sals3.cnf.write_logs = false;

        auto temp = sals3.localSearch( T0, alpha, max_mvt, 1, 1 );
        if( temp.size() < dfvs.size() ){
            dfvs = temp;
            return T0;
        }

        T0 -= 0.075;
//        T0 -= 0.015; // #TEST
//        T0 *= 0.9; // #TEST
    }

    return -1; // returning default initial temperature
}

void SALS3::moveRandomNodeToOrder() {
    // now moving some random node from DFVS to order.
    int v = getNodeForMove();
    int pos = getRandomInt(node_on_pos.size());
    insertIntoOrder(v, pos);
}

int SALS3::getNodeToRemoveFromOrder() {

    { // using RandomSelectionSet
        int v = prs_order.getRandomElement();

        if(v == -1){
            for(int i=0; i<N; i++) if(in_order[i] != -1){ v = i; break; }
        }

//        if(!prs_order.get(v) > 0){
//            DEBUG(V.size());
//            DEBUG(v);
//            DEBUG(getCurrentDFVS().size());
//            int in_order_cnt1 = 0;
//            for(int i=0; i<N; i++) if(in_order[i] != -1) in_order_cnt1++;
//            int in_order_cnt2 = 0;
//            for(int i=0; i<N; i++) if(prs_order.get(v) > 0) in_order_cnt2++;
//            DEBUG(in_order_cnt1);
//            DEBUG(in_order_cnt2);
//        }
//
//        assert(prs_order.get(v) > 0 );
//        assert(v != -1); // if this assertion fails, then it means that there are no nodes in order with backgoing arcs
//        assert( in_order[v] != -1 );
//        assert( backgoing_arcs_for_node[v] > 0 );
        return v;
    }

}



void SALS3::updateDFVSNodePrsProbabilities(int v) {

//    const int SCALE = 1 + log(N);
//    const int SCALE = 1 + sqrt(N);
    const int SCALE = N;

    for( int before = 0; before <= 1; before++ ) {
        LL val = INVERSE_PRS;

        if( in_order[v] != -1 ) val = 0;
        else {
            if (before) {
                if (delta[v].second <= 0) val *= SCALE; // if v is a good node, then probability of taking it is high
                else if (delta[v].second > 0) val /= score_fun_dfvs(delta[v].second);
            } else {
                if (delta[v].first <= 0) val *= SCALE;
                else if (delta[v].first > 0) val /= score_fun_dfvs(delta[v].first);
            }

//            val = ceil(1.0 * val / 1'000);
        }


//        prs_dfvs.set(2 * v + (before ? 1 : 0), score_fun_dfvs(val));
        prs_dfvs.set(2 * v + (before ? 1 : 0), val);
    }
}

void SALS3::setUseMetropolisSa(bool useMetropolisSa) {
    use_metropolis_sa = useMetropolisSa;
}





