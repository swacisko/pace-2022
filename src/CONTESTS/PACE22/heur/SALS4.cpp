//
// Created by sylwester on 3/25/22.
//

//
// Created by sylwester on 1/15/22.
//
#include <CONTESTS/PACE22/heur/SALS4.h>
#include <combinatorics/CombinatoricUtils.h>
#include <CONTESTS/PACE22/Utils.h>
#include <graphs/toposort/TopoSort.h>
#include <graphs/GraphInducer.h>
#include <CONTESTS/PACE22/heur/DFVSSolverH.h>
#include "graphs/GraphUtils.h"
#include "StandardUtils.h"
#include "CollectionOperators.h"

unsigned long long SALS4::x=123456789, SALS4::y=362436069, SALS4::z=521288629;

SALS4::SALS4(VVI V, VI dfvs, Config c)
        : prs_order(V.size()),
          rnd(0,1'000'000ll * 1'000'000ll), rnd_d(0,1),
          cnf(c)
{
    assert(GraphUtils::isSimple(V));
    this->V = V;
    revV = GraphUtils::reverseGraph(V);
    N = V.size();
    best_dfvs = dfvs;
    helper = VB(N,false);

    points.reserve(N);
    buckets = VVI(N);

    initialize(dfvs);

}

void SALS4::initialize(VI dfvs) {

    dfvs_size = dfvs.size();
    in_order = VI(N,-1);
    node_on_pos = VI(N * SPARSITY, -1 );
    delta = vector<vector<INS_PT>>(N);
    is_valid_delta = VB(N,false);

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
    for( int i=0; i<N; i++ ){
        backgoing_arcs_for_node[i] = 0;
    }

    prs_order = RandomSelectionSet(N);
    moveRandomNodeToOrder();

}


VI SALS4::localSearch(double T0, double alpha, int maxMvt, int maxFail, int max_iters) {
    if(cnf.write_logs) clog << "Starting SALS4 local search, T0: " << T0 << endl;
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

        if(cnf.write_logs) clog << "applications_done / nodes_checked = " << (1.0 * applications_done / nodes_checked) << endl;

        LL safety_counter = 0;
        LL cnt = 0;

        while( nbMvt < maxMvt ){

//            if( (safety_counter & 15) == 0 && cnf.sw.tle("main") ) return best_dfvs;
            if( (cnt++ & 255) == 0 && cnf.tle() ) return best_dfvs;

            int v = getNodeForMove();
            assert(in_order[v] == -1);

            if(!is_valid_delta[v]) updateNode(v);

            int insert_pos = -1;
            {
                int val = 1e9;
                for (auto ipt : delta[v]) {
                    if( ipt.u == -1 ){
                        insert_pos = getRandomInt(node_on_pos.size());
                        break;
                    }
                    double r = getRandomDouble();
                    int vvv = -1 + ipt.cnt;
//                    if (vvv <= 0 || exp(-1.0 * vvv / T) >= r) {
                    if (ipt.cnt <= 0 || exp(-1.0 * vvv / T) >= r) {
                        insert_pos = in_order[ipt.u];
                        if(!ipt.before) insert_pos++;
                        val = ipt.cnt;
                        break;
                    }
                }
            }

            nodes_checked++;


            if(insert_pos != -1) {
                applyMove(v, insert_pos);
                nbMvt++;
                applications_done++;

                if (total_backgoing_arcs == 0) {
                    last_improvement_T = T;
                    failure = false;

                    // we have no access to [dfvs], so we need to create the order anew.
                    best_dfvs.clear();
                    best_dfvs = getCurrentDFVS();
                    assert( Utils::isFVS(V,best_dfvs) );

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
        if(nbFail % 3 == 2) initialize(best_dfvs);
    }

    return best_dfvs;
}

int SALS4::getNodeForMove() {
    int v = getRandomInt(N);

    while( in_order[v] != -1 ) v = getRandomInt(N);
    assert(in_order[v] == -1);

    return v;
}

void SALS4::insertIntoOrder(int v, int pos) {
    dfvs_size--;

    assert(in_order[v] == -1);

    while(pos >= node_on_pos.size()) node_on_pos.push_back(-1);
    if( node_on_pos[pos] != -1 ) shiftRight(pos); // should make free space for node [v]

    node_on_pos[pos] = v;
    in_order[v] = pos;
    is_valid_delta[v] = false;

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

    for( int u : V[v] ) is_valid_delta[u] = false;
    for( int u : revV[v] ) is_valid_delta[u] = false;
}

void SALS4::removeFromOrder(int v) {
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

    is_valid_delta[v] = false;

    node_on_pos[in_order[v]] = -1;
    in_order[v] = -1;


    for( int u : V[v] ) is_valid_delta[u] = false;
    for( int u : revV[v] ) is_valid_delta[u] = false;

}

void SALS4::shiftRight(int pos) {
    int v = node_on_pos[pos];
    if( v == -1 ) return;
    if( pos+1 == node_on_pos.size() ) node_on_pos.push_back(-1);

    if( node_on_pos[pos+1] != -1 ) shiftRight(pos+1);

    node_on_pos[pos] = -1;
    node_on_pos[pos+1] = v;
    in_order[v] = pos+1;
}

void SALS4::shiftLeft(int pos) {
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


void SALS4::sparsifyOrder() {
    temp_order.clear();
    for( int i=0; i<N; i++ ) if( in_order[i] != -1 ) temp_order.emplace_back( i, in_order[i] );
    sort(ALL(temp_order), []( auto & a, auto & b ){
        return a.second < b.second;
    });

    for( auto & [v,ind] : temp_order) node_on_pos[ind] = in_order[v] = -1;
    int block_size = (0.9 * node_on_pos.size()) / (1 + temp_order.size() );

    for( int i=1; i<=temp_order.size(); i++ ){
        int v = temp_order[i-1].first;
        in_order[v] = (i-1) * block_size;
        node_on_pos[in_order[v]] = v;
    }
}


void SALS4::applyMove(int v, int pos) {
    int to_remove = getNodeToRemoveFromOrder();
    removeFromOrder(to_remove);
    insertIntoOrder(v, pos);
}

void SALS4::updateNode(int v) {
    delta[v].clear();

    points.clear();

    for( int d : revV[v] ) if(in_order[d] != -1) points.emplace_back( in_order[d], d, true );
    for( int d : V[v] ) if(in_order[d] != -1) points.emplace_back( in_order[d], d, false );

    sort(ALL(points), []( auto a, auto b ){ return a.pos < b.pos; });
    if(points.empty()){
        // #CAUTION! we insert node u with id -1 - this means that we can insert node anywhere
        delta[v].emplace_back( -1,0,false );
        is_valid_delta[v] = true;
        return;
    }


    int cnt = 0;
    for( int d : revV[v] ) if(in_order[d] != -1) cnt++;

    for( int i=0; i<points.size(); i++ ){
        delta[v].emplace_back( points[i].u, cnt, true );

        int p = i;
        while( p < points.size() && points[i].pos == points[p].pos ){
            if( points[p].in_neigh ) cnt--;
            else cnt++;
            p++;
        }

        i = p-1;
    }

    // this is the case of inserting node v after its last neighbor
    delta[v].emplace_back( points.back().u, cnt, false );

    constexpr bool use_fast_shuffle = true; // much faster when not using UniformIntGenerator
    if(use_fast_shuffle){
        for( int i=(int)delta[v].size()-1; i>0; i-- ){
            swap( delta[v][i], delta[v][ getRandomInt(i+1) ] );
        }
    }
    else StandardUtils::shuffle(delta[v]);

    constexpr bool use_bucket_sort = false; // much slower when using bucket sort here
    if(use_bucket_sort){
        int B = min(50, N-1);
        int BS = (node_on_pos.size()+B) / B;

        for( int i=0; i<delta[v].size(); i++ ){
            int pos = in_order[delta[v][i].u];
            if(!delta[v][i].before) pos++;

            int ind = (pos / BS);
            buckets[ind].push_back(i);
        }

        vector<INS_PT> temp; temp.reserve(delta[v].size());
        for( int i=0; i<=B; i++ ){
            if( buckets[i].size() > 1 ){
                stable_sort(ALL(buckets[i]), [&](int a, int b) { return delta[v][a].cnt < delta[v][b].cnt; });
            }
            for( int d : buckets[i] ) temp.push_back( delta[v][d] );
        }

        swap(delta[v], temp);
        for( int i=0; i<=B; i++ ) buckets[i].clear();
    }
    else sort(ALL(delta[v]), []( auto a, auto b ){ return a.cnt < b.cnt; }); // this can be sped up using bucket sort

    is_valid_delta[v] = true;
}

unsigned long long SALS4::xorshf96() {          //period 2^96-1
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

int SALS4::getRandomInt(int mod) {
    return xorshf96() % mod;
}

double SALS4::getRandomDouble() {
    return (double)(xorshf96() & ( (1ll<<31)-1 )) / ((1ll<<31)-1); // this is roughly number of the form   a / 2^31
}


VI SALS4::getCurrentDFVS() {
    VI dfvs;
    for( int i=0; i<N; i++ ) if( in_order[i] == -1 ) dfvs.push_back(i);
    if(dfvs.size() != dfvs_size){
        DEBUG(dfvs.size());
        DEBUG(dfvs_size);
        assert(dfvs.size() == dfvs_size);
    }
    return dfvs;
}

double SALS4::findInitialTemperature(VVI V, VI& dfvs, Config cnf, int max_mvt) {
    double T0 = 0.25; // original
//    double T0 = 0.28; // #TEST
    double alpha = 0.99; // anything here, it does not matter

    while( T0 >= 0.1 ){

        SALS4 sals4( V, dfvs, cnf );
        sals4.cnf.write_logs = true; // #TEST

        auto temp = sals4.localSearch( T0, alpha, max_mvt, 1, 1 );
        if( temp.size() < dfvs.size() ){
            dfvs = temp;
            return T0;
        }

        T0 -= 0.04;
    }

    return -1; // returning default initial temperature
}

void SALS4::moveRandomNodeToOrder() {
    // now moving some random node from DFVS to order.
    int v = getNodeForMove();
    int pos = getRandomInt(node_on_pos.size());
    insertIntoOrder(v, pos);
}

int SALS4::getNodeToRemoveFromOrder() {

    { // using RandomSelectionSet
        int v = prs_order.getRandomElement();
        while( v == -1 || in_order[v] == -1 ) v = getRandomInt(N);

//        assert(v != -1); // if this assertion fails, then it means that there are no nodes in order with backgoing arcs
//        assert( in_order[v] != -1 );
//        assert(prs_order.get(v) > 0 );
//        assert( backgoing_arcs_for_node[v] > 0 );
        return v;
    }

}




