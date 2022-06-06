//
// Created by sylwester on 1/15/22.
//
#include <CONTESTS/PACE22/heur/SALS2.h>
#include <combinatorics/CombinatoricUtils.h>
#include <CONTESTS/PACE22/Utils.h>
#include <graphs/toposort/TopoSort.h>
#include <graphs/GraphInducer.h>
#include <CONTESTS/PACE22/heur/DFVSSolverH.h>
#include "graphs/GraphUtils.h"
#include "StandardUtils.h"
#include "CollectionOperators.h"

//int SALS2::last_seed = 129103;
unsigned long long SALS2::x=123456789, SALS2::y=362436069, SALS2::z=521288629;

SALS2::SALS2(VVI V, VI dfvs, Config c, bool ufa)
:
rnd(0,1'000'000ll * 1'000'000ll), rnd_d(0,1),
cnf(c),
use_fast_updates(ufa)
{
//    rng.seed(last_seed++);
    assert(GraphUtils::isSimple(V));
    this->V = V;
    revV = GraphUtils::reverseGraph(V);
    N = V.size();
    best_dfvs = dfvs;

    helper = VB(N,false);

    initialize(dfvs);
//    {
//        dfvs_size = dfvs.size();
//        in_order = VI(N,-1);
//        node_on_pos = VI(N * SPARSITY, -1 );
//        delta = VPII(N);
//        is_valid_delta = VB(N,false);
//        insertion_point = VPII(N,{-1,-1});
//
//        VI order = CombinatoricUtils::getFullSetDifference(N, dfvs);
//        {
//            InducedGraph g = GraphInducer::induce(V, order);
//            TopoSort ts(g.V);
//            order = ts.sortTopologically();
//            for(int & d : order) d = g.nodes[d];
//        }
//        for( int i=0; i<order.size(); i++){
//            in_order[order[i]] = i;
//            node_on_pos[i] = order[i];
//        }
//        sparsifyOrder(); // make the in_order sparser
//
//        for(int i=0; i<N; i++) if( in_order[i] == -1 ) updateNode(i);
//    }

   preprocessRandoms();

}

void SALS2::initialize(VI dfvs) {

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
}

VI SALS2::localSearch(double T0, double alpha, int maxMvt, int maxFail, int max_iters) {
    if(cnf.write_logs) clog << "Starting SALS2 local search" << ( use_fast_updates ? " with fast updates, " : ", " ) << "T0: " << T0 << endl;
    T = T0;
    nbFail = 0;
    last_improvement_T = 2.0;
    int iters_done = 0;

    while( nbFail < maxFail ){
        if( iters_done >= max_iters ){
            // do not terminate when in find_init_temperature_mode
            if( find_init_temperature_mode && last_improvement_T != 2.0 );
            else{
//                DEBUG(find_init_temperature_mode);
//                DEBUG(last_improvement_T);
                break;
            }
        }
        iters_done++;
        if(cnf.write_logs){
            clog << endl << "Starting iteration #" << iters_done << endl;
            clog << "Main elapsed seconds: " << cnf.sw.getTime() / 1'000 << endl;
        }

        nbMvt = 0;
        bool failure = true;

        sparsifyOrder();

        LL safety_counter = 0;
        LL cnt = 0;

        while( nbMvt < maxMvt ){

//            if( (nbMvt & 15) == 0 && cnf.sw.tle("main") ) return best_dfvs;
            if( (cnt++ & 255) == 0 && cnf.tle() ) return best_dfvs;

            int v = getNodeForMove();
            assert(in_order[v] == -1);

            /**
             * If before is true, then v will be inserted jest before its first out-neighbors.
             * Otherwise it will be inserted just after its last in-neighbor
             */
            bool before = getRandomInt(2);

            int D_before = evaluateNode(v, before);
            int D_after = evaluateNode(v, !before);

            double r = getRandomDouble();
            bool cond1 = (D_before <= 0 || exp(-1.0 * D_before / T) >= r);
            bool cond2 = (D_after <= 0 || exp(-1.0 * D_after / T) >= r);

//            if(cond1 ) { // single check metropolis
            if(cond1 || cond2 ) { // double check metropolis
//                clog << "Applying move #" << iters_done++ << endl;

                if( !cond1 ) before = !before;

                applyMove(v, before);
                nbMvt++;

                if (dfvs_size < best_dfvs.size()) {
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
                }
            }else{
                // this ought to disable getting stuck into unbreaking solutions in small graphs
                safety_counter++;
                if(safety_counter > 1ll * N * maxMvt) break;
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
            const int steps_back = max_steps_without_improvement;
            // move [steps_back] temperature change steps back
            for( int i=0; i<steps_back; i++ ) T *= 1.0 / alpha;
            if(cnf.write_logs) clog << "************** Reverting temperature just above last improvement T: " << T << endl;
        }

        if(cnf.write_logs) clog << "Updating T: " << T << ", dfvs.size(): " << dfvs_size << ", best_dfvs.size(): "
             << best_dfvs.size() << endl;

//        if(nbFail & 1) initialize(getCurrentDFVS()); // #TEST #CAUTION! - will create an order anew in each iteration
        if(nbFail & 1) initialize(best_dfvs); // #TEST #CAUTION! - will create an order anew in each iteration from best_order

        preprocessRandoms();
    }

    return best_dfvs;
}

int SALS2::getNodeForMove() {
    int v = getRandomInt(N);

    while( in_order[v] != -1 ) v = getRandomInt(N);
    assert(in_order[v] == -1);

    return v;
}

void SALS2::insertIntoOrder(int v, int pos) {
    dfvs_size--;

    assert(in_order[v] == -1);

    if(pos == node_on_pos.size()){
//        if(write_logs) clog << "Pushing back an element in insertIntoOrder, T: " << T << endl;
        node_on_pos.push_back(-1);
    }

    if( node_on_pos[pos] != -1 ) shiftRight(pos); // should make free space for node [v]

    node_on_pos[pos] = v;
    in_order[v] = pos;

    if(use_fast_updates){ // updating deltas
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

}

void SALS2::removeFromOrder(int v) {
    dfvs_size++;
    assert(in_order[v] != -1);

    if(use_fast_updates){ // updating deltas

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
}



void SALS2::shiftRight(int pos) {
    int v = node_on_pos[pos];
    if( v == -1 ) return;
    if( pos+1 == node_on_pos.size() ){
//        if(write_logs) clog << "In shiftRight, pos+1 == node_on_pos.size(), should push_back!!!, T: " << T << endl;
        node_on_pos.push_back(-1);
    }

    if( node_on_pos[pos+1] != -1 ) shiftRight(pos+1);

    node_on_pos[pos] = -1;
    node_on_pos[pos+1] = v;
    in_order[v] = pos+1;
}

void SALS2::shiftLeft(int pos) {
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


void SALS2::sparsifyOrder() {
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

int SALS2::evaluateNode(int v, int before) {
    // UNCOMMENT
//    auto old_delta = delta[v];
//    auto old_ins_point = insertion_point[v];
//    bool valid = is_valid_delta[v];

    if( !is_valid_delta[v] )
        updateNode(v);

//    if( valid ){
//        assert(delta[v] == old_delta);
//        if(use_fast_updates){
//
//            if(insertion_point[v] != old_ins_point){
//                DEBUG(v);
//                DEBUG(insertion_point[v]);
//                DEBUG(old_ins_point);
//            }
//            assert(insertion_point[v] == old_ins_point);
//        }
//    }

    if(!before) return delta[v].first;
    else return delta[v].second;
}

void SALS2::applyMove(int v, int before) {
    int pos = -1;

    if(before){
        if(!use_fast_updates){
            int max_ind = node_on_pos.size()+1;
            pos = max_ind;
            for( int d : V[v] ) if( in_order[d] != -1 ) pos = min(pos, in_order[d]);
            if( pos == max_ind ) pos = getRandomInt(node_on_pos.size() );
        }else{
            int w = insertion_point[v].first;
            if(w == -1) pos = getRandomInt(node_on_pos.size() );
            else pos = in_order[w];
        }
    }else{
        if(!use_fast_updates){
            pos = -1;
            for( int d : revV[v] ) if( in_order[d] != -1 ) pos = max(pos, in_order[d]+1);
            if( pos == -1 ) pos = getRandomInt(node_on_pos.size() );
        }else{
            int w = insertion_point[v].second;
            if(w == -1) pos = getRandomInt(node_on_pos.size() );
            else pos = in_order[w]+1;
        }
    }

    conflicts.clear();

    for( int d : V[v] ) if( in_order[d] != -1 && in_order[d] < pos ) conflicts.push_back(d);
    for( int d : revV[v] ) if( in_order[d] != -1 && in_order[d] >= pos ) conflicts.push_back(d);

    assert( set<int>(ALL(conflicts)).size() == conflicts.size() );

    // removes from order all elements in [conflicts]
    for( int d : conflicts ) removeFromOrder(d);

    insertIntoOrder( v, pos );

//    if(use_fast_updates){ // #TEST
//        for(int i=0; i<N; i++) {
//            if( in_order[i] != -1 ) continue;
//            auto old_delta = delta[i];
//            auto old_ins_point = insertion_point[i];
//            bool valid = is_valid_delta[i];
//
//            updateNode(i);
//
//            if (valid) {
//                assert(delta[i] == old_delta);
//                assert(insertion_point[i] == old_ins_point);
//            }
//        }
//    }

    if(!use_fast_updates){
        conflicts.push_back(v);
        updateValidDelta(conflicts);
    }

}

void SALS2::updateNode(int v) {
    int m  = 1e9, M = -1;
    insertion_point[v] = {-1,-1};

    for( int d : V[v] ){
        if( in_order[d] != -1 ){
            if(!use_fast_updates) m = min( m, in_order[d] );
            else {
                if (in_order[d] < m) {
                    m = in_order[d];
                    insertion_point[v].first = d;
                }
            }
        }
    }
    for( int d : revV[v] ){
        if( in_order[d] != -1 ){
            if(!use_fast_updates) M = max( M, in_order[d] );
            else {
                if (in_order[d] > M) {
                    M = in_order[d];
                    insertion_point[v].second = d;
                }
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


void SALS2::updateValidDelta(VI &conflicts) {
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

unsigned long long SALS2::xorshf96() {          //period 2^96-1
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

int SALS2::getRandomInt(int mod) {
    return xorshf96() % mod;
//    return rnd.nextInt(mod);
//    return unif(rng) % mod;

    LL res = random_ints[random_int_index];
    random_int_index++;
    if(random_int_index == random_ints.size()) random_int_index = 0;
    return res % mod;
}

double SALS2::getRandomDouble() {
    return (double)(xorshf96() & ( (1ll<<31)-1 )) / ((1ll<<31)-1); // this is roughly number of the form   a / 2^31
//    return rnd_d.rand();
//    return unif_d(rng);

    double res = random_doubles[random_double_index];
    random_double_index++;
    if(random_double_index == random_doubles.size()) random_double_index = 0;
    return res;
}

void SALS2::preprocessRandoms() {
    return;

    { // creating lists of random numbers
        int MUL = 13;

        int E = GraphUtils::countEdges(V,true);

        random_ints = VLL(E * MUL);
        random_int_index = 0;
        for(int i=0; i<random_ints.size(); i++) random_ints[i] = rnd.rand();
//        for(int i=0; i<random_ints.size(); i++) random_ints[i] = unif(rng);

        random_doubles = vector<double>(E * MUL);
        random_double_index = 0;
        for(int i=0; i<random_doubles.size(); i++) random_doubles[i] = rnd_d.rand();
//        for(int i=0; i<random_doubles.size(); i++) random_doubles[i] = unif_d(rng);
    }
}

VI SALS2::getCurrentDFVS() {
    VI dfvs;
    for( int i=0; i<N; i++ ) if( in_order[i] == -1 ) dfvs.push_back(i);
    if(dfvs.size() != dfvs_size){
        DEBUG(dfvs.size());
        DEBUG(dfvs_size);
        assert(dfvs.size() == dfvs_size);
    }
    return dfvs;
}

double SALS2::findInitialTemperature(VVI V, VI& dfvs, Config c, int max_mvt, bool use_fast_updates) {
    double T0 = 0.33;
    double alpha = 0.99; // anything here, it does not matter

    int N = V.size();

    while( T0 >= 0.1 ){

        SALS2 sals2( V, dfvs, c, use_fast_updates );
        sals2.cnf.write_logs = false;

//        sals2.find_init_temperature_mode = true; // #TEST - this will not terminate after 1 iteration if improvement was done
//        auto temp = sals2.localSearch( T0, alpha, max_mvt, 10, 1 ); // #TEST

        auto temp = sals2.localSearch( T0, alpha, max_mvt, 1, 1 );
        if( temp.size() < dfvs.size() ){
            dfvs = temp;
//            clog << "Returning initial temperature T0: " << T0 << endl;
            return T0;
        }

        T0 -= 0.0075;
    }

    return -1; // returning default initial temperature
}







