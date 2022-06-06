//
// Created by sylwester on 4/8/22.
//



#include <CONTESTS/PACE22/heur/HittingSetLS.h>
#include <utils/RandomNumberGenerators.h>
#include <combinatorics/CombinatoricUtils.h>
#include <utils/StandardUtils.h>
#include <datastructures/RandomSelectionSet.h>


HittingSetLS::HittingSetLS(VVI sets, VI hs, Config cnf) :
    rnd(0ll, 1'000'000'000ll * 1'000'000'000ll)
{
    this->sets = sets;
    this->hs = hs;
    this->cnf = cnf;

    N = 0;
    for(VI C : sets) for(int d : C) N = max(N,d+1);

    S = sets.size();
    V = VVI(N); // V[i] is the list of ids of sets that contain node i

    {
        VI degs(N,0);
        for( int i=0; i<S; i++ ) for( int d : sets[i] ) degs[d]++;
        for(int i=0; i<N; i++) V[i].reserve(degs[i]);
    }

    for( int i=0; i<S; i++ ){
        for( int d : sets[i] ){
            V[d].push_back(i);
        }
    }

    hashes = VLL(N);
    rnd = UniformIntGenerator(0, 1'000'000'000ll * 1'000'000'000);
    for(int i=0; i<N; i++) hashes[i] = rnd.rand();
}


VI HittingSetLS::hsImprovementLS2( int iters, int lower_bound, int deviate_perm_frequency,
                                   int add_violation_frequency, int remove_violation_frequency ) {
    if(hs.empty()) return hs;

    constexpr bool use_visited_states = true;

    int uncovered = 0;
    in_hs = StandardUtils::toVB(N, hs);
    VI valid_hs;

    // covered_by[i] is the number of elements from current hitting set that belong to sets[i]
    covered_by = VI(S,0);

    for( int v : hs ) for( int d : V[v] ) covered_by[d]++;
    for( int i=0; i<S; i++ ) if( covered_by[i] == 0 ) uncovered++;



    unordered_set<LL> visited_states;
    LL current_state_hash = 0;
    for(int d : hs) current_state_hash ^= hashes[d];
    visited_states.insert(current_state_hash);

    VI perm = CombinatoricUtils::getRandomPermutation(N);



    VI covers_alone(N,0); // the number of sets that a node covers ALONE (without any other node)
    for( int i=0; i<S; i++ ){
        if(covered_by[i] == 1){
            for( int d : sets[i] ) if(in_hs[d]) covers_alone[d]++;
        }
    }
    auto cmp_hs = [&]( int a, int b ){
        if( covers_alone[a] != covers_alone[b] ) return covers_alone[a] < covers_alone[b];
        else return perm[a] > perm[b]; // anything here
    };
    set<int,decltype(cmp_hs)> hs_nodes(cmp_hs);
    for( int v : hs ) hs_nodes.insert(v);



    VI hits_uncovered(N,0); // number of uncovered sets that node hits
    for( int i=0; i<S; i++ ){
        if(covered_by[i] == 0){
            for(int d : sets[i]) hits_uncovered[d]++;
        }
    }
    auto cmp_nonhs = [&](int a, int b){
        if( hits_uncovered[a] != hits_uncovered[b] ) return hits_uncovered[a] > hits_uncovered[b];
        else return perm[a] < perm[b]; // anything here
    };
    set<int,decltype(cmp_nonhs)> nonhs_nodes(cmp_nonhs);
    for( int i=0; i<N; i++ ) if( !in_hs[i] ) nonhs_nodes.insert(i);


    // greedily select the best node to remove from current state - the one that uncovers smallest number of cycles
    auto getNodeToRemove = [&](){
        int res = -1;

        int cnt = 0;
        for( int u : hs_nodes ){
            if( use_visited_states && visited_states.count( current_state_hash ^ hashes[u] ) ) continue;
            res = u;
            break;
        }

        if(res == -1){
            do{
                res = rnd.nextInt(N);
            }while( !in_hs[res] );
        }
        return res;
    };


    // greedily select the best node to add to current hs - the one that will cover greatest number of uncovered sets
    auto getNodeToAdd = [&](){
        int res = -1, val = -1;

        for( int u : nonhs_nodes ){
            if( use_visited_states && visited_states.count( current_state_hash ^ hashes[u] ) ) continue;
            res = u;
            break;
        }

        if(res == -1) {
            do{
                res = rnd.nextInt(N);
            }while( in_hs[res] );
        }
        return res;
    };

    VB was(N,false);


    const int freq_remove = remove_violation_frequency;
    const int freq_add = add_violation_frequency;

    const bool use_continuous_perm_deviation = cnf.hsls_use_continuous_perm_deviation;

    ls_loop:
    for( int itr = 0; itr < iters; itr++ ){

//        if( (itr & 63) == 0 && cnf.sw.tle("main") ) break;
        if( (itr & 63) == 0 && cnf.tle() ) break;

        if(use_continuous_perm_deviation){
            int a = rnd.nextInt(N);
            int b = rnd.nextInt(N);
            while(b == a) b = rnd.nextInt(N);

            if(in_hs[a]) hs_nodes.erase(a);
            else nonhs_nodes.erase(a);

            if(in_hs[b]) hs_nodes.erase(b);
            else nonhs_nodes.erase(b);

            swap(perm[a], perm[b]);

            if(in_hs[a]) hs_nodes.insert(a);
            else nonhs_nodes.insert(a);

            if(in_hs[b]) hs_nodes.insert(b);
            else nonhs_nodes.insert(b);
        }

        {
            if (itr % deviate_perm_frequency == 0) {
                hs_nodes.clear();
                nonhs_nodes.clear();
                StandardUtils::shuffle(perm);
                for (int i = 0; i < N; i++) {
                    if (!in_hs[i]) nonhs_nodes.insert(i);
                    else hs_nodes.insert(i);
                }
            }
        }

        // this can be used to skip adding node back to HS if a valid has was found. This way
        // we may not terminate after HS was found, but try to improve it further
        bool add_back = (uncovered > 0 );

        //******************************//****************************** removing node

        { // remove node from current hs
            int v = getNodeToRemove();
            if( itr % freq_remove == 0 ) do{ v = rnd.nextInt(N); } while( !in_hs[v] );

            assert(in_hs[v]);
            current_state_hash ^= hashes[v];
            if(use_visited_states) visited_states.insert(current_state_hash);

            assert( hs_nodes.count(v) );
            hs_nodes.erase(v); // need to do that before reducing [covered_by]
            covers_alone[v] = 0;

            in_hs[v] = false;
            for( int d : V[v] ){
                covered_by[d]--;
                if(covered_by[d] == 0) uncovered++;
            }

            { // updating hs_nodes other than v (v was updated earlier)
                VI temp;
                for (int d : V[v]) {
                    if (covered_by[d] == 1) {
                        int u = -1;
                        for (int uu : sets[d]) if(in_hs[uu]) { u = uu; break; }
                        assert(u != -1);
                        assert(u != v);

                        if (!was[u]) {
                            assert(hs_nodes.count(u));
                            hs_nodes.erase(u);
                            was[u] = true;
                            temp.push_back(u);
                        }
                        covers_alone[u]++;
                    }
                }
                for (int u : temp) hs_nodes.insert(u);
                for(int d : temp) was[d] = false;
            }

            { // updating nonhs_nodes
                VI temp;

                for(int d : V[v]){
                    if(covered_by[d] == 0) {
                        for (int u : sets[d]) {
                            assert(!in_hs[u]);
                            if (!was[u] ) {
                                if(u != v) assert( nonhs_nodes.count(u) );
                                nonhs_nodes.erase(u);
                                was[u] = true;
                                temp.push_back(u);
                            }

                            hits_uncovered[u]++;
                        }
                    }
                }

                for(int u : temp) nonhs_nodes.insert(u);
                nonhs_nodes.insert(v);
                for(int d : temp) was[d] = false;
            }
        }

        //******************************//****************************** adding node back

        if(add_back){ // add node to current hs
            int v = getNodeToAdd();
            if( itr % freq_add == 0 ) do{ v = rnd.nextInt(N); } while( in_hs[v] );

            assert(!in_hs[v]);
            current_state_hash ^= hashes[v];
            if(use_visited_states) visited_states.insert(current_state_hash);

            { // update hs_nodes
                assert(nonhs_nodes.count(v));
                nonhs_nodes.erase(v);

                VI temp;
                for (int d : V[v]) {
                    if (covered_by[d] == 1) {
                        int u = -1;
                        for (int uu : sets[d]) if(in_hs[uu]) { u = uu; break; }
                        assert(u != -1);
                        assert(u != v);

                        if (!was[u]) {
                            hs_nodes.erase(u);
                            was[u] = true;
                            temp.push_back(u);
                        }
                        covers_alone[u]--;
                    }else if( covered_by[d] == 0 ){
                        covers_alone[v]++;
                    }
                }
                for (int u : temp) hs_nodes.insert(u);
                hs_nodes.insert(v);
                for(int d : temp) was[d] = false;
            }

            { // udate nonhs_nodes
                VI temp;
                for(int d : V[v]){
                    if( covered_by[d] == 0 ){
                        for( int u : sets[d] ){
                            if(!was[u]){
                                was[u] = true;
                                temp.push_back(u);
                                nonhs_nodes.erase(u);
                            }
                            hits_uncovered[u]--;
                        }
                    }
                }

                for( int u : temp ) if(u != v) nonhs_nodes.insert(u);
                for(int d : temp) was[d] = false;
            }

            in_hs[v] = true;
            for( int d : V[v] ){
                if(covered_by[d] == 0) uncovered--;
                covered_by[d]++;
            }
        }

        if(uncovered == 0){
            valid_hs = StandardUtils::toVI(in_hs);
            { // just an assertion if it is really valid hs
                fill(ALL(covered_by),0);
                for(int v : valid_hs) for( int d : V[v] ) covered_by[d]++;
                for( int i=0; i<S; i++ ) assert( covered_by[i] > 0 );
            }
//            clog << "Found valid hs of size: " << valid_hs.size() << endl;

            if( valid_hs.size() == lower_bound ) break;
        }
    }

    if( persistent_search_iterations_left > 0 && valid_hs.empty() ){
        persistent_search_iterations_left--;
        iters = iters * 1.25;
        lower_bound = hs.size()-1;
        if(cnf.write_logs){
            clog << "Persistent search, iters: " << iters << ", persistent_search_iterations_left"
                 << persistent_search_iterations_left << endl;
        }
        goto ls_loop;
    }

    return valid_hs;
}

