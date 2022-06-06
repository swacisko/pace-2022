//
// Created by sylwester on 8/27/19.
//

#include "combinatorics/CombinatoricUtils.h"
#include <utils/RandomNumberGenerators.h>


namespace CombinatoricUtils{


    VI getRandomPermutation(int N){
        UniformIntGenerator rnd(0,1e9);
        return getRandomPermutation(N, rnd.rand());
    }

    VI getRandomPermutation(int N, unsigned seed){
        VI perm(N);
        iota(ALL(perm),0);
        shuffle(ALL(perm), mt19937(seed) );
        return perm;
    }


    VI getRandomSequence( int U, int N ) {
        VI seq(N);
        UniformIntGenerator g( 0,U );
        for(int i=0; i<N; i++) seq[i] = g.rand();
        return seq;
    }


    VI getRandomSubset( int U, int L ){
        UniformIntGenerator rnd(0,1e9);
        return getRandomSubset(U,L,rnd.rand());
    }

    VI getRandomSubset( int U, int L, unsigned seed ){
        UniformIntGenerator gen(0,U, seed);

        VI res;
        if( L < U / 20 ){
            unordered_set<int> zb;
            while( zb.size() < L ) zb.insert( gen.rand() );
            res = VI(ALL(zb));
        }else{
            VI perm = getRandomPermutation(U+1, seed);
            if( L < perm.size() ) return VI( perm.begin(), perm.begin() + L );
            else return VI( ALL(perm) );
        }
        return res;
    }

    VI getFullSetDifference(int N, VI A) {
        VB helper(N,false);
        for(int d : A) helper[d] = true;
        VI res; res.reserve(N-A.size());
        for(int i=0; i<N; i++) if(!helper[i]) res.push_back(i);
        return res;
    }

}

