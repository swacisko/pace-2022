//
// Created by sylwester on 3/26/22.
//

#include "datastructures/RandomSelectionSet.h"

RandomSelectionSet::RandomSelectionSet(int n) : rnd(0, 1'000'000'000ll * 1'000'000'000ll) {
    n++;
    N = 1;
    while( N < n ) N <<= 1;

    probab = VLL(2*N,0);
    cover = VLL(2*N,0);
    max_el = VLL(2*N,0);
    cnt_num = VI(2*N,false);
}


int RandomSelectionSet::getRandomElement() {
    LL max_val = max_el[1];
    if( max_val == 0 ) return -1;
    LL r = rnd.nextInt(max_val);
    int ind = lowerBound(r);
    return ind;
}

void RandomSelectionSet::set(int id, LL val) {
    if( val == 0 ){
        if( probab[id] > 0 ) {
            int x = N + id;
            while (x >= 1) {
                cnt_num[x]--;
                x = par(x);
            }
        }
    }else if( probab[id] == 0 ){
        int x = N + id;
        while( x >= 1 ){
            cnt_num[x]++;
            x = par(x);
        }
    }

    LL last_val = probab[id];
    probab[id] = val;
    LL diff = val - last_val;

    add( id, diff );
}

void RandomSelectionSet::add(int beg, LL val) {
    beg += N;
    int max_beg_ind = 2*N-1;

    while( beg > 1 ){

        int p = par(beg);
        while( beg == L(p) ){
            beg = p;
            p = par(p);
            max_beg_ind >>= 1;
        }

        cover[beg] += val;
        max_el[beg] += val;

        {
            int x = beg;
            int p = par(beg);
            while( p >= 1 && x == R(p) ){
                max_el[p] = max_el[x] + cover[p];
                x = p;
                p = par(p);
            }
        }

        if( beg < max_beg_ind ) beg++;
        else break;

        beg = par(beg);
        max_beg_ind >>= 1;
    }
}


int RandomSelectionSet::lowerBound(LL val) {
    if( cnt_num[1] == 0 ) return -1; // all elements are 0
    if( get(N-1) < val ) return -1; // largest element is smaller than val

    int id = 1;

    while( id < N ){
        val -= cover[id];

        int l = L(id);
        if( cnt_num[l] > 0 && max_el[l] >= val ) id = l;
        else id = R(id);
    }

    return id - N;
}

LL RandomSelectionSet::get(int id) {
    LL res = 0;
    id += N;
    while( id >= 1 ){
        res += cover[id];
        id >>= 1;
    }
    return res;
}

