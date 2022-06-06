//
// Created by sylwester on 9/6/19.
//

#ifndef ALGORITHMSPROJECT_STANDARDUTILS_H
#define ALGORITHMSPROJECT_STANDARDUTILS_H

#include "Makros.h"
#include "RandomNumberGenerators.h"

namespace StandardUtils{

    extern void removeFromArrayPreserveOrderInplace(VI & A, VI & B, VB & helper );

    extern void removeFromArrayInplace(VI & A, VI B );

    extern VI toVI(VB& v);

    extern VB toVB(int N, VI& v);

    extern VVI partitionToLayers( VI & partition );

    extern VI layersToPartition( VVI & layers );

    vector<string> split( string & s, string pat = " " );

    template<class _T>
    void append( vector<_T> & v1, vector<_T> & v2 ){ v1.insert( v1.end(), ALL(v2) ); }

    template<class _T>
    vector<_T> slice( vector<_T> & v, int a, int b, int step=1 ){
        assert(step != 0);
        vector<_T> res;
        res.reserve( abs( (b-a) / step ) );
        if(a<=b){
            assert(step > 0);
            while(a<b){
                res.push_back( v[a] );
                a += step;
            }
        }
        else{
            assert(step < 0);
            while(a>b){
                res.push_back(v[a]);
                a += step;
            }
        }
        return res;
    }


    template<class _T, class _rnd>
    void shuffle( vector<_T> & V, _rnd rnd ){
        std::uniform_int_distribution<long long> unif( 0, 10ll * V.size() );

        for( int i=(int)V.size()-1; i>=0; i-- ){
            int ind = unif(rnd) % (i+1);
            if( ind != i ) swap( V[i], V[ind] );
        }
    }

    template<class _T>
    void shuffle( vector<_T> & V ){
        UniformIntGenerator rnd(0,1e9);
        shuffle(V, rnd.getRNG());
    }

    template<class _t, class _s>
    vector< pair<_t,_s> > zip( vector<_t> a, vector<_s> b ){
        vector<pair<_t,_s>> res(min(a.size(), b.size()));
        for( int i=0; i<min(a.size(), b.size()); i++ ) res[i] = {a[i], b[i]};
        return res;
    }

    template<class _t, class _s>
    pair<vector<_t>, vector<_s> > unzip( vector< pair<_t,_s> > &&a ){
        pair<vector<_t>, vector<_s> > res;
        res.first.reserve(a.size());
        res.second.reserve(a.size());
        for( int i=0; i<a.size(); i++ ){
            res.first.push_back(a[i].first);
            res.second.push_back(a[i].second);
        }
        return res;
    }

    template<class T>
    vector<pair<T,T>> product( vector<T> a, vector<T> b ){
        vector<pair<T,T>> res;
        res.reserve(a.size() * b.size());

        for( int i=0; i<a.size(); i++ ){
            for(int j=0; j<b.size(); j++){
                res.emplace_back( a[i], b[j] );
            }
        }
        return res;
    }

    template<class _T>
    void makeUnique( vector<_T> & v ){
        sort(ALL(v));
        v.resize( unique(ALL(v)) - v.begin() );
    }

    template<class _T>
    bool find( vector<_T> & v, _T el ){ return find(ALL(v),el) != v.end(); }

    template<class _T, class _S>
    VI setUnion( _T A, _S B, VB & helper ){
        VI res = A;
        for( auto a : A ) helper[a] = true;
        for( auto b : B ) if(!helper[b]){ helper[b] = true; res.push_back(b); }
        for( auto a : A ) helper[a] = false; for( auto a : B ) helper[a] = false;
        return res;
    }

    template<class _T, class _S>
    VI setIntersection( _T A, _S B, VB & helper ){
        VI res;
        for( auto a : A ) helper[a] = true;
        for( auto b : B ) if(helper[b]){ helper[b] = false; res.push_back(b); }
        for( auto a : A ) helper[a] = false; for( auto a : B ) helper[a] = false;
        return res;
    }

    template<class _T, class _S>
    VI setDifference( _T A, _S B, VB & helper ){
        VI res;
        for( auto b : B ) helper[b] = true;
        for( auto a : A ) if(!helper[a]){ helper[a] = true; res.push_back(a); }
        for( auto a : A ) helper[a] = false; for( auto a : B ) helper[a] = false;
        return res;
    }
}

#endif //ALGORITHMSPROJECT_STANDARDUTILS_H
