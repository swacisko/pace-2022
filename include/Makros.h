//
// Created by sylwester on 8/7/19.
//

#ifndef ALGORITHMSPROJECT_MAKROS_H
#define ALGORITHMSPROJECT_MAKROS_H



#include<cstdio>
#include<iostream>
#include<vector>
#include<string>
#include<map>
#include<complex>
#include<stack>
#include<list>
#include<bitset>
#include<set>
#include<unordered_set>
#include<unordered_map>
#include<iterator>
#include<cmath>
#include<queue>
#include<ctime>
#include<string.h>
#include<fstream>
#include<sstream>
#include<algorithm>
#include <numeric>
#include<chrono>
#include<random>
#include<functional>
#include<utility>
#include <assert.h>

using namespace std;

#define REP( x,y ) for( int x=0; x<(y); ++x )
#define FORD( x,y,z ) for( int x=y; x>=(z); --x )
#define FOR(x,b,e) for( int x = b; x <= (e); ++x )
#define SIZE(v) (int)v.size()
#define ALL(c) c.begin(),c.end()
#define VAR(v,n) __typeof(n) v=(n)
#define FOREACH(i,c) for( VAR(i,c.begin());i!=c.end();++i )
#define PB push_back
#define MP make_pair
#define ST first
#define ND second
#define WRITE( V ){ FOREACH(it,V) cerr << *it << ", "; cerr << endl; }

#define DEBUG_NOW_AND_HERE true
#define DEBUG(x) if( DEBUG_NOW_AND_HERE ) clog << #x << ": " << x << endl;

#define WRP(p) "(" << p.ST << "," << p.ND << ")"
#define WRITE_ALL(V,s,t) { cerr << s << endl;  REP( i,SIZE(V) ){ cerr  << i+t << " ---- ";  FOREACH(it,V[i]) cerr << *it+t << ", "; cerr << endl;     } }
#define ENDL(x) REP(crow,(x)) clog << endl;
#define ENDLS(x,c) REP(crow,(x)) clog << c << flush;


typedef long long LL;

typedef vector<int> VI;
typedef vector< VI > VVI;
typedef vector<VVI> VVVI;

typedef vector<double> VD;

typedef vector<bool> VB;
typedef vector< VB > VVB;

typedef pair<int,int> PII;
typedef vector<PII> VPII;
typedef vector<VPII> VVPII;

typedef vector<LL> VLL;


template<class _T, class _E>
ostream& operator<<( ostream& str, const pair<_T,_E> & pair){
    str << "(" << pair.first << "," << pair.second << ")";
    return str;
}

template<class _T>
void writeCollectionToStream(ostream& str, _T& col ){
    str << "{";
    int ile = 0;
    for(auto t : col){
        if(ile++ > 0) str << ", ";
        str << t;
    }
    str << "}";
}

template<class _P, class _Q, class _R>
void writeCollectionToStream(ostream& str, tuple<_P,_Q,_R> col ){
    str << "{" << get<0>(col) << ", " << get<1>(col) << ", " << get<2>(col) << "}";
}



template<class _T>
ostream& operator<<( ostream& str, vector<_T> v ){ writeCollectionToStream(str, v); return str; }

template<class _T>
ostream& operator<<( ostream& str, deque<_T> v ){ writeCollectionToStream(str, v); return str; }

template<class _P, class _Q, class _R>
ostream& operator<<( ostream& str, tuple<_P,_Q,_R> v ){ writeCollectionToStream(str, v); return str; }

template<class _T>
ostream& operator<<( ostream& str, set<_T> v ){ writeCollectionToStream(str, v); return str; }

template<class _T>
ostream& operator<<( ostream& str, unordered_set<_T> v ){ writeCollectionToStream(str, v); return str; }

template<class _T, class _H>
ostream& operator<<( ostream& str, unordered_set<_T, _H> v ){ writeCollectionToStream(str, v); return str; }

template<class _T, class _E>
ostream& operator<<( ostream& str, map<_T, _E> v ){ writeCollectionToStream(str, v); return str; }

template<class _T, class _E>
ostream& operator<<( ostream& str, unordered_map<_T, _E> v ){ writeCollectionToStream(str, v); return str; }

template<class _T, class _E, class _H>
ostream& operator<<( ostream& str, unordered_map<_T, _E, _H> v ){ writeCollectionToStream(str, v); return str; }

template<class _T>
struct Triple{
    Triple(){}
    Triple(_T first, _T second, _T third){ st = first; nd = second; rd = third; }
    _T st,nd,rd;
};

template<class _T>
struct Quadruple{
    Quadruple(){}
    Quadruple(_T first, _T second, _T third, _T foruth){  st = first; nd = second; rd = third; th = foruth; }
    _T st,nd,rd, th;
};




#endif //ALGORITHMSPROJECT_MAKROS_H
