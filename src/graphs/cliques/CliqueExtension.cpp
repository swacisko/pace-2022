//
// Created by sylwester on 3/25/20.
//

#include <graphs/cliques/CliqueExtension.h>

#include "graphs/cliques/CliqueExtension.h"
#include "graphs/cliques/CliqueUtils.h"
#include "graphs/GraphUtils.h"
#include "graphs/GraphInducer.h"
#include "graphs/GraphReader.h"


VI CliqueExtension::maximizeCliqueGreedy(VVI& V, VI clq){
    assert( CliqueUtils::isClique(V,clq) );

    VI A = clq;
    unordered_set<int> zbA(ALL(A));
    map<int,int> deg;
    for( int a : A ){
        for(int d : V[a]){
            if( !zbA.count(d) ) deg[d]++;
        }
    }

    VI B;
    for( PII p : deg ){
        if( p.second == A.size() ) B.push_back(p.first);
    }

    random_shuffle(ALL(B));

    VI inducer = A; inducer.insert( inducer.end(), ALL(B) );

    InducedGraph g = GraphInducer::induce( V, inducer );
    for( int & a : A ) a = g.perm[a];
    for( int & b : B ) b = g.perm[b];

    int n = g.V.size();
    VI neighCnt( n,0 );
    VB inA(n,false), inB(n,false);
    for(int a : A) inA[a] = true;
    for(int b : B) inB[b] = true;

    for( int p : B ){
        for(int d : g.V[p]){
            if( inB[d] ) neighCnt[d]++;
        }
    }

    VB Nx(n,false);

    auto comp =  [&neighCnt]( int a, int b ){ return neighCnt[a] < neighCnt[b]; };

    while( !B.empty() ){

        if( neighCnt[ *min_element( ALL(B), comp ) ] == (int)B.size()-1 ){
            A.insert(A.end(), ALL(B));
            break;
        }

        auto it = max_element( ALL(B), comp );
        int x = *it;

        A.push_back(x);
        inB[x] = false;
        inA[x] = true;

        for( int d : g.V[x] ) Nx[d] = true;
        for( int i=(int)B.size()-1; i>=0; i-- ){
            int b = B[i];
            if( !Nx[b] ){
                swap( B[i], B.back() );
                B.pop_back();
                inB[b] = false;

                for( int d : g.V[b] ) neighCnt[d]--;

            }else{
                neighCnt[b]--;
            }
        }
        for( int d : g.V[x] ) Nx[d] = false;
    }

    for(int& a : A) a = g.nodes[a];

    assert(CliqueUtils::isClique(V,A));

    return A;
}


VI CliqueExtension::findMaximalNodeCliqueExtension(VVI &V, bool sparse_check) {
    VI clq, bestClq;

    VI to_check(V.size());
    iota(ALL(to_check), 0);
    if(sparse_check){
        random_shuffle(ALL(to_check));
        to_check.resize(1 + pow(V.size(),0.66) );
        sort(ALL(to_check));
    }

    for(int i : to_check){

        if( !V[i].empty() ) clq = maximizeCliqueGreedy(V, {i} );
        else clq = {i};

        if( clq.size() > bestClq.size() ) bestClq = clq;
    }
    return bestClq;
}


