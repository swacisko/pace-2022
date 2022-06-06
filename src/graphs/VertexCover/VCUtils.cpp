//
// Created by sylwester on 8/8/19.
//

#include <graphs/VertexCover/VCUtils.h>
#include <graphs/VertexCover/kernelization/KernelizerVC.h>
#include <graphs/GraphUtils.h>
#include <utils/RandomNumberGenerators.h>
#include <graphs/VertexCover/approximation/LibMVC/fastvc.h>
#include <graphs/VertexCover/approximation/NuMVC/NuMVC.h>
#include <combinatorics/CombinatoricUtils.h>
#include "CollectionOperators.h"

bool VCUtils::isVertexCover( VVI & V, VI & vc ){
    VB inVC(V.size(),false);
    for(int p : vc) inVC[p] = true;
    for( int i=0; i<V.size(); i++ ){
        if(inVC[i]) continue;
        for(int d : V[i]){
            if( !inVC[d] ) return false;
        }
    }
    return true;
}


bool VCUtils::isIndependentSet(VVI &V, VI &S) {
    VI perm(V.size());
    iota(ALL(perm),0);
    VI S2 = S;
    sort(ALL(S2));
    VI vc;
    set_difference( ALL(perm), ALL(S2),  std::back_inserter(vc) );
    bool res = isVertexCover( V,vc );
    return res;
}

VI VCUtils::getMinCVUsingFastVC(VVI V, int milliseconds, VI init_vc) {
    if( !init_vc.empty() ){
        if( !isVertexCover(V, init_vc) ){
            DEBUG(GraphUtils::countEdges(V));
            DEBUG(init_vc.size());
            assert(false);
        }
    }

    int N = V.size();

    bool use_numvc = ( N <= 20'000 );
    bool use_fastvc = (!use_numvc);
    VI mis;

    VI kern_vc;
    const bool use_kernelization = true;
    if(use_kernelization) {
        KernelizerVC kern;
        auto[x, y] = kern.initialKernelization(V);
        kern_vc = x;
        int E = GraphUtils::countEdges(V);
        if(E == 0) return kern_vc;
    }

    if(use_fastvc){
        UniformIntGenerator rnd(0,1e6);
        int seed = rnd.rand();
        VPII edges = GraphUtils::getGraphEdges(V);
        libmvc::FastVC fastvc(edges, N, 0, std::chrono::milliseconds(milliseconds), false, seed, init_vc);
        fastvc.cover_LS();

        mis = fastvc.get_independent_set(false);
        assert(VCUtils::isIndependentSet(V,mis));
    }
    else{ // using NuMVC
        NuMVC numvc;
        mis = numvc.solve(V, 1.0 * milliseconds / 1000, init_vc); // CAUTION! NuMVC indices are from 1 to N
        for (int &d : mis) d--;
        assert(VCUtils::isIndependentSet(V, mis));
    }

    if(use_kernelization) {
        VI vc = CombinatoricUtils::getFullSetDifference(V.size(), mis);
        vc += kern_vc;
        mis = CombinatoricUtils::getFullSetDifference(V.size(), vc);
    }

    VI vc = CombinatoricUtils::getFullSetDifference(N,mis);
    return vc;
}


