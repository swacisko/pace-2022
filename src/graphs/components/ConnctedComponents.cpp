//
// Created by sylwester on 3/16/20.
//

#include <graphs/GraphUtils.h>
#include <graphs/components/ConnectedComponents.h>


void ConnectedComponents::dfs(VVI& V, int &num, int &par, VB& was, VVI & comps){
    was[num] = true;
    comps.back().push_back(num);

    for( int t : V[num] ){
        if( t != par && !was[t] ) dfs(V,t,num,was,comps);
    }

};

VVI ConnectedComponents::getConnectedComponents( VVI & V ){
    VVI comps;
    VB was(V.size(),false);

    for(int i=0; i<V.size(); i++){
        if(!was[i]){
            comps.push_back( getConnectedComponentForNode(V,i,was) );
        }
    }

    return comps;
}


VI ConnectedComponents::getConnectedComponentForNode(VVI &V, int v, VB &was) {
    VI comp;

    function< void(int,int) > dfsLambda = [&dfsLambda,&was,&V,&comp](int num, int par){
        if(was[num]) return;
        was[num] = true;
        comp.push_back(num);
        for( int t : V[num] ){
            if( t != par && !was[t] ) dfsLambda(t,num);
        }
    };

    dfsLambda(v,v);


    return comp;
}
