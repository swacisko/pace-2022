//
// Created by sylwester on 12/20/21.
//

#include <CONTESTS/PACE22/heur/AgentFlow.h>
#include <utils/RandomNumberGenerators.h>

AgentFlow::AgentFlow(VVI &v) : V(v) {
    N = v.size();
}

VD AgentFlow::flowContinuous(int iterations, int void_steps) {
    VD total(N,0);
    VD tokens(N,1);
    VD new_tokens(N,0);

    for(int i=0; i<iterations; i++){
        fill(ALL(new_tokens),0);

        for( int a = 0; a < N; a++ ){
            if(V[a].size() == 0) continue;
            double d = tokens[a] / V[a].size();
            for( int b : V[a] ) new_tokens[b] += d;
        }

        if( i >= void_steps ){
            for(int a=0; a<N; a++) total[a] += new_tokens[a];
        }

        swap(tokens,new_tokens);
    }

    return total;
}

VD AgentFlow::flowTokens(int iterations, int void_steps, int initial_tokens) {
    LL MAX_VAL = 1'000'000'000ll * 1'000'000;
    UniformIntGenerator rnd(0,MAX_VAL);

    VD total(N,0);
    VI tokens(N,initial_tokens);
    VI new_tokens(N,0);

    for(int i=0; i<iterations; i++){
        fill(ALL(new_tokens),0);

        for( int a = 0; a < N; a++ ){
            if(V[a].empty()) continue;
            for( int t = 0; t < tokens[a]; t++ ){
                int ind = rnd.nextInt( V[a].size() );
                int b = V[a][ind];
                new_tokens[b]++;
            }
        }

        if( i >= void_steps ){
            for(int a=0; a<N; a++) total[a] += new_tokens[a];
        }

        swap(tokens,new_tokens);
    }

    return total;
}

