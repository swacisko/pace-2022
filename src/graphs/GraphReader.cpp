//
// Created by sylwester on 8/8/19.
//

#include "graphs/GraphReader.h"


namespace GraphReader{


    VVI readGraphStandardEdges(istream &cin, bool directed ) {
        int N,M;
        cin >> N >> M;
        VVI V(N);
        for(int i=0; i<M; i++){
            int a,b;
            cin >> a >> b;

            clog << "CATUION - changed readGraphStandardEdges()" << endl;
            a++, b++; // #TEST!!

            V[a-1].push_back(b-1);
            if( !directed ) V[b-1].push_back(a-1);
        }

        return V;
    }


    VVI readGraphDIMACSWunweighed(istream &cin, bool edgeFoolowE ) {
        string s;
        VVI V;

        int N,M;
        int edges_read = 0;
        int cnt = 0;

        while( true ){
            getline(cin,s,char(10));

            if( s[0] == 'c' ){
                // nothing to do here, this is a comment
            }else if( s[0] == 'p' ){
                stringstream str(s);
                string nothingBox;
                str >> nothingBox >> nothingBox >> N >> M;

                V = VVI(N);
            }else{
                stringstream str(s);
                int a,b;
                char e;

                if(edgeFoolowE) str >> e >> a >> b;
                else str >> a >> b;

                edges_read++;

                // CAUTION - THIS SHOULD BE HERE, COMMENTED ONLY FOR TESTING
                a--;
                b--;

                V[a].push_back(b);
                V[b].push_back(a);

                if( edges_read == M ) break;
            }

            if(cnt++ > 1e9){
                ENDL(5);
                clog << "GraphReader endless loop, M = " << M <<", but only "
                << edges_read << " edges were read" << endl;
                ENDL(5);
                break;
            }
        }

        for( int i=0; i<N; i++ ){
            sort(ALL(V[i]));
            V[i].resize( unique(ALL(V[i])) - V[i].begin() );
        }

        return V;
    }

}