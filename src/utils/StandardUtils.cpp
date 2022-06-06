//
// Created by sylwester on 9/6/19.
//

#include <utils/StandardUtils.h>

namespace StandardUtils{

    VI toVI(VB &v) {
        VI res; res.reserve( max( 5, (int)v.size() / 2 ) );
        for(int i=0; i <v.size(); i++) if( v[i] ) res.push_back(i);
        return res;
    }

    VB toVB(int N, VI &v) {
        VB res(N,false);
        for(int d : v) res[d] = true;
        return res;
    }

    vector<string> split(string &s, string pat) {
        while( !s.empty() && s.back() == ' ' ) s.pop_back();
        if(s.empty()) return {};

        if( pat.size() == 1 ){
            vector<string> res;
            int p = 0, q = 1;
            while( p < s.size() ){
                while( q < s.size() && s[q] != pat[0] ) q++;
                res.push_back(s.substr(p, q-p));
                p = q;
                q++;
            }
            return res;
        }else{
            clog << "Split not implemented yet for patterns with length > 1" << endl;
            exit(18);
        }
    }

    VVI partitionToLayers(VI &partition) {
        int M = *max_element(ALL(partition));
        VVI v(M+1);
        for( int i=0; i<partition.size(); i++ ) if( partition[i] != -1 ) v[partition[i]].push_back(i);
        return v;
    }

    VI layersToPartition(VVI &layers) {
        int M = -1;
        for( int i=0; i<layers.size(); i++ ) M = max( M, *max_element( ALL(layers[i]) ) );

        VI part(M+1,-1);
        for( int i=0; i<layers.size(); i++ ) for(int t : layers[i]) part[t] = i;
        return part;
    }

    void removeFromArrayPreserveOrderInplace(VI &A, VI &B, VB &helper) {
        for(int d : B) helper[d] = true;
        int p = 0;
        for( int i=0; i<A.size(); i++ ){
            if( !helper[A[i]] ){
                A[p] = A[i];
                p++;
            }
        }
        A.resize(p);
        for(int d : B) helper[d] = false;
    }

    void removeFromArrayInplace(VI &A, VI B) {
        if(A.empty()) return;

        sort(ALL(A));
        sort(ALL(B));
        int p = (int)A.size()-1, q = (int)B.size()-1;

        while( p >= 0 && q >= 0 ){
            if( A[p] == B[q] ){
                swap(A[p], A.back());
                A.pop_back();
                p--;
                q--;
            }
            else if( A[p] > B[q] ) p--;
            else q--;
        }
    }

}