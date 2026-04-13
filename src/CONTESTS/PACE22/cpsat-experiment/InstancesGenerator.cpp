//
// Created by sylwe on 08/04/2026.
//

#include "InstancesGenerator.h"

#include "GraphUtils.h"
#include "IntGenerator.h"
#include "StandardUtils.h"

void InstancesGenerator::createHardcodedTests() {

}

vector<string> classes = {"er", "torus", "cyclic"};
// VI Ns = {5'000, 10'000, 20'000, 40'000};
// VD avg_outdegs = {2.5, 5, 10, 20};

// VI Ns = {250, 500, 1'000};
// VD avg_outdegs = {3, 5, 10, 20};

VI Ns = {500};
VD avg_outdegs = {5};

// int A = Ns.size() * avg_outdegs.size(); // 12 instances altogether
// int B = Ns.size() * avg_outdegs.size(); // 12 instances - no need for the 'super dense' one, neighborhood sizes: 4, 4+8=12, 4+8+12=24, 4+8+12+16=40, 4+8+12+16+20 = 60
// // int B = Ns.size() * (degs.size()-1); // 9 instances - no need for the 'super dense' one, neighborhood sizes: 4, 4+8=12, 4+8+12=24, 4+8+12+16=40, 4+8+12+16+20 = 60
// int C = Ns.size() * avg_outdegs.size(); // 12 instances
int total = classes.size() * Ns.size() * avg_outdegs.size();

void InstancesGenerator::createRandomTest(int test_id, ofstream &out_in) {
    test_id -= hardcoded_test_count;

    int A = classes.size();
    int B = Ns.size();
    int C = avg_outdegs.size();

    string class_name = classes[test_id / (B*C) % A];
    int N = Ns[test_id / C % B];
    double deg = avg_outdegs[test_id % C];

    clog << "Creating instance " << test_id << " of class " << class_name << " with N=" << N << " and avg_outdeg=" << deg << endl;
    string instance_name = class_name + "_" + to_string(N) + "_" + to_string((int)deg);
    clog << "\t\tinstance_name: " << instance_name << endl;

    input_files_to_rename[test_id + hardcoded_test_count] = instance_name;

    VVI V(N);

    if (class_name == "er") {
        set<PII> edges;
        int M = deg * N;
        IntGenerator rnd;

        while( edges.size() < M ){
            int a = rnd.nextInt(N);
            int b = rnd.nextInt(N);
            if(a != b) edges.insert({a,b});
        }

        for( auto [a,b] : edges ) V[a].push_back(b);
    }

    if (class_name == "torus") {

        auto randomNeighborhoodClosure = [&](VVI V, int D, int K){
            int N = V.size();
            VVI W(N);
            VB was(N);
            VI neigh;
            VI cur, prev;

            auto findNeighborhood = [&](int beg) {
                was[beg] = true;
                neigh.clear();
                cur.clear();
                prev.clear();
                prev.push_back(beg);

                for ( int d=1; d<=D; d++ ) {
                    cur.clear();
                    for ( int v : prev ) {
                        for ( int u : V[v] ) {
                            if ( !was[u] ) {
                                was[u] = true;
                                cur.push_back(u);
                            }
                        }
                    }
                    neigh += cur;
                    swap(prev,cur);
                }

                was[beg] = false;
                for (int d : neigh) was[d] = false;
            };

            IntGenerator rnd;
            for (int i=0; i<N; i++) {
                findNeighborhood(i);
                for(int d : neigh) assert(d != i);
                StandardUtils::shuffle(neigh, rnd);
                if ( neigh.size() > K ) neigh.resize(K);
                assert( neigh.size() == K );
                W[i] += neigh;
            }

            assert( GraphUtils::isSimple(W) );

            return W;
        };

        auto getGridTorus = [&](int rows, int columns) {

            VVI V( rows*columns );

            // adding edges to the right
            for( int i=0; i<columns-1; i++ ){
                for( int j=0; j<rows; j++ ){
                    int a = j*columns+i;
                    int b = a+1;
                    V[a].push_back(b);
                    V[b].push_back(a);
                }
            }
            for( int j=0; j<rows; j++ ){
                int a = j*columns + columns-1;
                int b = j*columns;
                V[a].push_back(b);
                V[b].push_back(a);
            }

            // adding edges down
            for( int i=0; i<rows-1; i++ ){
                for( int j=0; j<columns; j++ ){
                    int a = i*columns+j;
                    int b = a+columns;
                    V[a].push_back(b);
                    V[b].push_back(a);
                }
            }
            for( int j=0; j<columns; j++ ){
                int a = (rows-1)*columns+j;
                int b = columns;
                V[a].push_back(b);
                V[b].push_back(a);
            }
            return V;
        };

        int rows = sqrt(N);
        int columns = N / rows;
        V = getGridTorus(rows, columns);
        int max_distance = 5;
        V = randomNeighborhoodClosure(V,max_distance,ceil(deg));
    }

    if (class_name == "cyclic") {
        VI perm(N);
        iota(ALL(perm),0);
        IntGenerator rnd;
        StandardUtils::shuffle(perm,rnd);
        int offset = 0;
        int window_size = 3 * max(sqrt(N), deg);
        assert(window_size > deg);
        if (window_size < deg) window_size = deg;

        unordered_set<int> zb;
        for ( int i=0; i<N; i++ ) {
            zb.clear();
            while ( zb.size() < ceil(deg) ) zb.insert(rnd.nextInt(window_size));
            for ( int d : zb ) V[perm[i]].push_back( perm[ (i+1+offset+d) % perm.size() ] );
        }
    }

    N = V.size();
    assert( GraphUtils::isSimple(V) );
    VPII arcs = GraphUtils::getGraphEdges(V,true);
    if (arcs.size() > deg * N ) {
        StandardUtils::shuffle(arcs);
        arcs.resize(deg * N);
        sort(ALL(arcs));
    }


    // DEBUG(N); DEBUG(deg); DEBUG(arcs.size());
    assert(arcs.size() == N*deg);

    int M = arcs.size();
    out_in << N << " " << M << "\n";
    for (auto [a,b] : arcs) out_in << a << " " << b << "\n";
    out_in << flush;
}

void InstancesGenerator::createExemplarySolution(ifstream &in, ofstream &out) {
}

int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);


    // InstancesGenerator ig("dfvs-instances-small-334", total);
    // InstancesGenerator ig("dfvs-instances-representatives-small", total);
    InstancesGenerator ig("dfvs-instances-representatives-minimal", total);
    ig.threads = 1;
    ig.generate();


    return 0;
}


