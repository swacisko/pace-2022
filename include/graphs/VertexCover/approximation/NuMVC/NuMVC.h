//
// Created by sylwester on 11/29/19.
//

#ifndef ALGORITHMSPROJECT_NUMVC_H
#define ALGORITHMSPROJECT_NUMVC_H


/************************************************
** This is a local search solver for Minimum Vertex Cover.
************************************************/


/************************************************
** Date:	2011.7.1
** TSEWF (Two Stage Exchange and Weighting with Forgetting)
** Author: Shaowei Cai, shaoweicai.cs@gmail.com
**		   School of EECS, Peking University
**		   Beijing, China
**
** Date:    	2011.10.28
** Modify: Shaowei Cai
** use dynamic memory for v_adj[][] and v_edge[][], tidy codes.
************************************************/

/************************************************
** NuMVC Version 2011.11.7
************************************************/


/************************************************
** NuMVC Version 2015.8.20
** This version uses dynamic memory allocation for arrays, and
** implements the initialization function via heap (an advanced data
** structure). These two improvements make NuMVC can handle massive
** graphs well and much more efficient.
************************************************/



//#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <string.h>
//
//
#include <sys/times.h>
#include <cmath>
//#include <graphs/GraphUtils.h>
#include <graphs/GraphReader.h>

//using namespace std;


#define pop_stack_numvc(stack) stack[--stack ## _fill_pointer]
#define push_stack_numvc(item, stack) stack[stack ## _fill_pointer++] = item

/**
 * NuMVC solver for minimum vertex cover.
 *
 * CAUTION!
 * Use a new instance of NuMVC for each call of solve function!
 */
class NuMVC{

//public:

struct Edge {
    int v1;
    int v2;
};

tms start, finish;
double start_time;


/*parameters of algorithm*/
long long max_steps;            //step limit
double cutoff_time;            //time limit
long long step;
int optimal_size;            //terminate the algorithm before step limit if it finds a vertex cover of optimal_size

/*parameters of the instance*/
int v_num;//|V|: 1...v
int e_num;//|E|: 0...e-1

/*structures about edge*/
Edge *edge;
int *edge_weight;

/*structures about vertex*/
int *dscore;                        //dscore of v
long long *time_stamp;
int best_cov_v;        //the vertex of the highest dscore in C

//from vertex to it's edges and neighbors
int **v_edges;            //edges related to v, v_edges[i][k] means vertex v_i's k_th edge
int **v_adj;                //v_adj[v_i][k] = v_j(actually, that is v_i's k_th neighbor)
int *v_edge_count;        //amount of edges (neighbors) related to v


/* structures about solution */
//current candidate solution
int c_size;                        //cardinality of C
int *v_in_c;                        //a flag indicates whether a vertex is in C
int *remove_cand;                //remove candidates, an array consists of only vertices in C, not including tabu_remove
int *index_in_remove_cand;
int remove_cand_size;

//best solution found
int best_c_size;
int *best_v_in_c;                //a flag indicates whether a vertex is in best solution
double best_comp_time;
long best_step;


//uncovered edge stack
int *uncov_stack;                //store the uncov edge number
int uncov_stack_fill_pointer;
int *index_in_uncov_stack;        //which position is an edge in the uncov_stack


//CC and taboo
int *conf_change;
int tabu_remove = 0;

//smooth
int ave_weight = 1;
int delta_total_weight = 0;
int threshold;
float p_scale = 0.3;//w=w*p




int *pos_in_my_heap;

int *my_heap;
int my_heap_count;

void my_heap_swap(int a, int b) {
    int t = my_heap[a];

    my_heap[a] = my_heap[b];
    pos_in_my_heap[my_heap[a]] = a;

    my_heap[b] = t;
    pos_in_my_heap[my_heap[b]] = b;
}

bool my_heap_is_leaf(int pos) {
    if ((pos >= my_heap_count / 2) && (pos < my_heap_count)) return true;
    else return false;
}

int my_heap_left_child(int pos) {
    return (2 * pos + 1);
}

int my_heap_right_child(int pos) {
    return (2 * pos + 2);
}

int my_heap_parent(int pos) {
    return (pos - 1) / 2;
}

void my_heap_shiftdown(int pos) {
    while (!my_heap_is_leaf(pos)) {
        int j = my_heap_left_child(pos);
        int rc = my_heap_right_child(pos);
        if ((rc < my_heap_count) && (dscore[my_heap[rc]] > dscore[my_heap[j]]))
            j = rc;
        if (dscore[my_heap[pos]] > dscore[my_heap[j]]) return;
        my_heap_swap(pos, j);
        pos = j;
    }
}

void my_heap_insert(int v) {
    int curr = my_heap_count++;
    my_heap[curr] = v;
    pos_in_my_heap[v] = curr;

    while (curr != 0 && dscore[my_heap[curr]] > dscore[my_heap[my_heap_parent(curr)]]) {
        my_heap_swap(curr, my_heap_parent(curr));
        curr = my_heap_parent(curr);
    }
}

int my_heap_remove_first() {
    my_heap_swap(0, --my_heap_count);
    if (my_heap_count != 0) my_heap_shiftdown(0);
    return my_heap[my_heap_count];
}

int my_heap_remove(int pos) {
    if (pos == (my_heap_count - 1)) my_heap_count--;
    else {
        my_heap_swap(pos, --my_heap_count);
        while ((pos != 0) && (dscore[my_heap[pos]] > dscore[my_heap[my_heap_parent(pos)]])) {
            my_heap_swap(pos, my_heap_parent(pos));
            pos = my_heap_parent(pos);
        }
        if (my_heap_count != 0) my_heap_shiftdown(pos);
    }
    return my_heap[my_heap_count];
}



/* functions declaration */
    int build_instance(char *filename);
    void cover_rest_edges();


    void update_best_sol()
    {
        int i;

        for (i=1;i<=v_num;i++)
        {
            best_v_in_c[i] = v_in_c[i];
        }

        best_c_size = c_size;
        times(&finish);
        best_comp_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        best_comp_time = round(best_comp_time * 100)/100.0;
        best_step = step;

    }




    int build_instance(/*char *filename*/ istream & infile)
    {
        char line[1024];
        char tempstr1[10];
        char tempstr2[10];
        int  v,e;

        char	tmp;
        int		v1,v2;

//	stringstream infile(filename);
//    if(infile==NULL) return 0;

        /*** build problem data structures of the instance ***/
        infile.getline(line,1024);
        while (line[0] != 'p')
            infile.getline(line,1024);
        sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &v_num, &e_num);

        edge = new Edge [e_num];						//be initialized here
        edge_weight = new int [e_num];					//be initialized in init_sol()
        uncov_stack = new int [e_num];				//only need to initialized uncov_stack_fill_pointer, has been done in init_sol()
        index_in_uncov_stack = new int [e_num];     //the same as above
        dscore = new int [v_num + 1];					//be initialized in init_sol()
        time_stamp = new long long [v_num + 1];			//be initialized in init_sol()
        v_edges = new int* [v_num + 1];					//be initialized here
        v_adj = new int* [v_num + 1];                   //the same as above
        v_edge_count = new int [v_num + 1];				//be initialized here
        v_in_c = new int [v_num + 1];					//be initialized in init_sol()
        remove_cand = new int [v_num + 1];              //be initialized in reset_remove_cand() in init_sol()
        index_in_remove_cand = new int [v_num + 1];     //the same as above
        best_v_in_c = new int [v_num + 1];				//be initialized in update_best_sol() in init_sol()
        conf_change = new int [v_num + 1];				//be initializede int init_sol()

        my_heap = new int [v_num + 1];
        pos_in_my_heap = new int [v_num + 1];
        my_heap_count = 0;

        /* read edges and compute v_edge_count */
        for (v=1; v<=v_num; v++)
            v_edge_count[v] = 0;

        for (e=0; e<e_num; e++)
        {
            infile>>tmp>>v1>>v2;
            v_edge_count[v1]++;
            v_edge_count[v2]++;

            edge[e].v1 = v1;
            edge[e].v2 = v2;
        }

//	infile.close();

        /* build v_adj and v_edges arrays */
        for (v=1; v<=v_num; v++)
        {
            v_adj[v] = new int[v_edge_count[v]];
            v_edges[v] = new int[v_edge_count[v]];
        }

        int* v_edge_count_tmp = new int [v_num + 1];
        for(v=1; v<=v_num; v++)
            v_edge_count_tmp[v]=0;
        for (e=0; e<e_num; e++)
        {

            v1=edge[e].v1;
            v2=edge[e].v2;

            v_edges[v1][v_edge_count_tmp[v1]] = e;
            v_edges[v2][v_edge_count_tmp[v2]] = e;

            v_adj[v1][v_edge_count_tmp[v1]] = v2;
            v_adj[v2][v_edge_count_tmp[v2]] = v1;

            v_edge_count_tmp[v1]++;
            v_edge_count_tmp[v2]++;
        }
        delete[] v_edge_count_tmp;

        return 1;

    }


    void free_memory()
    {
        for (int v=1; v<=v_num; v++)
        {
            delete[] v_adj[v];
            delete[] v_edges[v];
        }
        delete[] conf_change;
        delete[] best_v_in_c;
        delete[] index_in_remove_cand;
        delete[] remove_cand;
        delete[] v_in_c;
        delete[] v_edge_count;
        delete[] v_adj;
        delete[] v_edges;
        delete[] time_stamp;
        delete[] dscore;
        delete[] index_in_uncov_stack;
        delete[] uncov_stack;
        delete[] edge_weight;
        delete[] edge;

        delete[] my_heap;
        delete[] pos_in_my_heap;
    }

    void reset_remove_cand()
    {
        int v,j;
        j=0;
        for (v=1;v<=v_num;v++)
        {
            if(v_in_c[v]==1)// && v!=tabu_remove)
            {
                remove_cand[j] = v;
                index_in_remove_cand[v]=j;
                j++;
            }
            else index_in_remove_cand[v]=0;
        }

        remove_cand_size = j;

    }




    void update_target_size()
    {
        c_size--;

        int max_improvement;
        int max_vertex;//vertex with the highest improvement in C

        max_improvement=-100000000;
        for (int v=1; v<=v_num; ++v)
        {
            if(v_in_c[v]==0)continue;
            if (dscore[v]>max_improvement)
            {
                max_improvement = dscore[v];
                max_vertex = v;
            }
        }
        remove(max_vertex);

        reset_remove_cand();
    }


    void update_best_cov_v()
    {
        int i,v;
        best_cov_v = remove_cand[0];
        for (i=1; i<remove_cand_size; ++i)
        {
            v = remove_cand[i];
            if(v==tabu_remove) continue;

            if( dscore[v] < dscore[best_cov_v])
                continue;
            else if( dscore[v]> dscore[best_cov_v] )
                best_cov_v = v;
            else if (time_stamp[v]<time_stamp[best_cov_v])
                best_cov_v = v;
        }
    }



    inline
    void uncover(int e)
    {
        index_in_uncov_stack[e] = uncov_stack_fill_pointer;
        push_stack_numvc(e,uncov_stack);
    }


    inline
    void cover(int e)
    {
        int index,last_uncov_edge;

        //since the edge is satisfied, its position can be reused to store the last_uncov_edge
        last_uncov_edge = pop_stack_numvc(uncov_stack);
        index = index_in_uncov_stack[e];
        uncov_stack[index] = last_uncov_edge;
        index_in_uncov_stack[last_uncov_edge] = index;
    }

    void init_sol( VI init_vc = VI() )
    {
        int i,v,e;

        /*** build solution data structures of the instance ***/
        //init vertex cover
        for (v=1; v<=v_num; v++)
        {
            v_in_c[v] = 0;
            dscore[v] = 0;

            conf_change[v] = 1;
            time_stamp[v]= 0; // to break ties
        }

        for (e=0; e<e_num; e++)
        {
            edge_weight[e] = 1;
            dscore[edge[e].v1]+=edge_weight[e];
            dscore[edge[e].v2]+=edge_weight[e];
        }

        //init uncovered edge stack and cover_vertrex_count_of_edge array
        uncov_stack_fill_pointer = 0;

        int best_vertex_improvement = 0;
        int best_count = 0;
        int *best_array = new int [v_num + 1];
        for (e=0; e<e_num; e++)
            uncover(e);

        for(v=1; v<=v_num; v++)
        {
            my_heap_insert(v);
        }
        //cout<<my_heap_count<<endl;

        if( init_vc.empty() ){

            for(i=0; uncov_stack_fill_pointer>0; )
            {
                int best_v = my_heap[0];
                //cout<<best_v<<" "<<dscore[best_v]<<""<<my_heap_count<<endl;
                //cout<<uncov_stack_fill_pointer<<endl;
                if(dscore[best_v]>0)
                {
                    add_init(best_v);
                    i++;
                }
            }

        }else{
            for( int d : init_vc ){
                add_init(d+1);
                i++;
            }
        }

        delete[] best_array;

        c_size = i;

        update_best_sol();

        times(&finish);
        int init_sol_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        init_sol_time = round(init_sol_time * 100)/100.0;
//		cout << "c initial solution size = " << c_size << endl;
//		cout << "c initial solution time = " << init_sol_time << endl;

        reset_remove_cand();

        update_best_cov_v();

    }



    void add(int v)
    {
        v_in_c[v] = 1;
        dscore[v] = -dscore[v];

        int i,e,n;

        int edge_count = v_edge_count[v];

        for (i=0; i<edge_count; ++i)
        {
            e = v_edges[v][i];// v's i'th edge
            n = v_adj[v][i];//v's i'th neighbor

            if (v_in_c[n]==0)//this adj isn't in cover set
            {
                dscore[n] -= edge_weight[e];
                conf_change[n] = 1;

                cover(e);
            }
            else
            {
                dscore[n] += edge_weight[e];
            }
        }

    }

    void add_init(int v)
    {
        int pos;
        pos = pos_in_my_heap[v];
        my_heap_remove(pos);

        v_in_c[v] = 1;
        dscore[v] = -dscore[v];

        int i,e,n;

        int edge_count = v_edge_count[v];

        for (i=0; i<edge_count; ++i)
        {
            e = v_edges[v][i];// v's i'th edge
            n = v_adj[v][i];//v's i'th neighbor

            if (v_in_c[n]==0)//this adj isn't in cover set
            {
                dscore[n] -= edge_weight[e];
                conf_change[n] = 1;

                pos = pos_in_my_heap[n];
                my_heap_remove(pos);
                my_heap_insert(n);

                cover(e);
            }
            else
            {
                dscore[n] += edge_weight[e];
            }
        }
    }

    void remove(int v)
    {
        v_in_c[v] = 0;
        dscore[v] = -dscore[v];
        conf_change[v] = 0;

        int i,e,n;

        int edge_count = v_edge_count[v];
        for (i=0; i<edge_count; ++i)
        {
            e = v_edges[v][i];
            n = v_adj[v][i];

            if (v_in_c[n]==0)//this adj isn't in cover set
            {
                dscore[n] += edge_weight[e];
                conf_change[n] = 1;

                uncover(e);
            }
            else
            {
                dscore[n] -=edge_weight[e];
            }
        }

    }


    void print_solution()
    {
        int mis_vertex_count=0;

        for (int i=1; i<=v_num; i++)
        {
            if (best_v_in_c[i]!=1)
                mis_vertex_count++;
        }

        if(mis_vertex_count+best_c_size!=v_num)
            cout<<"The size of independent set + the size of vertex cover is not equal to |V(G)|!"<<endl;

        cout<<"c Best found independent set size = "<<mis_vertex_count<<endl;
        cout<<"c The following output is the found independent set."<<endl;


        for (int i=1; i<=v_num; i++)
        {
            if (best_v_in_c[i]!=1)//output max independent set
                cout<<i<< " ";
        }
        cout<<endl;

    }

    VI get_solution_mis()
    {
        int mis_vertex_count=0;

        for (int i=1; i<=v_num; i++)
        {
            if (best_v_in_c[i]!=1)
                mis_vertex_count++;
        }


        VI res;
        for (int i=1; i<=v_num; i++)
        {
            if (best_v_in_c[i]!=1)//output max independent set
                res.push_back(i);
        }
//	cout<<endl;
        return res;
    }

//check whether the solution found is a proper solution
    int check_solution()
    {
        int e;

        for(e=0; e<e_num; ++e)
        {
            if(best_v_in_c[edge[e].v1]!=1 && best_v_in_c[edge[e].v2]!=1)
            {
                cout<<"uncovered edge "<<e<<endl;
                return 0;
            }
        }

        return 1;
    }




//******************************************************
//******************************************************
//******************************************************//******************************************************//******************************************************
//******************************************************





    int try_step = 10;

    void forget_edge_weights() {
        int v, e;
        int new_total_weight = 0;

        for (v = 1; v <= v_num; v++)
            dscore[v] = 0;

        //scale_ave=ave_weight*q_scale;
        for (e = 0; e < e_num; e++) {
            edge_weight[e] = edge_weight[e] * p_scale;

            new_total_weight += edge_weight[e];

            //update dscore
            if (v_in_c[edge[e].v1] + v_in_c[edge[e].v2] == 0) {
                dscore[edge[e].v1] += edge_weight[e];
                dscore[edge[e].v2] += edge_weight[e];
            } else if (v_in_c[edge[e].v1] + v_in_c[edge[e].v2] == 1) {
                if (v_in_c[edge[e].v1] == 1)dscore[edge[e].v1] -= edge_weight[e];
                else dscore[edge[e].v2] -= edge_weight[e];
            }
        }
        ave_weight = new_total_weight / e_num;

    }


    void update_edge_weight() {
        int i, e;
        for (i = 0; i < uncov_stack_fill_pointer; ++i) {
            e = uncov_stack[i];

            edge_weight[e] += 1;
            dscore[edge[e].v1] += 1;
            dscore[edge[e].v2] += 1;
        }


        delta_total_weight += uncov_stack_fill_pointer;
        if (delta_total_weight >= e_num) {
            ave_weight += 1;
            delta_total_weight -= e_num;
        }

        //smooth weights
        if (ave_weight >= threshold) {
            forget_edge_weights();
        }

    }


    void cover_LS(int seed) {
//        DEBUG("cover_LS");
        int best_add_v;
        int e, v1, v2;

        step = 1;

        UniformIntGenerator rnd(0, 1e9, seed);

        while (1)// wihin cutoff_time
            //while(step<=max_steps)
        {
            if (uncov_stack_fill_pointer == 0)//update best solution if needed
            {
                update_best_sol();

                if (c_size == optimal_size)
                    return;

                update_target_size();

                continue;
            }

            //if(step>=try_step)// wihin cutoff_time
            if (step % try_step == 0) {
                times(&finish);
                double elap_time = ((double)finish.tms_utime + (double)finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);


                if (elap_time >= cutoff_time)return;
            }


            update_best_cov_v();

            remove(best_cov_v);

            e = uncov_stack[rnd.rand() % uncov_stack_fill_pointer];
            v1 = edge[e].v1;
            v2 = edge[e].v2;

            if (conf_change[v1] == 0) best_add_v = v2;
            else if (conf_change[v2] == 0) best_add_v = v1;

            else {
                if (dscore[v1] > dscore[v2] || (dscore[v1] == dscore[v2] && time_stamp[v1] < time_stamp[v2]))
                    best_add_v = v1;
                else best_add_v = v2;
            }

            add(best_add_v);

            int index = index_in_remove_cand[best_cov_v];
            index_in_remove_cand[best_cov_v] = 0;

            remove_cand[index] = best_add_v;
            index_in_remove_cand[best_add_v] = index;

            time_stamp[best_add_v] = time_stamp[best_cov_v] = step;

            tabu_remove = best_add_v;

            update_edge_weight();

            step++;

        }

    }


public:
    /**
     *
     * @param V
     * @param cutoff max time in seconds
     * @param opt_size if we know the optimal size, than we can pass it here
     * @param seeed
     * @return independent set of V
     */
    VI solve( VVI & V, double cutoff, VI init_vc = VI(), int opt_size = 1, int seed = 171 ) {

        //float scale_threshold;

        stringstream infile;
        int edges = 0;
        for( auto & v : V ) edges += v.size();
        edges >>= 1;
        infile << "p td " << V.size() << " " << edges<< endl;
        for( int i=0; i<V.size(); i++ ){
            for( int d : V[i] ){
                if( d > i ) infile << "e " << i+1 << " " << d+1 << endl;
            }
        }

        if (build_instance(infile/*argv[1]*/) != 1) {
            cout << "can't open instance file in NuMVC" << endl;
            cerr << "can't open instance file in NuMVC" << endl;
            exit(1);
        }

        optimal_size = opt_size;
        cutoff_time = cutoff;

        threshold = (int) (0.5 * v_num);

        times(&start);
        start_time = start.tms_utime + start.tms_stime;

        init_sol(init_vc);

//        clog << "NuMVC initialized" << endl;

        if (c_size + uncov_stack_fill_pointer > optimal_size) {
//            cerr<<"c Start local search..."<<endl;
            cover_LS(seed);
        }

        VI res;
        //check solution
        if (check_solution() == 1) {
//			cout << "c Best found vertex cover size = " << best_c_size << endl;
//			cout << best_c_size << endl;
//			print_solution();
//			cout << "c searchSteps = " << best_step << endl;
//			cout << "c solveTime = " << best_comp_time << endl;
//
//            cout<<best_c_size<<' '<<best_comp_time<<' '<<best_step<<endl;

            res = get_solution_mis();

        } else {
            cout << "the solution is wrong." << endl;
            exit(1);
        }

        free_memory();

        return res;
    }


    static void test(){

        VVI V = GraphReader::readGraphDIMACSWunweighed(cin);
//        DEBUG(V.size());
//        DEBUG(V);

        NuMVC numvc;

        int seconds = 60;
        auto mis = numvc.solve(V,seconds);
        set<int> zb(ALL(mis));
        cout << "s vc " << V.size() << " " << V.size() - mis.size() << endl;
        for( int i=0; i<V.size(); i++ ){
            if( zb.count(i) == 0 ) cout << i+1 << endl;
        }


//        exit(1);

    }


};



#endif //ALGORITHMSPROJECT_NUMVC_H
