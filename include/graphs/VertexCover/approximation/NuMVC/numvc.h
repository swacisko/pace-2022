//
//
//#ifndef ALGORITHMSPROJECT_NuMVC_H
//#define ALGORITHMSPROJECT_NuMVC_H
//
//#include "graphs/VertexCover/approximation/NuMVC/my_heap.h"
//#include <graphs/GraphUtils.h>
//
//#include "Makros.h"
//
////
//namespace NuMVC {
//
///* functions declaration */
//int build_instance(char *filename);
//void init_sol( VI initMis  );
//void cover_LS();
//void add(int v);
//void add_init(int v);
//void remove(int v);
//void update_edge_weight();
//void cover_rest_edges();
//int check_solution();
////
////
////
////void update_best_sol()
////{
////	int i;
////
////	for (i=1;i<=v_num;i++)
////	{
////		best_v_in_c[i] = v_in_c[i];
////	}
////
////	best_c_size = c_size;
////	times(&finish);
////	best_comp_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
////	best_comp_time = round(best_comp_time * 100)/100.0;
////	best_step = step;
////
////}
////
////
////
////
////int build_instance(/*char *filename*/ istream & infile)
////{
////	char line[1024];
////	char tempstr1[10];
////	char tempstr2[10];
////	int  v,e;
////
////	char	tmp;
////	int		v1,v2;
////
//////	stringstream infile(filename);
//////    if(infile==NULL) return 0;
////
////	/*** build problem data structures of the instance ***/
////	infile.getline(line,1024);
////	while (line[0] != 'p')
////		infile.getline(line,1024);
////	sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &v_num, &e_num);
////
////	edge = new Edge [e_num];						//be initialized here
////	edge_weight = new int [e_num];					//be initialized in init_sol()
////	uncov_stack = new int [e_num];				//only need to initialized uncov_stack_fill_pointer, has been done in init_sol()
////	index_in_uncov_stack = new int [e_num];     //the same as above
////	dscore = new int [v_num + 1];					//be initialized in init_sol()
////	time_stamp = new long long [v_num + 1];			//be initialized in init_sol()
////	v_edges = new int* [v_num + 1];					//be initialized here
////	v_adj = new int* [v_num + 1];                   //the same as above
////	v_edge_count = new int [v_num + 1];				//be initialized here
////	v_in_c = new int [v_num + 1];					//be initialized in init_sol()
////	remove_cand = new int [v_num + 1];              //be initialized in reset_remove_cand() in init_sol()
////	index_in_remove_cand = new int [v_num + 1];     //the same as above
////	best_v_in_c = new int [v_num + 1];				//be initialized in update_best_sol() in init_sol()
////	conf_change = new int [v_num + 1];				//be initializede int init_sol()
////
////	my_heap = new int [v_num + 1];
////	pos_in_my_heap = new int [v_num + 1];
////	my_heap_count = 0;
////
////	/* read edges and compute v_edge_count */
////	for (v=1; v<=v_num; v++)
////		v_edge_count[v] = 0;
////
////	for (e=0; e<e_num; e++)
////	{
////		infile>>tmp>>v1>>v2;
////		v_edge_count[v1]++;
////		v_edge_count[v2]++;
////
////		edge[e].v1 = v1;
////		edge[e].v2 = v2;
////	}
////
//////	infile.close();
////
////	/* build v_adj and v_edges arrays */
////	for (v=1; v<=v_num; v++)
////	{
////		v_adj[v] = new int[v_edge_count[v]];
////		v_edges[v] = new int[v_edge_count[v]];
////	}
////
////	int* v_edge_count_tmp = new int [v_num + 1];
////	for(v=1; v<=v_num; v++)
////		v_edge_count_tmp[v]=0;
////	for (e=0; e<e_num; e++)
////	{
////
////		v1=edge[e].v1;
////		v2=edge[e].v2;
////
////		v_edges[v1][v_edge_count_tmp[v1]] = e;
////		v_edges[v2][v_edge_count_tmp[v2]] = e;
////
////		v_adj[v1][v_edge_count_tmp[v1]] = v2;
////		v_adj[v2][v_edge_count_tmp[v2]] = v1;
////
////		v_edge_count_tmp[v1]++;
////		v_edge_count_tmp[v2]++;
////	}
////	delete[] v_edge_count_tmp;
////
////	return 1;
////
////}
////
////
////void free_memory()
////{
////	for (int v=1; v<=v_num; v++)
////	{
////		delete[] v_adj[v];
////		delete[] v_edges[v];
////	}
////	delete[] conf_change;
////	delete[] best_v_in_c;
////	delete[] index_in_remove_cand;
////	delete[] remove_cand;
////	delete[] v_in_c;
////	delete[] v_edge_count;
////	delete[] v_adj;
////	delete[] v_edges;
////	delete[] time_stamp;
////	delete[] dscore;
////	delete[] index_in_uncov_stack;
////	delete[] uncov_stack;
////	delete[] edge_weight;
////	delete[] edge;
////
////	delete[] my_heap;
////	delete[] pos_in_my_heap;
////}
////
////void reset_remove_cand()
////{
////	int v,j;
////	j=0;
////	for (v=1;v<=v_num;v++)
////	{
////		if(v_in_c[v]==1)// && v!=tabu_remove)
////		{
////			remove_cand[j] = v;
////			index_in_remove_cand[v]=j;
////			j++;
////		}
////		else index_in_remove_cand[v]=0;
////	}
////
////	remove_cand_size = j;
////
////}
////
////
////
////
////void update_target_size()
////{
////	c_size--;
////
////	int max_improvement;
////	int max_vertex;//vertex with the highest improvement in C
////
////	max_improvement=-100000000;
////	for (int v=1; v<=v_num; ++v)
////	{
////		if(v_in_c[v]==0)continue;
////		if (dscore[v]>max_improvement)
////		{
////			max_improvement = dscore[v];
////			max_vertex = v;
////		}
////	}
////	remove(max_vertex);
////
////	reset_remove_cand();
////}
////
////
////
////
//////update the best vertex in C
////
////void update_best_cov_v()
////{
////	int i,v;
////	best_cov_v = remove_cand[0];
////	for (i=1; i<remove_cand_size; ++i)
////	{
////		v = remove_cand[i];
////		if(v==tabu_remove) continue;
////
////		if( dscore[v] < dscore[best_cov_v])
////			continue;
////		else if( dscore[v]> dscore[best_cov_v] )
////			best_cov_v = v;
////		else if (time_stamp[v]<time_stamp[best_cov_v])
////			best_cov_v = v;
////	}
////}
////
////
////
////inline
////void uncover(int e)
////{
////	index_in_uncov_stack[e] = uncov_stack_fill_pointer;
////	push(e,uncov_stack);
////}
////
////
////inline
////void cover(int e)
////{
////	int index,last_uncov_edge;
////
////	//since the edge is satisfied, its position can be reused to store the last_uncov_edge
////	last_uncov_edge = pop(uncov_stack);
////	index = index_in_uncov_stack[e];
////	uncov_stack[index] = last_uncov_edge;
////	index_in_uncov_stack[last_uncov_edge] = index;
////}
////
////
////
/////*
////void init_sol()
////{
////	int i,v,e;
////
////	// build solution data structures of the instance
////	//init vertex cover
////	for (v=1; v<=v_num; v++)
////	{
////		v_in_c[v] = 0;
////		dscore[v] = 0;
////
////		conf_change[v] = 1;
////		time_stamp[v]= 0; // to break ties
////	}
////
////	for (e=0; e<e_num; e++)
////	{
////		edge_weight[e] = 1;
////		dscore[edge[e].v1]+=edge_weight[e];
////		dscore[edge[e].v2]+=edge_weight[e];
////	}
////
////	//init uncovered edge stack and cover_vertrex_count_of_edge array
////	uncov_stack_fill_pointer = 0;
////
////	int best_vertex_improvement = 0;
////	int best_count = 0;
////	int *best_array = new int [v_num + 1];
////	for (e=0; e<e_num; e++)
////		uncover(e);
////
////
////
////	for (i=0; uncov_stack_fill_pointer>0; )
////	{
////		best_vertex_improvement = 0;
////		best_count = 0;
////		for (v=1; v<=v_num; ++v)
////		{
////			if(v_in_c[v]==1)continue;
////
////			if (dscore[v]>best_vertex_improvement)
////			{
////				best_vertex_improvement = dscore[v];
////				best_array[0] = v;
////				best_count = 1;
////			}
////			else if (dscore[v]==best_vertex_improvement)
////			{
////				best_array[best_count] = v;
////				best_count++;
////			}
////		}
////
////		if(best_count>0)
////		{
////			add(best_array[rand()%best_count]);
////			++i;
////		}
////	}
////	delete[] best_array;
////
////	c_size = i;
////
////	update_best_sol();
////
////	times(&finish);
////	int init_sol_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
////	init_sol_time = round(init_sol_time * 100)/100.0;
////	cout << "c initial solution size = " << c_size << endl;
////	cout << "c initial solution time = " << init_sol_time << endl;
////
////	reset_remove_cand();
////
////	update_best_cov_v();
////
////}
////*/
////
////void init_sol()
////{
////	int i,v,e;
////
////	/*** build solution data structures of the instance ***/
////	//init vertex cover
////	for (v=1; v<=v_num; v++)
////	{
////		v_in_c[v] = 0;
////		dscore[v] = 0;
////
////		conf_change[v] = 1;
////		time_stamp[v]= 0; // to break ties
////	}
////
////	for (e=0; e<e_num; e++)
////	{
////		edge_weight[e] = 1;
////		dscore[edge[e].v1]+=edge_weight[e];
////		dscore[edge[e].v2]+=edge_weight[e];
////	}
////
////	//init uncovered edge stack and cover_vertrex_count_of_edge array
////	uncov_stack_fill_pointer = 0;
////
////	int best_vertex_improvement = 0;
////	int best_count = 0;
////	int *best_array = new int [v_num + 1];
////	for (e=0; e<e_num; e++)
////		uncover(e);
////
////	for(v=1; v<=v_num; v++)
////	{
////		my_heap_insert(v);
////	}
////	//cout<<my_heap_count<<endl;
////
////	for(i=0; uncov_stack_fill_pointer>0; )
////	{
////		int best_v = my_heap[0];
////		//cout<<best_v<<" "<<dscore[best_v]<<""<<my_heap_count<<endl;
////		//cout<<uncov_stack_fill_pointer<<endl;
////		if(dscore[best_v]>0)
////		{
////			add_init(best_v);
////			i++;
////		}
////	}
////
////	/*
////	for (i=0; uncov_stack_fill_pointer>0; )
////	{
////		best_vertex_improvement = 0;
////		best_count = 0;
////		for (v=1; v<=v_num; ++v)
////		{
////			if(v_in_c[v]==1)continue;
////
////			if (dscore[v]>best_vertex_improvement)
////			{
////				best_vertex_improvement = dscore[v];
////				best_array[0] = v;
////				best_count = 1;
////			}
////			else if (dscore[v]==best_vertex_improvement)
////			{
////				best_array[best_count] = v;
////				best_count++;
////			}
////		}
////
////		if(best_count>0)
////		{
////			add(best_array[rand()%best_count]);
////			++i;
////		}
////	}
////	*/
////	delete[] best_array;
////
////	c_size = i;
////
////	update_best_sol();
////
////	times(&finish);
////	int init_sol_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
////	init_sol_time = round(init_sol_time * 100)/100.0;
////	cout << "c initial solution size = " << c_size << endl;
////	cout << "c initial solution time = " << init_sol_time << endl;
////
////	reset_remove_cand();
////
////	update_best_cov_v();
////
////}
////
////
////
////void add(int v)
////{
////	v_in_c[v] = 1;
////	dscore[v] = -dscore[v];
////
////	int i,e,n;
////
////	int edge_count = v_edge_count[v];
////
////	for (i=0; i<edge_count; ++i)
////	{
////		e = v_edges[v][i];// v's i'th edge
////		n = v_adj[v][i];//v's i'th neighbor
////
////		if (v_in_c[n]==0)//this adj isn't in cover set
////		{
////			dscore[n] -= edge_weight[e];
////			conf_change[n] = 1;
////
////			cover(e);
////		}
////		else
////		{
////			dscore[n] += edge_weight[e];
////		}
////	}
////
////}
////
////void add_init(int v)
////{
////	int pos;
////	pos = pos_in_my_heap[v];
////	my_heap_remove(pos);
////
////	v_in_c[v] = 1;
////	dscore[v] = -dscore[v];
////
////	int i,e,n;
////
////	int edge_count = v_edge_count[v];
////
////	for (i=0; i<edge_count; ++i)
////	{
////		e = v_edges[v][i];// v's i'th edge
////		n = v_adj[v][i];//v's i'th neighbor
////
////		if (v_in_c[n]==0)//this adj isn't in cover set
////		{
////			dscore[n] -= edge_weight[e];
////			conf_change[n] = 1;
////
////			pos = pos_in_my_heap[n];
////			my_heap_remove(pos);
////			my_heap_insert(n);
////
////			cover(e);
////		}
////		else
////		{
////			dscore[n] += edge_weight[e];
////		}
////	}
////}
////
////void remove(int v)
////{
////	v_in_c[v] = 0;
////	dscore[v] = -dscore[v];
////	conf_change[v] = 0;
////
////	int i,e,n;
////
////	int edge_count = v_edge_count[v];
////	for (i=0; i<edge_count; ++i)
////	{
////		e = v_edges[v][i];
////		n = v_adj[v][i];
////
////		if (v_in_c[n]==0)//this adj isn't in cover set
////		{
////			dscore[n] += edge_weight[e];
////			conf_change[n] = 1;
////
////			uncover(e);
////		}
////		else
////		{
////			dscore[n] -= edge_weight[e];
////		}
////	}
////
////}
////
/////*On solution*/
////
////void print_solution()
////{
////	int mis_vertex_count=0;
////
////	for (int i=1; i<=v_num; i++)
////	{
////		if (best_v_in_c[i]!=1)
////			mis_vertex_count++;
////	}
////
////	if(mis_vertex_count+best_c_size!=v_num)
////		cout<<"The size of independent set + the size of vertex cover is not equal to |V(G)|!"<<endl;
////
////	cout<<"c Best found independent set size = "<<mis_vertex_count<<endl;
////	cout<<"c The following output is the found independent set."<<endl;
////
////
////	for (int i=1; i<=v_num; i++)
////	{
////		if (best_v_in_c[i]!=1)//output max independent set
////			cout<<i<< " ";
////	}
////	cout<<endl;
////
////}
////
////VI get_solution_mis()
////{
////	int mis_vertex_count=0;
////
////	for (int i=1; i<=v_num; i++)
////	{
////		if (best_v_in_c[i]!=1)
////			mis_vertex_count++;
////	}
////
//////	if(mis_vertex_count+best_c_size!=v_num)
//////		cout<<"The size of independent set + the size of vertex cover is not equal to |V(G)|!"<<endl;
////
//////	cout<<"c Best found independent set size = "<<mis_vertex_count<<endl;
//////	cout<<"c The following output is the found independent set."<<endl;
////
////
////	VI res;
////	for (int i=1; i<=v_num; i++)
////	{
////		if (best_v_in_c[i]!=1)//output max independent set
////			res.push_back(i);
////	}
//////	cout<<endl;
////	return res;
////}
////
//////check whether the solution found is a proper solution
////int check_solution()
////{
////	int e;
////
////	for(e=0; e<e_num; ++e)
////	{
////		if(best_v_in_c[edge[e].v1]!=1 && best_v_in_c[edge[e].v2]!=1)
////		{
////			cout<<"uncovered edge "<<e<<endl;
////			return 0;
////		}
////	}
////
////	return 1;
////}
////
////
////
////
//////******************************************************
//////******************************************************
//////******************************************************//******************************************************//******************************************************
//////******************************************************
////
////
////
////
////
////	int try_step = 100;
////
////	void forget_edge_weights() {
////		int v, e;
////		int new_total_weight = 0;
////
////		for (v = 1; v <= v_num; v++)
////			dscore[v] = 0;
////
////		//scale_ave=ave_weight*q_scale;
////		for (e = 0; e < e_num; e++) {
////			edge_weight[e] = edge_weight[e] * p_scale;
////
////			new_total_weight += edge_weight[e];
////
////			//update dscore
////			if (v_in_c[edge[e].v1] + v_in_c[edge[e].v2] == 0) {
////				dscore[edge[e].v1] += edge_weight[e];
////				dscore[edge[e].v2] += edge_weight[e];
////			} else if (v_in_c[edge[e].v1] + v_in_c[edge[e].v2] == 1) {
////				if (v_in_c[edge[e].v1] == 1)dscore[edge[e].v1] -= edge_weight[e];
////				else dscore[edge[e].v2] -= edge_weight[e];
////			}
////		}
////		ave_weight = new_total_weight / e_num;
////
////	}
////
////
////	void update_edge_weight() {
////		int i, e;
////		for (i = 0; i < uncov_stack_fill_pointer; ++i) {
////			e = uncov_stack[i];
////
////			edge_weight[e] += 1;
////			dscore[edge[e].v1] += 1;
////			dscore[edge[e].v2] += 1;
////		}
////
////
////		delta_total_weight += uncov_stack_fill_pointer;
////		if (delta_total_weight >= e_num) {
////			ave_weight += 1;
////			delta_total_weight -= e_num;
////		}
////
////		//smooth weights
////		if (ave_weight >= threshold) {
////			forget_edge_weights();
////		}
////
////	}
////
////
////	void cover_LS() {
////		int best_add_v;
////		int e, v1, v2;
////
////		step = 1;
////
////		while (1)// wihin cutoff_time
////			//while(step<=max_steps)
////		{
////			if (uncov_stack_fill_pointer == 0)//update best solution if needed
////			{
////				update_best_sol();
////
////				if (c_size == optimal_size)
////					return;
////
////				update_target_size();
////
////				continue;
////			}
////
////			//if(step>=try_step)// wihin cutoff_time
////			if (step % try_step == 0) {
////				times(&finish);
////				double elap_time = (finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
////				if (elap_time >= cutoff_time)return;
////			}
////
////
////			update_best_cov_v();
////
////			remove(best_cov_v);
////
////			e = uncov_stack[rand() % uncov_stack_fill_pointer];
////			v1 = edge[e].v1;
////			v2 = edge[e].v2;
////
////			if (conf_change[v1] == 0) best_add_v = v2;
////			else if (conf_change[v2] == 0) best_add_v = v1;
////
////			else {
////				if (dscore[v1] > dscore[v2] || (dscore[v1] == dscore[v2] && time_stamp[v1] < time_stamp[v2]))
////					best_add_v = v1;
////				else best_add_v = v2;
////			}
////
////			add(best_add_v);
////
////			int index = index_in_remove_cand[best_cov_v];
////			index_in_remove_cand[best_cov_v] = 0;
////
////			remove_cand[index] = best_add_v;
////			index_in_remove_cand[best_add_v] = index;
////
////			time_stamp[best_add_v] = time_stamp[best_cov_v] = step;
////
////			tabu_remove = best_add_v;
////
////			update_edge_weight();
////
////			step++;
////
////		}
////
////	}
////
////	/**
////	 *
////	 * @param V
////	 * @param cutoff max time in seconds
////	 * @param opt_size if we know the optimal size, than we can pass it here
////	 * @param seeed
////	 * @return vertex cover of V
////	 */
////	VI solve(/*int argc, char *argv[]*/ VVI & V, int cutoff, int opt_size = 0, int seeed = 171 ) {
//////		if (argc != 5) {
//////			cout << "usage: " << endl;
//////			return 1;
//////		}
////		int seed;
////		//float scale_threshold;
////
////		stringstream infile;
////		infile << "p td " << V.size() << " " << GraphUtils::countEdges(V) << endl;
////		VPII edges = GraphUtils::getGraphEdges(V);
////		for( auto p : edges ) infile << "e " << p.first+1 << " " << p.second+1 << endl;
////
////		if (build_instance(infile/*argv[1]*/) != 1) {
////			cout << "can't open instance file in NuMVC" << endl;
////			cerr << "can't open instance file in NuMVC" << endl;
////			exit(1);
////			return VI();
////		}
////
//////		sscanf(argv[2], "%d", &optimal_size);//if you want to stop the algorithm only cutoff time is reached, set optimal_size to 0.
//////		sscanf(argv[3], "%d", &seed);
//////		sscanf(argv[4], "%d", &cutoff_time);
////
////		optimal_size = opt_size;
////		seed = seeed;
////		cutoff_time = cutoff;
////
////		threshold = (int) (0.5 * v_num);
////
////		srand(seed);
////
////		//cout<<seed<<' ';
//////		cout<<"c This is NuMVC, a local search solver for the Minimum Vertex Cover (and also Maximum Independent Set) problem."<<endl;
////
////		times(&start);
////		start_time = start.tms_utime + start.tms_stime;
////
////		init_sol();
////
////		if (c_size + uncov_stack_fill_pointer > optimal_size) {
////			//cout<<"c Start local search..."<<endl;
////			cover_LS();
////		}
////
////		VI res;
////		//check solution
////		if (check_solution() == 1) {
//////			cout << "c Best found vertex cover size = " << best_c_size << endl;
//////			cout << best_c_size << endl;
//////			print_solution();
//////			cout << "c searchSteps = " << best_step << endl;
//////			cout << "c solveTime = " << best_comp_time << endl;
////
////			//cout<<best_c_size<<' '<<best_comp_time<<' '<<best_step<<endl;
////
////			res = get_solution_mis();
////
////		} else {
////			cout << "the solution is wrong." << endl;
////			exit(1);
////			//print_solution();
////		}
////
////		free_memory();
////
//////		return 0;
////
////
////		return res;
////	}
////
////
////}
////
//
//
////namespace NuMVC{
//
//
//	void update_best_sol()
//	{
//		int i;
//
//		for (i=1;i<=v_num;i++)
//		{
//			best_v_in_c[i] = v_in_c[i];
//		}
//
//		best_c_size = c_size;
//		times(&finish);
//		best_comp_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
//		best_comp_time = round(best_comp_time * 100)/100.0;
//		best_step = step;
//
//	}
//
//
//
//
//	int build_instance(/*char *filename*/ istream & infile)
//	{
//		char line[1024];
//		char tempstr1[10];
//		char tempstr2[10];
//		int  v,e;
//
//		char	tmp;
//		int		v1,v2;
//
////	stringstream infile(filename);
////    if(infile==NULL) return 0;
//
//		/*** build problem data structures of the instance ***/
//		infile.getline(line,1024);
//		while (line[0] != 'p')
//			infile.getline(line,1024);
//		sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &v_num, &e_num);
//
//		edge = new Edge [e_num];						//be initialized here
//		edge_weight = new int [e_num];					//be initialized in init_sol()
//		uncov_stack = new int [e_num];				//only need to initialized uncov_stack_fill_pointer, has been done in init_sol()
//		index_in_uncov_stack = new int [e_num];     //the same as above
//		dscore = new int [v_num + 1];					//be initialized in init_sol()
//		time_stamp = new long long [v_num + 1];			//be initialized in init_sol()
//		v_edges = new int* [v_num + 1];					//be initialized here
//		v_adj = new int* [v_num + 1];                   //the same as above
//		v_edge_count = new int [v_num + 1];				//be initialized here
//		v_in_c = new int [v_num + 1];					//be initialized in init_sol()
//		remove_cand = new int [v_num + 1];              //be initialized in reset_remove_cand() in init_sol()
//		index_in_remove_cand = new int [v_num + 1];     //the same as above
//		best_v_in_c = new int [v_num + 1];				//be initialized in update_best_sol() in init_sol()
//		conf_change = new int [v_num + 1];				//be initializede int init_sol()
//
//		my_heap = new int [v_num + 1];
//		pos_in_my_heap = new int [v_num + 1];
//		my_heap_count = 0;
//
//		/* read edges and compute v_edge_count */
//		for (v=1; v<=v_num; v++)
//			v_edge_count[v] = 0;
//
//		for (e=0; e<e_num; e++)
//		{
//			infile>>tmp>>v1>>v2;
//			v_edge_count[v1]++;
//			v_edge_count[v2]++;
//
//			edge[e].v1 = v1;
//			edge[e].v2 = v2;
//		}
//
////	infile.close();
//
//		/* build v_adj and v_edges arrays */
//		for (v=1; v<=v_num; v++)
//		{
//			v_adj[v] = new int[v_edge_count[v]];
//			v_edges[v] = new int[v_edge_count[v]];
//		}
//
//		int* v_edge_count_tmp = new int [v_num + 1];
//		for(v=1; v<=v_num; v++)
//			v_edge_count_tmp[v]=0;
//		for (e=0; e<e_num; e++)
//		{
//
//			v1=edge[e].v1;
//			v2=edge[e].v2;
//
//			v_edges[v1][v_edge_count_tmp[v1]] = e;
//			v_edges[v2][v_edge_count_tmp[v2]] = e;
//
//			v_adj[v1][v_edge_count_tmp[v1]] = v2;
//			v_adj[v2][v_edge_count_tmp[v2]] = v1;
//
//			v_edge_count_tmp[v1]++;
//			v_edge_count_tmp[v2]++;
//		}
//		delete[] v_edge_count_tmp;
//
//		return 1;
//
//	}
//
//
//	void free_memory()
//	{
//		for (int v=1; v<=v_num; v++)
//		{
//			delete[] v_adj[v];
//			delete[] v_edges[v];
//		}
//		delete[] conf_change;
//		delete[] best_v_in_c;
//		delete[] index_in_remove_cand;
//		delete[] remove_cand;
//		delete[] v_in_c;
//		delete[] v_edge_count;
//		delete[] v_adj;
//		delete[] v_edges;
//		delete[] time_stamp;
//		delete[] dscore;
//		delete[] index_in_uncov_stack;
//		delete[] uncov_stack;
//		delete[] edge_weight;
//		delete[] edge;
//
//		delete[] my_heap;
//		delete[] pos_in_my_heap;
//	}
//
//	void reset_remove_cand()
//	{
//		int v,j;
//		j=0;
//		for (v=1;v<=v_num;v++)
//		{
//			if(v_in_c[v]==1)// && v!=tabu_remove)
//			{
//				remove_cand[j] = v;
//				index_in_remove_cand[v]=j;
//				j++;
//			}
//			else index_in_remove_cand[v]=0;
//		}
//
//		remove_cand_size = j;
//
//	}
//
//
//
//
//	void update_target_size()
//	{
//		c_size--;
//
//		int max_improvement;
//		int max_vertex;//vertex with the highest improvement in C
//
//		max_improvement=-100000000;
//		for (int v=1; v<=v_num; ++v)
//		{
//			if(v_in_c[v]==0)continue;
//			if (dscore[v]>max_improvement)
//			{
//				max_improvement = dscore[v];
//				max_vertex = v;
//			}
//		}
//		remove(max_vertex);
//
//		reset_remove_cand();
//	}
//
//
//
//
////update the best vertex in C
//
//	void update_best_cov_v()
//	{
//		int i,v;
//		best_cov_v = remove_cand[0];
//		for (i=1; i<remove_cand_size; ++i)
//		{
//			v = remove_cand[i];
//			if(v==tabu_remove) continue;
//
//			if( dscore[v] < dscore[best_cov_v])
//				continue;
//			else if( dscore[v]> dscore[best_cov_v] )
//				best_cov_v = v;
//			else if (time_stamp[v]<time_stamp[best_cov_v])
//				best_cov_v = v;
//		}
//	}
//
//
//
//	inline
//	void uncover(int e)
//	{
//		index_in_uncov_stack[e] = uncov_stack_fill_pointer;
//		push(e,uncov_stack);
//	}
//
//
//	inline
//	void cover(int e)
//	{
//		int index,last_uncov_edge;
//
//		//since the edge is satisfied, its position can be reused to store the last_uncov_edge
//		last_uncov_edge = remove(uncov_stack);
//		index = index_in_uncov_stack[e];
//		uncov_stack[index] = last_uncov_edge;
//		index_in_uncov_stack[last_uncov_edge] = index;
//	}
//
//
//
///*
//void init_sol()
//{
//	int i,v,e;
//
//	// build solution data structures of the instance
//	//init vertex cover
//	for (v=1; v<=v_num; v++)
//	{
//		v_in_c[v] = 0;
//		dscore[v] = 0;
//
//		conf_change[v] = 1;
//		time_stamp[v]= 0; // to break ties
//	}
//
//	for (e=0; e<e_num; e++)
//	{
//		edge_weight[e] = 1;
//		dscore[edge[e].v1]+=edge_weight[e];
//		dscore[edge[e].v2]+=edge_weight[e];
//	}
//
//	//init uncovered edge stack and cover_vertrex_count_of_edge array
//	uncov_stack_fill_pointer = 0;
//
//	int best_vertex_improvement = 0;
//	int best_count = 0;
//	int *best_array = new int [v_num + 1];
//	for (e=0; e<e_num; e++)
//		uncover(e);
//
//
//
//	for (i=0; uncov_stack_fill_pointer>0; )
//	{
//		best_vertex_improvement = 0;
//		best_count = 0;
//		for (v=1; v<=v_num; ++v)
//		{
//			if(v_in_c[v]==1)continue;
//
//			if (dscore[v]>best_vertex_improvement)
//			{
//				best_vertex_improvement = dscore[v];
//				best_array[0] = v;
//				best_count = 1;
//			}
//			else if (dscore[v]==best_vertex_improvement)
//			{
//				best_array[best_count] = v;
//				best_count++;
//			}
//		}
//
//		if(best_count>0)
//		{
//			add(best_array[rand()%best_count]);
//			++i;
//		}
//	}
//	delete[] best_array;
//
//	c_size = i;
//
//	update_best_sol();
//
//	times(&finish);
//	int init_sol_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
//	init_sol_time = round(init_sol_time * 100)/100.0;
//	cout << "c initial solution size = " << c_size << endl;
//	cout << "c initial solution time = " << init_sol_time << endl;
//
//	reset_remove_cand();
//
//	update_best_cov_v();
//
//}
//*/
//
//	void init_sol( VI initMis = VI() )
//	{
//		int i,v,e;
//
//		/*** build solution data structures of the instance ***/
//		//init vertex cover
//		for (v=1; v<=v_num; v++)
//		{
//			v_in_c[v] = 0;
//			dscore[v] = 0;
//
//			conf_change[v] = 1;
//			time_stamp[v]= 0; // to break ties
//		}
//
//		for (e=0; e<e_num; e++)
//		{
//			edge_weight[e] = 1;
//			dscore[edge[e].v1]+=edge_weight[e];
//			dscore[edge[e].v2]+=edge_weight[e];
//		}
//
//		//init uncovered edge stack and cover_vertrex_count_of_edge array
//		uncov_stack_fill_pointer = 0;
//
//		int best_vertex_improvement = 0;
//		int best_count = 0;
//		int *best_array = new int [v_num + 1];
//		for (e=0; e<e_num; e++)
//			uncover(e);
//
//			for(v=1; v<=v_num; v++)
//			{
//				my_heap_insert(v);
//			}
//		//cout<<my_heap_count<<endl;
//
//		if( initMis.empty() ){
//
//			for(i=0; uncov_stack_fill_pointer>0; )
//			{
//				int best_v = my_heap[0];
//				//cout<<best_v<<" "<<dscore[best_v]<<""<<my_heap_count<<endl;
//				//cout<<uncov_stack_fill_pointer<<endl;
//				if(dscore[best_v]>0)
//				{
//					add_init(best_v);
//					i++;
//				}
//			}
//
//		}else{
//			for( int d : initMis ){
//				add_init(d+1);
//				i++;
//			}
//		}
//
//
//		/*
//        for (i=0; uncov_stack_fill_pointer>0; )
//        {
//            best_vertex_improvement = 0;
//            best_count = 0;
//            for (v=1; v<=v_num; ++v)
//            {
//                if(v_in_c[v]==1)continue;
//
//                if (dscore[v]>best_vertex_improvement)
//                {
//                    best_vertex_improvement = dscore[v];
//                    best_array[0] = v;
//                    best_count = 1;
//                }
//                else if (dscore[v]==best_vertex_improvement)
//                {
//                    best_array[best_count] = v;
//                    best_count++;
//                }
//            }
//
//            if(best_count>0)
//            {
//                add(best_array[rand()%best_count]);
//                ++i;
//            }
//        }
//        */
//		delete[] best_array;
//
//		c_size = i;
//
//		update_best_sol();
//
//		times(&finish);
//		int init_sol_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
//		init_sol_time = round(init_sol_time * 100)/100.0;
////		cout << "c initial solution size = " << c_size << endl;
////		cout << "c initial solution time = " << init_sol_time << endl;
//
//		reset_remove_cand();
//
//		update_best_cov_v();
//
//	}
//
//
//
//	void add(int v)
//	{
//		v_in_c[v] = 1;
//		dscore[v] = -dscore[v];
//
//		int i,e,n;
//
//		int edge_count = v_edge_count[v];
//
//		for (i=0; i<edge_count; ++i)
//		{
//			e = v_edges[v][i];// v's i'th edge
//			n = v_adj[v][i];//v's i'th neighbor
//
//			if (v_in_c[n]==0)//this adj isn't in cover set
//			{
//				dscore[n] -= edge_weight[e];
//				conf_change[n] = 1;
//
//				cover(e);
//			}
//			else
//			{
//				dscore[n] += edge_weight[e];
//			}
//		}
//
//	}
//
//	void add_init(int v)
//	{
//		int pos;
//		pos = pos_in_my_heap[v];
//		my_heap_remove(pos);
//
//		v_in_c[v] = 1;
//		dscore[v] = -dscore[v];
//
//		int i,e,n;
//
//		int edge_count = v_edge_count[v];
//
//		for (i=0; i<edge_count; ++i)
//		{
//			e = v_edges[v][i];// v's i'th edge
//			n = v_adj[v][i];//v's i'th neighbor
//
//			if (v_in_c[n]==0)//this adj isn't in cover set
//			{
//				dscore[n] -= edge_weight[e];
//				conf_change[n] = 1;
//
//				pos = pos_in_my_heap[n];
//				my_heap_remove(pos);
//				my_heap_insert(n);
//
//				cover(e);
//			}
//			else
//			{
//				dscore[n] += edge_weight[e];
//			}
//		}
//	}
//
//	void remove(int v)
//	{
//		v_in_c[v] = 0;
//		dscore[v] = -dscore[v];
//		conf_change[v] = 0;
//
//		int i,e,n;
//
//		int edge_count = v_edge_count[v];
//		for (i=0; i<edge_count; ++i)
//		{
//			e = v_edges[v][i];
//			n = v_adj[v][i];
//
//			if (v_in_c[n]==0)//this adj isn't in cover set
//			{
//				dscore[n] += edge_weight[e];
//				conf_change[n] = 1;
//
//				uncover(e);
//			}
//			else
//			{
//				dscore[n] -=edge_weight[e];
//			}
//		}
//
//	}
//
///*On solution*/
//
//	void print_solution()
//	{
//		int mis_vertex_count=0;
//
//		for (int i=1; i<=v_num; i++)
//		{
//			if (best_v_in_c[i]!=1)
//				mis_vertex_count++;
//		}
//
//		if(mis_vertex_count+best_c_size!=v_num)
//			cout<<"The size of independent set + the size of vertex cover is not equal to |V(G)|!"<<endl;
//
//		cout<<"c Best found independent set size = "<<mis_vertex_count<<endl;
//		cout<<"c The following output is the found independent set."<<endl;
//
//
//		for (int i=1; i<=v_num; i++)
//		{
//			if (best_v_in_c[i]!=1)//output max independent set
//				cout<<i<< " ";
//		}
//		cout<<endl;
//
//	}
//
//	VI get_solution_mis()
//	{
//		int mis_vertex_count=0;
//
//		for (int i=1; i<=v_num; i++)
//		{
//			if (best_v_in_c[i]!=1)
//				mis_vertex_count++;
//		}
//
////	if(mis_vertex_count+best_c_size!=v_num)
////		cout<<"The size of independent set + the size of vertex cover is not equal to |V(G)|!"<<endl;
//
////	cout<<"c Best found independent set size = "<<mis_vertex_count<<endl;
////	cout<<"c The following output is the found independent set."<<endl;
//
//
//		VI res;
//		for (int i=1; i<=v_num; i++)
//		{
//			if (best_v_in_c[i]!=1)//output max independent set
//				res.push_back(i);
//		}
////	cout<<endl;
//		return res;
//	}
//
////check whether the solution found is a proper solution
//	int check_solution()
//	{
//		int e;
//
//		for(e=0; e<e_num; ++e)
//		{
//			if(best_v_in_c[edge[e].v1]!=1 && best_v_in_c[edge[e].v2]!=1)
//			{
//				cout<<"uncovered edge "<<e<<endl;
//				return 0;
//			}
//		}
//
//		return 1;
//	}
//
//
//
//
////******************************************************
////******************************************************
////******************************************************//******************************************************//******************************************************
////******************************************************
//
//
//
//
//
//	int try_step = 100;
//
//	void forget_edge_weights() {
//		int v, e;
//		int new_total_weight = 0;
//
//		for (v = 1; v <= v_num; v++)
//			dscore[v] = 0;
//
//		//scale_ave=ave_weight*q_scale;
//		for (e = 0; e < e_num; e++) {
//			edge_weight[e] = edge_weight[e] * p_scale;
//
//			new_total_weight += edge_weight[e];
//
//			//update dscore
//			if (v_in_c[edge[e].v1] + v_in_c[edge[e].v2] == 0) {
//				dscore[edge[e].v1] += edge_weight[e];
//				dscore[edge[e].v2] += edge_weight[e];
//			} else if (v_in_c[edge[e].v1] + v_in_c[edge[e].v2] == 1) {
//				if (v_in_c[edge[e].v1] == 1)dscore[edge[e].v1] -= edge_weight[e];
//				else dscore[edge[e].v2] -= edge_weight[e];
//			}
//		}
//		ave_weight = new_total_weight / e_num;
//
//	}
//
//
//	void update_edge_weight() {
//		int i, e;
//		for (i = 0; i < uncov_stack_fill_pointer; ++i) {
//			e = uncov_stack[i];
//
//			edge_weight[e] += 1;
//			dscore[edge[e].v1] += 1;
//			dscore[edge[e].v2] += 1;
//		}
//
//
//		delta_total_weight += uncov_stack_fill_pointer;
//		if (delta_total_weight >= e_num) {
//			ave_weight += 1;
//			delta_total_weight -= e_num;
//		}
//
//		//smooth weights
//		if (ave_weight >= threshold) {
//			forget_edge_weights();
//		}
//
//	}
//
//
//	void cover_LS() {
//		int best_add_v;
//		int e, v1, v2;
//
//		step = 1;
//
//		while (1)// wihin cutoff_time
//			//while(step<=max_steps)
//		{
//			if (uncov_stack_fill_pointer == 0)//update best solution if needed
//			{
//				update_best_sol();
//
//				if (c_size == optimal_size)
//					return;
//
//				update_target_size();
//
//				continue;
//			}
//
//			//if(step>=try_step)// wihin cutoff_time
//			if (step % try_step == 0) {
//				times(&finish);
//				double elap_time = ((double)finish.tms_utime + (double)finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
//
//
//				if (elap_time >= cutoff_time)return;
//			}
//
//
//			update_best_cov_v();
//
//			remove(best_cov_v);
//
//			e = uncov_stack[rand() % uncov_stack_fill_pointer];
//			v1 = edge[e].v1;
//			v2 = edge[e].v2;
//
//			if (conf_change[v1] == 0) best_add_v = v2;
//			else if (conf_change[v2] == 0) best_add_v = v1;
//
//			else {
//				if (dscore[v1] > dscore[v2] || (dscore[v1] == dscore[v2] && time_stamp[v1] < time_stamp[v2]))
//					best_add_v = v1;
//				else best_add_v = v2;
//			}
//
//			add(best_add_v);
//
//			int index = index_in_remove_cand[best_cov_v];
//			index_in_remove_cand[best_cov_v] = 0;
//
//			remove_cand[index] = best_add_v;
//			index_in_remove_cand[best_add_v] = index;
//
//			time_stamp[best_add_v] = time_stamp[best_cov_v] = step;
//
//			tabu_remove = best_add_v;
//
//			update_edge_weight();
//
//			step++;
//
//		}
//
//	}
//
//	/**
//	 *
//	 * @param V
//	 * @param cutoff max time in seconds
//	 * @param opt_size if we know the optimal size, than we can pass it here
//	 * @param seeed
//	 * @return vertex cover of V
//	 */
//	VI solve(/*int argc, char *argv[]*/ VVI & V, double cutoff, VI initMis = VI(), int opt_size = 0, int seeed = 171 ) {
////		if (argc != 5) {
////			cout << "usage: " << endl;
////			return 1;
////		}
//		int seed;
//		//float scale_threshold;
//
//		stringstream infile;
//		int edges = 0;
//		for( auto & v : V ) edges += v.size();
//		edges >>= 1;
//		infile << "p td " << V.size() << " " << edges<< endl;
////		for( auto p : edges ) infile << "e " << p.first+1 << " " << p.second+1 << endl;
//		for( int i=0; i<V.size(); i++ ){
//			for( int d : V[i] ){
//				if( d > i ) infile << "e " << i+1 << " " << d+1 << endl;
//			}
//		}
//
//		if (build_instance(infile/*argv[1]*/) != 1) {
//			cout << "can't open instance file in NuMVC" << endl;
//			cerr << "can't open instance file in NuMVC" << endl;
//			exit(1);
//			return VI();
//		}
//
////		sscanf(argv[2], "%d", &optimal_size);//if you want to stop the algorithm only cutoff time is reached, set optimal_size to 0.
////		sscanf(argv[3], "%d", &seed);
////		sscanf(argv[4], "%d", &cutoff_time);
//
//		optimal_size = opt_size;
//		seed = seeed;
//		cutoff_time = cutoff;
//
//		threshold = (int) (0.5 * v_num);
//
//		srand(seed);
//
//		//cout<<seed<<' ';
////		cout<<"c This is NuMVC, a local search solver for the Minimum Vertex Cover (and also Maximum Independent Set) problem."<<endl;
//
//		times(&start);
//		start_time = start.tms_utime + start.tms_stime;
//
//		init_sol(initMis);
//
//		if (c_size + uncov_stack_fill_pointer > optimal_size) {
//			//cout<<"c Start local search..."<<endl;
//			cover_LS();
//		}
//
//		VI res;
//		//check solution
//		if (check_solution() == 1) {
////			cout << "c Best found vertex cover size = " << best_c_size << endl;
////			cout << best_c_size << endl;
////			print_solution();
////			cout << "c searchSteps = " << best_step << endl;
////			cout << "c solveTime = " << best_comp_time << endl;
//
//			//cout<<best_c_size<<' '<<best_comp_time<<' '<<best_step<<endl;
//
//			res = get_solution_mis();
//
//		} else {
//			cout << "the solution is wrong." << endl;
//			exit(1);
//			//print_solution();
//		}
//
//		free_memory();
//
////		return 0;
//
//
//		return res;
//	}
//
//
//
//
//
//};
//
//
//#endif //ALGORITHMSPROJECT_NuMVC_H