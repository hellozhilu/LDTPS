#include"func_state.h"
#include"global_variables.h"



//Depth-First Search (DFS)
void dfs(int v, int *sol, int *visited)
{
	visited[v] = 1;
	Compo_sol[v] = Compo_cnt;
	Compo_len++;
//	printf("v: %d, visited: %d\n", v + 1, visited[v]);

	for (int i = 0; i < Num_v; i++)
	{
		if (Edge[v][i] > 0 && sol[i] > 0)
		{
//			printf("v1: %d, visited: %d\n", i + 1, visited[i]);
			if (!visited[i])
				dfs(i, sol, visited);
		}
	}
//	printf("\n");
}


//calculate connected Compos of the current solution IUC (MPC)
void calculate_connect_Compo(int *sol)
{
	int *visited = new int[Num_v];
	memset(visited, 0, sizeof(int) * Num_v);
	Compo_cnt = 0;
	Compo_sol = new int[Num_v];
	for (int i = 0; i < Num_v; i++)
		Compo_sol[i] = -1;

	for (int i = 0; i < Num_v; i++)
	{
		if (!visited[i] && sol[i])
		{
			Compo_len = 0;
			dfs(i, sol, visited);
			if (Compo_len > 1)
				Compo_cnt++;
			else
			{
				visited[i] = 0;
				Compo_sol[i] = -1;
				continue;
			}
		}
	}

	printf("\nFinding totally %d Compos\n", Compo_cnt);
#ifdef DEBUG
	for (int i = 0; i < Compo_cnt; i++)
	{
		printf("%d Compo: [", i);
		for (int j = 0; j < Num_v; j++)
			if (Compo_sol[j] == i)
				printf("%d ", j + 1);
		printf("]\n");
	}
	printf("\n");
#endif

	//greedy algorithm for finding clusters
	int **cur_sol = new int*[Compo_cnt];
	for (int i = 0; i < Compo_cnt; i++)
		cur_sol[i] = new int[Num_v];
	int *cur_sol_len = new int[Compo_cnt];
	memset(cur_sol_len, 0, sizeof(int) * Compo_cnt);
	int *sol_temp = new int[Num_v];
	int *deg_temp = new int[Compo_cnt];
	memset(sol_temp, 0, sizeof(int) * Num_v);
	memset(deg_temp, 0, sizeof(int) * Compo_cnt);

	for (int i = 0; i < Compo_cnt; i++)
		for (int j = 0; j < Num_v; j++)
			if (Compo_sol[j] == i)
				cur_sol[i][cur_sol_len[i]++] = j;

	for (int i = 0; i < Num_v; i++)
	{
		if (!visited[i])
		{
			int max_deg = -1;
			int candidate_deg_c[MAXCAN], candidate_len = 0;
			memset(deg_temp, 0, sizeof(int) * Compo_cnt);

			for (int x = 0; x < Compo_cnt; x++)
			{
				for (int y = 0; y < cur_sol_len[x]; y++)
				{
					int elem = cur_sol[x][y];
					if (Edge[elem][i] > 0)
						deg_temp[x]++;
				}
				if (deg_temp[x] == cur_sol_len[x] && deg_temp[x] > max_deg)
				{
					max_deg = deg_temp[x];
					candidate_len = 0;
					candidate_deg_c[candidate_len++] = x;
				}
				else if (deg_temp[x] == cur_sol_len[x] && deg_temp[x] == max_deg)
					candidate_deg_c[candidate_len++] = x;
			}

			if (candidate_len > 0)
			{
				int sel_c = candidate_deg_c[rand() % candidate_len];
				cur_sol[sel_c][cur_sol_len[sel_c]++] = i;
				visited[i] = 1;
			}
		}
	} //for i

	int count1 = 0;
	int clust = 0;
	printf("After greedy, finding totally %d Compos\n", Compo_cnt);
	for (int i = 0; i < Compo_cnt; i++)
	{
		count1 += cur_sol_len[i];
		if (cur_sol_len[i] > 1)
		{
			printf("%d Compo: [", clust);
			for (int j = 0; j < cur_sol_len[i]; j++)
			{
				sol_temp[cur_sol[i][j]] = clust + 1;
				printf("%d ", cur_sol[i][j] + 1);
			}
			clust++;
			printf("]\n");
		}
	}
	printf("total count: %d\n\n", count1);

	//output result for cliques label
	FILE *fp1, *fp2;
	char buff[MAXCAN];
	char *graph_name = basename(Instance_name);

	sprintf(buff, "./%s_split_test.csv", graph_name);
	fp1 = fopen(buff, "w+");
	if (fp1 == NULL)
		exit(-1);

	fprintf(fp1, "Id,Modularity Class\n");
	for (int i = 0; i < Num_v; i++)
		fprintf(fp1, "%d,%d\n", i + 1, sol_temp[i]);
	fclose(fp1);

	//assign rest vertices to form clusters
	for (int i = 0; i < Num_v; i++)
	{
		if (!visited[i])
		{
			int max_deg = -1;
			int candidate_deg_c[MAXCAN], candidate_len = 0;
			memset(deg_temp, 0, sizeof(int) * Compo_cnt);

			for (int x = 0; x < Compo_cnt; x++)
			{
				for (int y = 0; y < cur_sol_len[x]; y++)
				{
					int elem = cur_sol[x][y];
					if (Edge[elem][i] > 0)
						deg_temp[x]++;
				}
				if (deg_temp[x] > max_deg)
				{
					max_deg = deg_temp[x];
					candidate_len = 0;
					candidate_deg_c[candidate_len++] = x;
				}
				else if (deg_temp[x] == max_deg)
					candidate_deg_c[candidate_len++] = x;
			}

			if (candidate_len > 0)
			{
				int sel_c = candidate_deg_c[rand() % candidate_len];
				cur_sol[sel_c][cur_sol_len[sel_c]++] = i;
				visited[i] = 1;
			}
		}
	} //for i

	int count2 = 0;
	printf("After clustering, finding totally %d clusters\n", Compo_cnt);
	for (int i = 0; i < Compo_cnt; i++)
	{
		count2 += cur_sol_len[i];
		printf("%d Compo: [", i);
		for (int j = 0; j < cur_sol_len[i]; j++)
		{
			sol_temp[cur_sol[i][j]] = i;
			printf("%d ", cur_sol[i][j] + 1);
		}
		printf("]\n");
	}
	printf("total count: %d\n\n", count2);

	//output result for clusters
	sprintf(buff, "./%s_split_test_all.csv", graph_name);
	fp2 = fopen(buff, "w+");
	if (fp2 == NULL)
		exit(-1);

	fprintf(fp2, "Id,Modularity Class\n");
	for (int i = 0; i < Num_v; i++)
		fprintf(fp2, "%d,%d\n", i + 1, sol_temp[i]);
	fclose(fp2);

	//release memory
	delete[] Compo_sol; Compo_sol = NULL;
	for (int i = 0; i < Compo_cnt; i++)
	{
		delete[] cur_sol[i]; cur_sol[i] = NULL;
	}
	delete[] cur_sol; cur_sol = NULL;
	delete[] visited; visited = NULL;
	delete[] cur_sol_len; cur_sol_len = NULL;
	delete[] sol_temp; sol_temp = NULL;
	delete[] deg_temp; deg_temp = NULL;
}


//time complexity: O(V)
void build_prob_matrix()
{
	for (int i = 0; i < Num_v; i++)
		for (int j = 0; j < 2; j++)
			Prob_matrix[i][j] = 0.5;
}


//time complexity: O(V)
void update_prob_matrix(int *tmp_sol, int *cur_sol)
{
	for (int i = 0; i < Num_v; i++)
	{
		if (cur_sol[i] == tmp_sol[i])  //if vertex i stays in its original set
		{
			int k = tmp_sol[i];
			double pk = Prob_matrix[i][k];

			Prob_matrix[i][k] = Alpha + (1 - Alpha) * pk;  //reward
			Prob_matrix[i][1 - k] = (1 - Alpha) * pk;
		}
		else //if vertex i makes a move
		{
			int m = tmp_sol[i];
			int n = cur_sol[i];
			double pm = Prob_matrix[i][m];
			double pn = Prob_matrix[i][n];

			Prob_matrix[i][m] = (1 - Gamma) * (1 - Beta) * pm;    //penalize
			Prob_matrix[i][n] = Gamma + (1 - Gamma) * Beta
								+ (1 - Gamma) * (1 - Beta) * pn;  //compensate
		}
	}
}


#ifdef DEBUG
void smooth_prob_matrix_tmp()
{
	for (int i = 0; i < Num_v; i++)
		for (int j = 0; j < 2; j++)
			if (Prob_matrix[i][j] > P0 + PRECISION)
				Prob_matrix[i][j] *= P;
}
#endif


#ifdef DEBUG
//time complexity: O(V)
void smooth_prob_matrix()
{
	for (int i = 0; i < Num_v; i++)
	{
		double max_prob = -MAXNUM;
		int sel = 0;
		for (int j = 0; j < 2; j++)
		{
			if (Prob_matrix[i][j] > max_prob + PRECISION)
			{
				max_prob = Prob_matrix[i][j];
				sel = j;
			}
		}

		if (max_prob > 0.995 + PRECISION)
		{
			for (int k = 0; k < 2; k++)
			{
				if (fabs(Prob_matrix[i][k] - max_prob) <= PRECISION)
					Prob_matrix[i][k] *= 0.5;
				else
					Prob_matrix[i][k] = (1 - 0.5) * max_prob + Prob_matrix[i][k];
			}
		}
	}
}
#endif


//time complexity: O(V * K_fixed)
void build_delta_matrix(int *sol, int help_len, int &ff)
{
	memset(Delta_matrix[0], 0, sizeof(int) * Num_v);
	memset(Delta_matrix[1], 0, sizeof(int) * Num_v);

	for (int i = 0; i < Num_v; i++)
	{
		for (int j = 0; j < help_len; j++)
		{
			int vj = Help_sol[j];
			for (int k = 0; k < help_len; k++)
			{
				int vk = Help_sol[k];
				if (vk != vj && vj != i && i != vk && Edge[vk][vj] && Edge[vj][i] && !Edge[i][vk])
					Delta_matrix[0][i]++;
			}
		}

		for (int j = 0; j < help_len; j++)
		{
			int vj = Help_sol[j];
			for (int k = 0; k < help_len; k++)
			{
				int vk = Help_sol[k];
				if (vk != i && i != vj && vj != vk && Edge[vk][i] && Edge[i][vj] && !Edge[vj][vk])
					Delta_matrix[1][i]++;
			}
		}
	}

	ff = 0;
	for (int i = 0; i < Num_v; i++)
	{
		if (sol[i])
		{
			ff += Delta_matrix[0][i];
			ff += Delta_matrix[1][i];
		}
	}
	ff /= 4;

#ifdef DEBUG
	ff = 0;
	for (int i = 0; i < Num_v; i++)
	{
		if (sol[i])
			ff += Delta_matrix[0][i]; //or ff += Delta_matrix[1][i];
	}
	ff /= 2;

	for (int i = 0; i < Num_v; i++)
		printf("i: %d, sol: %d, delta_0: %d, delta_1: %d\n",
				i, sol[i], Delta_matrix[0][i], Delta_matrix[1][i]);
	printf("in build func, ff: %d\n", ff);
#endif
}


#ifdef DEBUG
void build_delta_matrix_tmp(int *sol, int help_len, int &ff)
{
	memset(Delta_matrix_tmp[0], 0, sizeof(int) * Num_v);
	memset(Delta_matrix_tmp[1], 0, sizeof(int) * Num_v);

	for (int i = 0; i < Num_v; i++)
	{
		for (int j = 0; j < help_len; j++)
		{
			int vj = Help_sol[j];
			for (int k = 0; k < help_len; k++)
			{
				int vk = Help_sol[k];
				if (vk != vj && vj != i && i != vk && Edge[vk][vj] && Edge[vj][i] && !Edge[i][vk])
					Delta_matrix_tmp[0][i]++;
			}
		}

		for (int j = 0; j < help_len; j++)
		{
			int vj = Help_sol[j];
			for (int k = 0; k < help_len; k++)
			{
				int vk = Help_sol[k];
				if (vk != i && i != vj && vj != vk && Edge[vk][i] && Edge[i][vj] && !Edge[vj][vk])
					Delta_matrix_tmp[1][i]++;
			}
		}

	}
	ff = 0;
	for (int i = 0; i < Num_v; i++)
	{
		if (sol[i])
		{
			ff += Delta_matrix_tmp[0][i];
			ff += Delta_matrix_tmp[1][i];
		}
	}
	ff /= 4;

//	for (int i = 0; i < Num_v; i++)
//		printf("i: %d, sol: %d, delta_0: %d, delta_1: %d\n",
//				i, sol[i], Delta_matrix_tmp[0][i], Delta_matrix_tmp[1][i]);
//	printf("in build_tmp func, ff: %d\n", ff);
}


void compare_delta(int f1, int f2)
{	
	for (int i = 0; i < Num_v; i++)
	{
		if (Delta_matrix[0][i] != Delta_matrix_tmp[0][i] || Delta_matrix[1][i] != Delta_matrix_tmp[1][i])
		{
			printf("### Error, delta_matrix, in node: %d, delta_0: %d %d, delta_1: %d %d ###\n",
					i, Delta_matrix[0][i], Delta_matrix_tmp[0][i], Delta_matrix[1][i], Delta_matrix_tmp[1][i]);
			exit(-1);
		}
	}
	if (f1 != f2)
	{
		printf("### Error, f1 != f2, f1: %d, f2: %d ###\n", f1, f2);
		exit(-1);
	}
}
#endif


//vertex u is out, vertex v is in
//time complexity: O(V * K_fixed)
void update_delta_matrix(int u, int v, int help_len)
{
	for (int i = 0; i < Num_v; i++)
	{
		for (int k = 0; k < help_len; k++)
		{
			int vk = Help_sol[k];

			//update Delta_matrix[0][i]
			if (vk != i && Edge[i][u] && Edge[u][vk] && !Edge[vk][i]) //triple <i, u, vk>
				Delta_matrix[0][i]--;
			if (u != i && Edge[i][vk] && Edge[vk][u] && !Edge[u][i])  //triple <i, vk, u>
				Delta_matrix[0][i]--;			
			if (vk != i && vk != u && Edge[i][v] && Edge[v][vk] && !Edge[vk][i]) //triple <i, v, vk>
				Delta_matrix[0][i]++;
			if (v != i && vk != u && Edge[i][vk] && Edge[vk][v] && !Edge[v][i])  //triple <i, vk, v>
				Delta_matrix[0][i]++;		
		
			//update Delta_matrix[1][i]
			if (vk != u && Edge[u][i] && Edge[i][vk] && !Edge[vk][u]) //triple <u, i, vk>
				Delta_matrix[1][i]--;
			if (u != vk && Edge[vk][i] && Edge[i][u] && !Edge[u][vk]) //triple <vk, i, u>
				Delta_matrix[1][i]--;
			if (vk != v && vk != u && Edge[v][i] && Edge[i][vk] && !Edge[vk][v]) //triple <v, i, vk>
				Delta_matrix[1][i]++;
			if (vk != v && vk != u && Edge[vk][i] && Edge[i][v] && !Edge[v][vk]) //triple <vk, i, v>
				Delta_matrix[1][i]++;
		}				
	}

#ifdef DEBUG	
	for (int i = 0; i < Num_v; i++)			
		printf("i: %d, delta_0: %d, delta_1: %d\n", i, Delta_matrix[0][i], Delta_matrix[1][i]);
	printf("in update_delta_matrix func\n");
#endif
}


//greedy create a clique, return its size as the initial value of k
//time complexity: O(V)
int initial_k()
{
	int *flag_in_c = new int[Num_v];
	int *candidate_v = new int[Num_v];
	int *clique = new int[Num_v];
	memset(flag_in_c, 0, sizeof(int) * Num_v);
	int clique_size = 0;
	int rand_v = rand() % Num_v;
	clique[clique_size++] = rand_v;
	flag_in_c[rand_v] = 1;
	int edge_clique = 0;

	while (1)
	{
		int max_deg = -1;
		int candidate_len = 0;
		for (int i = 0; i < Num_v; i++)
		{
			int ele1 = i;
			if (!flag_in_c[ele1])
			{
				int deg = 0;
				for (int j = 0; j < clique_size; j++)
				{
					int ele2 = clique[j];
					if (Edge[ele1][ele2])
						deg++;
				}
				if (deg > max_deg)
				{
					max_deg = deg;
					candidate_len = 0;
					candidate_v[candidate_len++] = ele1;
				}
				else if (deg == max_deg)
					candidate_v[candidate_len++] = ele1;
			}
		}
		int necessary_num_edge = clique_size * (clique_size + 1) / 2; //because we add a new vertex in a candidate clique
		if (edge_clique + max_deg >= necessary_num_edge) //form a clique, or ==
		{
			int sel_v = candidate_v[rand() % candidate_len];
			clique[clique_size++] = sel_v;
			flag_in_c[sel_v] = 1;
			edge_clique += max_deg;
		}
		else
			break;
	}

#ifdef DEBUG
	printf("clique_size: %d ", clique_size);
	for (int i = 0; i < clique_size; i++)
	{
		int vi = clique[i];
		for (int j = i + 1; j < clique_size; j++)
		{
			int vj = clique[j];
			if (!Edge[vi][vj])
			{
				printf("not clique, vi: %d, vj: %d\n", vi, vj);
				exit(-1);
			}
		}
	}
#endif

	if (clique_size > 0)
	{
		Best_k = clique_size;
		memcpy(Best_sol_one_run, flag_in_c, sizeof(int) * Num_v);
		Run_time = (double) (clock() - Start_time) / CLOCKS_PER_SEC;
	}

	//free
	delete[] flag_in_c; flag_in_c = NULL;
	delete[] candidate_v; candidate_v = NULL;
	delete[] clique; clique = NULL;

	return clique_size;
}


void random_selection(int *sol)
{
	int len = 0;
	memset(sol, 0, sizeof(int) * Num_v);
	while (len < K_fixed)
	{
		int rx = rand() % Num_v;
		while (sol[rx])
			rx = rand() % Num_v;
		sol[rx] = 1;
		len++;
	}
}


void greedy_selection(int *sol)
{
#ifdef DEBUG
	for (int i = 0; i < Num_v; i++)
	{
		for (int j = 0; j < 2; j++)
			printf("%.2f ", Prob_matrix[i][j]);
		printf("\n");
	}
#endif

	int len = 0;
	memset(sol, 0, sizeof(int) * Num_v);
	while (len < K_fixed)
	{
		int rx = rand() % Num_v;
		while (sol[rx])
			rx = rand() % Num_v;

		if (Prob_matrix[rx][1] > Prob_matrix[rx][0] + PRECISION)
		{
			sol[rx] = 1;
			len++;
		}
		else if (fabs(Prob_matrix[rx][1] - Prob_matrix[rx][0]) <= PRECISION)
		{
			if (rand() % 100 < 50)
			{
				sol[rx] = 1;
				len++;
			}
			else
				sol[rx] = 0;
		}
		else  //in case of endless loop
		{
			if (rand() % 100 < 50)
			{
				sol[rx] = 1;
				len++;
			}
		}
	} //while
}


//time complexity: O(K_fixed)
void hybrid_initial_sol(int *sol)
{
	if (((double) (rand() % 101) / 100.0) < 0.2)
		random_selection(sol);
	else
		greedy_selection(sol);
}


//swap neighborhood based descent local search
//time complexity: O(V^2)
int descent_local_search(int *sol)
{
	int help_len = 0;
	for (int i = 0; i < Num_v; i++)
	{
		if (sol[i] > 0)
		{
			Help_sol[help_len] = i;
			Address_help_sol[i] = help_len;
			help_len++;
		}
	}

	int f_sol;
	build_delta_matrix(sol, help_len, f_sol);
	memcpy(Local_best_sol, sol, sizeof(int) * Num_v);
	int f_local_best = f_sol;
	int flag_improve = 1;
//	printf("in descent local search func, f_sol: %d\n", f_sol);

	while (flag_improve && (double) (clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	{	
		flag_improve = 0;
		for (int i = 0; i < Num_v; i++)
		{
			if (sol[i])
			{
				for (int j = 0; j < Num_v; j++)
				{
					if (!sol[j])
					{
						double delta = (Delta_matrix[0][j] + Delta_matrix[1][j] / 2) - (Delta_matrix[0][i] + Delta_matrix[1][i] / 2);
						//calculate the number of open triangles related to j and i
						int dev = 0;
						for (int m = 0; m < help_len; m++)  //triple <j, i, vm>
						{
							int vm = Help_sol[m];
							if (Edge[j][i] && Edge[i][vm] && !Edge[vm][j] && vm != i)
								dev++;
						}
						for (int m = 0; m < help_len; m++)  //triple <j, vm, i>
						{
							int vm = Help_sol[m];
							if (Edge[j][vm] && Edge[vm][i] && !Edge[i][j] && vm != i)
								dev++;
						}
						for (int m = 0; m < help_len; m++)  //triple <i, j, vm>
						{
							int vm = Help_sol[m];
							if (Edge[i][j] && Edge[j][vm] && !Edge[vm][i] && vm != i)
								dev++;
						}
						delta -= dev;  //the complexity of calculating move gain is bounded by O(K_fixed)

						if (delta < 0)
						{
							update_delta_matrix(i, j, help_len);
							sol[i] = 0;
							sol[j] = 1;
							f_sol += delta;
							int pos = Address_help_sol[i];
							Help_sol[pos] = j;
							Address_help_sol[j] = pos;
							Freq[i]++;
							Freq[j]++;
							flag_improve = 1;
							if (f_sol < f_local_best)
							{
								memcpy(Local_best_sol, sol, sizeof(int) * Num_v);
								f_local_best = f_sol;
							}
							if (f_sol == 0)
							{
								memcpy(sol, Local_best_sol, sizeof(int) * Num_v);
								return f_sol;
							}
							break;
						}					
					}
				} //for j
			}
		} //for i
	} //while

	memcpy(sol, Local_best_sol, sizeof(int) * Num_v);
	f_sol = f_local_best;
	return f_sol;
}


//ascending order
int compare_ascending(const void *a, const void *b)
{
	return (*(int*) a - *(int*) b);
}


//swap neighborhood based tabu search
//time complexity: O(V * K_fixed)
int tabu_search(int *sol)
{
	int *data_ins = new int[Num_v];
	int *data_outs = new int[Num_v];

	int tabu_swap_best_delta, swap_best_delta, delta;
	int tabu_swap_best_num, swap_best_num;  
	int swap_best_u[MAXCAN], swap_best_v[MAXCAN];
	int tabu_swap_best_u[MAXCAN], tabu_swap_best_v[MAXCAN];

	int non_improve = 0;
	int iter = 0;
	memset(Tabu_list[0], 0, sizeof(int) * Num_v);
	memset(Tabu_list[1], 0, sizeof(int) * Num_v);
	int help_len = 0;
	for (int i = 0; i < Num_v; i++)
	{
		if (sol[i] > 0)
		{
			Help_sol[help_len] = i;
			Address_help_sol[i] = help_len;
			help_len++;
		}
	}
	int f_sol;
	build_delta_matrix(sol, help_len, f_sol);
	memcpy(Local_best_sol, sol, sizeof(int) * Num_v);
	int f_local_best = f_sol;

	while (non_improve < Tabu_depth && (double) (clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	{
		swap_best_num = 0;
		tabu_swap_best_num = 0;
		swap_best_delta = MAXNUM;
		tabu_swap_best_delta = MAXNUM;
		int len_ins = 0, len_outs = 0;

		//prepare for parametric constraint neighborhood
		for (int i = 0; i < Num_v; i++)
		{
			if (sol[i])
				data_ins[len_ins++] = Delta_matrix[0][i] + Delta_matrix[1][i] / 2;
			else
				data_outs[len_outs++] = Delta_matrix[0][i] + Delta_matrix[1][i] / 2;
		}
		qsort(data_ins, len_ins, sizeof(int), compare_ascending);
		qsort(data_outs, len_outs, sizeof(int), compare_ascending);

		int pos_ins = int((1 - Cn_rho) * len_ins);  //TODO important
		int pos_outs = int(Cn_rho * len_outs);
		if (pos_ins == len_ins)
			pos_ins = len_ins - 1;
		if (pos_outs == len_outs)
			pos_outs = len_outs - 1;
		int delta_pos_ins = data_ins[pos_ins];
		int delta_pos_outs = data_outs[pos_outs];

#ifdef DEBUG
		for (int i = 0; i < len_ins; i++)
			printf("%d ", data_ins[i]);
		printf("\n");
		for (int i = 0; i < len_outs; i++)
			printf("%d ", data_outs[i]);
		printf("\nlen_ins: %d, len_outs: %d\n", len_ins, len_outs);
#endif

		//parametric constraint neighborhood, neighborhood size: Cn_rho^2 * K_fixed * (V - K_fixed)
		for (int i = 0; i < Num_v; i++)
		{
			if (sol[i] && Delta_matrix[0][i] + Delta_matrix[1][i] / 2 >= delta_pos_ins)
			{				
				for (int j = 0; j < Num_v; j++)
				{
					if (!sol[j] && Delta_matrix[0][j] + Delta_matrix[1][j] / 2 <= delta_pos_outs)
					{
						delta = (Delta_matrix[0][j] + Delta_matrix[1][j] / 2) - (Delta_matrix[0][i] + Delta_matrix[1][i] / 2);
						//calculate the number of open triangles related to j and i
						int dev = 0;
						for (int m = 0; m < help_len; m++)  //triple <j, i, vm>
						{
							int vm = Help_sol[m];
							if (Edge[j][i] && Edge[i][vm] && !Edge[vm][j] && vm != i)
								dev++;
						}
						for (int m = 0; m < help_len; m++)  //triple <j, vm, i>
						{
							int vm = Help_sol[m];
							if (Edge[j][vm] && Edge[vm][i] && !Edge[i][j] && vm != i)
								dev++;
						}
						for (int m = 0; m < help_len; m++)  //triple <i, j, vm>
						{
							int vm = Help_sol[m];
							if (Edge[i][j] && Edge[j][vm] && !Edge[vm][i] && vm != i)
								dev++;
						}
						delta -= dev;  //the complexity of calculating move gain is bounded by O(K_fixed)

						if (Tabu_list[0][i] <= iter && Tabu_list[1][j] <= iter)
						{
							if (delta < swap_best_delta)
							{
								swap_best_num = 0;
								swap_best_u[swap_best_num] = i;
								swap_best_v[swap_best_num] = j;
								swap_best_delta = delta;
								swap_best_num++;
							}
							else if (delta == swap_best_delta && swap_best_num < MAXCAN)
							{
								swap_best_u[swap_best_num] = i;
								swap_best_v[swap_best_num] = j;
								swap_best_num++;
							}
						}
						else if (Tabu_list[0][i] > iter || Tabu_list[1][j] > iter) //note ||
						{
							if (delta < tabu_swap_best_delta)
							{
								tabu_swap_best_num = 0;
								tabu_swap_best_u[tabu_swap_best_num] = i;
								tabu_swap_best_v[tabu_swap_best_num] = j;
								tabu_swap_best_delta = delta;
								tabu_swap_best_num++;
							}
							else if (delta == tabu_swap_best_delta && tabu_swap_best_num < MAXCAN)
							{
								tabu_swap_best_u[tabu_swap_best_num] = i;
								tabu_swap_best_v[tabu_swap_best_num] = j;
								tabu_swap_best_num++;
							}
						}
					}
				}
			}
		}
		int u = -1;
		int v = -1;
		int f_move = 0;
		if ((tabu_swap_best_num > 0 && tabu_swap_best_delta < swap_best_delta && f_sol + tabu_swap_best_delta < f_local_best)
			|| (swap_best_num == 0 && tabu_swap_best_num > 0))  //aspiration criterion
		{
			int rx = rand() % tabu_swap_best_num;
			u = tabu_swap_best_u[rx];
			v = tabu_swap_best_v[rx];
			f_move = tabu_swap_best_delta;
		}
		else
		{
			if (swap_best_num > 0)
			{
				int rx = rand() % swap_best_num;
				u = swap_best_u[rx];
				v = swap_best_v[rx];
				f_move = swap_best_delta;
			}
		}
		if (u != -1 && v != -1)  //perform the swap move (u, v)
		{
			update_delta_matrix(u, v, help_len);			
			sol[u] = 0;
			sol[v] = 1;
			f_sol += f_move;						
			int pos = Address_help_sol[u];
			Help_sol[pos] = v;
			Address_help_sol[v] = pos;
			Tabu_list[1][u] = iter + Tabu_tenure;
			Tabu_list[0][v] = iter + Tabu_tenure /*+ (Num_v - K_fixed) / K_fixed*/;
			Freq[u]++;
			Freq[v]++;
		}

		iter++;
		if (f_sol < f_local_best)
		{
			memcpy(Local_best_sol, sol, sizeof(int) * Num_v);
			f_local_best = f_sol;
			non_improve = 0;
		}
		else
			non_improve++;

#ifdef DEBUG
		int tmp_ff;
		build_delta_matrix_tmp(sol, help_len, tmp_ff);
		compare_delta(f_sol, tmp_ff);
#endif
//		if (iter % 1000 == 0)
//			printf("iter: %d, K_fixed: %d, f_sol: %d, f_move: %d, u: %d, v: %d\n",
//					iter, K_fixed, f_sol, f_move, u, v);
		if (f_sol == 0)
			break;
	}	

	memcpy(sol, Local_best_sol, sizeof(int) * Num_v);
	f_sol = f_local_best;
	delete[] data_ins; data_ins = NULL;
	delete[] data_outs; data_outs = NULL;

	return f_sol;
}


//destroy process (frequency remove) + repair procedure (random add)
//time complexity: O(V) + O(K_fixed)
void freq_driven_perturabtion(int *sol)
{
	//calculate average freq value (with respect to subset D)
	int sum_freq = 0;
	double avg_freq;
	for (int i = 0; i < Num_v; i++)
		if (sol[i])
			sum_freq += Freq[i];
	avg_freq = (double) sum_freq / (double) K_fixed;

	//destroy process: remove vertices in D with small freq value,
	//freq value < average freq value
	int len = K_fixed;
	for (int i = 0; i < Num_v; i++)
	{
		if (sol[i] && (double)Freq[i] < avg_freq - PRECISION)
		{
			sol[i] = 0;
			len--;
		}
		if (len == 0)
			break;
	}
	
	//repair procedure: random add vertices to D
	while (len < K_fixed)
	{
		int rx = rand() % Num_v;
		while (sol[rx])
			rx = rand() % Num_v;
		sol[rx] = 1;
		len++;
	}
}


//time complexity: O(V)
void fresh_freq()
{
	for (int i = 0; i < Num_v; i++)
		if (Freq[i] > K_fixed)
			Freq[i] = 0;
}


#ifdef DEBUG
int greedy_expansion(int *sol)
{
	int *flag_in_c = new int[Num_v];
	int *candidate_v = new int[Num_v];
	memset(flag_in_c, 0, sizeof(int) * Num_v);

	int help_len = 0;
	for (int i = 0; i < Num_v; i++)
	{
		if (sol[i] > 0)
		{
			Help_sol[help_len] = i;
			Address_help_sol[i] = help_len;
			help_len++;
			flag_in_c[i] = 1;
		}
	}
	int edge_clique = help_len * (help_len - 1) / 2;

	while (1)
	{
		int max_deg = -1;
		int candidate_len = 0;
		for (int i = 0; i < Num_v; i++)
		{
			int ele1 = i;
			if (!flag_in_c[ele1])
			{
				int deg = 0;
				for (int j = 0; j < help_len; j++)
				{
					int ele2 = Help_sol[j];
					if (Edge[ele1][ele2])
						deg++;
				}
				if (deg > max_deg)
				{
					max_deg = deg;
					candidate_len = 0;
					candidate_v[candidate_len++] = ele1;
				}
				else if (deg == max_deg)
					candidate_v[candidate_len++] = ele1;
			}
		}
		int necessary_num_edge = help_len * (help_len + 1) / 2; //because we add a new vertex in a candidate clique
		if (edge_clique + max_deg >= necessary_num_edge) //form a clique, or ==
		{
			int sel_v = candidate_v[rand() % candidate_len];
			Help_sol[help_len++] = sel_v;
			flag_in_c[sel_v] = 1;
			edge_clique += max_deg;
		}
		else
			break;
	}

#ifdef DEBUG
	printf("help_len: \n", help_len);
	for (int i = 0; i < help_len; i++)
	{
		int vi = Help_sol[i];
		for (int j = i + 1; j < help_len; j++)
		{
			int vj = Help_sol[j];
			if (!Edge[vi][vj])
			{
				printf("not clique, vi: %d, vj: %d\n", vi, vj);
				exit(-1);
			}
		}
	}
#endif

	if (help_len == K_opt)  //find the optimal solution
	{
		Best_k = K_fixed;
		memcpy(Best_sol_one_run, flag_in_c, sizeof(int) * Num_v);
		Run_time = (double) (clock() - Start_time) / CLOCKS_PER_SEC;
	}

	delete[] flag_in_c; flag_in_c = NULL;
	delete[] candidate_v; candidate_v = NULL;

	return help_len;
}
#endif


//learning driven three-phase search (LDTPS)
//total time complexity: O(V^2 + V * K_fixed)
void local_search_ldtps()
{
	//obtained initial k value
	K_fixed = initial_k();  //greedy create a clique, return its size as the initial value of K_fixed
	printf(", initial K_fixed: %d, Run_time: %.2f\n", K_fixed, Run_time);
	if (K_fixed == K_opt || K_fixed == Num_v) //find the optimal solution
		return;

	build_prob_matrix();  //RL: initial probability matrix

	//TODO initial solution is NOT necessary a clique, but with the size of K_fixed
	K_fixed++;
	hybrid_initial_sol(Cur_sol);  //RL: group selection

	memset(Temp_sol, 0, sizeof(int) * Num_v);
	memset(Freq, 0, sizeof(int) * Num_v);

	while ((double) (clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	{
		memcpy(Temp_sol, Cur_sol, sizeof(int) * Num_v);
		int ff_des = descent_local_search(Cur_sol);
		update_prob_matrix(Temp_sol, Cur_sol);  //RL: probability updating

		int ff_ts = -1;
		if (ff_des > 0)
		{
//			printf("ff_des: %d\n", ff_des);
			memcpy(Temp_sol, Cur_sol, sizeof(int) * Num_v);
			ff_ts = tabu_search(Cur_sol);
			update_prob_matrix(Temp_sol, Cur_sol);  //RL: probability updating
		}

		//the number of open triangles = 0,
		//stop and store: Best_sol_one_run[] and K_fixed,
		//then K_fixed++
		if (ff_des == 0 || ff_ts == 0)
		{
#ifdef DEBUG
			//greedy expansion
			int tmp_k = K_fixed;
			int cur_k = greedy_expansion(Cur_sol);
			if (cur_k > K_fixed)
				K_fixed = cur_k;
			printf("variation: %d\n", tmp_k - K_fixed);
#endif

			Best_k = K_fixed;
			memcpy(Best_sol_one_run, Cur_sol, sizeof(int) * Num_v);
			Run_time = (double) (clock() - Start_time) / CLOCKS_PER_SEC;

			if (Best_k == K_opt || K_fixed == Num_v)  //find the optimal solution
				break;

			K_fixed++;
			hybrid_initial_sol(Cur_sol);  //RL: group selection
		}
		else		
			freq_driven_perturabtion(Cur_sol);  //frequency driven perturbation

//		smooth_prob_matrix();  //RL: smoothing probability matrix
		fresh_freq();  //freshing freq array

		printf("next K_fixed: %d, Best_k: %d, Run_time: %.2f, ff_des: %d, ff_ts: %d\n",
				K_fixed, Best_k, Run_time, ff_des, ff_ts);
	}

	//update global best
	if (Best_k > G_best_k)
	{
		G_best_k = Best_k;
		memcpy(Global_best_sol, Best_sol_one_run, sizeof(int) * Num_v);
	}
}
