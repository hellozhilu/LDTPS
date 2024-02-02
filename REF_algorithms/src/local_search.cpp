#include "func_state.h"
#include "global_variables.h"


int calculate_sol(int *sol)
{
	int f_obj = 0;
	int len = 0;
	for (int i = 0; i < Num_v; i++)
	{
		if (sol[i] > 0)
		{
			len++;
			for (int j = i + 1; j < Num_v; j++)
			{
				if (sol[j] > 0)
				{
					for (int k = j + 1; k < Num_v; k++)
					{
						if (sol[k] > 0)
						{
							if (Edge[i][j] && Edge[j][k] && !Edge[k][i])
								f_obj++;
							if (Edge[i][j] && !Edge[j][k] && Edge[k][i])
								f_obj++;
							if (!Edge[i][j] && Edge[j][k] && Edge[k][i])
								f_obj++;
						}
					} //for k
				}
			} //for j
		}
	} //for i
	if (len != K_fixed)
	{
		printf("in calculate_sol func, an error is detected, len != K_fixed, && len: %d, K_fixed: %d\n", len, K_fixed);
		exit(-1);
	}

	return f_obj;
}


void build_delta_matrix(int *sol, int &ff)
{
	memset(Delta_matrix[0], 0, sizeof(int) * Num_v);
	memset(Delta_matrix[1], 0, sizeof(int) * Num_v);

	for (int i = 0; i < Num_v; i++)
	{
		for (int j = 0; j < Help_len; j++)
		{
			int vj = Help_sol[j];
			for (int k = 0; k < Help_len; k++)
			{
				int vk = Help_sol[k];
				if (vk != vj && vj != i && i != vk && Edge[vk][vj] && Edge[vj][i] && !Edge[i][vk])
					Delta_matrix[0][i]++;
			}
		}

		for (int j = 0; j < Help_len; j++)
		{
			int vj = Help_sol[j];
			for (int k = 0; k < Help_len; k++)
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
		if (sol[i] > 0)
		{
			ff += Delta_matrix[0][i];
			ff += Delta_matrix[1][i];
		}
	}
	ff /= 4;

#ifdef DEBUG
	ff = 0;
	for (int i = 0; i < Num_v; i++)
		if (sol[i] > 0)
			ff += Delta_matrix[0][i];  //or ff += Delta_matrix[1][i];
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
		if (sol[i] > 0)
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
void update_delta_matrix(int u, int v)
{
	for (int i = 0; i < Num_v; i++)
	{
		for (int k = 0; k < Help_len; k++)
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
			if (vk != u && Edge[u][i] && Edge[i][vk] && !Edge[vk][u])  //triple <u, i, vk>
				Delta_matrix[1][i] --;
			if (u != vk && Edge[vk][i] && Edge[i][u] && !Edge[u][vk])  //triple <vk, i, u>
				Delta_matrix[1][i]--;
			if (vk != v && vk != u && Edge[v][i] && Edge[i][vk] && !Edge[vk][v])  //triple <v, i, vk>
				Delta_matrix[1][i]++;
			if (vk != v && vk != u && Edge[vk][i] && Edge[i][v] && !Edge[v][vk])  //triple <vk, i, v>
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
int initial_k()
{
	int *flag_in_c = new int[Num_v];
	int *candidate_v = new int[Num_v];
	int *clique = new int[Num_v];
	memset(clique, 0, sizeof(int) * Num_v);
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
		int necessary_num_edge = clique_size * (clique_size + 1) / 2;  //because we add a new vertex in a candidate clique
		if (edge_clique + max_deg >= necessary_num_edge)  //form a clique, or ==
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
		Run_time = (double)(clock() - Start_time) / CLOCKS_PER_SEC;
	}

	//release memory
	delete[] flag_in_c; flag_in_c = NULL;
	delete[] candidate_v; candidate_v = NULL;
	delete[] clique; clique = NULL;

	return clique_size;
}


void random_initial_func(int *sol)
{
	int len = 0;
	memset(sol, 0, sizeof(int) * Num_v);
	while (len < K_fixed)
	{
		int rx = rand() % Num_v;
		while (sol[rx] > 0)
			rx = rand() % Num_v;
		sol[rx] = 1;
		len++;
	}
}


int annealing_search(double init_temp, int *sol)
{
	int frozen_cnt = 0;
	int accp_cnt = 0;
	int saiter = 0;
	double T = init_temp;
	int L = 2 * Num_v * Saiter_factor;

	Help_len = 0;
	for (int i = 0; i < Num_v; i++)
	{
		if (sol[i] > 0)
		{
			Help_sol[Help_len] = i;
			Address_help_sol[i] = Help_len;
			Help_len++;
		}
	}
	int f_sol;
	build_delta_matrix(sol, f_sol);
	memcpy(Local_best_sol, sol, sizeof(int) * Num_v);
	int f_local_best = f_sol;

	while (1)
	{
		int i = rand() % Num_v;
		while (sol[i] == 0)
			i = rand() % Num_v;

		int j = rand() % Num_v;
		while (sol[j] == 1)
			j = rand() % Num_v;

		double delta = (Delta_matrix[0][j] + Delta_matrix[1][j] / 2) - (Delta_matrix[0][i] + Delta_matrix[1][i] / 2);

		//calculate num. of the open triangles related to j and i
		int dev = 0;
		for (int m = 0; m < Help_len; m++) //triple <j, i, vm>
		{
			int vm = Help_sol[m];
			if (Edge[j][i] && Edge[i][vm] && !Edge[vm][j] && vm != i)
				dev++;
		}
		for (int m = 0; m < Help_len; m++) //triple <j, vm, i>
		{
			int vm = Help_sol[m];
			if (Edge[j][vm] && Edge[vm][i] && !Edge[i][j] && vm != i)
				dev++;
		}
		for (int m = 0; m < Help_len; m++) //triple <i, j, vm>
		{
			int vm = Help_sol[m];
			if (Edge[i][j] && Edge[j][vm] && !Edge[vm][i] && vm != i)
				dev++;
		}
		delta -= dev;  //the complexity of calculating move gain is bounded by O(K_fixed)

		if (delta < 0)
		{
			update_delta_matrix(i, j);
			sol[i] = 0;
			sol[j] = 1;
			f_sol += delta;
			int pos = Address_help_sol[i];
			Help_sol[pos] = j;
			Address_help_sol[j] = pos;
			accp_cnt++;
		}
		else
		{
			double prob = exp(-(double) delta / T);
			if (rand() % 1000 < prob * 1000)
			{
				update_delta_matrix(i, j);
				sol[i] = 0;
				sol[j] = 1;
				f_sol += delta;
				int pos = Address_help_sol[i];
				Help_sol[pos] = j;
				Address_help_sol[j] = pos;
				accp_cnt++;
			}
		}
//		printf("frozen_cnt: %d, saiter: %d, K_fixed: %d, f_sol: %d\n",
//				frozen_cnt, saiter, K_fixed, f_sol);

		// update
		if (f_sol < f_local_best)
		{
			memcpy(Local_best_sol, sol, sizeof(int) * Num_v);
			f_local_best = f_sol;
		}
		if (f_sol == 0)
		{
			memcpy(sol, Local_best_sol, sizeof(int) * Num_v);  //return the best local optimum
			return f_sol;
		}

		saiter++;
		if (saiter > L)
		{
//			printf("frozen_cnt: %d, accp_cnt/L: %.6f\n", frozen_cnt, (double) accp_cnt / L);
			saiter = 0;
			T = T * Cool_ratio;  // cool down
			if (accp_cnt * 100 < L)  // = (double) accp_cnt / L
				frozen_cnt++;
			if (frozen_cnt > 5)
				break;

			accp_cnt = 0;
			saiter = 0;
		}

		if ((double) (clock() - Start_time) / CLOCKS_PER_SEC > Time_limit)
			break;
	}

	memcpy(sol, Local_best_sol, sizeof(int) * Num_v);  //return the best local optimum
	f_sol = f_local_best;
	return f_sol;
}


//restart simulated annealing algorithm (RSA)
void local_search_rsa()
{
	//obtained initial k value
	K_fixed = initial_k();  //greedy create a clique, return its size as the initial value of K_fixed
	printf(", initial K_fixed: %d, Run_time: %.2f\n", K_fixed, Run_time);
	if (K_fixed == K_opt)  //find the optimal solution
		return;
	K_fixed++;

	//initial solution is NOT necessary a clique (but with the size of K_fixed)
	random_initial_func(Cur_sol);

	while ((double) (clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	{
		int ff_des = annealing_search(Init_temp, Cur_sol);

		//the number of open triangles = 0
		//stop and store: Best_sol_one_run[] and K_fixed,
		//then K_fixed++
		if (ff_des == 0)
		{
			Best_k = K_fixed;
			memcpy(Best_sol_one_run, Cur_sol, sizeof(int) * Num_v);
			Run_time = (double) (clock() - Start_time) / CLOCKS_PER_SEC;
			if (Best_k == K_opt)
				break;
			K_fixed++;
		}

		//restart
		random_initial_func(Cur_sol);

		printf("next K_fixed: %d, Best_k: %d, Run_time: %.2f, ff_des: %d\n",
				K_fixed, Best_k, Run_time, ff_des);
	}

	//update global best
	if (Best_k > G_best_k)
	{
		G_best_k = Best_k;
		memcpy(Global_best_sol, Best_sol_one_run, sizeof(int) * Num_v);
	}
}
