//============================================================================
// Name        : maximum_iuc_refs.cpp
// Author      : Zhi Lu
// Version     : 30 April, 2022 / revised until October 2023
// Copyright   : Business School, University of Shanghai for Science & Technology (USST), Shanghai, China
// Description : The maximum independent union of cliques problem (IUC), or the complementory MPC problem
//               restart simulated annealing search (RSA): 3 paras
//               genetic algorithm (GA): 2 paras
//============================================================================

#include "func_state.h"
#include "global_variables.h"

char *Instance_name;	//instance name
int Num_v;				//number of vertices
int Num_e;              //number of edges
double Density;         //density of graph
int K_fixed;            //size of IUC
int K_opt;              //optimal size of IUC (available for some instances)

//std::bitset<MAXLENGTH> Triple_set;	//triple set <i, j, k>, <i, j ,k>=1 if i connects j and j connects k
int **Edge;				//adjacent matrix between each pair of vertices
int *Degree;            //degree of each vertex
int **Delta_matrix;		//incremental matrix, Delta_matrix[0][i] records the number of open triangles that vertex i is in one end, i.e, <i, j, k> or <k, j, i> or <i, k, j> or <j, k, i>
						//Delta_matrix[1][j] records the number of open triangles that vertex j is in the middle, i.e, <i, j, k> or <k, j, i>
//int **Delta_matrix_tmp;
int *Help_sol;			//contains the vertices with value of 1
int *Address_help_sol;
int Help_len;
int **Similarity;       //similarity matrix for each pair of solution in population
int **Population;       //solutions of population
int *Pool_cost;         //fitness function value of population

int *Cur_sol;
int *Best_sol;
int *Local_best_sol;
int *Best_sol_one_run;
int *Global_best_sol;
int f_best_sol;
int Best_k;				//corresponds to the Best_sol solution, i.e., Best_sol_one_run, preserves the best value of k
int G_best_k;           //corresponds to the Global_best_sol, storing the best value of k over 20 runs
int Worst_pool_idx;

int Runs;
double Time_limit, Start_time, Run_time;

//parameters
double Init_temp;       //SA initial temperature
int Saiter_factor;      //SA iteration factor
double Cool_ratio;      //SA cooling ratio
int Pool_size;          //GA population size
double Mutation_ratio;  //GA mutation ratio


int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		printf("show usage: maximum_iuc_refs.exe input_file optimal_k");
		exit(0);
	}
	Instance_name = argv[1];
	K_opt = atoi(argv[2]);

	//parameters
	Time_limit = 3600.0;
	Runs = 20;
	Init_temp = 1.0;
	Saiter_factor = 8;
	Cool_ratio = 0.96;
	Pool_size = 50;
	Mutation_ratio = 0.1;

	printf("PARAMETERS: \n");
	printf("Exe_name: %s, Instance_name: %s, K_opt: %d\n", argv[0], Instance_name, K_opt);
	printf("Time_limit: %.1f(s), Runs: %d\n", Time_limit, Runs);
	printf("RSA: Init_temp: %.2f, Saiter_factor: %d, Cool_ratio: %.2f\n",
			Init_temp, Saiter_factor, Cool_ratio);
	printf("GA: Pool_size: %d, Mutation_ratio: %.2f\n\n", Pool_size, Mutation_ratio);

	srand(unsigned(time(NULL)));

	read_instance();
	allocate_memory();

	G_best_k = -MAXNUM;
	int ff_total[Runs];
	double time_total[Runs];
	for (int i = 0; i < Runs; i++)
	{
		Start_time = clock();

		printf("Runs: %d", i);
		local_search_rsa();   //restart simulated annealing algorithm (RSA)
//		genetic_algorithm();  //genetic algorithm (GA)
		printf("Runs: %d, K_opt: %d, Best_k: %d, Run_time: %.2f\n\n", i, K_opt, Best_k, Run_time);

		out_results_one_run(Instance_name, i, Best_k, Run_time, Best_sol_one_run);  //output result per run
		out_results_one_run_tmp(Instance_name, Best_k);  //output solution per run
		proof(Best_sol_one_run, Best_k);
		ff_total[i] = Best_k;
		time_total[i] = Run_time;
	}

	printf("FINISHED: Instance_name: %s, K_opt: %d, G_best_k: %d\n\n", argv[1], K_opt, G_best_k);
	out_total_results(Instance_name, ff_total, time_total);  //output total result
	proof(Global_best_sol, G_best_k);
	release_memory();
	return 0;
}
