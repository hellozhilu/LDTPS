//============================================================================
// Name        : maximum_iuc_ldtps.cpp
// Author      : Zhi Lu
// Version     : 30 April, 2022
// Copyright   : Business School, University of Shanghai for Science & Technology (USST), Shanghai, China
// Description : Learning driven three-phase search for the maximum independent union of cliques problem
//============================================================================

#include "func_state.h"
#include "global_variables.h"


char *Instance_name;	//instance name
int **Edge;				//adjacent matrix between each pair of vertices
int *Degree;            //degree of each vertex
int Num_v;				//number of vertices
int Num_e;              //number of edges
double Density;         //density of graph
int K_fixed;            //size of IUC
int K_opt;              //optimal size of IUC (available for some instances)

//std::bitset<MAXLENGTH> Triple_set;	//triple set <i, j, k>, <i, j ,k>=1 if i connects j and j connects k
int **Delta_matrix;		//incremental matrix, Delta_matrix[0][i] records the number of open triangles that vertex i is in one end, i.e, <i, j, k> or <k, j, i> or <i, k, j> or <j, k, i>
						//Delta_matrix[1][j] records the number of open triangles that vertex j is in the middle, i.e, <i, j, k> or <k, j, i>
//int **Delta_matrix_tmp;
int **Tabu_list;
double **Prob_matrix;   //probability matrix in Reinforcement learning, size: V * 2
int *Help_sol;			//contains the vertices with value of 1
int *Address_help_sol;
int *Freq;
int *Compo_sol;
int Compo_cnt;
int Compo_len;

int *Temp_sol;
int *Cur_sol;
int *Local_best_sol;
int *Best_sol_one_run;
int *Global_best_sol;
int Best_k;				//corresponds to the Best_sol_one_run, storing the best value of k per run
int G_best_k;           //corresponds to the Global_best_sol, storing the best value of k over 20 runs

int Runs;
double Time_limit, Start_time, Run_time;

//6 parameters
//for neighborhood (1 parameter)
double Cn_rho;			//para: parametric constrained neighborhood, Cn_rho^2 * k * (n - k)
//for tabu search (2 parameters)
int Tabu_depth;         //para: maximum non-improvement iterations of TS
int Tabu_tenure;        //para: tabu tenure
//for reinforcement learning (3 parameters)
double Alpha;           //para: reward factor for correct subset
double Beta;            //para: penalization factor for incorrect subset
double Gamma;           //para: compensation factor for expected subset


int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		printf("show usage: maximum_iuc.exe input_file optimal_k");
		exit(-1);
	}
	Instance_name = argv[1];
	K_opt = atoi(argv[2]);

	srand(unsigned(time(NULL)));

	//parameters (6)
	Time_limit = 3600.0;
	Runs = 20;
	Cn_rho = 0.1;
	Tabu_depth = 5000;
	Tabu_tenure = 85;
	Alpha = 0.4;
	Beta = 0.1;
	Gamma = 0.5;

	printf("PARAMETERS: \n");
	printf("Exe_name: %s, Instance_name: %s, K_opt: %d\n", argv[0], Instance_name, K_opt);
	printf("Time_limit: %.2f(s), Runs: %d\n", Time_limit, Runs);
	printf("LS: Cn_rho: %.2f, Tabu_depth: %d, Tabu_tenure: %d\n", Cn_rho, Tabu_depth, Tabu_tenure);
	printf("RL: Alpha: %.2f, Beta: %.2f, Gamma: %.2f\n\n", Alpha, Beta, Gamma);

	read_instance();
	allocate_memory();

	G_best_k = -MAXNUM;
	int ff_total[Runs];
	double time_total[Runs];
	for (int i = 0; i < Runs; i++)
	{
		Start_time = clock();

		printf("Runs: %d", i);
		local_search_ldtps();
//		calculate_connect_component(Best_sol_one_run);  //TODO
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
