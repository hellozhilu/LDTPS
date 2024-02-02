#ifndef _GLOBAL_VARIABLES_H
#define _GLOBAL_VARIABLES_H

#include <iostream>
#include <fstream>
#include <cstdlib>
//#include <bitset>
#include <string.h>
#include <math.h>
#include <time.h>

using namespace std;

#define MAXCAN 1024
#define MAXNUM 99999999
#define PRECISION 1.0e-6
#define TRUE 1
#define FALSE 0
#define TRISTATE -1
//#define MAXLENGTH (4000000000)

extern char *Instance_name;		//instance name
extern int **Edge;				//adjacent matrix between each pair of vertices
//extern int **Adjacent;          //adjacent list
extern int *Degree;             //degree of each vertex
extern int Num_v;				//number of vertices
extern int Num_e;               //number of edges
extern double Density;          //density of graph
extern int K_fixed;             //size of IUC
extern int K_opt;               //optimal size of IUC (available for some instances)

//extern std::bitset<MAXLENGTH> Triple_set;	//triple set <i, j, k>, <i, j ,k>=1 if i connects j and j connects k
extern int **Delta_matrix;		//incremental matrix, Delta_matrix[0][i] records the number of open triangles that vertex i is in one end, i.e, <i, j, k> or <k, j, i> or <i, k, j> or <j, k, i>
								//Delta_matrix[1][j] records the number of open triangles that vertex j is in the middle, i.e, <i, j, k> or <k, j, i>
//extern int **Delta_matrix_tmp;
extern int **Tabu_list;
extern double **Prob_matrix;   //probability matrix in Reinforcement learning, size: V * 2
extern int *Help_sol;			//contains the vertices with value of 1
extern int *Address_help_sol;
extern int *Freq;
extern int *Compo_sol;
extern int Compo_cnt;
extern int Compo_len;

extern int *Temp_sol;
extern int *Cur_sol;
extern int *Local_best_sol;
extern int *Best_sol_one_run;
extern int *Global_best_sol;
extern int Best_k;				//corresponds to the Best_sol_one_run, storing the best value of k per run
extern int G_best_k;            //corresponds to the Global_best_sol, storing the best value of k over 20 runs

extern int Runs;
extern double Time_limit, Start_time, Run_time;

//6 parameters
//for neighborhood (1 parameter)
extern double Cn_rho;		   //para: parametric constrained neighborhood, i.e., Cn_rho^2 * k * (n - k)
//for tabu search (2 parameters)
extern int Tabu_depth;         //para: maximum non-improvement iterations of TS
extern int Tabu_tenure;        //para: tabu tenure
//for reinforcement learning (3 parameters)
extern double Alpha;           //para: reward factor for correct subset
extern double Beta;            //para: penalization factor for incorrect subset
extern double Gamma;           //para: compensation factor for expected subset

#endif
