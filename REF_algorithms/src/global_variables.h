#ifndef _GLOBAL_VARIABLES_H
#define _GLOBAL_VARIABLES_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string.h>
//#include <bitset>
#include <math.h>

using namespace std;

#define MAXTRY 500
#define MAXCAN 1024
#define MAXNUM 99999999
#define PRECISION 1.0e-6
//#define MAXLENGTH (4000000000)

extern char *Instance_name;		//instance name
extern int Num_v;				//number of vertices
extern int Num_e;				//number of edges
extern double Density;          //density of graph
extern int K_fixed;             //size of IUC
extern int K_opt;               //optimal size of IUC (available for some instances)

//extern std::bitset<MAXLENGTH> Triple_set;	//triple set <i, j, k>, <i, j ,k>=1 if i connects j and j connects k
extern int **Edge;				//adjacent matrix between each pair of vertices
extern int *Degree;             //degree of each vertex
extern int **Delta_matrix;		//incremental matrix, Delta_matrix[0][i] records the number of open triangles that vertex i is in one end, i.e, <i, j, k> or <k, j, i> or <i, k, j> or <j, k, i>
								//Delta_matrix[1][j] records the number of open triangles that vertex j is in the middle, i.e, <i, j, k> or <k, j, i>
//extern int **Delta_matrix_tmp;
extern int *Help_sol;			//contains the vertices with value of 1
extern int *Address_help_sol;
extern int Help_len;
extern int **Similarity;        //similarity matrix for each pair of solution in population
extern int **Population;        //solutions of population
extern int *Pool_cost;          //fitness function value of population

extern int *Cur_sol;
extern int *Best_sol;
extern int *Local_best_sol;
extern int *Best_sol_one_run;
extern int *Global_best_sol;
extern int f_best_sol;
extern int Best_k;				//corresponds to the Best_sol solution, preserves the best value of k
extern int G_best_k;            //corresponds to the Global_best_sol, storing the best value of k over 20 runs
extern int Worst_pool_idx;

extern int Runs;
extern double Time_limit, Start_time, Run_time;

//parameters
extern double Init_temp;        //SA initial temperature
extern int Saiter_factor;       //SA iteration factor
extern double Cool_ratio;       //SA cooling ratio
extern int Pool_size;           //GA population size
extern double Mutation_ratio;   //GA mutation ratio


#endif
