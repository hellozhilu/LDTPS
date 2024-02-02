#ifndef _FUNC_STATE_H
#define _FUNC_STATE_H


void read_instance();
void allocate_memory();
void release_memory();
void out_results_one_run_tmp(char *instance_name, int ff);
void out_results_one_run(char *instance_name, int runs, int ff, double run_time, int *sol);
void out_total_results(char *instance_name, int *cost_total, double *time_total);
void proof(int *sol, int ff);

int calculate_sol(int *sol);
void build_delta_matrix(int *sol, int &ff);
//void build_delta_matrix_tmp(int *sol, int help_len, int &ff);
//void compare_delta(int f1, int f2);
void update_delta_matrix(int u, int v);
int initial_k();
void random_initial_func(int *sol);
int annealing_search(double init_temp, int *sol);
void local_search_rsa();

void print_current_population(int pop_cnt);
void append_to_population(int pop_cnt, int ff_des, int *sol);
int check_if_different(int idx, int *sol);
int build_initial_population(int *sol);
int uniform_crossover(int *sol, int *parent1, int *parent2);
int mutation();
int evaluate_offspring(int pop_size, int f_sol, int *sol);
void genetic_algorithm();


#endif
