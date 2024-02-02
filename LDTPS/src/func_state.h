#ifndef _FUNC_STATE_H
#define _FUNC_STATE_H


void read_instance();
void allocate_memory();
void release_memory();
void out_results_one_run_tmp(char *instance_name, int ff);
void out_results_one_run(char *instance_name, int runs, int ff, double run_time, int *sol);
void out_total_results(char *instance_name, int *cost_total, double *time_total);
void proof(int *sol, int ff);

void dfs(int v, int *sol, int *visited);
void calculate_connect_component(int *sol);
void build_prob_matrix();
void update_prob_matrix(int *tmp_sol, int *cur_sol);
//void smooth_prob_matrix_tmp()
//void smooth_prob_matrix();
void build_delta_matrix(int *sol, int help_len, int &ff);
//void build_delta_matrix_tmp(int *sol, int help_len, int &ff);
//void compare_delta(int f1, int f2);
void update_delta_matrix(int u, int v, int help_len);
int initial_k();
void random_selection(int *sol);
void greedy_selection(int *sol);
void hybrid_initial_sol(int *sol);
int descent_local_search(int *sol);
int compare_ascending(const void *a, const void *b);
int tabu_search(int *sol);
void freq_driven_perturabtion(int *sol);
void fresh_freq();
//int greedy_expansion(int *sol);
void local_search_ldtps();


#endif
