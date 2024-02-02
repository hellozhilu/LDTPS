/*
 * genetic_algorithm.cpp
 *
 *  Created on: Jul 13, 2023
 *      Author: zhi
 */
#include"func_state.h"
#include"global_variables.h"


void print_current_population(int pop_size)
{
	int obj_best = MAXNUM;
	int obj_worst = -MAXNUM;
	double obj_avg = 0;
	for (int i = 1; i <= pop_size; i++)
	{
		if (Pool_cost[i] < obj_best)
			obj_best = Pool_cost[i];
		if (Pool_cost[i] > obj_worst)
			obj_worst = Pool_cost[i];
		obj_avg += Pool_cost[i];
	}
	obj_avg /= (double) pop_size;

	printf("f_best_sol: %d, f_worst: %d, f_avg: %.2f\n", obj_best, obj_worst, obj_avg);
	printf("Population = [");
	for (int i = 1; i <= pop_size; i++)
		printf("%d, ", Pool_cost[i]);
	printf("]\n\n");
}


void append_to_population(int pop_cnt, int ff_des, int *sol)
{
	memcpy(Population[pop_cnt], sol, sizeof(int) * Num_v);
	Pool_cost[pop_cnt] = ff_des;

	if (pop_cnt == 1 || ff_des > Pool_cost[0])
	{
		memcpy(Population[0], sol, sizeof(int) * Num_v);
		Pool_cost[0] = ff_des;
		Worst_pool_idx = pop_cnt;
	}
}


// count == partition number means two solution are identical
// return 1: two solutions are different
// return 0: two solutions are identical
int check_if_different(int idx, int *sol)
{
	int count = 0;
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			Similarity[i][j] = 0;

	for (int i = 0; i < Num_v; i++)
	{
		int x = sol[i];
		int y = Population[idx][i];
		Similarity[x][y]++;
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (Similarity[i][j] > 0)
				count++;
			if (count > 2)  //count > 2 means two solutions are different
				return 1;
		}
	}
	return 0;
}


int build_initial_population(int *sol)
{
	int pop_cnt = 0;
	int ntrials = 0;
	int append, differ;

	while (pop_cnt < Pool_size)
	{
		random_initial_func(sol);
		int ff_sol = calculate_sol(sol);

		if (ff_sol == 0)
		{
			Best_k = K_fixed;
			memcpy(Best_sol_one_run, sol, sizeof(int) * Num_v);
			Run_time = (double) (clock() - Start_time) / CLOCKS_PER_SEC;
			if (Best_k == K_opt)
				break;
			K_fixed++;

			pop_cnt = 0;
			continue;
		}

		ntrials++;

		if (pop_cnt == 0)
			append = 1;
		else
		{
			differ = 1;
			for (int i = 1; i <= pop_cnt; i++)
			{
				differ = check_if_different(i, sol);
				if (differ == 0)
					break;
			}
			if (differ > 0)
				append = 1;
			else
				append = 0;
		}

		// add solution to population
		if (append > 0)
		{
			pop_cnt++;
			append_to_population(pop_cnt, ff_sol, sol);
		}
		if (ntrials >= MAXTRY)
			break;
	}

	return pop_cnt;
}


int uniform_crossover(int *sol, int *parent1, int *parent2)
{
	//only 2 parents are supported
	int len = 0;
	for (int i = 0; i < Num_v; i++)
	{
		if (rand() % 100 < 50)
			sol[i] = parent1[i];
		else
			sol[i] = parent2[i];
	}
	for (int i = 0; i < Num_v; i++)
		if (sol[i] > 0)
			len++;

	//repair procedure
	// if len < K_fixed, randomly add vertices to D
	while (len < K_fixed)
	{
		int rx = rand() % Num_v;
		while (sol[rx] == 1)
			rx = rand() % Num_v;
		sol[rx] = 1;
		len++;
	}
	// if len > K_fixed, randomly remove vertices from D
	while (len > K_fixed)
	{
		int rx = rand() % Num_v;
		while (sol[rx] == 0)
			rx = rand() % Num_v;
		sol[rx] = 0;
		len--;
	}

	int f_sol = calculate_sol(sol);
	return f_sol;
}


int mutation(int *sol)
{
	int len = K_fixed;
	for (int i = 0; i < Num_v; i++)
	{
		if (rand() % 100 < Mutation_ratio * 100)
		{
			if (sol[i] > 0)
			{
				sol[i] = 0;
				len--;
			}
			else
			{
				sol[i] = 1;
				len++;
			}
		}
	}

	//repair procedure
	// if len < K_fixed, randomly add vertices to D
	while (len < K_fixed)
	{
		int rx = rand() % Num_v;
		while (sol[rx] == 1)
			rx = rand() % Num_v;
		sol[rx] = 1;
		len++;
	}
	// if len > K_fixed, randomly remove vertices from D
	while (len > K_fixed)
	{
		int rx = rand() % Num_v;
		while (sol[rx] == 0)
			rx = rand() % Num_v;
		sol[rx] = 0;
		len--;
	}

	int f_sol = calculate_sol(sol);
	return f_sol;
}


int evaluate_offspring(int pop_size, int f_sol, int *sol)
{
	int idx;
	int include, differ, replaced;
	double max_cost = -MAXNUM;

	include = 0;
	differ = 1;
	replaced = 0;
	for (int i = 1; i <= Pool_size; i++)
	{
		differ = check_if_different(i, sol);
		if (differ == 0)  //differ = 0 means two solutions are identical
			return pop_size;
	}
	if (differ > 0)  //differ = 1 means two solutions are different
	{
		if (pop_size < Pool_size || f_sol <= Pool_cost[0])
			include = 1;
		else
			include = 0;
	}

	if (include == 0)
		return pop_size;

	if (pop_size < Pool_size)  //insert offspring into the population
	{
		pop_size++;
		idx = pop_size;
	}
	else  //replace worst individual of the population
	{
		idx = Worst_pool_idx;
		replaced = 1;
	}
	memcpy(Population[idx], sol, sizeof(int) * Num_v);  //include offspring into population
	Pool_cost[idx] = f_sol;

	//recalculate the worst solution of the NEW population
	if (replaced > 0)
	{
		max_cost = Pool_cost[1];
		idx = 1;
		for (int i = 2; i <= pop_size; i++)
		{
			if (Pool_cost[i] > max_cost)
			{
				max_cost = Pool_cost[i];
				idx = i;
			}
		}
		memcpy(Population[0], Population[idx], sizeof(int) * Num_v);
		Pool_cost[0] = max_cost;
		Worst_pool_idx = idx;
	}

	return pop_size;
}


//genetic algorithm (GA)
void genetic_algorithm()
{
	int generation = 1;
	memset(Best_sol, 0, sizeof(int) * Num_v);
	Pool_cost[0] = -MAXNUM;

	//obtained initial k value
	K_fixed = initial_k();  //greedy create a clique, return its size as the initial value of K_fixed
	printf(", initial K_fixed: %d, Run_time: %.2f\n", K_fixed, Run_time);
	if (K_fixed == K_opt)  //find the optimal solution
		return;
	K_fixed++;

	//build initial pool with K_fixed
	int pop_size = build_initial_population(Cur_sol);
	if (K_fixed == K_opt)  //find the optimal solution
		return;
	if (pop_size < 2)
	{
		printf("### Initial population is too small to proceed: Pool_size: %d ###\n", Pool_size);
		exit(-1);
	}

	while ((double) (clock() - Start_time) / CLOCKS_PER_SEC < Time_limit)
	{
		f_best_sol = MAXNUM;

		/* parent selection and crossover operator */
		int ff_cross = -1;
		int idx1 = rand() % Pool_size;
		int idx2 = rand() % Pool_size;
		while (idx1 == idx2)
			idx2 = rand() % Pool_size;
		ff_cross = uniform_crossover(Cur_sol, Population[idx1], Population[idx2]);
		if (ff_cross > 0 && ff_cross < f_best_sol)
		{
			f_best_sol = ff_cross;
			memcpy(Best_sol, Cur_sol, sizeof(int) * Num_v);
		}

		/* mutation operator */
		int ff_mutation = -1;
		if (ff_cross > 0)
			ff_mutation = mutation(Cur_sol);
		if (ff_mutation > 0 && ff_mutation < f_best_sol)
		{
			f_best_sol = ff_mutation;
			memcpy(Best_sol, Cur_sol, sizeof(int) * Num_v);
		}

		//the number of open triangles = 0
		//stop and store: Best_sol_one_run[] and K_fixed,
		//then K_fixed++
		if (ff_cross == 0 || ff_mutation == 0)
		{
			printf("In generation: %d, next K_fixed: %d, Best_k: %d, Run_time: %.2f\n", generation, K_fixed, Best_k, Run_time);

			Best_k = K_fixed;
			memcpy(Best_sol_one_run, Cur_sol, sizeof(int) * Num_v);
			Run_time = (double) (clock() - Start_time) / CLOCKS_PER_SEC;
			if (Best_k == K_opt)
				break;
			K_fixed++;
		}

		/* pool updating */
		//build new population according to current K_fixed value
		if (ff_cross == 0 || ff_mutation == 0)
		{
			pop_size = build_initial_population(Cur_sol);
			if (pop_size < 2)
			{
				printf("### Current population is too small to proceed: Pool_size: %d ###\n", Pool_size);
				exit(-1);
			}
		}
		else
			pop_size = evaluate_offspring(pop_size, f_best_sol, Cur_sol);  //update current population

//		print_current_population(pop_size);

		generation++;
	}

	//update global best
	if (Best_k > G_best_k)
	{
		G_best_k = Best_k;
		memcpy(Global_best_sol, Best_sol_one_run, sizeof(int) * Num_v);
	}
}
