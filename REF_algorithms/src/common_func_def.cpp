#include"func_state.h"
#include"global_variables.h"


//read instance
void read_instance()
{
	ifstream FIC;
	FIC.open(Instance_name);
	if (FIC.fail())
	{
		printf("### Error open, Instance_name: %s ###\n", Instance_name);
		exit(-1);
	}
	char StrReading[MAXCAN];
	FIC >> StrReading;
	int max_edg = 0;
	while (!FIC.eof())
	{
		char bidon[100];
		if (strcmp(StrReading, "p") == 0)
		{
			FIC >> bidon >> Num_v >> Num_e;
			Density = (double) 2 * Num_e / (Num_v * (Num_v - 1));

			Edge = new int*[Num_v];
			for (int i = 0; i < Num_v; i++)
				Edge[i] = new int[Num_v];
			for (int i = 0; i < Num_v; i++)
				memset(Edge[i], 0, sizeof(int) * Num_v);

		}
		if (strcmp(StrReading, "e") == 0)
		{
			int x1, x2;
			FIC >> x1 >> x2;
			x1--; x2--;
			if (x1 < 0 || x2 < 0 || x1 >= Num_v || x2 >= Num_v)
			{
				printf("### Error of node x1: %d, x2: %d ###\n", x1, x2);
				exit(-1);
			}
			Edge[x1][x2] = Edge[x2][x1] = 1;				
			max_edg++;
		}
		FIC >> StrReading;
	}
	if (max_edg != Num_e)
	{
		printf("### Error max_edge != nb_edge, Num_e: %d, max_edge: %d ###\n", Num_e, max_edg);
		exit(-1);
	}
	FIC.close();

#ifdef DEBUG
	//TODO MPC: The MPC number was found by computing the IUC number of the complement graph
	printf("running MPC\n");
	for (int i = 0; i < Num_v; i++)
	{
		for (int j = 0; j < Num_v; j++)
		{
			if (Edge[i][j] > 0)
				Edge[i][j] = 0;
			else
				Edge[i][j] = 1;

			if (i == j)
				Edge[i][j] = 0;
		}
	}
//	nb_edg = (double) Num_v * (Num_v - 1) / 2.0 - nb_edg;
#endif

	printf("running IUC\n");
	printf("Instance_name: %s, Num_v: %d, Num_e: %d, Density: %.2f\n\n",
			Instance_name, Num_v, Num_e, Density);

#ifdef DEBUG
	printf("Num_v: %d, Num_e: %d\n", Num_v, Num_e);
	for (int i = 0; i < Num_v; i++)
		for (int j = 0; j < Num_v; j++)
			if (Edge[i][j] > 0)
				printf("e %d %d\n", i + 1, j + 1);

	printf("Degree[]: [")
	for (int i = 0; i < Num_v; i++)
		printf("%d ", Degree[i]);
	printf("]\n");
#endif
}


void allocate_memory()
{
	Delta_matrix = new int*[2];
//	Delta_matrix_tmp = new int*[2];
	for (int i = 0; i < 2; i++)
	{
		Delta_matrix[i] = new int[Num_v];
//		Delta_matrix_tmp[i] = new int[Num_v];
	}

	Cur_sol = new int[Num_v];
	Best_sol = new int[Num_v];
	Local_best_sol = new int[Num_v];
	Global_best_sol = new int[Num_v];
	Best_sol_one_run = new int[Num_v];
	Help_sol = new int[Num_v];
	Address_help_sol = new int[Num_v];

	Similarity = new int*[2];
	for (int i = 0; i < 2; i++)
		Similarity[i] = new int[2];

	Population = new int*[Pool_size + 1];
	for (int i = 0; i <= Pool_size; i++)
		Population[i] = new int[Num_v];
	Pool_cost = new int[Pool_size + 1];
}


void release_memory()
{
	for (int i = 0; i < Num_v; i++)
	{
		delete[] Edge[i]; Edge[i] = NULL;
	}	
	delete[] Edge; Edge = NULL;
	delete[] Degree; Degree = NULL;

	for (int i = 0; i < 2; i++)
	{
		delete[] Delta_matrix[i]; Delta_matrix[i] = NULL;
//		delete[] Delta_matrix_tmp[i]; Delta_matrix_tmp[i] = NULL;
		delete[] Similarity[i]; Similarity[i] = NULL;
	}
	delete[] Delta_matrix; Delta_matrix = NULL;
//	delete[] Delta_matrix_tmp; Delta_matrix_tmp = NULL;
	delete[] Similarity; Similarity = NULL;

	delete[] Cur_sol; Cur_sol = NULL;
	delete[] Best_sol; Best_sol = NULL;
	delete[] Local_best_sol; Local_best_sol = NULL;
	delete[] Best_sol_one_run; Best_sol_one_run = NULL;
	delete[] Global_best_sol; Global_best_sol = NULL;
	delete[] Help_sol; Help_sol = NULL;
	delete[] Address_help_sol; Address_help_sol = NULL;

	for (int i = 0; i < Pool_size + 1; i++)
	{
		delete[] Population[i]; Population[i] = NULL;
	}
	delete[] Population; Population = NULL;
	delete[] Pool_cost; Pool_cost = NULL;

	printf("Success release all!!!\n");
}


void out_results_one_run_tmp(char *instance_name, int ff)
{
	FILE *fp;
	char buff[MAXCAN];
	char *graph_name = basename(instance_name);

	sprintf(buff, "./REF_algorithms/output_dir_temp/%s_rec_temp.txt", graph_name);
	fp = fopen(buff, "a+");
	if (fp == NULL)
		exit(-1);

	fprintf(fp, "%d ", ff);
	fclose(fp);
}


void out_results_one_run(char *instance_name, int runs, int ff, double run_time, int *sol)
{
	FILE *fp;
	char buff[MAXCAN];
	char *graph_name = basename(instance_name);

	sprintf(buff, "./REF_algorithms/output_dir/%s_rec.txt", graph_name);
	fp = fopen(buff, "a+");
	if (fp == NULL)
		exit(-1);

	fprintf(fp, "%d %s %d %d %f\n", runs, graph_name, K_opt, ff, run_time);
	for (int i = 0; i < Num_v; i++)
		if (sol[i])
			fprintf(fp, "%d ", i + 1);
	fprintf(fp, "\n");
	fclose(fp);
}


void out_total_results(char *instance_name, int *cost_total, double *time_total)
{
	int ff_best = -MAXNUM;
	int ff_worst = MAXNUM;
	int hit = 0;
	double ff_avg = 0.0;
	double avg_time = 0.0;
	double std = 0.0;

	for (int i = 0; i < Runs; i++)
	{
		ff_avg += cost_total[i];
		avg_time += time_total[i];
	}
	ff_avg /= Runs;
	avg_time /= Runs;
	for (int i = 0; i < Runs; i++)
	{
		if (cost_total[i] > ff_best)
			ff_best = cost_total[i];
		if (cost_total[i] < ff_worst)
			ff_worst = cost_total[i];
	}
	for (int i = 0; i < Runs; i++)
	{
		if (cost_total[i] == ff_best)
			hit++;
	}
	for (int i = 0; i < Runs; i++)
		std += pow((cost_total[i] - ff_avg), 2);
	std /= Runs;
	std = sqrt(std);

	//output
	FILE *fp;
	char buff[MAXCAN];
	sprintf(buff, "./REF_algorithms/Sol_stat_total.txt");
	fp = fopen(buff, "a+");
	if (fp == NULL)
		exit(-1);

	fprintf(fp, "%s %d %d %.2f %d %d %f %d %d %f %f\n",
			instance_name, Num_v, Num_e, Density, K_opt, ff_best, ff_avg, ff_worst, hit, avg_time, std);
	fclose(fp);
}


void proof(int *sol, int ff)
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
					}
				}
			}
		}
	}
	if (f_obj > 0)
		printf("in proof func, an error is detected, f_obj != 0, && f_obj = %d\n", f_obj);
	if (len != ff)
		printf("in proof func, an error is detected, len != ff, && len = %d, ff = %d\n", len, ff);
}

