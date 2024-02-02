# LDTPS
Source code for the article "Learning driven three-phase search for the maximum independent union of cliques"

1. For the source code, you can specify the parameters to match your needs when you execute the code. The code was tested on a computer under the Linux operating system. If you have any questions feel free to contact me (Zhi Lu: zhilusix@gmail.com).     

2. The proposed LDTPS algorithm was used to solve both the maximum independent union of cliques (IUC) problem and maximum multi-partite clique (MPC) problem. Specifically, for the maximum independent union of cliques (IUC) problem,
   ```
   g++ ./LDTPS/src/main.cpp ./LDTPS/src/common_func_def.cpp ./LDTPS/src/local_search.cpp -o ./LDTPS/IUC -O3
   ```
   ```
   ./LDTPS/IUC ./Instances/Test_Set_I/brock200_2.clq 15
   ```
   Among them,  
   ```
   IUC                  //binary code for the maximum IUC problem
   ./Instances/Test_Set_I/brock200_2.clq //input instance file brock200_2.clq
   15                   //the proven maximum IUC number for instance brock200_2.clq
   ```
and for the equivalent maximum MPC problem was solved by computing the IUC number of the complement graph $\overline{G}$. One just need uncomment the lines 56-73 in the `./LDTPS/src/common_func_def.cpp` 
```
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
//	Num_e = (double) Num_v * (Num_v - 1) / 2.0 - Num_e;
#endif
```
 
2. Please make sure that the following paper is cited if you use the code in your research.    
   Lu, Z., Gao, J., Hao, J. K., Yang, P., & Zhou, L. (2024). Learning driven three-phase search for the maximum independent union of cliques problem. Computers & Operations Research, 106549.

3. The source code is distributed for academic purposes only.    
   If you wish to use it for commercial applications, please contact the authors (Zhi Lu: zhilusix@gmail.com, Jin-Kao Hao: jin-kao.hao@univ-angers.fr).

