# LDTPS
Source codes for the article "Learning driven three-phase search for the maximum independent union of cliques"

1. For all the source codes, you can specify the parameters to match your needs when you execute the codes. The codes were tested on a computer under the Linux operating system. If you have any questions feel free to contact me (Zhi Lu: zhilusix@gmail.com).     

2. The proposed LDTPS algorithm is used to solve both the maximum independent union of cliques (IUC) problem and the maximum multi-partite clique (MPC) problem.
   Specifically, the maximum IUC problem is tested as follows,
   ```bash
   g++ ./LDTPS/src/main.cpp ./LDTPS/src/common_func_def.cpp ./LDTPS/src/local_search.cpp -o ./LDTPS/IUC -O3
   ```
   ```bash
   ./LDTPS/IUC ./Instances/Test_Set_I/brock200_2.clq 15
   ```
   among them,  
   ```
   IUC                  //LDTPS binary code for the maximum IUC problem
   ./Instances/Test_Set_I/brock200_2.clq //input instance file brock200_2.clq
   15                   //the proven maximum IUC number for the instance brock200_2.clq
   ```

   The equivalent maximum MPC problem is then solved by computing the IUC number of the complement graph $\overline{G}$. You need to uncomment lines 57-71 in the `./LDTPS/src/common_func_def.cpp` (see below) and test it as before,
   ```C++
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
   ```

3. The source codes of the compared restart simulated annealing (RSA) algorithm and the genetic algorithm (GA) are also available in the directory `./REF_algorithms`. One of the RSA algorithm is tested as follows,
   ```bash
   g++ ./REF_algorithms/src/main.cpp ./REF_algorithms/src/common_func_def.cpp ./REF_algorithms/src/local_search.cpp ./REF_algorithms/src/genetic_algorithm.cpp -o ./REF_algorithms/RSA -O3
   ```
   ```bash
   ./REF_algorithms/RSA ./Instances/Test_Set_I/brock200_2.clq 15
   ```
   among them,  
   ```bash
   RSA                  //RSA binary code for the maximum IUC problem, and replace it with GA if needed
   ./Instances/Test_Set_I/brock200_2.clq //input instance file brock200_2.clq
   15                   //the proven maximum IUC number for the instance brock200_2.clq
   ```
 
5. Please make sure that the following paper is cited if you use the codes in your research.    
   Lu, Z., Gao, J., Hao, J. K., Yang, P., & Zhou, L. (2024). Learning driven three-phase search for the maximum independent union of cliques problem. Computers & Operations Research, 106549.

3. The source codes are distributed for academic purposes only.    
   If you wish to use it for commercial applications, please contact the authors (Zhi Lu: zhilusix@gmail.com, Jin-Kao Hao: jin-kao.hao@univ-angers.fr).

