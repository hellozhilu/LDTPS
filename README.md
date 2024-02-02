# LDTPS
Source code for the article "Learning driven three-phase search for the maximum independent union of cliques"

1. For the source code, you can specify the parameters to match your needs when you execute the code. The code was tested on a computer under Linux operating system. If you have any questions feel free to contact me (Zhi Lu: zhilusix@gmail.com).     
   
   ```
   g++ ./LDTPS/main.cpp ./LDTPS/common_func_def.cpp ./LDTPS/local_search.cpp -o ./LDTPS/IUC -O3
   ```
   ```
   ./LDTPS/IUC.exe ./Instances/Test_Set_I/brock200_2.clq 15
   ```
   Among them,  
   ```
   IUC.exe              //binary code for the maximum IUC problem
   ./Instances/Test_Set_I/brock200_2.clq //input instance file
   15                   //the proven maximum IUC number for instance brock200_2.clq
   ```
  
2. Please make sure that the following paper is cited if you use the code in your research.    
   Lu, Z., Gao, J., Hao, J. K., Yang, P., & Zhou, L. (2024). Learning driven three-phase search for the maximum independent union of cliques problem. Computers & Operations Research, 106549.

3. The source code is distributed for academic purposes only.    
   If you wish to use it for commercial applications, please contact the authors (Zhi Lu: zhilusix@gmail.com).  

