# LDTPS
Source code for the article "Learning driven three-phase search for the maximum independent union of cliques"

1. For the source code, you can specify the parameters to match your needs when you execute the code. The code was tested on a computer under Linux operating system. If you have any questions feel free to contact me (Zhi Lu: zhilusix@gmail.com).     
   
   ```
   g++ ./src/main.cpp ./src/common_func_def.cpp ./src/local_search.cpp -o ./IUC -O3
   ```
   ```
   IUC.exe ./Instances/Test_Set_I/brock200_1.clq 21
   ```
   Among them,  
   ```
   IUC.exe              //binary code
   ./Instances/Test_Set_I/brock200_1.clq //input instance file
   21                   //best-known result for current instance (brock200_1.clq) for the maximum IUC problem
   ```
  
2. Please make sure that the following paper is cited if you use the code in your research.    
   Lu, Z., Gao, J., Hao, J. K., Yang, P., & Zhou, L. (2024). Learning driven three-phase search for the maximum independent union of cliques problem. Computers & Operations Research, 106549.

3. The source code is distributed for academic purposes only.    
   If you wish to use it for commercial applications, please contact the authors (Zhi Lu: zhilusix@gmail.com).  

