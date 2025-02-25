![screenshot](long_cylindrical_pore.png)

# CGDFT_3D: 3D Coarse-Grained Lattice Gas Density Functional Theory

This is a the 3D version. Parallelized with MPI. To run, ensure that proper number of processors are assigned in the input file. The visualization shows gas-to-liquid phase transition in a long cylindrical pore. 

## Building 

```
mpic++ -std=c++0x code.cpp -o code
mpiexec -np NumProc code 
```



