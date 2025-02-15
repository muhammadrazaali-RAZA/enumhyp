# enumhyp
Source code used for the performance measurements in [Bläsius, Thomas; Friedrich, Tobias; Lischeid, Julius; Meeks, Kitty; Schirneck, Martin. *Efficiently Enumerating Hitting Sets of Hypergraphs Arising in Data Profiling.* Algorithm Engineering and Experiments (ALENEX) 2019](https://hpi.de/friedrich/research/enumdat). This tool generates unique column combination (UCC) hypergraphs from CSVs and enumerates all minimal hitting sets of hypergraphs using the enumeration algorithm described in the paper above. Use commit `8e6dfdba80a16d9d5eb1060561d9dbd747d564ef` to reproduce results.

**New in this fork**: A parallel MPI implementation for distributed enumeration of minimal hitting sets. See the [`parallel-enumhyp`](https://github.com/goodefroi/enumhyp/tree/parallel-enumhyp) branch for details.

## Documentation
The detailed documentation is attached. It discusses the parallelization of hitting set enumeration algorithms for data profiling, specifically for discovering minimal unique column combinations in databases. The report explores the efficiency of a parallel approach using MPI and compares the performance against a sequential method. The methodology involves a master-slave architecture to distribute the workload efficiently, and experimental results show significant speedups on complex datasets.

## Dependencies
- **Parallel MPI version**:  
  [OpenMPI](https://www.open-mpi.org/) or equivalent MPI library.
  
  [Boost](https://www.boost.org/).

---

## Build using CMake

### Parallel Version (MPI)
Navigate to the `bin` subdirectory and use MPI compilers:
```bash
cd enumhyp/bin
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc ..
make
```

### Enumerating Hitting Sets

#### Parallel Execution (MPI)
Use `mpirun` with `-np` to specify the number of nodes (e.g., 16 nodes):
```bash
GRAPH_FILE="table.graph"
OUTPUT_FILE="results_parallel.out"
NODES=16

mpirun -np $NODES ./enumhyp enumerate -I legacy "$GRAPH_FILE" -o "$OUTPUT_FILE"
```

## Hypergraph Files
Graphs are saved as plain text files:
1. The **first line** contains the number of vertices.
2. **Subsequent lines** list edges as comma-separated vertex indices.

Example:
```
5          # 5 vertices
0,1,2      # Edge 1: vertices 0,1,2
0,1,2,3,4  # Edge 2: vertices 0,1,2,3,4
3,4        # Edge 3: vertices 3,4
1,2,3      # Edge 4: vertices 1,2,3
```

---

## Notes for Parallel Execution
- **Cluster Compatibility**: Tested on the ALMA cluster (16 nodes). Adjust `NODES` to match your infrastructure.
- **Overhead**: Smaller datasets may see slower performance due to MPI communication costs.
- **Reproducibility**: The original sequential results use commit `8e6dfdba`. The MPI version is available in the [`parallel-enumhyp`](https://github.com/goodefroi/enumhyp/tree/parallel-enumhyp) branch.

For questions or issues with the parallel implementation, contact [Muhammad Raza Ali](mailto:muhammadrazaali.raza@gmail.com).
