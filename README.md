# Content
The repository contains a main function to solve a Laplace equation. Boundary conditions are Dirichlet conditions (can be homogeneus or heterogeneus).

# Requirements
Linking static libraries from `PACS/Examples` is needed: 
- `pacs`;
- `muparser`;
- `Mesh1D`.
 
The code relays also on other stuffs from `PACS/Examples`:
 - `Matrix` to handle the data structure;
 - `partitioner.hpp` to scatter and gather data independently from their dimensions;
 - `GetPot` to pass parameters from the user.
  
In order to do so is sufficient changing the enviromental variable `PACS_ROOT` in `Makefile`.

# Data structure
For solution `U` and the two grids (horizontal and vertical) has been used a row-wise stored matrix of double.
It is necessary to use two matrices for the grid since is not possible to communicate through MPI a `std::pair<double,double>`.

# Running the code
Compiling:
~~~
make
~~~

Running in parallel:
~~~
mpirun -np 4 ./main
~~~
where 4 can be substitute with the number of processes wanted.

Documentation: 
~~~
make doc
~~~
create a directory `doc`, in which there are two subfolders (`html` and `latex`): in `html`, opening the file `index.html`, will open on the browser the code documentation constructed with Doxygen.

Cleaning from object files, executables, file solution and documentation folder:
~~~
make distclean
~~~

# Code description
The solver is constructed using hybrid parallelization: `MPI` to handle the majority of the algorithm, and `OMP` to further parallelize it internally.
Solution and grid are scattered through the process. Communication is done to switch rows needed from other processes. Then the updating is performed, trivially.
Static variables and `MPI_Wtime()` are used to measure performances.
The convergence is check using OMP reduction and `MPI_Allreduce` with `MPI_LAND`.
Then everything is saved in the file `Solution`, that can be opened through `ParaView`.

Parameters are passed through a `GetPot` file, so once the code has been compiled multiple setting of the experiment can be run.
Function and boundary conditions are pased through `muparser` encoding for the string (`pi` has to be passed as `_pi`).

# Weak points
I apologize in advance for the code, and the challenge in general, since it has been "solved" horribly by me. I am not taking it as an excuse but I could dedicate to it barely 4 days since I was returning home from Erasmus experience, having had there the exam session in the past weeks. I wanted to try the same to do it, I know that is not respecting minimally course's expectations. I was barely able to write a first draft code. I am sorry if correcting it will make you lose time.

Generally:
 - I have no idea of the code efficiency from parallel point of view. I think is quite bad, especially in solution updating.
 - It lacks all the requirements of the challenge about results show and their discussion.
 - `ParaView` has to be opened manually.
 - Sometimes is necessary to make it run multiple times to make it works.
