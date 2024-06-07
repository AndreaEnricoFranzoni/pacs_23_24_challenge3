#include <omp.h>
#include <mpi.h>
#include <vector>
#include <functional>
#include <iostream>

#include "partitioner.hpp"
#include "GetPot"
#include "readParameters.hpp"
#include "muparser_fun_2D.hpp"

#include "Jacobi.hpp"




void
printHelp()
{
  std::cout
    << "USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"
    << std::endl;
  std::cout << "-h this help" << std::endl;
  std::cout << "-v verbose output" << std::endl;
}

/*!
* Static variables to measure time of parallel code
*/
static double c_start, c_diff;
#define tic() c_start = MPI_Wtime();
#define toc(x)                                       \
  {                                                  \
    c_diff = MPI_Wtime() - c_start;                  \
    std::cout << x << c_diff << " [s]" << std::endl; \
  }





int main(int argc,char** argv)
{
    using namespace std;
    int status(0);

    //reading parameters
    GetPot cl(argc, argv);
    if(cl.search(2, "-h", "--help"))
    {
      printHelp();
      return 0;
    }
    bool verbose = cl.search(1, "-v");

    // Get file with parameter values
    string filename = cl.follow("parameters.pot", "-p");

    parameters param=readParameters(filename, verbose);
    //parameters
    const auto &[f,boundary,n,eps,itermax] = param;


    //Parsing the function and the boundary conditions
    //wrapper for the function
    MuparserFun2 f_pars(f);
    function<double(double,double)> fun = f_pars;       
    //wrappers for boundary conditions
    MuparserFun2 boundary_pars(boundary);
    function<double(double,double)> bc = boundary_pars;

    //solver for the solution
    APSC::JACOBI_SOLVER::Jacobi solver(n,bc);

    //condition if you are not actually partitioning the grid: already returning
    if (n==1)
    {
      return status;
    }
    
    //starting hybrid parallel part
    tic();
    int provided;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_SINGLE,&provided);

    //processes and number of processses
    int rank, size;
    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    MPI_Comm_rank(mpi_comm, &rank);
    MPI_Comm_size(mpi_comm, &size);

    //to scatter data structure
    std::vector<int> counts(size);
    std::vector<int> displacements(size);

    if (rank == 0) 
    { 
        //partitioning of the two matrices: solution and grid
        apsc::MatrixPartitioner matrix_part(n,n,size);
        auto cd = apsc::counts_and_displacements(matrix_part);
        counts = cd[0];
        displacements = cd[1];        
    }

    //Chunks sizes scattered across the ranks
    int local_size;
    MPI_Scatter(counts.data(),1,MPI_INT, &local_size,1,MPI_INT,0,mpi_comm);

    //each process gets its chunk of U and grid (grid is needed for function evaluation)                  
    apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> local_U(local_size/n,n);
    apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> local_grid_x(local_size/n,n);
    apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> local_grid_y(local_size/n,n); 

    //scattering U
    MPI_Scatterv(solver.U().data(),counts.data(),displacements.data(),MPI_DOUBLE,
                 local_U.data(),local_size, MPI_DOUBLE,0,mpi_comm);

    //scattering the grids
    MPI_Scatterv(solver.grid_x().data(),counts.data(),displacements.data(),MPI_DOUBLE,
                 local_grid_x.data(),local_size,MPI_DOUBLE,0,mpi_comm);
    MPI_Scatterv(solver.grid_y().data(),counts.data(),displacements.data(),MPI_DOUBLE,
                 local_grid_y.data(),local_size,MPI_DOUBLE,0,mpi_comm);


    //constants useful during computations
    constexpr double const_factor = 0.25;
    const double h = solver.h();
    const double h_squared = std::pow(h,2);
    //to check itermax
    std::size_t count = 0;
    bool check_local = false;
    bool converged =false;


    //synchronizing after scattering
    MPI_Barrier(mpi_comm);

    // partner to which I send the last row sent (next process)
    // only rank from 0 to size-2 will do it actually
    int partner_send_last_row = (rank+1)%size; // 0 --> 1, 1 --> 2.. p->0
    // partner from which I receive the last row received (previous process)
    // only rank from 1 to size-1 will do it actually
    int partner_receive_last_row = (rank+size-1)%size;

    //partner to which I send the first row sent (previous process)
    //only rank from 1 to size-1 will do it actually
    int partner_send_first_row = (rank+size-1)%size;
    //partner from which I receive the first row received (next process)
    //only rank from 0 to size-2 will do it actually
    int partner_receive_first_row = (rank+1)%size;

    //Jacobi algo
    while (count < itermax and converged == false)
    { 
      //what to send to previous (my first row)
      apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> send_to_prev(1,n);
      //what to receive from next (its first row)
      apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> recv_from_next(1,n);
      //what to send to next (my last row)
      apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> send_to_next(1,n);
      //what to receive from previous (its last row)
      apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> recv_from_prev(1,n);


      //number of row in each process 
      int n_rows_local = local_size/static_cast<int>(n);
      //new matrix to make update (to be created at each iteration in order to check the convergence)
      apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> u_new(n_rows_local,n);

      //boundaries in all the local solutions
      if (rank==0)  //copying the first row
      {
#pragma omp parallel for num_threads(2)
        for (size_t j = 0; j < n; ++j)
        {
          u_new(0,j) = local_U(0,j);
        }      
      }
      if (rank==size-1) //copying the last row
      {
#pragma omp parallel for num_threads(2)
        for (size_t j = 0; j < n; ++j)
        {
          u_new(n_rows_local-1,j) = local_U(n_rows_local-1,j);
        }    
      }

      //copying the first and last element of each row
#pragma omp parallel for num_threads(2)
        for (int i = 0; i < n_rows_local; ++i)
        {
          u_new(i,0) = local_U(i,0);
          u_new(i,n-1) = local_U(i,n-1);
        }   
      
      
      //preparing what to send

      //rank 0 sends only to next one its last row
      if (rank==0)
      { 
#pragma omp parallel for num_threads(2)
        for (std::size_t j=0; j<n; ++j)
        {
          send_to_next(0,j) = local_U(n_rows_local-1,j);
        }
      }
      //rank size-1 sends only to previous one its first row
      else if (rank==(size-1))
      {
#pragma omp parallel for num_threads(2)
        for (std::size_t j=0; j<n; ++j)
        {
          send_to_prev(0,j) = local_U(0,j);
        }
      }
      //all the other processes send first and last
      else
      {
        //int last_row = (local_size/static_cast<int>(n))-1;
#pragma omp parallel for num_threads(2)
        for (std::size_t j=0; j<n; ++j)
        {
          send_to_prev(0,j) = local_U(0,j);
          send_to_next(0,j) = local_U(n_rows_local-1,j);
        }        
      }

      //in order to have all the processes having detected what to send: synchronizing
      MPI_Barrier(mpi_comm);

      
      //communication: using MPI_Sendrecv to avoid deadlocks
      MPI_Status status;
      //communicating first row to previous     
      MPI_Sendrecv(send_to_prev.data(),n,MPI_DOUBLE,partner_send_first_row,0,
                   recv_from_next.data(),n,MPI_DOUBLE,partner_receive_first_row,0,mpi_comm,&status);
      
      MPI_Barrier(mpi_comm);
      
      //communicating last row to next
      MPI_Sendrecv(send_to_next.data(),n,MPI_DOUBLE,partner_send_last_row,0,
                   recv_from_prev.data(),n,MPI_DOUBLE,partner_receive_last_row,0,mpi_comm,&status);
      
      MPI_Barrier(mpi_comm);

      //UPDATING the solution
      if (rank==0)
      { 
        //no need for updating the first row (boundary conditions)

        //updating central rows
        for (int i = 1; i < n_rows_local-1; ++i)
        {
          for (size_t j = 1; j < n-1; ++j)
          {
              u_new(i,j) = const_factor*( local_U(i-1,j) + local_U(i+1,j) + local_U(i,j-1) + local_U(i,j+1) + h_squared*fun(local_grid_x(i,j),local_grid_y(i,j)) ); 
          }
        }
        //updating the last row with data from the other process
        for (size_t j = 1; j < n-1; ++j)
        {
              u_new(n_rows_local-1,j) = const_factor*( local_U(n_rows_local-2,j) + recv_from_next(0,j) + local_U(n_rows_local-1,j-1) + local_U(n_rows_local-1,j+1) + h_squared*fun(local_grid_x(n_rows_local-1,j),local_grid_y(n_rows_local-1,j)) ); 
        }
      }
      else if (rank==(size-1))
      { 
        //no need for updating last row: boundary conditions

        //first row: takes from previous process
        for (size_t j = 1; j < n-1; ++j)
        {
              u_new(0,j) = const_factor*( recv_from_prev(0,j) + local_U(1,j) + local_U(0,j-1) + local_U(0,j+1) + h_squared*fun(local_grid_x(0,j),local_grid_y(0,j)) ); 
        }
        //central rows: all data are in the corrispective process
        for (int i = 1; i < n_rows_local-1; ++i)
        {
          for (size_t j = 1; j < n-1; ++j)
          {
              u_new(i,j) = const_factor*( local_U(i-1,j) + local_U(i+1,j) + local_U(i,j-1) + local_U(i,j+1) + h_squared*fun(local_grid_x(i,j),local_grid_y(i,j)) ); 
          }
        }
      }
      else
      { 
        //first row: takes from previous process
        for (size_t j = 1; j < n-1; ++j)
        {
              u_new(0,j) = const_factor*( recv_from_prev(0,j) + local_U(1,j) + local_U(0,j-1) + local_U(0,j+1) + h_squared*fun(local_grid_x(0,j),local_grid_y(0,j)) ); 
        }
        //central rows: all data are in the corrispective process
        for (int i = 1; i < n_rows_local-1; ++i)
        {
          for (size_t j = 1; j < n-1; ++j)
          {
              u_new(i,j) = const_factor*( local_U(i-1,j) + local_U(i+1,j) + local_U(i,j-1) + local_U(i,j+1) + h_squared*fun(local_grid_x(i,j),local_grid_y(i,j)) ); 
          }
        }
        //last row: data from next process are needed
        for (size_t j = 1; j < n-1; ++j)
        {
            u_new(n_rows_local-1,j) = const_factor*( local_U(n_rows_local-2,j) + recv_from_next(0,j) + local_U(n_rows_local-1,j-1) + local_U(n_rows_local-1,j+1) + h_squared*fun(local_grid_x(n_rows_local-1,j),local_grid_y(n_rows_local-1,j)) ); 
        }        
      }

      MPI_Barrier(mpi_comm);


      //CHECKING CONVERGENCE
      double sum;

      //for each process: checking sum of squared in two consecutive iterations
#pragma omp parallel for shared(local_U,u_new) reduction(+:sum)
      for (int i=0;i<local_size;++i) sum +=std::pow(local_U.data()[i] - u_new.data()[i],2);

      MPI_Barrier(mpi_comm);

      check_local = std::sqrt(h*sum) < eps;

      //updating solution actually
      std::swap(local_U,u_new);

      MPI_Barrier(mpi_comm);    

      //freeing memory
      u_new.clear();
      send_to_prev.clear();
      recv_from_next.clear();
      send_to_next.clear();
      recv_from_prev.clear();

      //checking the convergence in each process
      MPI_Allreduce(&check_local,&converged,sizeof(check_local),
                    MPI_BYTE,MPI_LAND,mpi_comm);

      count ++;
    }
    
    //gathering the solution
    MPI_Gatherv(local_U.data(),local_size,MPI_DOUBLE,solver.U().data(),counts.data(),
                displacements.data(),MPI_DOUBLE,0,mpi_comm);

    MPI_Finalize();
    toc("Jacobi algorithm finisehd: grid scattered, convergence reached in all the processes, and then gathered: ");

    //file to be opened in ParaView
    solver.fileVTK();

    return status;
}