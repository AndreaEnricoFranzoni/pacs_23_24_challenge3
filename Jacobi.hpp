#ifndef PACS_CHALLENGE3_JACOBI_SOLVER_HPP
#define PACS_CHALLENGE3_JACOBI_SOLVER_HPP

#include <omp.h>
#include <functional>
#include <string>
#include "mesh.hpp"
#include "Matrix.hpp"
#include "writeVTK.hpp"




namespace APSC::JACOBI_SOLVER
{    
    //boundaries for the domain
    constexpr double a = 0.0;
    constexpr double b = 1.0;

    class Jacobi
    {
        public:

        //constructor
        Jacobi(std::size_t N = 1, std::function<double(double,double)> boundary_condition = [](double i1, double i2){return 0.0;})
            :   m_n{N},
                m_omega(a,b), 
                m_grid_1D(m_omega,N-1),
                m_U(N,N),
                m_grid_x(N,N),
                m_grid_y(N,N)
            {   
                //mesh spacing
                N > 1 ? (m_h = 1/(N-1)) : (m_h = 1);

                //grid and boundary conditions init, relaying on simple parallel code
#pragma omp parallel for num_threads(2)
                for (std::size_t i=0; i<N; ++i)
                {   
                    double tmp = this->m_grid_1D[i];
                    this->m_U(0,i) = boundary_condition(0.0,tmp);
                    this->m_U(i,0) = boundary_condition(tmp,0.0);
                    this->m_U(N-1,i) = boundary_condition(1.0,tmp);
                    this->m_U(i,N-1) = boundary_condition(tmp,1.0);

                    for (size_t j = 0; j < N; ++j)
                    {
                        this->m_grid_x(i,j) = this->m_grid_1D[j];
                        this->m_grid_y(i,j) = tmp;
                    }              
                }
            }


        /*!
        * Getter for number of points in the grid for each dimension
        * @return the private member m_n: read-only
        */
       auto n() const {return m_n;};

       /*!
        * Getter for mesh spacing
        * @return the private member m_h: read-only
        */
       auto h() const {return m_h;};
        
        /*!
        * Getter for 1-D grid. 
        * @return the private member m_grid_h: read-only
        */
        auto grid_1D() const {return m_grid_1D;};

        /*!
        * Getter for the x-coordinates of the grid. 
        * @return the private member m_grid_x: read-only. It is only necessary to access the elements of it
        */
        auto grid_x() const {return m_grid_x;};

        /*!
        * Getter for the y-coordinates of the grid. 
        * @return the private member m_grid_x: read-only. It is only necessary to access the elements of it
        */
        auto grid_y() const {return m_grid_y;};

        /*!
        * Getter for the solution. 
        * @return the private member m_U: read-only
        */
        auto U() const {return m_U;};

        /*!
        * Setter for the solution. 
        * @return the private member m_m, modifying it (necessary to update the solution)
        */
        auto &U() {return m_U;};

        /*!
        * To create a VTK file
        */
       void fileVTK() const
       { const std::string filename="Solution";
        generateVTKFile(filename,this->m_U,this->m_n,this->m_n, this->m_h,this->m_h);
       }
        
        

        private:
        
        std::size_t m_n = 1;                                                                                            //number of nodes in the grid for each dimension
        std::size_t m_h = 1;                                                                                            //mesh spacing
        Geometry::Domain1D m_omega;                                                                                     //domain
        Geometry::Mesh1D m_grid_1D;                                                                                     //1-D grid
        apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> m_grid_x;                           //x-grid (n x n matrix)
        apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> m_grid_y;                           //y-grid (n x n matrix)       
        apsc::LinearAlgebra::Matrix<double,apsc::LinearAlgebra::ORDERING::ROWMAJOR> m_U;                                //storing the results (n x n matrix)

    };
}


#endif