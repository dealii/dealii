// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// test FE_Q_DG0 (modified step-22)

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>
#include <sstream>

namespace Step22
{
  using namespace dealii;

  template <int dim>
  struct InnerPreconditioner;

  template <>
  struct InnerPreconditioner<2>
  {
    typedef SparseDirectUMFPACK type;
  };

  template <>
  struct InnerPreconditioner<3>
  {
    typedef SparseILU<double> type;
  };

  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem (const unsigned int degree, FESystem<dim> &fe_);
    void run ();

  private:
    void setup_dofs ();
    void assemble_system ();
    void solve ();
    void output_results (const unsigned int refinement_cycle);

    void divergence_velocity(const BlockVector<double> &calc_solution,
                             Vector<double> &output_vector,
                             const Quadrature<dim> &quadrature, bool norm);

    const unsigned int   degree;

    Triangulation<dim>   triangulation;
    FESystem<dim>        &fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     constraints;

    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;

    std_cxx11::shared_ptr<typename InnerPreconditioner<dim>::type>
    A_preconditioner;

    ConvergenceTable convergence_table;
  };

  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution () : Function<dim>(dim+1) {}

    /*virtual*/ double value (const Point<dim>   &p,
                              const unsigned int  component) const;

    /*virtual*/
    Tensor<1,dim> gradient (const Point< dim> &p,
                            const unsigned int component) const;

    /*virtual*/
    double laplacian (const Point<dim>   &p,
                      const unsigned int  component) const;
  };


  template <int dim>
  double
  ExactSolution<dim>::value (const Point<dim>  &p,
                             const unsigned int component) const
  {
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));

    double x = p[0];
    double y = p[1];

    switch (component)
      {
        //velocity
      case 0:
        return 2*(x-1)*(x-1)*x*x*(y-1)*y*(2*y-1);
        break;
      case 1:
        return -2*(y-1)*(y-1)*y*y*(x-1)*x*(2*x-1);
        break;
        //pressure
      case 2:
        //discontinuous Boffi
        return y*(1-y)*exp((x-.5)*(x-.5))-.5+(x<.5);
        //discontinuous simple
        //return (x<.5)-.5;
        //continuous
        //return .5*x*x-1/6;

      }
    ExcNotImplemented();

    return 0;
  }

  template <int dim>
  Tensor<1,dim>
  ExactSolution<dim>::gradient (const Point<dim> &p,
                                const unsigned int component) const
  {
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));

    double x = p[0];
    double y = p[1];

    Tensor<1,dim> gradient;

    switch (component)
      {
        //velocity
      case 0:
        gradient[0]=x*(x*(x*y*(y*(16*y-24)+8)+y*((36-24*y)*y-12))
                       +y*(y*(8*y-12)+4));
        gradient[1]=x*x*(x*(x*(y*(12*y-12)+2)+(24-24*y)*y-4)+y*(12*y-12)+2);
        break;
      case 1:
        gradient[0]=x*(x*((24-12*y)*y-12)*y*y+(y*(12*y-24)+12)*y*y)
                    +((4-2*y)*y-2)*y*y;
        gradient[1]=x*(x*(x*y*((24-16*y)*y-8)+y*(y*(24*y-36)+12))
                       +y*((12-8*y)*y-4));
        break;
        //pressure
      case 2:
        //discontinuous Boffi
        gradient[0]=-exp((x-.5)*(x-.5))*(2*x-1)*(y-1)*y;
        gradient[1]=-exp((x-.5)*(x-.5))*(2*y-1);
        //discontinuous simple
        //gradient[0]=0;
        //gradient[1]=0;
        //continuous
        //gradient[0]=x;
        //gradient[1]=0;
      }

    return gradient;
  }


  template <int dim>
  double
  ExactSolution<dim>::laplacian (const Point<dim>  &p,
                                 const unsigned int component) const
  {
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));

    double x = p[0];
    double y = p[1];

    switch (component)
      {
        //velocity
      case 0:
        return x*(x*(x*(x*(24*y-12)-48*y+24)
                     +y*(y*(48*y-72)+48)-12)
                  +y*((72-48*y)*y-24))
               +y*(y*(8*y-12)+4);

      case 1:
        return x*(x*(x*((48-48*y)*y-8)
                     +y*(72*y-72)+12)
                  +y*(y*((48-24*y)*y-48)+24)-4)
               +(y*(12*y-24)+12)*y*y;
      }
    ExcNotImplemented();
    return 0;
  }

  template <int dim>
  class JumpFunction : public Function<dim>
  {
  public:
    JumpFunction () : Function<dim>(1) {}

    double jump (const Point< dim> &p,
                 const Point<dim> &normal) const;
  };

  template <int dim>
  double
  JumpFunction<dim>::jump (const Point<dim>  &p,
                           const Point<dim> &normal) const
  {
    double x = p[0];
    double y = p[1];
    //discontinuous
    if (std::abs(x-.5)>1e-10)
      return 0;
    if (normal[0]>0)
      return -1;
    return 1;
    //continuous
    //return 0;
  }

  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>(dim+1) {}

    /*virtual*/ double value (const Point<dim>   &p,
                              const unsigned int  component) const;

    const ExactSolution<dim> solution;
  };

  template <int dim>
  double
  RightHandSide<dim>::value (const Point<dim> &p,
                             const unsigned int component) const
  {
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));
    if (component==dim)
      return 0;
    //grad p -laplace u
    return solution.gradient(p,dim)[component]
           -solution.laplacian(p,component);
  }


  template <class Matrix, class Preconditioner>
  class InverseMatrix : public Subscriptor
  {
  public:
    InverseMatrix (const Matrix         &m,
                   const Preconditioner &preconditioner);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const Matrix> matrix;
    const SmartPointer<const Preconditioner> preconditioner;
  };


  template <class Matrix, class Preconditioner>
  InverseMatrix<Matrix,Preconditioner>::InverseMatrix (const Matrix &m,
                                                       const Preconditioner &preconditioner)
    :
    matrix (&m),
    preconditioner (&preconditioner)
  {}

  template <class Matrix, class Preconditioner>
  void InverseMatrix<Matrix,Preconditioner>::vmult (Vector<double>       &dst,
                                                    const Vector<double> &src) const
  {
    SolverControl solver_control (src.size(), 1e-6*src.l2_norm());
    SolverCG<>    cg (solver_control);

    dst = 0;

    cg.solve (*matrix, dst, src, *preconditioner);

    /*std::cout << "  "
    << solver_control.last_step()
    << " inner CG iterations for pressure"
    << std::endl;*/
  }

  template <class Preconditioner>
  class SchurComplement : public Subscriptor
  {
  public:
    SchurComplement (const BlockSparseMatrix<double> &system_matrix,
                     const InverseMatrix<SparseMatrix<double>, Preconditioner>
                     &A_inverse);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>,
          Preconditioner> > A_inverse;

    mutable Vector<double> tmp1, tmp2;
  };



  template <class Preconditioner>
  SchurComplement<Preconditioner>::
  SchurComplement (const BlockSparseMatrix<double> &system_matrix,
                   const InverseMatrix<SparseMatrix<double>,Preconditioner>
                   &A_inverse)
    :
    system_matrix (&system_matrix),
    A_inverse (&A_inverse),
    tmp1 (system_matrix.block(0,0).m()),
    tmp2 (system_matrix.block(0,0).m())
  {}


  template <class Preconditioner>
  void SchurComplement<Preconditioner>::vmult (Vector<double>       &dst,
                                               const Vector<double> &src) const
  {
    system_matrix->block(0,1).vmult (tmp1, src);
    A_inverse->vmult (tmp2, tmp1);
    system_matrix->block(1,0).vmult (dst, tmp2);
  }

  template <int dim>
  StokesProblem<dim>::StokesProblem (const unsigned int degree, FESystem<dim> &fe_)
    :
    degree (degree),
    triangulation (Triangulation<dim>::maximum_smoothing),
    fe (fe_),
    dof_handler (triangulation)
  {}

  template <int dim>
  void StokesProblem<dim>::setup_dofs ()
  {
    A_preconditioner.reset ();
    system_matrix.clear ();

    dof_handler.distribute_dofs (fe);
    DoFRenumbering::Cuthill_McKee (dof_handler);

    std::vector<unsigned int> block_component (dim+1,0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise (dof_handler, block_component);

    {
      constraints.clear ();
      std::vector<bool> component_mask (dim+1, true);
      component_mask[dim] = false;

      DoFTools::make_hanging_node_constraints (dof_handler,
                                               constraints);

      VectorTools::interpolate_boundary_values (dof_handler,
                                                0,
                                                ExactSolution<dim>(),
                                                constraints,
                                                component_mask);

      /*std::vector<bool> boundary_dofs (dof_handler.n_dofs(), false);

      std::vector<bool>boundary_mask (dim+1, false);
      boundary_mask[dim]=true;

      DoFTools::extract_boundary_dofs (dof_handler,boundary_mask,boundary_dofs);

      const unsigned int first_boundary_dof
              = std::distance (boundary_dofs.begin(),std::find
                              (boundary_dofs.begin(),boundary_dofs.end(),true));

      constraints.add_line (first_boundary_dof);
      for (unsigned int i=first_boundary_dof+1; i<dof_handler.n_dofs();++i)
        if (boundary_dofs[i] == true)
          constraints.add_entry (first_boundary_dof, i, -1.);*/

      /*const unsigned int dofs_per_cell = fe.dofs_per_cell;

      typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

      cell->get_dof_indices(local_dof_indices);
      unsigned int first_disc_dof=local_dof_indices[dofs_per_cell-1];
      constraints.add_line (first_disc_dof);

      for (++cell; cell!=endc; ++cell)
      {
        cell->get_dof_indices(local_dof_indices);
        if(cell->at_boundary())
          constraints.add_entry (first_disc_dof,
                                 local_dof_indices[dofs_per_cell-1],-1.);
      }*/
    }

    constraints.close ();

    std::vector<types::global_dof_index> dofs_per_block (2);
    DoFTools::count_dofs_per_block (dof_handler, dofs_per_block,
                                    block_component);
    const unsigned int n_u = dofs_per_block[0],
                       n_p = dofs_per_block[1];

    deallog << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')'
            << std::endl;

    {
      BlockCompressedSimpleSparsityPattern csp (2,2);

      csp.block(0,0).reinit (n_u, n_u);
      csp.block(1,0).reinit (n_p, n_u);
      csp.block(0,1).reinit (n_u, n_p);
      csp.block(1,1).reinit (n_p, n_p);

      csp.collect_sizes();

      DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
      sparsity_pattern.copy_from (csp);
    }

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (2);
    solution.block(0).reinit (n_u);
    solution.block(1).reinit (n_p);
    solution.collect_sizes ();

    system_rhs.reinit (2);
    system_rhs.block(0).reinit (n_u);
    system_rhs.block(1).reinit (n_p);
    system_rhs.collect_sizes ();
  }

  template <int dim>
  void StokesProblem<dim>::assemble_system ()
  {
    system_matrix=0;
    system_rhs=0;

    QGauss<dim>   quadrature_formula(degree+2);
    QGauss<dim-1>   quadrature_face(degree+2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |
                             update_quadrature_points  |
                             update_JxW_values |
                             update_gradients);

    FEFaceValues<dim> fe_v_face (fe, quadrature_face,
                                 update_values | update_quadrature_points |
                                 update_JxW_values | update_normal_vectors);

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;

    const unsigned int   n_q_points      = quadrature_formula.size();
    const unsigned int   n_q_face      = quadrature_face.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const RightHandSide<dim>          right_hand_side;
    std::vector<Vector<double> >      rhs_values (n_q_points,
                                                  Vector<double>(dim+1));

    const JumpFunction<dim> jumpfunction;

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
    std::vector<double>                  div_phi_u   (dofs_per_cell);
    std::vector<double>                  phi_p       (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        local_matrix = 0;
        local_rhs = 0;

        right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                          rhs_values);
        for (unsigned int q=0; q<n_q_points; ++q)
          {
            for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k,
                                   q);
                div_phi_u[k]     = fe_values[velocities].divergence (k, q);
                phi_p[k]         = fe_values[pressure].value (k, q);
              }

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<=i; ++j)
                  {
                    local_matrix(i,j) += (2*symgrad_phi_u[i] * symgrad_phi_u[j]
                                          - div_phi_u[i] * phi_p[j]
                                          - phi_p[i] * div_phi_u[j]
                                          + phi_p[i] * phi_p[j])
                                         * fe_values.JxW(q);

                  }

                const unsigned int component_i =
                  fe.system_to_component_index(i).first;
                local_rhs(i) += fe_values.shape_value(i,q) *
                                rhs_values[q](component_i) *
                                fe_values.JxW(q);
              }
          }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=i+1; j<dofs_per_cell; ++j)
            local_matrix(i,j) = local_matrix(j,i);

        for (unsigned int face_no=0;
             face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
          {
            typename DoFHandler<dim>::face_iterator face= cell->face(face_no);
            if (face->at_boundary()==false)
              {
                typename DoFHandler<dim>::cell_iterator neighbor=
                  cell->neighbor(face_no);
                if (neighbor->index() > cell->index())
                  {
                    fe_v_face.reinit (cell, face_no);

                    const std::vector<Point<dim> > &normals =
                      fe_v_face.get_normal_vectors ();
                    const std::vector<Point<dim> > &quad_points=
                      fe_v_face.get_quadrature_points();

                    for (unsigned int q=0; q<n_q_face; ++q)
                      {
                        double jump=jumpfunction.jump(quad_points[q],normals[q]);
                        for (unsigned int i=0; i<dofs_per_cell; ++i)
                          {
                            const unsigned int component_i =
                              fe.system_to_component_index(i).first;
                            if (component_i<dim)
                              local_rhs(i)+=fe_v_face.shape_value(i,q)
                                            *jump
                                            *normals[q][component_i]
                                            *fe_v_face.JxW(q);
                          }
                      }
                  }
              }
          }

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (local_matrix, local_rhs,
                                                local_dof_indices,
                                                system_matrix, system_rhs);
      }

    A_preconditioner
      = std_cxx11::shared_ptr<typename InnerPreconditioner<dim>::type>(new
          typename InnerPreconditioner<dim>::type());
    A_preconditioner->initialize (system_matrix.block(0,0),
                                  typename
                                  InnerPreconditioner<dim>::type::AdditionalData());

  }

  template <int dim>
  void StokesProblem<dim>::solve ()
  {
    const InverseMatrix<SparseMatrix<double>,
          typename InnerPreconditioner<dim>::type>
          A_inverse (system_matrix.block(0,0), *A_preconditioner);
    Vector<double> tmp (solution.block(0).size());

    {
      Vector<double> schur_rhs (solution.block(1).size());
      A_inverse.vmult (tmp, system_rhs.block(0));
      system_matrix.block(1,0).vmult (schur_rhs, tmp);
      schur_rhs -= system_rhs.block(1);

      SchurComplement<typename InnerPreconditioner<dim>::type>
      schur_complement (system_matrix, A_inverse);

      SolverControl solver_control (solution.block(1).size(),
                                    1e-8*schur_rhs.l2_norm());
      SolverCG<>    cg (solver_control);

      SparseILU<double> preconditioner;
      preconditioner.initialize (system_matrix.block(1,1),
                                 SparseILU<double>::AdditionalData());

      InverseMatrix<SparseMatrix<double>,SparseILU<double> >
      m_inverse (system_matrix.block(1,1), preconditioner);

      cg.solve (schur_complement, solution.block(1), schur_rhs,
                m_inverse);

      constraints.distribute (solution);

    }

    {
      system_matrix.block(0,1).vmult (tmp, solution.block(1));
      tmp *= -1;
      tmp += system_rhs.block(0);

      A_inverse.vmult (solution.block(0), tmp);

      constraints.distribute (solution);
    }
  }

  template <int dim>
  void
  StokesProblem<dim>::output_results (const unsigned int refinement_cycle)
  {
    const ComponentSelectFunction<dim> pressure_mask (dim, dim+1);
    const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
                                                     dim+1);

    ExactSolution<dim> exactsolution;

    const unsigned int n_active_cells=triangulation.n_active_cells();

    Vector<double> difference_per_cell (n_active_cells);

    QGauss<dim> quadrature (degree+3);

    VectorTools::integrate_difference (dof_handler, solution, exactsolution,
                                       difference_per_cell, quadrature,
                                       VectorTools::L2_norm,&velocity_mask);
    const double L2_error_velocity = difference_per_cell.l2_norm();

    VectorTools::integrate_difference (dof_handler, solution, exactsolution,
                                       difference_per_cell, quadrature,
                                       VectorTools::H1_seminorm,&velocity_mask);
    const double H1_error_velocity = difference_per_cell.l2_norm();

    VectorTools::integrate_difference (dof_handler, solution, exactsolution,
                                       difference_per_cell, quadrature,
                                       VectorTools::Linfty_norm,&velocity_mask);
    const double Linfty_error_velocity = difference_per_cell.linfty_norm();

    divergence_velocity(solution, difference_per_cell, quadrature, true);

    const double L2_div_velocity =
      sqrt(difference_per_cell.mean_value()*n_active_cells);

    divergence_velocity(solution, difference_per_cell, quadrature, false);
//    std::cout<<"maximum divergence per cell: "
//             <<difference_per_cell.linfty_norm()<<std::endl;

    //int_\Omega (f(x)-c)^2 dx is minimized for c=1/|\Omega|\int_\Omega f(x)dx.
    //That gives \int_\Omega f(x)^2 dx - 1/|\Omega|(\int_\Omega f(x) dx)^2
    //=\int_\Omega f(x)^2 dx - |\Omega|c^2=\int_\Omega f(x)^2 dx - c^2
    VectorTools::integrate_difference (dof_handler, solution, exactsolution,
                                       difference_per_cell, quadrature,
                                       VectorTools::L2_norm, &pressure_mask);
    /*std::cout<<"Maximal difference per cell:"
             <<difference_per_cell.linfty_norm()<<std::endl;*/

    double L2_error_pressure = difference_per_cell.l2_norm() *
                               difference_per_cell.l2_norm();

//    std::cout<<"l2 difference "<<L2_error_pressure<<std::endl;

    VectorTools::integrate_difference (dof_handler, solution, exactsolution,
                                       difference_per_cell, quadrature,
                                       VectorTools::mean,&pressure_mask);
    double integral = difference_per_cell.mean_value()*n_active_cells;
    //std::cout<<"mean difference "<<integral *integral<<std::endl;
    L2_error_pressure -= integral*integral;

//    std::cout<<"Pressure error squared: "<<L2_error_pressure<<std::endl;
    L2_error_pressure = (L2_error_pressure>0)?sqrt(L2_error_pressure):0;

    /*Vector<double> difference(dim+1);
    const Point<dim> point(.25,0);
    VectorTools::point_difference (dof_handler, solution, exactsolution,
                                    difference, point);
    std::cout<<"difference at 0.25,0 "<<difference[dim]<<std::endl;*/

    VectorTools::integrate_difference (dof_handler, solution, exactsolution,
                                       difference_per_cell, quadrature,
                                       VectorTools::H1_seminorm,&pressure_mask);
    const double H1_error_pressure = difference_per_cell.l2_norm();

    VectorTools::integrate_difference (dof_handler, solution, exactsolution,
                                       difference_per_cell, quadrature,
                                       VectorTools::Linfty_norm,&pressure_mask);
    const double Linfty_error_pressure = difference_per_cell.linfty_norm();

    convergence_table.add_value("cycle", refinement_cycle);
    convergence_table.add_value("cells", n_active_cells);
    convergence_table.add_value("dofs", dof_handler.n_dofs());
    convergence_table.add_value("L2_v", L2_error_velocity);
    convergence_table.add_value("H1_v", H1_error_velocity);
    convergence_table.add_value("Linfty_v", Linfty_error_velocity);
    convergence_table.add_value("L2_div_v", L2_div_velocity);
    convergence_table.add_value("L2_p", L2_error_pressure);
    convergence_table.add_value("H1_p", H1_error_pressure);
    convergence_table.add_value("Linfty_p", Linfty_error_pressure);

    std::vector<std::string> solution_names (dim, "velocity");
    solution_names.push_back ("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, solution_names,
                              DataOut<dim>::type_dof_data,
                              data_component_interpretation);
    data_out.build_patches ();

  }

  template <int dim>
  void StokesProblem<dim>::run ()
  {
    Assert(dim==2, ExcNotImplemented());
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global (1);

    for (unsigned int refinement_cycle = 0; refinement_cycle<5;
         ++refinement_cycle)
      {
        deallog << "Refinement cycle " << refinement_cycle << std::endl;

        if (refinement_cycle > 0)
          triangulation.refine_global(1);

        setup_dofs ();

        assemble_system ();

        solve ();

        output_results (refinement_cycle);

        deallog << std::endl;
      }

    convergence_table.set_precision("L2_v", 7);
    convergence_table.set_precision("H1_v", 7);
    convergence_table.set_precision("Linfty_v", 7);
    convergence_table.set_precision("L2_div_v", 7);
    convergence_table.set_precision("L2_p", 7);
    convergence_table.set_precision("H1_p", 7);
    convergence_table.set_precision("Linfty_p", 7);

    convergence_table.set_scientific("L2_v", true);
    convergence_table.set_scientific("H1_v", true);
    convergence_table.set_scientific("Linfty_v", true);
    convergence_table.set_scientific("L2_div_v", true);
    convergence_table.set_scientific("L2_p", true);
    convergence_table.set_scientific("H1_p", true);
    convergence_table.set_scientific("Linfty_p", true);

    convergence_table.set_tex_caption("cells", "\\# cells");
    convergence_table.set_tex_caption("dofs", "\\# dofs");
    convergence_table.set_tex_caption("L2_v", "$L^2$-error velocity");
    convergence_table.set_tex_caption("H1_v", "$H^1$-error velocity");
    convergence_table.set_tex_caption("Linfty_v",
                                      "$L^\\infty$-error velocity");
    convergence_table.set_tex_caption("L2_div_v",
                                      "$L^2$-error divergence velocity");
    convergence_table.set_tex_caption("L2_p", "$L^2$-error pressure");
    convergence_table.set_tex_caption("H1_p", "$H^1$-error pressure");
    convergence_table.set_tex_caption("Linfty_p",
                                      "$L^\\infty$-error pressure");

    convergence_table.set_tex_format("cells", "r");
    convergence_table.set_tex_format("dofs", "r");

    //std::cout << std::endl;
    //convergence_table.write_text(std::cout);


    convergence_table.add_column_to_supercolumn("cycle", "n cells");
    convergence_table.add_column_to_supercolumn("cells", "n cells");

    std::vector<std::string> new_order;
    new_order.push_back("n cells");
    new_order.push_back("L2_v");
    new_order.push_back("H1_v");
    new_order.push_back("L2_div_v");
    new_order.push_back("L2_p");
    new_order.push_back("H1_p");
    convergence_table.set_column_order (new_order);

    convergence_table.evaluate_convergence_rates
    ("L2_v",ConvergenceTable::reduction_rate_log2);
    convergence_table.evaluate_convergence_rates
    ("H1_v",ConvergenceTable::reduction_rate_log2);
    convergence_table.evaluate_convergence_rates
    ("L2_div_v",ConvergenceTable::reduction_rate_log2);
    convergence_table.evaluate_convergence_rates
    ("L2_p",ConvergenceTable::reduction_rate_log2);
    convergence_table.evaluate_convergence_rates
    ("H1_p",ConvergenceTable::reduction_rate_log2);

    convergence_table.write_text(deallog.get_file_stream());
  }

  //squared l^2 norm of divergence of velocity
  template<int dim>
  void StokesProblem<dim>::divergence_velocity
  (const BlockVector<double> &calc_solution,
   Vector<double>         &output_vector,
   const Quadrature<dim> &quadrature, bool norm)
  {
    output_vector = 0;

    FEValues<dim> fe_v(fe, quadrature, update_gradients | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature.size();

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    for (; cell!=endc; ++cell)
      {
        fe_v.reinit(cell);
        cell->get_dof_indices (local_dof_indices);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            double div = 0;
            for (unsigned int i=0; i < dofs_per_cell-1; ++i)
              {
                double tmp=0;
                for (unsigned int d = 0; d < dim; ++d)
                  tmp+= fe_v.shape_grad_component(i,q,d)[d];

                div+=tmp*calc_solution(local_dof_indices[i]);
              }
            if (norm)
              output_vector(cell->index()) += div * div * fe_v.JxW(q); //L^2-Norm
            else
              output_vector(cell->index()) += div * fe_v.JxW(q); //Integral
          }
      }
  }
}

int main ()
{
  using namespace dealii;
  using namespace Step22;

  initlog();

  deallog.depth_file (1);

  unsigned int degree;
  const unsigned int dim=2;
  {
    degree=1;
    FESystem<dim> fe(FE_Q<dim>(degree+1), dim, FE_Q<dim>(degree), 1);
    deallog << fe.get_name() << ":" << std::endl;
    StokesProblem<2> flow_problem(degree, fe);
    flow_problem.run ();
  }
  {
    degree=1;
    FESystem<2> fe(FE_Q<dim>(degree+1), dim, FE_Q_DG0<dim>(degree), 1);
    deallog << fe.get_name() << ":" << std::endl;
    StokesProblem<2> flow_problem(degree, fe);
    flow_problem.run ();
  }
  {
    degree=2;
    FESystem<2> fe(FE_Q<dim>(degree+1), dim, FE_Q_DG0<dim>(degree), 1);
    deallog << fe.get_name() << ":" << std::endl;
    StokesProblem<2> flow_problem(degree, fe);
    flow_problem.run ();
  }

  return 0;
}
