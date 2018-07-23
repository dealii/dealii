/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Guido Kanschat, Texas A&M University, 2009
 */


// The first few files have already been covered in previous examples and will
// thus not be further commented on:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q1.h>
// Here the discontinuous finite elements are defined. They are used in the
// same way as all other finite elements, though -- as you have seen in
// previous tutorial programs -- there isn't much user interaction with finite
// element classes at all: they are passed to <code>DoFHandler</code> and
// <code>FEValues</code> objects, and that is about it.
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
// We are going to use the simplest possible solver, called Richardson
// iteration, that represents a simple defect correction. This, in combination
// with a block SSOR preconditioner (defined in precondition_block.h), that
// uses the special block matrix structure of system matrices arising from DG
// discretizations.
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/precondition.h>
// We are going to use gradients as refinement indicator.
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/meshworker/mesh_loop.h>

// Like in all programs, we finish this section by including the needed C++
// headers and declaring we want to use objects in the dealii namespace
// without prefix.
#include <iostream>
#include <fstream>

#include <deal.II/fe/fe_facet.h>

namespace Step12
{
  using namespace dealii;

  template <int dim>
  class Solution:  public Function<dim>
  {
  public:
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double> &values,
                             const unsigned int component=0) const
    {
      for (unsigned int i=0; i<values.size(); ++i)
        values[i]=ref.value(points[i]);
    }
  private:
    Functions::LSingularityFunction ref;
  };


  template <int dim>
  class Viscosity:  public Function<dim>
  {
  public:
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double> &values,
                             const unsigned int component=0) const
    {
      for (unsigned int i=0; i<values.size(); ++i)
        values[i]=1.0;
    }
  };


  template <int dim>
  class RHS:  public Function<dim>
  {
  public:
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double> &values,
                             const unsigned int component=0) const;
  };

  template <int dim>
  void RHS<dim>::value_list(const std::vector<Point<dim> > &points,
                            std::vector<double> &values,
                            const unsigned int) const
  {
    Assert(values.size()==points.size(),
           ExcDimensionMismatch(values.size(),points.size()));

    for (unsigned int i=0; i<values.size(); ++i)
      values[i]=0;
  }


  template <int dim>
  class Beta
  {
  public:
    Beta () {}
    void value_list (const std::vector<Point<dim> > &points,
                     std::vector<Point<dim> > &values) const;
  };

  template <int dim>
  void Beta<dim>::value_list(const std::vector<Point<dim> > &points,
                             std::vector<Point<dim> > &values) const
  {
    Assert(values.size()==points.size(),
           ExcDimensionMismatch(values.size(),points.size()));

    for (unsigned int i=0; i<points.size(); ++i)
      {
        const Point<dim> &p=points[i];
        Point<dim> &beta=values[i];

        beta(0) = -p(1);
        beta(1) = p(0);
        beta /= std::sqrt(beta.square());
      }
  }

  // @sect3{Equation data}
  //
  // First, we define a class describing the inhomogeneous boundary
  // data. Since only its values are used, we implement value_list(), but
  // leave all other functions of Function undefined.
  template <int dim>
  class BoundaryValues:  public Function<dim>
  {
  public:
    BoundaryValues () {};
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double> &values,
                             const unsigned int component=0) const;
  };

  // Given the flow direction, the inflow boundary of the unit square
  // $[0,1]^2$ are the right and the lower boundaries. We prescribe
  // discontinuous boundary values 1 and 0 on the x-axis and value 0 on the
  // right boundary. The values of this function on the outflow boundaries
  // will not be used within the DG scheme.
  template <int dim>
  void BoundaryValues<dim>::value_list(const std::vector<Point<dim> > &points,
                                       std::vector<double> &values,
                                       const unsigned int) const
  {
    Assert(values.size()==points.size(),
           ExcDimensionMismatch(values.size(),points.size()));

    for (unsigned int i=0; i<values.size(); ++i)
      {
//        values[i] = sin(2*numbers::PI*points[i](0));
        if (points[i](0)<0.5)
          values[i]=1.;
        else
          values[i]=0.;
      }
  }
  // @sect3{The AdvectionProblem class}
  //
  // After this preparations, we proceed with the main class of this program,
  // called AdvectionProblem. It is basically the main class of step-6. We do
  // not have a ConstraintMatrix, because there are no hanging node
  // constraints in DG discretizations.

  // Major differences will only come up in the implementation of the assemble
  // functions, since here, we not only need to cover the flux integrals over
  // faces, we also use the MeshWorker interface to simplify the loops
  // involved.
  template <int dim>
  class AdvectionProblem
  {
  public:
    AdvectionProblem ();
    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void solve (Vector<double> &solution);
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;

    ConstraintMatrix constraints;

    // Furthermore we want to use DG elements of degree 1 (but this is only
    // specified in the constructor). If you want to use a DG method of a
    // different degree the whole program stays the same, only replace 1 in
    // the constructor by the desired polynomial degree.
    FE_DGQ<dim>          fe;
    DoFHandler<dim>      dof_handler;

    // The next four members represent the linear system to be
    // solved. <code>system_matrix</code> and <code>right_hand_side</code> are
    // generated by <code>assemble_system()</code>, the <code>solution</code>
    // is computed in <code>solve()</code>. The <code>sparsity_pattern</code>
    // is used to determine the location of nonzero elements in
    // <code>system_matrix</code>.
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;

    // Finally, we have to provide functions that assemble the cell, boundary,
    // and inner face terms. Within the MeshWorker framework, the loop over
    // all cells and much of the setup of operations will be done outside this
    // class, so all we have to provide are these three operations. They will
    // then work on intermediate objects for which first, we here define
    // typedefs to the info objects handed to the local integration functions
    // in order to make our life easier below.
    typedef MeshWorker::DoFInfo<dim> DoFInfo;
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;

  };


  // We start with the constructor. The 1 in the constructor call of
  // <code>fe</code> is the polynomial degree.
  template <int dim>
  AdvectionProblem<dim>::AdvectionProblem ()
    :
    mapping (),
    fe (3),
    dof_handler (triangulation)
  {}


  template <int dim>
  void AdvectionProblem<dim>::setup_system ()
  {
    // In the function that sets up the usual finite element data structures,
    // we first need to distribute the DoFs.
    dof_handler.distribute_dofs (fe);

    // We start by generating the sparsity pattern. To this end, we first fill
    // an intermediate object of type DynamicSparsityPattern with the
    // couplings appearing in the system. After building the pattern, this
    // object is copied to <code>sparsity_pattern</code> and can be discarded.

    // To build the sparsity pattern for DG discretizations, we can call the
    // function analogue to DoFTools::make_sparsity_pattern, which is called
    // DoFTools::make_flux_sparsity_pattern:
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    // Finally, we set up the structure of all components of the linear
    // system.
    system_matrix.reinit (sparsity_pattern);
    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
  }

  // @sect4{The assemble_system function}



  // Here we see the major difference to assembling by hand. Instead of
  // writing loops over cells and faces, we leave all this to the MeshWorker
  // framework. In order to do so, we just have to define local integration
  // functions and use one of the classes in namespace MeshWorker::Assembler
  // to build the global system.
  template <int dim>
  void AdvectionProblem<dim>::assemble_system ()
  {

    typedef decltype(dof_handler.begin_active()) Iterator;
    const Beta<dim> beta_function;
    const RHS<dim> rhs_function;
    const Viscosity<dim> viscosity_function;
    const Solution<dim> boundary_function;

    auto cell_worker = [&] (const Iterator &cell, ScratchData<dim> &scratch_data, CopyData &copy_data)
    {
      scratch_data.fe_values.reinit (cell);
      const unsigned int dofs_per_cell   = scratch_data.fe_values.get_fe().dofs_per_cell;
      const unsigned int n_q_points      = scratch_data.fe_values.get_quadrature().size();

      copy_data.cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      copy_data.cell_rhs.reinit (dofs_per_cell);

      copy_data.local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices (copy_data.local_dof_indices);



      const FEValues<dim> &fe_v = scratch_data.fe_values;
      const std::vector<double> &JxW = fe_v.get_JxW_values ();

      std::vector<double> nu (fe_v.n_quadrature_points);
      viscosity_function.value_list (fe_v.get_quadrature_points(), nu);
      std::vector<double> rhs (fe_v.n_quadrature_points);
      rhs_function.value_list (fe_v.get_quadrature_points(), rhs);

      for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
        for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
              copy_data.cell_matrix(i,j) +=
                  // nu \nabla u \nabla v
                  nu[point]
                  * fe_v.shape_grad(i,point)
                  * fe_v.shape_grad(j,point)
                  * JxW[point];

            copy_data.cell_rhs(i) += rhs[point] * fe_v.shape_value(i,point) * JxW[point];
          }
    };

    auto boundary_worker = [&] (const Iterator &cell, const unsigned int &face_no, ScratchData<dim> &scratch_data, CopyData &copy_data)
    {
        scratch_data.fe_facet_values.reinit(cell, face_no);
        const FEFaceValuesBase<dim> &fe_face =
          scratch_data.fe_facet_values.get_fe_values();

      const auto &q_points = fe_face.get_quadrature_points();

      const std::vector<double> &JxW = fe_face.get_JxW_values ();
      const std::vector<Tensor<1,dim> > &normals = fe_face.get_normal_vectors ();

      std::vector<double> nu (fe_face.n_quadrature_points);
      viscosity_function.value_list (fe_face.get_quadrature_points(), nu);

      std::vector<double> g(q_points.size());
      boundary_function.value_list (q_points, g);

      const double degree = std::max(1.0, static_cast<double> (fe_face.get_fe().degree));
      const double extent1 = cell->extent_in_direction(GeometryInfo<dim>::unit_normal_direction[face_no]);
      const double penalty = 4.0 * degree * (degree+1.0) / (extent1);

      for (unsigned int point=0; point<q_points.size(); ++point)
        {
          for (unsigned int i=0; i<fe_face.dofs_per_cell; ++i)
            for (unsigned int j=0; j<fe_face.dofs_per_cell; ++j)
              copy_data.cell_matrix(i,j) +=
                  (
                    // - nu (\nabla u . n) v
                    - nu[point]
                    * (fe_face.shape_grad(j,point) * normals[point])
                    * fe_face.shape_value(i,point)

                    // - nu u (\nabla v . n)  // NIPG: use +
                    - nu[point]
                    * fe_face.shape_value(j,point)
                    * (fe_face.shape_grad(i,point) * normals[point])

                    // + nu * penalty u v
                    + nu[point]
                    * penalty
                    * fe_face.shape_value(j,point)
                    * fe_face.shape_value(i,point)
                    ) * JxW[point];

          for (unsigned int i=0; i<fe_face.dofs_per_cell; ++i)
            copy_data.cell_rhs(i) +=
                (
                  // -nu g (\nabla v . n) // NIPG: use +
                  - nu[point]
                  * g[point]
                  * (fe_face.shape_grad(i,point) * normals[point])

                  // +nu penalty g v
                  + nu[point]
                  * penalty
                  * g[point]
                  * fe_face.shape_value(i,point)
                  ) * JxW[point];

        }
    };

    auto face_worker = [&]
                       (const Iterator &cell, const unsigned int &f, const unsigned int &sf,
                        const Iterator &ncell, const unsigned int &nf, const unsigned int &nsf,
                        ScratchData<dim> &scratch_data, CopyData &copy_data)
    {
      FEFacetValues<dim> &fe_facet = scratch_data.fe_facet_values;
      fe_facet.reinit(cell, f, sf, ncell, nf, nsf);
      const auto &q_points = fe_facet.get_fe_values().get_quadrature_points();

      copy_data.face_data.emplace_back();
      CopyDataFace &copy_data_face = copy_data.face_data.back();

      const unsigned int n_dofs        = fe_facet.n_facet_dofs();
      copy_data_face.joint_dof_indices = fe_facet.get_facet_dof_indices();

      copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);

      const std::vector<double> &        JxW = fe_facet.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals =
        fe_facet.get_normal_vectors();

      std::vector<double> nu (q_points.size());
      viscosity_function.value_list (q_points, nu);

      const double degree = std::max(1.0, static_cast<double>(fe_facet.get_fe_values().get_fe().degree));
      const double extent1 = cell->extent_in_direction(GeometryInfo<dim>::unit_normal_direction[f])
          * (cell->has_children() ? 2.0 : 1.0);
      const double extent2 = ncell->extent_in_direction(GeometryInfo<dim>::unit_normal_direction[nf])
          * (ncell->has_children() ? 2.0 : 1.0);
      const double penalty =4.0 * degree * (degree+1.0) * 0.5 * (1.0/extent1 + 1.0/extent2);

      for (unsigned int point=0; point<q_points.size(); ++point)
        {
          for (unsigned int i=0; i<n_dofs; ++i)
            for (unsigned int j=0; j<n_dofs; ++j)
              copy_data_face.cell_matrix(i,j) +=
                  (
                  // - nu {\nabla u}.n [v] (consistency)
                  - nu[point]
                  * (fe_facet.scalar().gradient_avg(j, point) * normals[point])
                  * fe_facet.scalar().jump(i, point)

                    // - nu [u] {\nabla v}.n  (symmetry) // NIPG: use +
                    - nu[point]
                    * fe_facet.scalar().jump(j, point)
                    * (fe_facet.scalar().gradient_avg(i, point) * normals[point])

                    // nu sigma [u] [v] (penalty)
                    + nu[point] * penalty
                    * fe_facet.scalar().jump(j, point)
                    * fe_facet.scalar().jump(i, point)

                  ) * JxW[point];
        }

    };

    auto copier = [&] (const CopyData &c)
    {
        copy(c, constraints, system_matrix, system_rhs);
      };

    const unsigned int n_gauss_points = dof_handler.get_fe().degree+1;

    ScratchData<dim> scratch_data(mapping, fe, n_gauss_points);
    CopyData cd;
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          cd,
                          MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces | MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker,
                          face_worker);


  }



  // @sect3{All the rest}
  //
  // For this simple problem we use the simplest possible solver, called
  // Richardson iteration, that represents a simple defect correction. This,
  // in combination with a block SSOR preconditioner, that uses the special
  // block matrix structure of system matrices arising from DG
  // discretizations. The size of these blocks are the number of DoFs per
  // cell. Here, we use a SSOR preconditioning as we have not renumbered the
  // DoFs according to the flow field. If the DoFs are renumbered in the
  // downstream direction of the flow, then a block Gauss-Seidel
  // preconditioner (see the PreconditionBlockSOR class with relaxation=1)
  // does a much better job.
  template <int dim>
  void AdvectionProblem<dim>::solve (Vector<double> &solution)
  {
    SolverControl           solver_control (1000, 1e-12);
    SolverGMRES<>      solver (solver_control);

    // Here we create the preconditioner,
    // then assign the matrix to it and set the right block size:
    PreconditionBlockSSOR<SparseMatrix<double> > preconditioner;
    preconditioner.initialize(system_matrix, fe.dofs_per_cell);


    //PreconditionSSOR<SparseMatrix<double> > preconditioner;
    //preconditioner.initialize(system_matrix, fe.dofs_per_cell);

    // After these preparations we are ready to start the linear solver.
    solver.solve (system_matrix, solution, system_rhs,
                  preconditioner);
  }


  // We refine the grid according to a very simple refinement criterion,
  // namely an approximation to the gradient of the solution. As here we
  // consider the DG(1) method (i.e. we use piecewise bilinear shape
  // functions) we could simply compute the gradients on each cell. But we do
  // not want to base our refinement indicator on the gradients on each cell
  // only, but want to base them also on jumps of the discontinuous solution
  // function over faces between neighboring cells. The simplest way of doing
  // that is to compute approximative gradients by difference quotients
  // including the cell under consideration and its neighbors. This is done by
  // the <code>DerivativeApproximation</code> class that computes the
  // approximate gradients in a way similar to the
  // <code>GradientEstimation</code> described in step-9 of this tutorial. In
  // fact, the <code>DerivativeApproximation</code> class was developed
  // following the <code>GradientEstimation</code> class of step-9. Relating
  // to the discussion in step-9, here we consider $h^{1+d/2}|\nabla_h
  // u_h|$. Furthermore we note that we do not consider approximate second
  // derivatives because solutions to the linear advection equation are in
  // general not in $H^2$ but in $H^1$ (to be more precise, in $H^1_\beta$)
  // only.
  template <int dim>
  void AdvectionProblem<dim>::refine_grid ()
  {
    // The <code>DerivativeApproximation</code> class computes the gradients
    // to float precision. This is sufficient as they are approximate and
    // serve as refinement indicators only.
    Vector<float> gradient_indicator (triangulation.n_active_cells());

    // Now the approximate gradients are computed
    DerivativeApproximation::approximate_gradient (mapping,
                                                   dof_handler,
                                                   solution,
                                                   gradient_indicator);

    // and they are cell-wise scaled by the factor $h^{1+d/2}$
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
      gradient_indicator(cell_no)*=std::pow(cell->diameter(), 1+1.0*dim/2);

    // Finally they serve as refinement indicator.
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     gradient_indicator,
                                                     0.3, 0.1);

    triangulation.execute_coarsening_and_refinement ();
  }


  // The output of this program consists of eps-files of the adaptively
  // refined grids and the numerical solutions given in gnuplot format. This
  // was covered in previous examples and will not be further commented on.
  template <int dim>
  void AdvectionProblem<dim>::output_results (const unsigned int cycle) const
  {
    // estimate
    {
      Vector<float> difference_per_cell (triangulation.n_active_cells());
      VectorTools::integrate_difference (dof_handler,
                                         solution,
                                         Solution<dim>(),
                                         difference_per_cell,
                                         QGauss<dim>(fe.degree+2),
                                         VectorTools::L2_norm);

      const double L2_error = VectorTools::compute_global_error(triangulation,
                                       difference_per_cell,
                                       VectorTools::L2_norm);
      std::cout << "L2: " << L2_error << std::endl;

    }


    // Write the grid in eps format.
    std::string filename = "grid-";
    filename += ('0' + cycle);
    Assert (cycle < 10, ExcInternalError());

    filename += ".eps";
    deallog << "Writing grid to <" << filename << ">" << std::endl;
    std::ofstream eps_output (filename.c_str());

    GridOut grid_out;
    grid_out.write_eps (triangulation, eps_output);

    // Output of the solution in gnuplot format.
    filename = "sol-";
    filename += ('0' + cycle);
    Assert (cycle < 10, ExcInternalError());

    filename += ".vtu";
    deallog << "Writing solution to <" << filename << ">" << std::endl;
    std::ofstream gnuplot_output (filename.c_str());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "u");

    data_out.build_patches ();

    data_out.write_vtu(gnuplot_output);
  }


  // The following <code>run</code> function is similar to previous examples.
  template <int dim>
  void AdvectionProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<6; ++cycle)
      {
        deallog << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          {
            //GridGenerator::hyper_cube (triangulation);
            GridGenerator::hyper_L (triangulation);

            triangulation.refine_global (2);
          }
        else
          refine_grid ();


        deallog << "Number of active cells:       "
                << triangulation.n_active_cells()
                << std::endl;

        setup_system ();

        deallog << "Number of degrees of freedom: "
                << dof_handler.n_dofs()
                << std::endl;

        assemble_system ();
        solve (solution);

        output_results (cycle);
      }
  }
}


// The following <code>main</code> function is similar to previous examples as
// well, and need not be commented on.
int main ()
{
  dealii::deallog.depth_console(99);
  try
    {
      Step12::AdvectionProblem<2> dgmethod;
      dgmethod.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };

  return 0;
}
