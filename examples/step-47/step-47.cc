/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Natasha Sharma, University of Texas at El Paso,
 *          Guido Kanschat, University of Heidelberg
 *          Timo Heister, Clemson University
 *          Wolfgang Bangerth, Colorado State University
 *          Zhuroan Wang, Colorado State University
 */


// @sect3{Include files}

// The first few include files have already been used in the previous
// example, so we will not explain their meaning here again. The principal
// structure of the program is very similar to that of, for example, step-4
// and so we include many of the same header files.

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

// The two most interesting header files will be these two:
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/meshworker/mesh_loop.h>
// The first of these is responsible for providing the class FEInterfaceValues
// that can be used to evaluate quantities such as the jump or average
// of shape functions (or their gradients) across interfaces between cells.
// This class will be quite useful in evaluating the penalty terms that appear
// in the C0IP formulation.


#include <fstream>
#include <iostream>
#include <cmath>


namespace Step47
{
  using namespace dealii;


  // In the following namespace, let us define the exact solution against
  // which we will compare the numerically computed one. It has the form
  // $u(x,y) = \sin(\pi x) \sin(\pi y)$ (only the 2d case is implemented),
  // and the namespace also contains a class that corresponds to the right
  // hand side that produces this solution.
  namespace ExactSolution
  {
    using numbers::PI;

    template <int dim>
    class Solution : public Function<dim>
    {
    public:
      static_assert(dim == 2, "Only dim==2 is implemented.");

      virtual double value(const Point<dim> &p,
                           const unsigned int /*component*/ = 0) const override
      {
        return std::sin(PI * p[0]) * std::sin(PI * p[1]);
      }

      virtual Tensor<1, dim>
      gradient(const Point<dim> &p,
               const unsigned int /*component*/ = 0) const override
      {
        Tensor<1, dim> r;
        r[0] = PI * std::cos(PI * p[0]) * std::sin(PI * p[1]);
        r[1] = PI * std::cos(PI * p[1]) * std::sin(PI * p[0]);
        return r;
      }

      virtual void
      hessian_list(const std::vector<Point<dim>> &       points,
                   std::vector<SymmetricTensor<2, dim>> &hessians,
                   const unsigned int /*component*/ = 0) const override
      {
        for (unsigned i = 0; i < points.size(); ++i)
          {
            const double x = points[i][0];
            const double y = points[i][1];

            hessians[i][0][0] = -PI * PI * std::sin(PI * x) * std::sin(PI * y);
            hessians[i][0][1] = PI * PI * std::cos(PI * x) * std::cos(PI * y);
            hessians[i][1][1] = -PI * PI * std::sin(PI * x) * std::sin(PI * y);
          }
      }
    };


    template <int dim>
    class RightHandSide : public Function<dim>
    {
    public:
      static_assert(dim == 2, "Only dim==2 is implemented");

      virtual double value(const Point<dim> &p,
                           const unsigned int /*component*/ = 0) const override

      {
        return 4 * std::pow(PI, 4.0) * std::sin(PI * p[0]) *
               std::sin(PI * p[1]);
      }
    };
  } // namespace ExactSolution



  // @sect3{The main class}
  //
  // The following is the principal class of this tutorial program. It has
  // the structure of many of the other tutorial programs and there should
  // really be nothing particularly surprising about its contents or
  // the constructor that follows it.
  template <int dim>
  class BiharmonicProblem
  {
  public:
    BiharmonicProblem(const unsigned int fe_degree);

    void run();

  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void compute_errors();
    void output_results(const unsigned int iteration) const;

    Triangulation<dim> triangulation;

    MappingQ<dim> mapping;

    FE_Q<dim>                 fe;
    DoFHandler<dim>           dof_handler;
    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
  };



  template <int dim>
  BiharmonicProblem<dim>::BiharmonicProblem(const unsigned int fe_degree)
    : mapping(1)
    , fe(fe_degree)
    , dof_handler(triangulation)
  {}



  // Next up are the functions that create the initial mesh (a once refined
  // unit square) and set up the constraints, vectors, and matrices on
  // each mesh. Again, both of these are essentially unchanged from many
  // previous tutorial programs.
  template <int dim>
  void BiharmonicProblem<dim>::make_grid()
  {
    GridGenerator::hyper_cube(triangulation, 0., 1.);
    triangulation.refine_global(1);

    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Total number of cells: " << triangulation.n_cells()
              << std::endl;
  }



  template <int dim>
  void BiharmonicProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             ExactSolution::Solution<dim>(),
                                             constraints);
    constraints.close();


    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints, true);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }



  // @sect4{Assembling the linear system}
  //
  // The following pieces of code are more interesting. They all relate to the
  // assembly of the linear system. While assembling the cell-interior terms
  // is not of great difficulty -- that works in essence like the assembly
  // of the corresponding terms of the Laplace equation, and you have seen
  // how this works in step-4 or step-6, for example -- the difficulty
  // is with the penalty terms in the formulation. These require the evaluation
  // of gradients of shape functions at interfaces of cells. At the least,
  // one would therefore need to use two FEFaceValues objects, but if one of the
  // two sides is adaptively refined, then one actually needs an FEFaceValues
  // and one FESubfaceValues objects; one also needs to keep track which
  // shape functions live where, and finally we need to ensure that every
  // face is visited only once. All of this is a substantial overhead to the
  // logic we really want to implement (namely the penalty terms in the
  // bilinear form). As a consequence, we will make use of the
  // FEInterfaceValues class -- a helper class in deal.II that allows us
  // to abstract away the two FEFaceValues or FESubfaceValues objects and
  // directly access what we really care about: jumps, averages, etc.
  //
  // But this doesn't yet solve our problem of having to keep track of
  // which faces we have already visited when we loop over all cells and
  // all of their faces. To make this process simpler, we use the
  // MeshWorker::mesh_loop() function that provides a simple interface
  // for this task: Based on the ideas outlined in the WorkStream
  // namespace documentation, MeshWorker::mesh_loop() requires three
  // functions that do work on cells, interior faces, and boundary
  // faces. These functions work on scratch objects for intermediate
  // results, and then copy the result of their computations into
  // copy data objects from where a copier function copies them into
  // the global matrix and right hand side objects.
  //
  // The following structures then provide the scratch and copy objects
  // that are necessary for this approach. You may look up the WorkStream
  // namespace as well as the
  // @ref threads "Parallel computing with multiple processors"
  // module for more information on how they typically work.
  template <int dim>
  struct ScratchData
  {
    ScratchData(const Mapping<dim> &      mapping,
                const FiniteElement<dim> &fe,
                const unsigned int        quadrature_degree,
                const UpdateFlags         update_flags,
                const UpdateFlags         interface_update_flags)
      : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags)
      , fe_interface_values(mapping,
                            fe,
                            QGauss<dim - 1>(quadrature_degree),
                            interface_update_flags)
    {}


    ScratchData(const ScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags())
      , fe_interface_values(scratch_data.fe_values.get_mapping(),
                            scratch_data.fe_values.get_fe(),
                            scratch_data.fe_interface_values.get_quadrature(),
                            scratch_data.fe_interface_values.get_update_flags())
    {}

    FEValues<dim>          fe_values;
    FEInterfaceValues<dim> fe_interface_values;
  };



  struct CopyData
  {
    CopyData(const unsigned int dofs_per_cell)
      : cell_matrix(dofs_per_cell, dofs_per_cell)
      , cell_rhs(dofs_per_cell)
      , local_dof_indices(dofs_per_cell)
    {}


    CopyData(const CopyData &) = default;


    struct FaceData
    {
      FullMatrix<double>                   cell_matrix;
      std::vector<types::global_dof_index> joint_dof_indices;
    };

    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<FaceData>                face_data;
  };



  // The more interesting part is where we actually assemble the linear system.
  // Fundamentally, this function has five parts:
  // - The definition of the `cell_worker` lambda function, a small
  //   function that is defined within the `assemble_system()`
  //   function and that will be responsible for computing the local
  //   integrals on an individual cell. It will work on a copy of the
  //   `ScratchData` class and put its results into the corresponding
  //   `CopyData` object.
  // - The definition of the `face_worker` lambda function that does
  //   the integration of all terms that live on the interfaces between
  //   cells.
  // - The definition of the `boundary_worker` function that does the
  //   same but for cell faces located on the boundary of the domain.
  // - The definition of the `copier` function that is responsible
  //   for copying all of the data the previous three functions have
  //   put into copy objects for a single cell, into the global matrix
  //   and right hand side.
  //
  // The fifth part is the one where we bring all of this together.
  //
  // Let us go through each of these pieces necessary for the assembly
  // in turns.
  template <int dim>
  void BiharmonicProblem<dim>::assemble_system()
  {
    using Iterator = typename DoFHandler<dim>::active_cell_iterator;

    // The first piece is the `cell_worker` that does the assembly
    // on the cell interiors. It is a (lambda) function that takes
    // a cell (input), a scratch object, and a copy object (output)
    // as arguments. It looks like the assembly functions of many
    // other of the tutorial programs, or at least the body of the
    // loop over all cells.
    //
    // The terms we integrate here are the cell contribution
    // @f{align*}{
    //    A^K_{ij} = \int_K \nabla^2\varphi_i(x) : \nabla^2\varphi_j(x) dx
    // @f}
    // to the global matrix, and
    // @f{align*}{
    //    f^K_i = \int_K \varphi_i(x) f(x) dx
    // @f}
    // to the right hand side vector.
    //
    // We use the same technique as used in the assembly of step-22
    // to accelerate the function: Instead of calling
    // `fe_values.shape_hessian(i, qpoint)` in the innermost loop,
    // we instead create a variable `hessian_i` that evaluates this
    // value once in the loop over `i` and re-use the so-evaluated
    // value in the loop over `j`. For symmetry, we do the same with a
    // variable `hessian_j`, although it is indeed only used once and
    // we could have left the call to `fe_values.shape_hessian(j,qpoint)`
    // in the instruction that computes the scalar product between
    // the two terms.
    auto cell_worker = [&](const Iterator &  cell,
                           ScratchData<dim> &scratch_data,
                           CopyData &        copy_data) {
      copy_data.cell_matrix = 0;
      copy_data.cell_rhs    = 0;

      FEValues<dim> &fe_values = scratch_data.fe_values;
      fe_values.reinit(cell);

      cell->get_dof_indices(copy_data.local_dof_indices);

      const ExactSolution::RightHandSide<dim> right_hand_side;

      const unsigned int dofs_per_cell =
        scratch_data.fe_values.get_fe().dofs_per_cell;

      for (unsigned int qpoint = 0; qpoint < fe_values.n_quadrature_points;
           ++qpoint)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const Tensor<2, dim> hessian_i =
                fe_values.shape_hessian(i, qpoint);

              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  const Tensor<2, dim> hessian_j =
                    fe_values.shape_hessian(j, qpoint);

                  copy_data.cell_matrix(i, j) +=
                    scalar_product(hessian_i,   // nabla^2 phi_i(x)
                                   hessian_j) * // nabla^2 phi_j(x)
                    fe_values.JxW(qpoint);      // dx
                }

              copy_data.cell_rhs(i) +=
                fe_values.shape_value(i, qpoint) * // phi_i(x)
                right_hand_side.value(
                  fe_values.quadrature_point(qpoint)) * // f(x)
                fe_values.JxW(qpoint);                  // dx
            }
        }
    };


    // The next building block is the one that assembles penalty terms on each
    // of the interior faces of the mesh. As described in the documentation of
    // MeshWorker::mesh_loop(), this function receives arguments that denote
    // a cell and its neighboring cell, as well as (for each of the two
    // cells) the face (and potentially sub-face) we have to integrate
    // over. Again, we also get a scratch object, and a copy object
    // for putting the results in.
    //
    // The function has three parts itself. At the top, we initialize
    // the FEInterfaceValues object and create a new `CopyData::FaceData`
    // object to store our input in. This gets pushed to the end of the
    // `copy_data.face_data` variable. We need to do this because
    // the number of faces (or subfaces) over which we integrate for a
    // given cell differs from cell to cell, and the sizes of these
    // matrices also differ, depending on what degrees of freedom
    // are adjacent to the face or subface. As discussed in the documentation
    // of MeshWorker::mesh_loop(), the copy object is reset every time a new
    // cell is visited, so that what we push to the end of
    // `copy_data.face_data()` is really all that the later `copier` function
    // gets to see when it copies the contributions of each cell to the global
    // matrix and right hand side objects.
    auto face_worker = [&](const Iterator &    cell,
                           const unsigned int &f,
                           const unsigned int &sf,
                           const Iterator &    ncell,
                           const unsigned int &nf,
                           const unsigned int &nsf,
                           ScratchData<dim> &  scratch_data,
                           CopyData &          copy_data) {
      FEInterfaceValues<dim> &fe_interface_values =
        scratch_data.fe_interface_values;
      fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf);

      copy_data.face_data.emplace_back();
      CopyData::FaceData &copy_data_face = copy_data.face_data.back();

      copy_data_face.joint_dof_indices =
        fe_interface_values.get_interface_dof_indices();

      const unsigned int n_interface_dofs =
        fe_interface_values.n_current_interface_dofs();
      copy_data_face.cell_matrix.reinit(n_interface_dofs, n_interface_dofs);

      // The second part deals with determining what the penalty
      // parameter should be. By looking at the units of the various
      // terms in the bilinear form, it is clear that the penalty has
      // to have the form $\frac{\gamma}{h_K}$ (i.e., one over length
      // scale), but it is not a priori obvious how one should choose
      // the dimension-less number $\gamma$. From the discontinuous
      // Galerkin theory for the Laplace equation, one might
      // conjecture that the right choice is $\gamma=p(p+1)$ is the
      // right choice, where $p$ is the polynomial degree of the
      // finite element used. We will discuss this choice in a bit
      // more detail in the results section of this program.
      //
      // In the formula above, $h_K$ is the size of cell $K$. But this
      // is not quite so straightforward either: If one uses highly
      // stretched cells, then a more involved theory says that $h$
      // should be replaced by the diameter of cell $K$ normal to the
      // direction of the edge in question.  It turns out that there
      // is a function in deal.II for that. Secondly, $h_K$ may be
      // different when viewed from the two different sides of a face.
      //
      // To stay on the safe side, we take the maximum of the two values.
      // We will note that it is possible that this computation has to be
      // further adjusted if one were to use hanging nodes resulting from
      // adaptive mesh refinement.
      const unsigned int p = fe.degree;
      const double       gamma_over_h =
        std::max((1.0 * p * (p + 1) /
                  cell->extent_in_direction(
                    GeometryInfo<dim>::unit_normal_direction[f])),
                 (1.0 * p * (p + 1) /
                  ncell->extent_in_direction(
                    GeometryInfo<dim>::unit_normal_direction[nf])));

      // Finally, and as usual, we loop over the quadrature points and
      // indices `i` and `j` to add up the contributions of this face
      // or sub-face. These are then stored in the
      // `copy_data.face_data` object created above. As for the cell
      // worker, we pull the evaluation of averages and jumps out of
      // the loops if possible, introducing local variables that store
      // these results. The assembly then only needs to use these
      // local variables in the innermost loop. Regarding the concrete
      // formula this code implements, recall that the interface terms
      // of the bilinear form were as follows:
      // @f{align*}{
      //  -\sum_{e \in \mathbb{F}} \int_{e}
      //  \jump{ \frac{\partial v_h}{\partial \mathbf n}}
      //  \average{\frac{\partial^2 u_h}{\partial \mathbf n^2}} \ ds
      // -\sum_{e \in \mathbb{F}} \int_{e}
      // \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}}
      // \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds
      // + \sum_{e \in \mathbb{F}}
      // \frac{\gamma}{h_e}
      // \int_e
      // \jump{\frac{\partial v_h}{\partial \mathbf n}}
      // \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds.
      // @f}
      for (unsigned int qpoint = 0;
           qpoint < fe_interface_values.n_quadrature_points;
           ++qpoint)
        {
          const auto &n = fe_interface_values.normal(qpoint);

          for (unsigned int i = 0; i < n_interface_dofs; ++i)
            {
              const double av_hessian_i_dot_n_dot_n =
                (fe_interface_values.average_hessian(i, qpoint) * n * n);
              const double jump_grad_i_dot_n =
                (fe_interface_values.jump_gradient(i, qpoint) * n);

              for (unsigned int j = 0; j < n_interface_dofs; ++j)
                {
                  const double av_hessian_j_dot_n_dot_n =
                    (fe_interface_values.average_hessian(j, qpoint) * n * n);
                  const double jump_grad_j_dot_n =
                    (fe_interface_values.jump_gradient(j, qpoint) * n);

                  copy_data_face.cell_matrix(i, j) +=
                    (-av_hessian_i_dot_n_dot_n       // - {grad^2 v n n }
                       * jump_grad_j_dot_n           // [grad u n]
                     - av_hessian_j_dot_n_dot_n      // - {grad^2 u n n }
                         * jump_grad_i_dot_n         // [grad v n]
                     +                               // +
                     gamma_over_h *                  // gamma/h
                       jump_grad_i_dot_n *           // [grad v n]
                       jump_grad_j_dot_n) *          // [grad u n]
                    fe_interface_values.JxW(qpoint); // dx
                }
            }
        }
    };


    // The third piece is to do the same kind of assembly for faces that
    // are at the boundary. The idea is the same as above, of course,
    // with only the difference that there are now penalty terms that
    // also go into the right hand side.
    //
    // As before, the first part of the function simply sets up some
    // helper objects:
    auto boundary_worker = [&](const Iterator &    cell,
                               const unsigned int &face_no,
                               ScratchData<dim> &  scratch_data,
                               CopyData &          copy_data) {
      FEInterfaceValues<dim> &fe_interface_values =
        scratch_data.fe_interface_values;
      fe_interface_values.reinit(cell, face_no);
      const auto &q_points = fe_interface_values.get_quadrature_points();

      copy_data.face_data.emplace_back();
      CopyData::FaceData &copy_data_face = copy_data.face_data.back();

      const unsigned int n_dofs =
        fe_interface_values.n_current_interface_dofs();
      copy_data_face.joint_dof_indices =
        fe_interface_values.get_interface_dof_indices();

      copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);

      const std::vector<double> &JxW = fe_interface_values.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals =
        fe_interface_values.get_normal_vectors();


      const ExactSolution::Solution<dim> exact_solution;
      std::vector<Tensor<1, dim>>        exact_gradients(q_points.size());
      exact_solution.gradient_list(q_points, exact_gradients);


      // Positively, because we now only deal with one cell adjacent to the
      // face (as we are on the boundary), the computation of the penalty
      // factor $\gamma$ is substantially simpler:
      const unsigned int p = fe.degree;
      const double       gamma_over_h =
        (1.0 * p * (p + 1) /
         cell->extent_in_direction(
           GeometryInfo<dim>::unit_normal_direction[face_no]));

      // The third piece is the assembly of terms. This is now
      // slightly more involved since these contains both terms for
      // the matrix and for the right hand side. The former is exactly
      // the same as for the interior faces stated above if one just
      // defines the jump and average appropriately (which is what the
      // FEInterfaceValues class does). The latter requires us to
      // evaluate the boundary conditions $j(\mathbf x)$, which in the
      // current case (where we know the exact solution) we compute
      // from $j(\mathbf x) = \frac{\partial u(\mathbf x)}{\partial
      // {\mathbf n}}$. The term to be added to the right hand side
      // vector is then
      // $\frac{\gamma}{h_e}\int_e
      // \jump{\frac{\partial v_h}{\partial \mathbf n}} j \ ds$.
      for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
        {
          const auto &n = normals[qpoint];

          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              const double av_hessian_i_dot_n_dot_n =
                (fe_interface_values.average_hessian(i, qpoint) * n * n);
              const double jump_grad_i_dot_n =
                (fe_interface_values.jump_gradient(i, qpoint) * n);

              for (unsigned int j = 0; j < n_dofs; ++j)
                {
                  const double av_hessian_j_dot_n_dot_n =
                    (fe_interface_values.average_hessian(j, qpoint) * n * n);
                  const double jump_grad_j_dot_n =
                    (fe_interface_values.jump_gradient(j, qpoint) * n);

                  copy_data_face.cell_matrix(i, j) +=
                    (-av_hessian_i_dot_n_dot_n  // - {grad^2 v n n}
                       * jump_grad_j_dot_n      //   [grad u n]
                                                //
                     - av_hessian_j_dot_n_dot_n // - {grad^2 u n n}
                         * jump_grad_i_dot_n    //   [grad v n]
                                                //
                     + gamma_over_h             //  gamma/h
                         * jump_grad_i_dot_n    // [grad v n]
                         * jump_grad_j_dot_n    // [grad u n]
                     ) *
                    JxW[qpoint]; // dx
                }

              copy_data.cell_rhs(i) +=
                (-av_hessian_i_dot_n_dot_n *       // - {grad^2 v n n }
                   (exact_gradients[qpoint] * n)   //   (grad u_exact . n)
                 +                                 // +
                 gamma_over_h                      //  gamma/h
                   * jump_grad_i_dot_n             // [grad v n]
                   * (exact_gradients[qpoint] * n) // (grad u_exact . n)
                 ) *
                JxW[qpoint]; // dx
            }
        }
    };

    // Part 4 is a small function that copies the data produced by the
    // cell, interior, and boundary face assemblers above into the
    // global matrix and right hand side vector. There really is not
    // very much to do here: We distribute the cell matrix and right
    // hand side contributions as we have done in almost all of the
    // other tutorial programs using the constraints objects. We then
    // also have to do the same for the face matrix contributions
    // that have gained content for the faces (interior and boundary)
    // and that the `face_worker` and `boundary_worker` have added
    // to the `copy_data.face_data` array.
    auto copier = [&](const CopyData &copy_data) {
      constraints.distribute_local_to_global(copy_data.cell_matrix,
                                             copy_data.cell_rhs,
                                             copy_data.local_dof_indices,
                                             system_matrix,
                                             system_rhs);

      for (auto &cdf : copy_data.face_data)
        {
          constraints.distribute_local_to_global(cdf.cell_matrix,
                                                 cdf.joint_dof_indices,
                                                 system_matrix);
        }
    };


    // Having set all of this up, what remains is to just create a scratch
    // and copy data object and call the MeshWorker::mesh_loop() function
    // that then goes over all cells and faces, calls the respective workers
    // on them, and then the copier function that puts things into the
    // global matrix and right hand side. As an additional benefit,
    // MeshWorker::mesh_loop() does all of this in parallel, using
    // as many processor cores as your machine happens to have.
    const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1;
    ScratchData<dim>   scratch_data(mapping,
                                  fe,
                                  n_gauss_points,
                                  update_values | update_gradients |
                                    update_hessians | update_quadrature_points |
                                    update_JxW_values,
                                  update_values | update_gradients |
                                    update_hessians | update_quadrature_points |
                                    update_JxW_values | update_normal_vectors);
    CopyData           copy_data(dof_handler.get_fe().dofs_per_cell);
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker,
                          face_worker);
  }



  // @sect4{Solving the linear system and postprocessing}
  //
  // The show is essentially over at this point: The remaining functions are
  // not overly interesting or novel. The first one simply uses a direct
  // solver to solve the linear system (see also step-29):
  template <int dim>
  void BiharmonicProblem<dim>::solve()
  {
    std::cout << "   Solving system..." << std::endl;

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);

    constraints.distribute(solution);
  }



  // The next function evaluates the error between the computed solution
  // and the exact solution (which is known here because we have chosen
  // the right hand side and boundary values in a way so that we know
  // the corresponding solution). In the first two code blocks below,
  // we compute the error in the $L_2$ norm and the $H^1$ semi-norm.
  template <int dim>
  void BiharmonicProblem<dim>::compute_errors()
  {
    {
      Vector<float> norm_per_cell(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        ExactSolution::Solution<dim>(),
                                        norm_per_cell,
                                        QGauss<dim>(fe.degree + 2),
                                        VectorTools::L2_norm);
      const double error_norm =
        VectorTools::compute_global_error(triangulation,
                                          norm_per_cell,
                                          VectorTools::L2_norm);
      std::cout << "   Error in the L2 norm       :     " << error_norm
                << std::endl;
    }

    {
      Vector<float> norm_per_cell(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        ExactSolution::Solution<dim>(),
                                        norm_per_cell,
                                        QGauss<dim>(fe.degree + 2),
                                        VectorTools::H1_seminorm);
      const double error_norm =
        VectorTools::compute_global_error(triangulation,
                                          norm_per_cell,
                                          VectorTools::H1_seminorm);
      std::cout << "   Error in the H1 seminorm       : " << error_norm
                << std::endl;
    }

    // Now also compute an approximation to the $H^2$ seminorm error. The actual
    // $H^2$ seminorm would require us to integrate second derivatives of the
    // solution $u_h$, but given the Lagrange shape functions we use, $u_h$ of
    // course has kinks at the interfaces between cells, and consequently second
    // derivatives are singular at interfaces. As a consequence, we really only
    // integrate over the interior of cells and ignore the interface
    // contributions. This is *not* an equivalent norm to the energy norm for
    // the problem, but still gives us an idea of how fast the error converges.
    //
    // We note that one could address this issue by defining a norm that
    // is equivalent to the energy norm. This would involve adding up not
    // only the integrals over cell interiors as we do below, but also adding
    // penalty terms for the jump of the derivative of $u_h$ across interfaces,
    // with an appropriate scaling of the two kinds of terms. We will leave
    // this for later work.
    {
      const QGauss<dim>            quadrature_formula(fe.degree + 2);
      ExactSolution::Solution<dim> exact_solution;
      Vector<double> error_per_cell(triangulation.n_active_cells());

      FEValues<dim> fe_values(mapping,
                              fe,
                              quadrature_formula,
                              update_values | update_hessians |
                                update_quadrature_points | update_JxW_values);

      FEValuesExtractors::Scalar scalar(0);
      const unsigned int         n_q_points = quadrature_formula.size();

      std::vector<SymmetricTensor<2, dim>> exact_hessians(n_q_points);
      std::vector<Tensor<2, dim>>          hessians(n_q_points);
      for (auto &cell : dof_handler.active_cell_iterators())
        {
          fe_values.reinit(cell);
          fe_values[scalar].get_function_hessians(solution, hessians);
          exact_solution.hessian_list(fe_values.get_quadrature_points(),
                                      exact_hessians);

          double local_error = 0;
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              local_error +=
                ((exact_hessians[q_point] - hessians[q_point]).norm_square() *
                 fe_values.JxW(q_point));
            }
          error_per_cell[cell->active_cell_index()] = std::sqrt(local_error);
        }

      const double error_norm = error_per_cell.l2_norm();
      std::cout << "   Error in the broken H2 seminorm: " << error_norm
                << std::endl;
    }
  }



  // Equally uninteresting is the function that generates graphical output.
  // It looks exactly like the one in step-6, for example.
  template <int dim>
  void
  BiharmonicProblem<dim>::output_results(const unsigned int iteration) const
  {
    std::cout << "   Writing graphical output..." << std::endl;

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();

    std::ofstream output_vtu(
      ("output_" + Utilities::int_to_string(iteration, 6) + ".vtu").c_str());
    data_out.write_vtu(output_vtu);
  }



  // The same is true for the `run()` function: Just like in previous
  // programs.
  template <int dim>
  void BiharmonicProblem<dim>::run()
  {
    make_grid();

    const unsigned int n_cycles = 4;
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        std::cout << "Cycle: " << cycle << " of " << n_cycles << std::endl;

        triangulation.refine_global(1);
        setup_system();

        assemble_system();
        solve();

        output_results(cycle);

        compute_errors();
        std::cout << std::endl;
      }
  }
} // namespace Step47



// @sect3{The main() function}
//
// Finally for the `main()` function. There is, again, not very much to see
// here: It looks like the ones in previous tutorial programs. There
// is a variable that allows selecting the polynomial degree of the element
// we want to use for solving the equation. Because the C0IP formulation
// we use requires the element degree to be at least two, we check with
// an assertion that whatever one sets for the polynomial degree actually
// makes sense.
int main()
{
  try
    {
      using namespace dealii;
      using namespace Step47;

      const unsigned int fe_degree = 2;
      Assert(fe_degree >= 2,
             ExcMessage("The C0IP formulation for the biharmonic problem "
                        "only works if one uses elements of polynomial "
                        "degree at least 2."));

      BiharmonicProblem<2> biharmonic_problem(fe_degree);
      biharmonic_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
