/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 - 2019 by the deal.II authors
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

 *      Author: Zhuoran Wang, Colorado State University, 2018
 */

// @sect3{Include files}
// This program is based on step-7, step-20 and step-51,
// so most of the following header files are familiar. We
// need the following:
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <fstream>
#include <iostream>


// Our first step, as always, is to put everything related to this tutorial
// program into its own namespace:
namespace Step61
{
  using namespace dealii;

  // @sect3{The WGDarcyEquation class template}

  // This is the main class of this program. We will solve for the numerical
  // pressure in the interior and on faces using the weak Galerkin (WG) method,
  // and calculate the $L_2$ error of pressure. In the post-processing step, we
  // will also calculate $L_2$-errors of the velocity and flux.
  //
  // The structure of the class is not fundamentally different from that of
  // previous tutorial programs, so there is little need to comment on the
  // details.
  template <int dim>
  class WGDarcyEquation
  {
  public:
    WGDarcyEquation(const unsigned int degree);
    void run();

  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void compute_velocity_errors();
    void compute_pressure_error();
    void output_results() const;

    Triangulation<dim> triangulation;

    FESystem<dim>   fe;
    DoFHandler<dim> dof_handler;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
  };



  // @sect3{Right hand side, boundary values, and exact solution}

  // Next, we define the coefficient matrix $\mathbf{K}$ (here, the
  // identity matrix), Dirichlet boundary conditions, the right-hand
  // side $f = 2\pi^2 \sin(\pi x) \sin(\pi y)$, and the exact solution
  // that corresponds to these choices for $K$ and $f$, namely $p =
  // \sin(\pi x) \sin(\pi y)$.
  template <int dim>
  class Coefficient : public TensorFunction<2, dim>
  {
  public:
    Coefficient()
      : TensorFunction<2, dim>()
    {}

    virtual void value_list(const std::vector<Point<dim>> &points,
                            std::vector<Tensor<2, dim>> &values) const override;
  };



  template <int dim>
  void Coefficient<dim>::value_list(const std::vector<Point<dim>> &points,
                                    std::vector<Tensor<2, dim>> &  values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));
    for (unsigned int p = 0; p < points.size(); ++p)
      values[p] = unit_symmetric_tensor<dim>();
  }



  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues()
      : Function<dim>(2)
    {}

    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };



  template <int dim>
  double BoundaryValues<dim>::value(const Point<dim> & /*p*/,
                                    const unsigned int /*component*/) const
  {
    return 0;
  }



  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };



  template <int dim>
  double RightHandSide<dim>::value(const Point<dim> &p,
                                   const unsigned int /*component*/) const
  {
    return (2 * numbers::PI * numbers::PI * std::sin(numbers::PI * p[0]) *
            std::sin(numbers::PI * p[1]));
  }



  // The class that implements the exact pressure solution has an
  // oddity in that we implement it as a vector-valued one with two
  // components. (We say that it has two components in the constructor
  // where we call the constructor of the base Function class.) In the
  // `value()` function, we do not test for the value of the
  // `component` argument, which implies that we return the same value
  // for both components of the vector-valued function. We do this
  // because we describe the finite element in use in this program as
  // a vector-valued system that contains the interior and the
  // interface pressures, and when we compute errors, we will want to
  // use the same pressure solution to test both of these components.
  template <int dim>
  class ExactPressure : public Function<dim>
  {
  public:
    ExactPressure()
      : Function<dim>(2)
    {}

    virtual double value(const Point<dim> & p,
                         const unsigned int component) const override;
  };



  template <int dim>
  double ExactPressure<dim>::value(const Point<dim> &p,
                                   const unsigned int /*component*/) const
  {
    return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]);
  }



  template <int dim>
  class ExactVelocity : public TensorFunction<1, dim>
  {
  public:
    ExactVelocity()
      : TensorFunction<1, dim>()
    {}

    virtual Tensor<1, dim> value(const Point<dim> &p) const override;
  };



  template <int dim>
  Tensor<1, dim> ExactVelocity<dim>::value(const Point<dim> &p) const
  {
    Tensor<1, dim> return_value;
    return_value[0] = -numbers::PI * std::cos(numbers::PI * p[0]) *
                      std::sin(numbers::PI * p[1]);
    return_value[1] = -numbers::PI * std::sin(numbers::PI * p[0]) *
                      std::cos(numbers::PI * p[1]);
    return return_value;
  }



  // @sect3{WGDarcyEquation class implementation}

  // @sect4{WGDarcyEquation::WGDarcyEquation}

  // In this constructor, we create a finite element space for vector valued
  // functions, which will here include the ones used for the interior and
  // interface pressures, $p^\circ$ and $p^\partial$.
  template <int dim>
  WGDarcyEquation<dim>::WGDarcyEquation(const unsigned int degree)
    : fe(FE_DGQ<dim>(degree), 1, FE_FaceQ<dim>(degree), 1)
    , dof_handler(triangulation)

  {}



  // @sect4{WGDarcyEquation::make_grid}

  // We generate a mesh on the unit square domain and refine it.
  template <int dim>
  void WGDarcyEquation<dim>::make_grid()
  {
    GridGenerator::hyper_cube(triangulation, 0, 1);
    triangulation.refine_global(2);

    std::cout << "   Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "   Total number of cells: " << triangulation.n_cells()
              << std::endl;
  }



  // @sect4{WGDarcyEquation::setup_system}

  // After we have created the mesh above, we distribute degrees of
  // freedom and resize matrices and vectors. The only piece of
  // interest in this function is how we interpolate the boundary
  // values for the pressure. Since the pressure consists of interior
  // and interface components, we need to make sure that we only
  // interpolate onto that component of the vector-valued solution
  // space that corresponds to the interface pressures (as these are
  // the only ones that are defined on the boundary of the domain). We
  // do this via a component mask object for only the interface
  // pressures.
  template <int dim>
  void WGDarcyEquation<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    std::cout << "   Number of pressure degrees of freedom: "
              << dof_handler.n_dofs() << std::endl;

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    {
      constraints.clear();
      const FEValuesExtractors::Scalar interface_pressure(1);
      const ComponentMask              interface_pressure_mask =
        fe.component_mask(interface_pressure);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               BoundaryValues<dim>(),
                                               constraints,
                                               interface_pressure_mask);
      constraints.close();
    }


    // In the bilinear form, there is no integration term over faces
    // between two neighboring cells, so we can just use
    // <code>DoFTools::make_sparsity_pattern</code> to calculate the sparse
    // matrix.
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
  }



  // @sect4{WGDarcyEquation::assemble_system}

  // This function is more interesting. As detailed in the
  // introduction, the assembly of the linear system requires us to
  // evaluate the weak gradient of the shape functions, which is an
  // element in the Raviart-Thomas space. As a consequence, we need to
  // define a Raviart-Thomas finite element object, and have FEValues
  // objects that evaluate it at quadrature points. We then need to
  // compute the matrix $C^K$ on every cell $K$, for which we need the
  // matrices $M^K$ and $G^K$ mentioned in the introduction.
  //
  // A point that may not be obvious is that in all previous tutorial
  // programs, we have always called FEValues::reinit() with a cell
  // iterator from a DoFHandler. This is so that one can call
  // functions such as FEValuesBase::get_function_values() that
  // extract the values of a finite element function (represented by a
  // vector of DoF values) on the quadrature points of a cell. For
  // this operation to work, one needs to know which vector elements
  // correspond to the degrees of freedom on a given cell -- i.e.,
  // exactly the kind of information and operation provided by the
  // DoFHandler class.
  //
  // On the other hand, we don't have such a DoFHandler object for the
  // Raviart-Thomas space in this program. In fact, we don't even have
  // an element that can represent the "broken" Raviart-Thomas space
  // we really want to use here (i.e., the restriction of the
  // Raviart-Thomas shape functions to individual cells, without the
  // need for any kind of continuity across cell interfaces). We solve
  // this conundrum by using the fact that one can call
  // FEValues::reinit() also with cell iterators into Triangulation
  // objects (rather than DoFHandler objects). In this case, FEValues
  // can of course only provide us with information that only
  // references information of cells, rather than degrees of freedom
  // enumerated on these cells. So we can't use
  // FEValuesBase::get_function_values(), but we can use
  // FEValues::shape_value() to obtain the values of shape functions
  // at quadrature points on the current cell. It is this kind of
  // functionality we will make use of below.
  //
  // Given this introduction, the following declarations should be
  // pretty obvious:
  template <int dim>
  void WGDarcyEquation<dim>::assemble_system()
  {
    const FE_RaviartThomas<dim> fe_rt(fe.base_element(0).degree);

    const QGauss<dim>     quadrature_formula(fe_rt.degree + 1);
    const QGauss<dim - 1> face_quadrature_formula(fe_rt.degree + 1);

    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values);
    FEFaceValues<dim> fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                       update_quadrature_points |
                                       update_JxW_values);

    FEValues<dim>     fe_values_rt(fe_rt,
                               quadrature_formula,
                               update_values | update_gradients |
                                 update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values_rt(fe_rt,
                                        face_quadrature_formula,
                                        update_values | update_normal_vectors |
                                          update_quadrature_points |
                                          update_JxW_values);

    const unsigned int dofs_per_cell    = fe.dofs_per_cell;
    const unsigned int dofs_per_cell_rt = fe_rt.dofs_per_cell;

    const unsigned int n_q_points    = fe_values.get_quadrature().size();
    const unsigned int n_q_points_rt = fe_values_rt.get_quadrature().size();

    const unsigned int n_face_q_points = fe_face_values.get_quadrature().size();

    const RightHandSide<dim> right_hand_side;
    std::vector<double>      right_hand_side_values(n_q_points);

    const Coefficient<dim>      coefficient;
    std::vector<Tensor<2, dim>> coefficient_values(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


    // Next, let us declare the various cell matrices discussed in the
    // introduction:
    FullMatrix<double> cell_matrix_M(dofs_per_cell_rt, dofs_per_cell_rt);
    FullMatrix<double> cell_matrix_G(dofs_per_cell_rt, dofs_per_cell);
    FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_rt);
    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    Vector<double>     cell_solution(dofs_per_cell);

    // We need <code>FEValuesExtractors</code> to access the @p interior and
    // @p face component of the shape functions.
    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar interior(0);
    const FEValuesExtractors::Scalar face(1);

    // This finally gets us in position to loop over all cells. On
    // each cell, we will first calculate the various cell matrices
    // used to construct the local matrix -- as they depend on the
    // cell in question, they need to be re-computed on each cell. We
    // need shape functions for the Raviart-Thomas space as well, for
    // which we need to create first an iterator to the cell of the
    // triangulation, which we can obtain by assignment from the cell
    // pointing into the DoFHandler.
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        const typename Triangulation<dim>::active_cell_iterator cell_rt = cell;
        fe_values_rt.reinit(cell_rt);

        right_hand_side.value_list(fe_values.get_quadrature_points(),
                                   right_hand_side_values);
        coefficient.value_list(fe_values.get_quadrature_points(),
                               coefficient_values);

        // The first cell matrix we will compute is the mass matrix
        // for the Raviart-Thomas space.  Hence, we need to loop over
        // all the quadrature points for the velocity FEValues object.
        cell_matrix_M = 0;
        for (unsigned int q = 0; q < n_q_points_rt; ++q)
          for (unsigned int i = 0; i < dofs_per_cell_rt; ++i)
            {
              const Tensor<1, dim> v_i = fe_values_rt[velocities].value(i, q);
              for (unsigned int k = 0; k < dofs_per_cell_rt; ++k)
                {
                  const Tensor<1, dim> v_k =
                    fe_values_rt[velocities].value(k, q);
                  cell_matrix_M(i, k) += (v_i * v_k * fe_values_rt.JxW(q));
                }
            }
        // Next we take the inverse of this matrix by using
        // FullMatrix::gauss_jordan(). It will be used to calculate
        // the coefficient matrix $C^K$ later. It is worth recalling
        // later that `cell_matrix_M` actually contains the *inverse*
        // of $M^K$ after this call.
        cell_matrix_M.gauss_jordan();

        // From the introduction, we know that the right hand side
        // $G^K$ of the equation that defines $C^K$ is the difference
        // between a face integral and a cell integral. Here, we
        // approximate the negative of the contribution in the
        // interior. Each component of this matrix is the integral of
        // a product between a basis function of the polynomial space
        // and the divergence of a basis function of the
        // Raviart-Thomas space. These basis functions are defined in
        // the interior.
        cell_matrix_G = 0;
        for (unsigned int q = 0; q < n_q_points; ++q)
          for (unsigned int i = 0; i < dofs_per_cell_rt; ++i)
            {
              const double div_v_i = fe_values_rt[velocities].divergence(i, q);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  const double phi_j_interior = fe_values[interior].value(j, q);

                  cell_matrix_G(i, j) -=
                    (div_v_i * phi_j_interior * fe_values.JxW(q));
                }
            }


        // Next, we approximate the integral on faces by quadrature.
        // Each component is the integral of a product between a basis function
        // of the polynomial space and the dot product of a basis function of
        // the Raviart-Thomas space and the normal vector. So we loop over all
        // the faces of the element and obtain the normal vector.
        for (unsigned int face_n = 0;
             face_n < GeometryInfo<dim>::faces_per_cell;
             ++face_n)
          {
            fe_face_values.reinit(cell, face_n);
            fe_face_values_rt.reinit(cell_rt, face_n);

            for (unsigned int q = 0; q < n_face_q_points; ++q)
              {
                const Tensor<1, dim> normal = fe_face_values.normal_vector(q);

                for (unsigned int i = 0; i < dofs_per_cell_rt; ++i)
                  {
                    const Tensor<1, dim> v_i =
                      fe_face_values_rt[velocities].value(i, q);
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      {
                        const double phi_j_face =
                          fe_face_values[face].value(j, q);

                        cell_matrix_G(i, j) +=
                          ((v_i * normal) * phi_j_face * fe_face_values.JxW(q));
                      }
                  }
              }
          }

        // @p cell_matrix_C is then the matrix product between the
        // transpose of $G^K$ and the inverse of the mass matrix
        // (where this inverse is stored in @p cell_matrix_M):
        cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M);

        // Finally we can compute the local matrix $A^K$.  Element
        // $A^K_{ij}$ is given by $\int_{E} \sum_{k,l} C_{ik} C_{jl}
        // (\mathbf{K} \mathbf{v}_k) \cdot \mathbf{v}_l
        // \mathrm{d}x$. We have calculated the coefficients $C$ in
        // the previous step, and so obtain the following after
        // suitably re-arranging the loops:
        local_matrix = 0;
        for (unsigned int q = 0; q < n_q_points_rt; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell_rt; ++k)
              {
                const Tensor<1, dim> v_k = fe_values_rt[velocities].value(k, q);
                for (unsigned int l = 0; l < dofs_per_cell_rt; ++l)
                  {
                    const Tensor<1, dim> v_l =
                      fe_values_rt[velocities].value(l, q);

                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        local_matrix(i, j) +=
                          (coefficient_values[q] * cell_matrix_C[i][k] * v_k) *
                          cell_matrix_C[j][l] * v_l * fe_values_rt.JxW(q);
                  }
              }
          }

        // Next, we calculate the right hand side, $\int_{K} f q \mathrm{d}x$:
        cell_rhs = 0;
        for (unsigned int q = 0; q < n_q_points; ++q)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) += (fe_values[interior].value(i, q) *
                              right_hand_side_values[q] * fe_values.JxW(q));
            }

        // The last step is to distribute components of the local
        // matrix into the system matrix and transfer components of
        // the cell right hand side into the system right hand side:
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          local_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }



  // @sect4{WGDarcyEquation<dim>::solve}

  // This step is rather trivial and the same as in many previous
  // tutorial programs:
  template <int dim>
  void WGDarcyEquation<dim>::solve()
  {
    SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
    SolverCG<>    solver(solver_control);
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
    constraints.distribute(solution);
  }



  // @sect4{WGDarcyEquation<dim>::compute_pressure_error}

  // This part is to calculate the $L_2$ error of the pressure.  We
  // define a vector that holds the norm of the error on each cell.
  // Next, we use VectorTool::integrate_difference() to compute the
  // error in the $L_2$ norm on each cell. However, we really only
  // care about the error in the interior component of the solution
  // vector (we can't even evaluate the interface pressures at the
  // quadrature points because these are all located in the interior
  // of cells) and consequently have to use a weight function that
  // ensures that the interface component of the solution variable is
  // ignored. This is done by using the ComponentSelectFunction whose
  // arguments indicate which component we want to select (component
  // zero, i.e., the interior pressures) and how many components there
  // are in total (two).
  template <int dim>
  void WGDarcyEquation<dim>::compute_pressure_error()
  {
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    const ComponentSelectFunction<dim> select_interior_pressure(0, 2);
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      ExactPressure<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 2),
                                      VectorTools::L2_norm,
                                      &select_interior_pressure);

    const double L2_error = difference_per_cell.l2_norm();
    std::cout << "L2_error_pressure " << L2_error << std::endl;
  }



  // @sect4{WGDarcyEquation<dim>::compute_velocity_errors}

  // In this function, we evaluate $L_2$ errors for the velocity on
  // each cell, and $L_2$ errors for the flux on faces.

  // We are going to evaluate velocities on each cell and calculate
  // the difference between numerical and exact velocities. The
  // velocity is defined as $\mathbf{u}_h = \mathbf{Q}_h \left(
  // -\mathbf{K}\nabla_{w,d}p_h \right)$, which requires us to compute
  // many of the same terms as in the assembly of the system matrix.
  // There are also the matrices $E^K,D^K$ we need to assemble (see
  // the introduction) but they really just follow the same kind of
  // pattern.
  //
  // Computing the same matrices here as we have already done in the
  // `assemble_system()` function is of course wasteful in terms of
  // CPU time. Likewise, we copy some of the code from there to this
  // function, and this is also generally a poor idea. A better
  // implementation might provide for a function that encapsulates
  // this duplicated code. One could also think of using the classic
  // trade-off between computing efficiency and memory efficiency to
  // only compute the $C^K$ matrices once per cell during the
  // assembly, storing them somewhere on the side, and re-using them
  // here. (This is what step-51 does, for example, where the
  // `assemble_system()` function takes an argument that determines
  // whether the local matrices are recomputed, and a similar approach
  // -- maybe with storing local matrices elsewhere -- could be
  // adapted for the current program.)
  template <int dim>
  void WGDarcyEquation<dim>::compute_velocity_errors()
  {
    const FE_RaviartThomas<dim> fe_rt(fe.base_element(0).degree);

    const QGauss<dim>     quadrature_formula(fe_rt.degree + 1);
    const QGauss<dim - 1> face_quadrature_formula(fe_rt.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values);

    FEFaceValues<dim> fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                       update_quadrature_points |
                                       update_JxW_values);

    FEValues<dim> fe_values_rt(fe_rt,
                               quadrature_formula,
                               update_values | update_gradients |
                                 update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values_rt(fe_rt,
                                        face_quadrature_formula,
                                        update_values | update_normal_vectors |
                                          update_quadrature_points |
                                          update_JxW_values);

    const unsigned int dofs_per_cell    = fe.dofs_per_cell;
    const unsigned int dofs_per_cell_rt = fe_rt.dofs_per_cell;

    const unsigned int n_q_points    = fe_values.get_quadrature().size();
    const unsigned int n_q_points_rt = fe_values_rt.get_quadrature().size();

    const unsigned int n_face_q_points = fe_face_values.get_quadrature().size();
    const unsigned int n_face_q_points_rt =
      fe_face_values_rt.get_quadrature().size();


    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    FullMatrix<double> cell_matrix_M(dofs_per_cell_rt, dofs_per_cell_rt);
    FullMatrix<double> cell_matrix_G(dofs_per_cell_rt, dofs_per_cell);
    FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_rt);

    FullMatrix<double> cell_matrix_D(dofs_per_cell_rt, dofs_per_cell_rt);
    FullMatrix<double> cell_matrix_E(dofs_per_cell_rt, dofs_per_cell_rt);

    Vector<double> cell_solution(dofs_per_cell);
    Vector<double> cell_velocity(dofs_per_cell_rt);

    double L2_err_velocity_cell_sqr_global = 0;
    double L2_err_flux_sqr                 = 0;

    const Coefficient<dim>      coefficient;
    std::vector<Tensor<2, dim>> coefficient_values(n_q_points_rt);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);
    const FEValuesExtractors::Scalar interior(0);
    const FEValuesExtractors::Scalar face(1);

    const ExactVelocity<dim> exact_velocity;

    // In the loop over all cells, we will calculate $L_2$ errors of velocity
    // and flux.

    // First, we calculate the $L_2$ velocity error.
    // In the introduction, we explained how to calculate the numerical velocity
    // on the cell. We need the pressure solution values on each cell,
    // coefficients of the Gram matrix and coefficients of the $L_2$ projection.
    // We have already calculated the global solution, so we will extract the
    // cell solution from the global solution. The coefficients of the Gram
    // matrix have been calculated when we assembled the system matrix for the
    // pressures. We will do the same way here. For the coefficients of the
    // projection, we do matrix multiplication, i.e., the inverse of the Gram
    // matrix times the matrix with $(\mathbf{K} \mathbf{w}, \mathbf{w})$ as
    // components. Then, we multiply all these coefficients and call them beta.
    // The numerical velocity is the product of beta and the basis functions of
    // the Raviart-Thomas space.
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        const typename Triangulation<dim>::active_cell_iterator cell_rt = cell;
        fe_values_rt.reinit(cell_rt);

        coefficient.value_list(fe_values_rt.get_quadrature_points(),
                               coefficient_values);

        // The component of this <code>cell_matrix_E</code> is the integral of
        // $(\mathbf{K} \mathbf{w}, \mathbf{w})$. <code>cell_matrix_M</code> is
        // the Gram matrix.
        cell_matrix_M = 0;
        cell_matrix_E = 0;
        for (unsigned int q = 0; q < n_q_points_rt; ++q)
          for (unsigned int i = 0; i < dofs_per_cell_rt; ++i)
            {
              const Tensor<1, dim> v_i = fe_values_rt[velocities].value(i, q);
              for (unsigned int k = 0; k < dofs_per_cell_rt; ++k)
                {
                  const Tensor<1, dim> v_k =
                    fe_values_rt[velocities].value(k, q);

                  cell_matrix_E(i, k) +=
                    (coefficient_values[q] * v_i * v_k * fe_values_rt.JxW(q));

                  cell_matrix_M(i, k) += (v_i * v_k * fe_values_rt.JxW(q));
                }
            }

        // To compute the matrix $D$ mentioned in the introduction, we
        // then need to evaluate $D=M^{-1}E$ as explained in the
        // introduction:
        cell_matrix_M.gauss_jordan();
        cell_matrix_M.mmult(cell_matrix_D, cell_matrix_E);

        // Then we also need, again, to compute the matrix $C$ that is
        // used to evaluate the weak discrete gradient. This is the
        // exact same code as used in the assembly of the system
        // matrix, so we just copy it from there:
        cell_matrix_G = 0;
        for (unsigned int q = 0; q < n_q_points; ++q)
          for (unsigned int i = 0; i < dofs_per_cell_rt; ++i)
            {
              const double div_v_i = fe_values_rt[velocities].divergence(i, q);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  const double phi_j_interior = fe_values[interior].value(j, q);

                  cell_matrix_G(i, j) -=
                    (div_v_i * phi_j_interior * fe_values.JxW(q));
                }
            }

        for (unsigned int face_n = 0;
             face_n < GeometryInfo<dim>::faces_per_cell;
             ++face_n)
          {
            fe_face_values.reinit(cell, face_n);
            fe_face_values_rt.reinit(cell_rt, face_n);

            for (unsigned int q = 0; q < n_face_q_points; ++q)
              {
                const Tensor<1, dim> normal = fe_face_values.normal_vector(q);

                for (unsigned int i = 0; i < dofs_per_cell_rt; ++i)
                  {
                    const Tensor<1, dim> v_i =
                      fe_face_values_rt[velocities].value(i, q);
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      {
                        const double phi_j_face =
                          fe_face_values[face].value(j, q);

                        cell_matrix_G(i, j) +=
                          ((v_i * normal) * phi_j_face * fe_face_values.JxW(q));
                      }
                  }
              }
          }
        cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M);

        // Finally, we need to extract the pressure unknowns that
        // correspond to the current cell:
        cell->get_dof_values(solution, cell_solution);

        // We are now in a position to compute the local velocity
        // unknowns (with respect to the Raviart-Thomas space we are
        // projecting the term $-\mathbf K \nabla_{w,d} p_h$ into):
        cell_velocity = 0;
        for (unsigned int k = 0; k < dofs_per_cell_rt; ++k)
          for (unsigned int j = 0; j < dofs_per_cell_rt; ++j)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_velocity(k) +=
                -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j));

        // Now, we can calculate the numerical velocity at each quadrature point
        // and compute the $L_2$ error on each cell.
        double L2_err_velocity_cell_sqr_local = 0;
        for (unsigned int q = 0; q < n_q_points_rt; ++q)
          {
            Tensor<1, dim> velocity;
            for (unsigned int k = 0; k < dofs_per_cell_rt; ++k)
              {
                const Tensor<1, dim> phi_k_u =
                  fe_values_rt[velocities].value(k, q);
                velocity += cell_velocity(k) * phi_k_u;
              }

            const Tensor<1, dim> true_velocity =
              exact_velocity.value(fe_values_rt.quadrature_point(q));

            L2_err_velocity_cell_sqr_local +=
              ((velocity - true_velocity) * (velocity - true_velocity) *
               fe_values_rt.JxW(q));
          }
        L2_err_velocity_cell_sqr_global += L2_err_velocity_cell_sqr_local;

        // For reconstructing the flux we need the size of cells and
        // faces.  Since fluxes are calculated on faces, we have the
        // loop over all four faces of each cell. To calculate the
        // face velocity, we use the coefficients `cell_velocity` we
        // have computed previously. Then, we calculate the squared
        // velocity error in normal direction. Finally, we calculate
        // the $L_2$ flux error on the cell and add it to the global
        // error.
        const double cell_area = cell->measure();
        for (unsigned int face_n = 0;
             face_n < GeometryInfo<dim>::faces_per_cell;
             ++face_n)
          {
            const double face_length = cell->face(face_n)->measure();
            fe_face_values.reinit(cell, face_n);
            fe_face_values_rt.reinit(cell_rt, face_n);

            double L2_err_flux_face_sqr_local = 0;
            for (unsigned int q = 0; q < n_face_q_points_rt; ++q)
              {
                Tensor<1, dim> velocity;
                for (unsigned int k = 0; k < dofs_per_cell_rt; ++k)
                  {
                    const Tensor<1, dim> phi_k_u =
                      fe_face_values_rt[velocities].value(k, q);
                    velocity += cell_velocity(k) * phi_k_u;
                  }
                const Tensor<1, dim> true_velocity =
                  exact_velocity.value(fe_face_values_rt.quadrature_point(q));

                const Tensor<1, dim> normal = fe_face_values.normal_vector(q);

                L2_err_flux_face_sqr_local +=
                  ((velocity * normal - true_velocity * normal) *
                   (velocity * normal - true_velocity * normal) *
                   fe_face_values_rt.JxW(q));
              }
            const double err_flux_each_face =
              L2_err_flux_face_sqr_local / (face_length) * (cell_area);
            L2_err_flux_sqr += err_flux_each_face;
          }
      }

    // After adding up errors over all cells and faces, we take the
    // square root and get the $L_2$ errors of velocity and
    // flux. These we output to screen.
    const double L2_err_velocity_cell =
      std::sqrt(L2_err_velocity_cell_sqr_global);
    const double L2_err_flux_face = std::sqrt(L2_err_flux_sqr);

    std::cout << "L2_error_vel:  " << L2_err_velocity_cell << std::endl
              << "L2_error_flux: " << L2_err_flux_face << std::endl;
  }


  // @sect4{WGDarcyEquation::output_results}

  // We have two sets of results to output: the interior solution and
  // the skeleton solution. We use <code>DataOut</code> to visualize
  // interior results. The graphical output for the skeleton results
  // is done by using the DataOutFaces class.
  //
  // In both of the output files, both the interior and the face
  // variables are stored. For the interface output, the output file
  // simply contains the interpolation of the interior pressures onto
  // the faces, but because it is undefined which of the two interior
  // pressure variables you get from the two adjacent cells, it is
  // best to ignore the interior pressure in the interface output
  // file. Conversely, for the cell interior output file, it is of
  // course impossible to show any interface pressures $p^\partial$,
  // because these are only available on interfaces and not cell
  // interiors. Consequently, you will see them shown as an invalid
  // value (such as an infinity).
  template <int dim>
  void WGDarcyEquation<dim>::output_results() const
  {
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "Pressure_Interior");
      data_out.build_patches(fe.degree);
      std::ofstream output("Pressure_Interior.vtu");
      data_out.write_vtu(output);
    }

    {
      DataOutFaces<dim> data_out_faces(false);
      data_out_faces.attach_dof_handler(dof_handler);
      data_out_faces.add_data_vector(solution, "Pressure_Face");
      data_out_faces.build_patches(fe.degree);
      std::ofstream face_output("Pressure_Face.vtu");
      data_out_faces.write_vtu(face_output);
    }
  }


  // @sect4{WGDarcyEquation::run}

  // This is the final function of the main class. It calls the other functions
  // of our class.
  template <int dim>
  void WGDarcyEquation<dim>::run()
  {
    std::cout << "Solving problem in " << dim << " space dimensions."
              << std::endl;
    make_grid();
    setup_system();
    assemble_system();
    solve();
    compute_pressure_error();
    compute_velocity_errors();
    output_results();
  }

} // namespace Step61


// @sect3{The <code>main</code> function}

// This is the main function. We can change the dimension here to run in 3d.
int main()
{
  try
    {
      dealii::deallog.depth_console(2);
      Step61::WGDarcyEquation<2> wg_darcy(0);
      wg_darcy.run();
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
