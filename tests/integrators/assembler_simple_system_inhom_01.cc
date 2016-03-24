// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2016 by the deal.II authors
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

/*
 * Test SystemSimple with mixed homogeneous and inhomogeneous constraints.
 * Compare SystemSimple with matrix and rhs computed from
 * MatrixSimple, ResidualSimple with homogeneous constraints
 * and global constraints.condense(matrix,rhs) operation to
 * incorporate the inhomogeneous constraints.
 * Tests are run for
 * -> FE_Q<2>(1), FE_DGQ<2>(1)
 * -> mesh: adaptive with hanging nodes
 * -> homogeneous hanging nodes constraints
 * -> inhomogeneous boundary constraints plus constraint of all degrees
 *    at the origin
 * -> matrix is the (nonsymmetric) advection matrix plus the jump of the
 *    normal flux times the average of the values along inner edges.
 * -> vector is standard rhs plus neumann boundary integral.
 *
 * Output: norm of the difference of the matrices and vectors
 */

#include "../tests.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/loop.h>

#include <iostream>
#include <fstream>

#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/base/tensor_function.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/numerics/data_out.h>


using namespace dealii;


/*
 * Right hand side
 */
template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
                                  const unsigned int component) const
{
  Assert (component == 0, ExcNotImplemented());

  double val = 1; // f = 1
  return val;
}


/*
 * Boundary values
 */
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
                                   const unsigned int component) const
{
  Assert (component == 0, ExcNotImplemented());

  double val = 1; // u_D = 1
  return val;
}


/*
 * Flux boundary values
 */
template <int dim>
class FluxBoundaryValues : public Function<dim>
{
public:
  FluxBoundaryValues () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <int dim>
double FluxBoundaryValues<dim>::value (const Point<dim> &p,
                                       const unsigned int component) const
{
  double val = 1; // g = 1
  return val;
}


/*
 * Advection coefficient
 */
template <int dim>
class Advection : public TensorFunction<1,dim>
{
public:
  Advection () : TensorFunction<1,dim> () {}

  virtual Tensor<1,dim> value (const Point<dim> &p) const;
};

template <int dim>
Tensor<1,dim>
Advection<dim>::value (const Point<dim> &p) const
{
  Point<dim> value; // beta = (1,0)^t
  value[0] = 1;
  value[1] = 0;

  return value;
}


/*
 * System integrator
 */
template <int dim>
class SystemIntegrator : public MeshWorker::LocalIntegrator<dim>
{
public:
  void cell(MeshWorker::DoFInfo<dim> &dinfo,
            typename MeshWorker::IntegrationInfo<dim> &info) const;
  void boundary(MeshWorker::DoFInfo<dim> &dinfo,
                typename MeshWorker::IntegrationInfo<dim> &info) const;
  void face(MeshWorker::DoFInfo<dim> &dinfo1,
            MeshWorker::DoFInfo<dim> &dinfo2,
            typename MeshWorker::IntegrationInfo<dim> &info1,
            typename MeshWorker::IntegrationInfo<dim> &info2) const;
};

template <int dim>
void SystemIntegrator<dim>::cell(
  MeshWorker::DoFInfo<dim> &dinfo,
  typename MeshWorker::IntegrationInfo<dim> &info) const
{
  // Matrix
  const FEValuesBase<dim> &fe = info.fe_values();

  FullMatrix<double> &local_matrix = dinfo.matrix(0).matrix;
  const unsigned int n_points = fe.n_quadrature_points;
  const unsigned int n_dofs = fe.dofs_per_cell;
  const Advection<dim> advection;

  for (unsigned int k = 0; k < n_points; ++k)
    {
      const Tensor<1,dim> advection_directions = advection.value(fe.quadrature_point(k));
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
          local_matrix(i, j) += ( fe.shape_grad(i,k)* fe.shape_grad(j,k)
                                  + fe.shape_value(i,k)*(advection_directions*fe.shape_grad(j,k))
                                ) * fe.JxW(k);
    }

  // RHS
  Vector<double> &b = dinfo.vector(0).block(0);
  const RightHandSide<dim>  right_hand_side;

  for (unsigned int k=0; k<n_points; ++k)
    {
      const double fval = right_hand_side.value(fe.quadrature_point(k));
      for (unsigned int i=0; i<n_dofs; ++i)
        b(i) += fe.JxW(k)*(fe.shape_value(i,k) * fval);
    }
}


template <int dim>
void SystemIntegrator<dim>::boundary(
  MeshWorker::DoFInfo<dim> &dinfo,
  typename MeshWorker::IntegrationInfo<dim> &info) const
{
  const FEValuesBase<dim> &fe = info.fe_values();
  const unsigned int n_points = fe.n_quadrature_points;
  const unsigned int n_dofs = fe.dofs_per_cell;
  Vector<double> &b = dinfo.vector(0).block(0);
  const FluxBoundaryValues<dim>  flux_bd;

  for (unsigned int k=0; k<n_points; ++k)
    {
      const double flux_bd_value = flux_bd.value(fe.quadrature_point(k));
      for (unsigned int i=0; i<n_dofs; ++i)
        b(i) += fe.JxW(k)*(fe.shape_value(i,k) * flux_bd_value);
    }
}


template <int dim>
void SystemIntegrator<dim>::face(
  MeshWorker::DoFInfo<dim> &dinfo1,
  MeshWorker::DoFInfo<dim> &dinfo2,
  typename MeshWorker::IntegrationInfo<dim> &info1,
  typename MeshWorker::IntegrationInfo<dim> &info2) const
{
  FullMatrix<double> &A11 = dinfo1.matrix(0,false).matrix;
  FullMatrix<double> &A12 = dinfo1.matrix(0,true).matrix;
  FullMatrix<double> &A21 = dinfo2.matrix(0,true).matrix;
  FullMatrix<double> &A22 = dinfo2.matrix(0,false).matrix;
  const FEValuesBase<dim> &fe1 = info1.fe_values(0);
  const FEValuesBase<dim> &fe2 = info2.fe_values(0);
  const unsigned int n_points = fe1.n_quadrature_points;
  const unsigned int n_dofs = fe1.dofs_per_cell;
  double h_e;
  if (dim==3) h_e = std::sqrt(dinfo1.face->measure());
  else h_e = dinfo1.face->measure();

  for (unsigned int k=0; k<n_points; ++k)
    {
      const double dx = fe1.JxW(k);
      const Point<dim> &n = (Point<dim>)fe1.normal_vector(k);
      for (unsigned int i=0; i<n_dofs; ++i)
        {
          // average
          const double value1i = 0.5*fe1.shape_value(i,k);
          const double value2i = 0.5*fe2.shape_value(i,k);

          for (unsigned int j=0; j<n_dofs; ++j)
            {
              // normal jump
              const double grad1j = -n * fe1.shape_grad(j,k);
              const double grad2j = n * fe2.shape_grad(j,k);

              A11(i,j) += fe1.JxW(k)*( value1i*grad1j );
              A12(i,j) += fe1.JxW(k)*( value1i*grad2j );
              A21(i,j) += fe1.JxW(k)*( value2i*grad1j );
              A22(i,j) += fe1.JxW(k)*( value2i*grad2j );
            }
        }
    }
}


/*
 * Matrix integrator
 */
template <int dim>
class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim>
{
public:
  void cell(MeshWorker::DoFInfo<dim> &dinfo,
            typename MeshWorker::IntegrationInfo<dim> &info) const;
  void boundary(MeshWorker::DoFInfo<dim> &dinfo,
                typename MeshWorker::IntegrationInfo<dim> &info) const;
  void face(MeshWorker::DoFInfo<dim> &dinfo1,
            MeshWorker::DoFInfo<dim> &dinfo2,
            typename MeshWorker::IntegrationInfo<dim> &info1,
            typename MeshWorker::IntegrationInfo<dim> &info2) const;
};

template <int dim>
void MatrixIntegrator<dim>::cell(
  MeshWorker::DoFInfo<dim> &dinfo,
  typename MeshWorker::IntegrationInfo<dim> &info) const
{
  const FEValuesBase<dim> &fe = info.fe_values();
  FullMatrix<double> &local_matrix = dinfo.matrix(0).matrix;

  const unsigned int n_points = fe.n_quadrature_points;
  const unsigned int n_dofs = fe.dofs_per_cell;
  const Advection<dim> advection;

  for (unsigned int k = 0; k < n_points; ++k)
    {
      const Tensor<1,dim> advection_directions = advection.value(fe.quadrature_point(k));
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
          local_matrix(i, j) += ( fe.shape_grad(i,k)* fe.shape_grad(j,k)
                                  + fe.shape_value(i,k)*(advection_directions*fe.shape_grad(j,k))
                                ) * fe.JxW(k);
    }
}


template <int dim>
void MatrixIntegrator<dim>::boundary(
  MeshWorker::DoFInfo<dim> &dinfo,
  typename MeshWorker::IntegrationInfo<dim> &info) const
{
}


template <int dim>
void MatrixIntegrator<dim>::face(
  MeshWorker::DoFInfo<dim> &dinfo1,
  MeshWorker::DoFInfo<dim> &dinfo2,
  typename MeshWorker::IntegrationInfo<dim> &info1,
  typename MeshWorker::IntegrationInfo<dim> &info2) const
{
  FullMatrix<double> &A11 = dinfo1.matrix(0,false).matrix;
  FullMatrix<double> &A12 = dinfo1.matrix(0,true).matrix;
  FullMatrix<double> &A21 = dinfo2.matrix(0,true).matrix;
  FullMatrix<double> &A22 = dinfo2.matrix(0,false).matrix;
  const FEValuesBase<dim> &fe1 = info1.fe_values(0);
  const FEValuesBase<dim> &fe2 = info2.fe_values(0);
  const unsigned int n_points = fe1.n_quadrature_points;
  const unsigned int n_dofs = fe1.dofs_per_cell;
  double h_e;
  if (dim==3) h_e = std::sqrt(dinfo1.face->measure());
  else h_e = dinfo1.face->measure();

  for (unsigned int k=0; k<n_points; ++k)
    {
      const double dx = fe1.JxW(k);
      const Point<dim> &n = (Point<dim>)fe1.normal_vector(k);
      for (unsigned int i=0; i<n_dofs; ++i)
        {
          // average
          const double value1i = 0.5*fe1.shape_value(i,k);
          const double value2i = 0.5*fe2.shape_value(i,k);

          for (unsigned int j=0; j<n_dofs; ++j)
            {
              // normal jump
              const double grad1j = -n * fe1.shape_grad(j,k);
              const double grad2j = n * fe2.shape_grad(j,k);

              A11(i,j) += fe1.JxW(k)*( value1i*grad1j );
              A12(i,j) += fe1.JxW(k)*( value1i*grad2j );
              A21(i,j) += fe1.JxW(k)*( value2i*grad1j );
              A22(i,j) += fe1.JxW(k)*( value2i*grad2j );
            }
        }
    }
}


/*
 * Vector integrator
 */
template <int dim>
class RHSIntegrator : public MeshWorker::LocalIntegrator<dim>
{
public:
  void cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
  void boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
  void face(MeshWorker::DoFInfo<dim> &dinfo1,
            MeshWorker::DoFInfo<dim> &dinfo2,
            typename MeshWorker::IntegrationInfo<dim> &info1,
            typename MeshWorker::IntegrationInfo<dim> &info2) const;
};


template <int dim>
void RHSIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
{
  const FEValuesBase<dim> &fe = info.fe_values();
  const unsigned int n_points = fe.n_quadrature_points;
  const unsigned int n_dofs = fe.dofs_per_cell;
  Vector<double> &b = dinfo.vector(0).block(0);
  const RightHandSide<dim>  right_hand_side;

  for (unsigned int k=0; k<n_points; ++k)
    {
      const double fval = right_hand_side.value(fe.quadrature_point(k));
      for (unsigned int i=0; i<n_dofs; ++i)
        b(i) += fe.JxW(k)*(fe.shape_value(i,k) * fval);
    }
}


template <int dim>
void RHSIntegrator<dim>::boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
{
  const FEValuesBase<dim> &fe = info.fe_values();
  const unsigned int n_points = fe.n_quadrature_points;
  const unsigned int n_dofs = fe.dofs_per_cell;
  Vector<double> &b = dinfo.vector(0).block(0);
  const FluxBoundaryValues<dim>  flux_bd;

  for (unsigned int k=0; k<n_points; ++k)
    {
      const double flux_bd_value = flux_bd.value(fe.quadrature_point(k));
      for (unsigned int i=0; i<n_dofs; ++i)
        b(i) += fe.JxW(k)*(fe.shape_value(i,k) * flux_bd_value);
    }
}


template <int dim>
void RHSIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &,
                              MeshWorker::DoFInfo<dim> &,
                              typename MeshWorker::IntegrationInfo<dim> &,
                              typename MeshWorker::IntegrationInfo<dim> &) const
{}


/*
 * Main class
 */
template<int dim>
class MeshWorkerConstraintMatrixTest
{
public:
  MeshWorkerConstraintMatrixTest(const FiniteElement<dim> &fe);
  ~MeshWorkerConstraintMatrixTest();

  void run();

private:
  void setup_system();
  void createInhomConstraints();
  void assemble_system_MeshWorker();
  void assemble_MeshWorker();

  Triangulation<dim> triangulation;
  const MappingQ1<dim> mapping;
  DoFHandler<dim> dof_handler;
  const FiniteElement<dim> &fe;

  SparsityPattern sparsity_pattern;
  SparseMatrix<double> system_matrix;
  SparseMatrix<double> matrix;
  SparseMatrix<double> error_matrix;

  Vector<double> system_rhs;
  Vector<double> rhs;
  Vector<double> error_rhs;

  ConstraintMatrix constraintsInhom;
  ConstraintMatrix constraints;
  ConstraintMatrix constraintsAll;
};


template<int dim>
MeshWorkerConstraintMatrixTest<dim>::MeshWorkerConstraintMatrixTest(const FiniteElement<dim> &fe) :
  mapping(),
  dof_handler(triangulation),
  fe(fe)
{
  GridGenerator::hyper_cube(this->triangulation, -1, 1);

  //refine with hanging node
  this->triangulation.refine_global (1);
  this->triangulation.begin_active()->set_refine_flag ();
  this->triangulation.execute_coarsening_and_refinement ();
  this->triangulation.refine_global (1);
}


template <int dim>
MeshWorkerConstraintMatrixTest<dim>::~MeshWorkerConstraintMatrixTest ()
{
  dof_handler.clear ();
}


template<int dim>
void MeshWorkerConstraintMatrixTest<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  system_rhs.reinit(dof_handler.n_dofs());
  rhs.reinit(dof_handler.n_dofs());
  error_rhs.reinit(dof_handler.n_dofs());

  constraints.clear();
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  constraints.close();

  DynamicSparsityPattern c_sparsity(dof_handler.n_dofs());

  std::string str("FE_DGQ");
  std::string fe_name = fe.get_name();
  std::size_t found = fe_name.find(str);
  if ( found!=std::string::npos)
    DoFTools::make_flux_sparsity_pattern(dof_handler, c_sparsity, constraints,
                                         false);
  else
    DoFTools::make_sparsity_pattern(dof_handler, c_sparsity, constraints,
                                    false);

  sparsity_pattern.copy_from(c_sparsity);

  system_matrix.reinit(sparsity_pattern);
  matrix.reinit(sparsity_pattern);
  error_matrix.reinit(sparsity_pattern);
}

/*
 * Assemble with SystemSimple
 */
template<int dim>
void MeshWorkerConstraintMatrixTest<dim>::assemble_system_MeshWorker()
{

  MeshWorker::IntegrationInfoBox<dim> info_box;

  const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1;
  info_box.initialize_gauss_quadrature(n_gauss_points, n_gauss_points,
                                       n_gauss_points);

  UpdateFlags update_flags = update_quadrature_points |
                             update_values |
                             update_gradients;

  info_box.add_update_flags_all(update_flags);
  info_box.initialize(fe, mapping);

  MeshWorker::DoFInfo<dim> dof_info(dof_handler);
  MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double> > assembler;

  assembler.initialize(system_matrix, system_rhs);
  assembler.initialize(constraintsAll);

  SystemIntegrator<dim> matrix_integrator;
  MeshWorker::integration_loop<dim, dim> (
    dof_handler.begin_active(), dof_handler.end(),
    dof_info, info_box, matrix_integrator, assembler);

}

/*
 * Assemble matrix and vector seperately
 */
template<int dim>
void MeshWorkerConstraintMatrixTest<dim>::assemble_MeshWorker()
{

  MeshWorker::IntegrationInfoBox<dim> info_box;

  const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1;
  info_box.initialize_gauss_quadrature(n_gauss_points, n_gauss_points,
                                       n_gauss_points);

  UpdateFlags update_flags = update_quadrature_points |
                             update_values |
                             update_gradients;

  info_box.add_update_flags_all(update_flags);
  info_box.initialize(fe, mapping);

  MeshWorker::DoFInfo<dim> dof_info(dof_handler);
  MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> > assemblerMatrix;
  assemblerMatrix.initialize(constraints);
  assemblerMatrix.initialize(matrix);
  MeshWorker::Assembler::ResidualSimple<Vector<double> > assemblerVector;
  assemblerVector.initialize(constraints);
  AnyData data;
  data.add<Vector<double>*>(&rhs, "RHS");
  assemblerVector.initialize(data);

  MatrixIntegrator<dim> matrix_integrator;
  MeshWorker::integration_loop<dim, dim> (
    dof_handler.begin_active(), dof_handler.end(),
    dof_info, info_box, matrix_integrator, assemblerMatrix);
  RHSIntegrator<dim> vector_integrator;
  MeshWorker::integration_loop<dim, dim> (
    dof_handler.begin_active(), dof_handler.end(),
    dof_info, info_box, vector_integrator, assemblerVector);
}


template <int dim>
void MeshWorkerConstraintMatrixTest<dim>::createInhomConstraints()
{
  this->constraintsInhom.clear();
  // boundary constraints
  VectorTools::interpolate_boundary_values (this->dof_handler,
                                            0,
                                            BoundaryValues<dim>(),
                                            this->constraintsInhom);
  // constraint of the origin
  std::map<types::global_dof_index, Point<dim> > support_points;
  DoFTools::map_dofs_to_support_points(this->mapping,this->dof_handler,support_points);
  for (unsigned int dof_index=0; dof_index < this->dof_handler.n_dofs(); ++dof_index)
    {
      if (!this->constraints.is_constrained(dof_index) && !this->constraintsInhom.is_constrained(dof_index))
        {
          const Point<dim> p = support_points[dof_index];

          if  (std::sqrt(p.square()) < 1e-6)
            {
              this->constraintsInhom.add_line (dof_index);
              this->constraintsInhom.set_inhomogeneity (dof_index, 2);
            }
        }
    }
  this->constraintsInhom.close();
}


template<int dim>
void MeshWorkerConstraintMatrixTest<dim>::run()
{
  setup_system();
  createInhomConstraints();
  constraintsAll.clear();
  constraintsAll.merge(constraints);
  constraintsAll.merge(constraintsInhom);
  constraintsAll.close();

  assemble_system_MeshWorker();

  assemble_MeshWorker();
  constraintsInhom.condense(matrix,rhs);

  // make diagonal entries and constraints component equal
  constraintsAll.distribute(system_rhs);
  constraintsAll.distribute(rhs);
  for (unsigned int i=0; i<error_matrix.m(); ++i)
    if (constraintsAll.is_constrained(i))
      {
        system_matrix.diag_element(i) = 1;
        matrix.diag_element(i) = 1;
      }

  // evaluate difference
  error_matrix.copy_from(system_matrix);
  error_matrix.add(-1.0,matrix);
  error_rhs = system_rhs;
  error_rhs -= rhs;
  deallog << "difference rhs "<< error_rhs.l2_norm() << std::endl;
  deallog << "difference matrix "<< error_matrix.frobenius_norm() << std::endl;

}


int main()
{
  initlog();
  deallog.threshold_double(1.e-10);

  FE_Q<2> fe(1);
  deallog.push(fe.get_name());
  MeshWorkerConstraintMatrixTest<2> test(fe);
  test.run();
  deallog.pop ();

  FE_DGQ<2> fe2(1);
  deallog.push(fe2.get_name());
  MeshWorkerConstraintMatrixTest<2> test2(fe2);
  test2.run();
  deallog.pop ();

}

