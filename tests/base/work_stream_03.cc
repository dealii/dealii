//----------------------------  work_stream_03.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  work_stream_03.cc  ---------------------------


// Moritz originally implemented thread local scratch objects for
// WorkStream in r24748 but it led to failures in the testsuite. what
// exactly went on was a mystery and this test is a first step in
// figuring out what happens by running a simplified version of one of
// the failing tests (deal.II/project_q_01) multiple times and
// verifying that it indeed works

#include "../tests.h"
#include <deal.II/base/work_stream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_abf.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_dgp_nonparametric.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <vector>
#include <iomanip>

using namespace dealii;


char logname[] = "work_stream_03/output";



template <int dim>
class F :  public Function<dim>
{
  public:
    F ()
		    :
		    q(1)
      {}

    virtual double value (const Point<dim> &p,
			  const unsigned int component) const
      {
	Assert (component == 0, ExcInternalError());
	double val = 0;
	for (unsigned int d=0; d<dim; ++d)
	  for (unsigned int i=0; i<=q; ++i)
	    val += (d+1)*(i+1)*std::pow (p[d], 1.*i);
	return val;
      }

  private:
    const unsigned int q;
};


Triangulation<2>     triangulation;
FE_Q<2> fe(3);
DoFHandler<2>        dof_handler(triangulation);
ConstraintMatrix constraints;

namespace
{
  template <int dim,
           int spacedim>
  struct Scratch
  {
    Scratch (const ::dealii::hp::FECollection<dim,spacedim> &fe,
             const UpdateFlags         update_flags,
             const Function<spacedim> *coefficient,
             const Function<spacedim> *rhs_function,
             const ::dealii::hp::QCollection<dim> &quadrature,
             const ::dealii::hp::MappingCollection<dim,spacedim> &mapping)
      :
      fe_collection (fe),
      quadrature_collection (quadrature),
      mapping_collection (mapping),
      x_fe_values (mapping_collection,
                   fe_collection,
                   quadrature_collection,
                   update_flags),
      coefficient_values(quadrature_collection.max_n_quadrature_points()),
      coefficient_vector_values (quadrature_collection.max_n_quadrature_points(),
                                 dealii::Vector<double> (fe_collection.n_components())),
      rhs_values(quadrature_collection.max_n_quadrature_points()),
      rhs_vector_values(quadrature_collection.max_n_quadrature_points(),
                        dealii::Vector<double> (fe_collection.n_components())),
      coefficient (coefficient),
      rhs_function (rhs_function),
      update_flags (update_flags)
    {}

    Scratch (const Scratch &data)
      :
      fe_collection (data.fe_collection),
      quadrature_collection (data.quadrature_collection),
      mapping_collection (data.mapping_collection),
      x_fe_values (mapping_collection,
                   fe_collection,
                   quadrature_collection,
                   data.update_flags),
      coefficient_values (data.coefficient_values),
      coefficient_vector_values (data.coefficient_vector_values),
      rhs_values (data.rhs_values),
      rhs_vector_values (data.rhs_vector_values),
      coefficient (data.coefficient),
      rhs_function (data.rhs_function),
      update_flags (data.update_flags)
    {}

    const ::dealii::hp::FECollection<dim,spacedim>      &fe_collection;
    const ::dealii::hp::QCollection<dim>                &quadrature_collection;
    const ::dealii::hp::MappingCollection<dim,spacedim> &mapping_collection;

    ::dealii::hp::FEValues<dim,spacedim> x_fe_values;

    std::vector<double>                  coefficient_values;
    std::vector<dealii::Vector<double> > coefficient_vector_values;
    std::vector<double>                  rhs_values;
    std::vector<dealii::Vector<double> > rhs_vector_values;

    std::vector<double> old_JxW;

    const Function<spacedim>   *coefficient;
    const Function<spacedim>   *rhs_function;

    const UpdateFlags update_flags;
  };


  struct CopyData
  {
    std::vector<types::global_dof_index> dof_indices;
    FullMatrix<double>        cell_matrix;
    dealii::Vector<double>    cell_rhs;
    const ConstraintMatrix   *constraints;
  };
}

template<int dim, int spacedim, typename CellIterator>
  void
  mass_assembler(const CellIterator &cell, Scratch<dim, spacedim> &data,
      CopyData &copy_data)
  {
    data.x_fe_values.reinit(cell);
    const FEValues<dim, spacedim> &fe_values =
        data.x_fe_values.get_present_fe_values();

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell, n_q_points =
        fe_values.n_quadrature_points;
    const FiniteElement<dim, spacedim> &fe = fe_values.get_fe();
    const unsigned int n_components = fe.n_components();

      {
        copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        copy_data.cell_rhs.reinit(dofs_per_cell);

        copy_data.dof_indices.resize(dofs_per_cell);
      }

      Assert (copy_data.cell_matrix.frobenius_norm() == 0, ExcInternalError());
    cell->get_dof_indices(copy_data.dof_indices);

    data.rhs_values.resize(n_q_points);
    data.rhs_function->value_list(fe_values.get_quadrature_points(),
        data.rhs_values);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const unsigned int component_i = fe.system_to_component_index(i).first;
        const double *phi_i = &fe_values.shape_value(i, 0);

        double add_data = 0;
        for (unsigned int point = 0; point < n_q_points; ++point)
          add_data += phi_i[point] * fe_values.JxW(point)
              * data.rhs_values[point];
        copy_data.cell_rhs(i) = add_data;
      }
  }


template <typename MatrixType,
         typename VectorType>
void copy_local_to_global (const CopyData &data,
                           MatrixType *matrix,
                           VectorType *right_hand_side)
{
  const unsigned int dofs_per_cell = data.dof_indices.size();
  Assert (data.cell_matrix.frobenius_norm() == 0, ExcInternalError());

  Assert (matrix->frobenius_norm() == 0, ExcInternalError());
    data.constraints->distribute_local_to_global(data.cell_matrix,
                                                 data.cell_rhs,
                                                 data.dof_indices,
                                                 *matrix, *right_hand_side);
    Assert (matrix->frobenius_norm() == 0, ExcInternalError());
//Q: why does this write anything into the matrix???
}



template <int dim, typename number, int spacedim>
void create_mass_matrix (const Mapping<dim,spacedim>       &mapping,
                         const DoFHandler<dim,spacedim>    &dof,
                         const Quadrature<dim>    &q,
                         SparseMatrix<number>     &matrix,
                         const Function<spacedim>      &rhs,
                         Vector<double>           &rhs_vector,
                         const Function<spacedim> *const coefficient,
                         const ConstraintMatrix   &constraints)
{
  Assert (matrix.m() == dof.n_dofs(),
          ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
          ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  hp::FECollection<dim,spacedim>      fe_collection (dof.get_fe());
  hp::QCollection<dim>                q_collection (q);
  hp::MappingCollection<dim,spacedim> mapping_collection (mapping);
  Scratch<dim, spacedim>
  assembler_data (fe_collection,
                  update_values |
                  update_JxW_values | update_quadrature_points,
                  coefficient, &rhs,
                  q_collection, mapping_collection);
  CopyData copy_data;
  copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
                                assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.constraints = &constraints;

  dealii::WorkStream::run (dof.begin_active(),
                   static_cast<typename DoFHandler<dim>::active_cell_iterator>(dof.end()),
                   &mass_assembler<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator>,
                   std_cxx1x::bind(&
                                   copy_local_to_global<SparseMatrix<number>, Vector<double> >,
                                   std_cxx1x::_1, &matrix, &rhs_vector),
                   assembler_data,
                   copy_data,
                   2*multithread_info.n_default_threads,
                   1);
}


template <int dim>
void do_project (const unsigned int        p)
{
  SparsityPattern sparsity;
  {
    CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints,
                                     false);

    sparsity.copy_from (csp);
  }

  SparseMatrix<double> mass_matrix (sparsity);
  Vector<double> tmp (mass_matrix.n());

  const Function<2>* dummy = 0;
  create_mass_matrix(StaticMappingQ1 < dim > ::mapping,
      dof_handler, QGauss < dim > (5), mass_matrix, F<dim>(), tmp, dummy,
      constraints);
  std::ostringstream x;
  x.precision(18);
  x << "Check1: " << mass_matrix.frobenius_norm() << " " << tmp.l1_norm() << std::endl;
  std::cout << x.str();
}




template <int dim>
void test ()
{
  Threads::TaskGroup<> g;
  for (unsigned int p=1; p<48*3; ++p)
    g += Threads::new_task (&do_project<dim>, 3);
  g.join_all ();
}


int main ()
{
  std::ofstream logfile(logname);
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  GridGenerator::hyper_cube (triangulation, 0, 1);
  triangulation.refine_global (1);
  dof_handler.distribute_dofs (fe);
  constraints.close ();

  test<2>();

  deallog << "OK" << std::endl;
}

