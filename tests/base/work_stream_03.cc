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
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
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
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/hp/fe_values.h>

#include <fstream>
#include <vector>

char logname[] = "work_stream_03/output";



template <int dim>
class F :  public Function<dim>
{
  public:
    virtual double value (const Point<dim> &p,
                          const unsigned int component) const
      {
        Assert ((component == 0) && (this->n_components == 1),
                ExcInternalError());
        double val = 0;
        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int i=0; i<=1; ++i)
            val += (d+1)*(i+1)*std::pow (p[d], 1.*i);
        return val;
      }
};

namespace
{
  template <int dim,
           int spacedim>
  struct Scratch
  {
    Scratch (const ::dealii::hp::FECollection<dim,spacedim> &fe,
             const UpdateFlags         update_flags,
             const ::dealii::hp::QCollection<dim> &quadrature)
      :
      fe_collection (fe),
      quadrature_collection (quadrature),
      x_fe_values (fe_collection,
                   quadrature_collection,
                   update_flags),
      rhs_values(quadrature_collection.max_n_quadrature_points()),
      update_flags (update_flags)
    {}

    Scratch (const Scratch &data)
      :
      fe_collection (data.fe_collection),
      quadrature_collection (data.quadrature_collection),
      x_fe_values (fe_collection,
                   quadrature_collection,
                   data.update_flags),
      rhs_values (data.rhs_values),
      update_flags (data.update_flags)
    {}

    const ::dealii::hp::FECollection<dim,spacedim>      &fe_collection;
    const ::dealii::hp::QCollection<dim>                &quadrature_collection;

    ::dealii::hp::FEValues<dim,spacedim> x_fe_values;

    std::vector<double>                  rhs_values;

    const UpdateFlags update_flags;
  };


  struct CopyData
  {
    Vector<double>    cell_rhs;
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
        copy_data.cell_rhs.reinit(dofs_per_cell);
      }

    data.rhs_values.resize(n_q_points);
    F<dim>().value_list(fe_values.get_quadrature_points(),
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


template <typename VectorType>
void copy_local_to_global (const CopyData &data,
                           VectorType *right_hand_side)
{
//  std::cout << data.cell_rhs.l2_norm() << ' ';
  right_hand_side->push_back (data.cell_rhs.l2_norm());
}



template <int dim, int spacedim>
void create_mass_matrix (const DoFHandler<dim,spacedim>    &dof,
                         const Quadrature<dim>    &q,
                         std::vector<double>           &rhs_vector)
{
  hp::FECollection<dim,spacedim>      fe_collection (dof.get_fe());
  hp::QCollection<dim>                q_collection (q);
  Scratch<dim, spacedim>
  assembler_data (fe_collection,
                  update_values |
                  update_JxW_values | update_quadrature_points,
                  q_collection);
  CopyData copy_data;
  copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());

  dealii::WorkStream::run (dof.begin_active(),
                   static_cast<typename DoFHandler<dim>::active_cell_iterator>(dof.end()),
                   &mass_assembler<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator>,
                   std_cxx1x::bind(&
                                   copy_local_to_global<std::vector<double> >,
                                   std_cxx1x::_1, &rhs_vector),
                   assembler_data,
                   copy_data,
                   2*multithread_info.n_default_threads,
                   1);
}


template <int dim>
void do_project (const FiniteElement<dim> &fe)
{
  Triangulation<dim>     triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (2);

  DoFHandler<dim>        dof_handler(triangulation);
  dof_handler.distribute_dofs (fe);

  std::vector<double> tmp;

  create_mass_matrix (dof_handler, QGauss<dim>(3),
      tmp);

  double sum=0;
  for (unsigned int i=0; i<tmp.size(); ++i)
    sum += std::fabs(tmp[i]);
  printf ("\nCheck: %5.13f\n", sum);
}



template <int dim>
void test_no_hanging_nodes (const FiniteElement<dim> &fe)
{
  for (unsigned int i=0; i<12; ++i)
    do_project (fe);
}


int main ()
{
  std::ofstream logfile(logname);
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_no_hanging_nodes (FE_Q<3>(1));
}
