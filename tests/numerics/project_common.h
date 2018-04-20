// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2017 by the deal.II authors
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


// common framework to check whether an element of polynomial order p can
// represent functions of order q

#include "../tests.h"
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
#include <deal.II/grid/manifold_lib.h>
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

#include <fstream>
#include <vector>


template <int dim>
void test ();


template <int dim>
class F :  public Function<dim>
{
public:
  F (const unsigned int q,
     const unsigned int n_components)
    :
    Function<dim>(n_components),
    q(q)
  {}

  virtual double value (const Point<dim> &p,
                        const unsigned int component) const
  {
    Assert ((component == 0) && (this->n_components == 1),
            ExcInternalError());
    double val = 0;
    for (unsigned int d=0; d<dim; ++d)
      for (unsigned int i=0; i<=q; ++i)
        val += (d+1)*(i+1)*std::pow (p[d], 1.*i);
    return val;
  }


  virtual void vector_value (const Point<dim> &p,
                             Vector<double>   &v) const
  {
    for (unsigned int c=0; c<v.size(); ++c)
      {
        v(c) = 0;
        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int i=0; i<=q; ++i)
            v(c) += (d+1)*(i+1)*std::pow (p[d], 1.*i)+c;
      }
  }

private:
  const unsigned int q;
};


template <int dim>
void do_project (const Triangulation<dim> &triangulation,
                 const FiniteElement<dim> &fe,
                 const unsigned int        p,
                 const unsigned int        order_difference)
{
  DoFHandler<dim>        dof_handler(triangulation);
  dof_handler.distribute_dofs (fe);

  deallog << "n_dofs=" << dof_handler.n_dofs() << std::endl;

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);
  constraints.close ();

  Vector<double> projection (dof_handler.n_dofs());
  Vector<float>  error (triangulation.n_active_cells());
  for (unsigned int q=0; q<=p+2-order_difference; ++q)
    {
      // project the function
      VectorTools::project (dof_handler,
                            constraints,
                            QGauss<dim>(p+2),
                            F<dim> (q, fe.n_components()),
                            projection);
      // just to make sure it doesn't get
      // forgotten: handle hanging node
      // constraints
      constraints.distribute (projection);

      // then compute the interpolation error
      VectorTools::integrate_difference (dof_handler,
                                         projection,
                                         F<dim> (q, fe.n_components()),
                                         error,
                                         QGauss<dim>(std::max(p,q)+1),
                                         VectorTools::L2_norm);
      deallog << fe.get_name() << ", P_" << q
              << ", rel. error=" << error.l2_norm() / projection.l2_norm()
              << std::endl;

      if (q<=p-order_difference)
        if (error.l2_norm() > 1e-10*projection.l2_norm())
          deallog << "Projection failed with relative error "
                  << error.l2_norm() / projection.l2_norm()
                  << std::endl;
    }
}



// check the given element of polynomial order p. the last parameter, if
// given, denotes a gap in convergence order; for example, the Nedelec element
// of polynomial degree p has normal components of degree p-1 and therefore
// can only represent polynomials of degree p-1 exactly. the gap is then 1.
template <int dim>
void test_no_hanging_nodes (const FiniteElement<dim> &fe,
                            const unsigned int        p,
                            const unsigned int        order_difference = 0)
{
  Triangulation<dim>     triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (3);

  do_project (triangulation, fe, p, order_difference);
}



// same test as above, but this time with a mesh that has hanging nodes
template <int dim>
void test_with_hanging_nodes (const FiniteElement<dim> &fe,
                              const unsigned int        p,
                              const unsigned int        order_difference = 0)
{
  Triangulation<dim>     triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (1);
  triangulation.begin_active()->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();
  triangulation.refine_global (1);

  do_project (triangulation, fe, p, order_difference);
}



// test with a 3d grid that has cells with face_orientation==false and hanging
// nodes. this trips up all sorts of pieces of code, for example there was a
// crash when computing hanging node constraints on such faces (see
// bits/face_orientation_crash), and it triggers all sorts of other
// assumptions that may be hidden in places
//
// the mesh we use is the 7 cells of the hyperball mesh in 3d, with each of
// the cells refined in turn. that then makes 7 meshes with 14 active cells
// each. this also cycles through all possibilities of coarser or finer cell
// having face_orientation==false
template <int dim>
void test_with_wrong_face_orientation (const FiniteElement<dim> &fe,
                                       const unsigned int        p,
                                       const unsigned int        order_difference = 0)
{
  if (dim != 3)
    return;

  for (unsigned int i=0; i<7; ++i)
    {
      Triangulation<dim>     triangulation;
      GridGenerator::hyper_ball (triangulation);
      triangulation.set_manifold(0);
      typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active();
      std::advance (cell, i);
      cell->set_refine_flag ();
      triangulation.execute_coarsening_and_refinement ();

      do_project (triangulation, fe, p, order_difference);
    }
}




// test with a 2d mesh that forms a square but subdivides it into 3
// elements. this tests the case of the sign_change thingy in
// fe_poly_tensor.cc
template <int dim>
void test_with_2d_deformed_mesh (const FiniteElement<dim> &fe,
                                 const unsigned int        p,
                                 const unsigned int        order_difference = 0)
{
  if (dim != 2)
    return;

  std::vector<Point<dim> > points_glob;
  std::vector<Point<dim> > points;

  points_glob.push_back (Point<dim> (0.0, 0.0));
  points_glob.push_back (Point<dim> (1.0, 0.0));
  points_glob.push_back (Point<dim> (1.0, 0.5));
  points_glob.push_back (Point<dim> (1.0, 1.0));
  points_glob.push_back (Point<dim> (0.6, 0.5));
  points_glob.push_back (Point<dim> (0.5, 1.0));
  points_glob.push_back (Point<dim> (0.0, 1.0));

  // Prepare cell data
  std::vector<CellData<dim> > cells (3);

  cells[0].vertices[0] = 0;
  cells[0].vertices[1] = 1;
  cells[0].vertices[2] = 4;
  cells[0].vertices[3] = 2;
  cells[0].material_id = 0;

  cells[1].vertices[0] = 4;
  cells[1].vertices[1] = 2;
  cells[1].vertices[2] = 5;
  cells[1].vertices[3] = 3;
  cells[1].material_id = 0;

  cells[2].vertices[0] = 0;
  cells[2].vertices[1] = 4;
  cells[2].vertices[2] = 6;
  cells[2].vertices[3] = 5;
  cells[2].material_id = 0;

  Triangulation<dim>     triangulation;
  triangulation.create_triangulation (points_glob, cells, SubCellData());

  do_project (triangulation, fe, p, order_difference);
}



// same as test_with_2d_deformed_mesh, but refine each element in turn. this
// makes sure we also check the sign_change thingy for refined cells
template <int dim>
void test_with_2d_deformed_refined_mesh (const FiniteElement<dim> &fe,
                                         const unsigned int        p,
                                         const unsigned int        order_difference = 0)
{
  if (dim != 2)
    return;

  for (unsigned int i=0; i<3; ++i)
    {
      std::vector<Point<dim> > points_glob;
      std::vector<Point<dim> > points;

      points_glob.push_back (Point<dim> (0.0, 0.0));
      points_glob.push_back (Point<dim> (1.0, 0.0));
      points_glob.push_back (Point<dim> (1.0, 0.5));
      points_glob.push_back (Point<dim> (1.0, 1.0));
      points_glob.push_back (Point<dim> (0.6, 0.5));
      points_glob.push_back (Point<dim> (0.5, 1.0));
      points_glob.push_back (Point<dim> (0.0, 1.0));

      // Prepare cell data
      std::vector<CellData<dim> > cells (3);

      cells[0].vertices[0] = 0;
      cells[0].vertices[1] = 1;
      cells[0].vertices[2] = 4;
      cells[0].vertices[3] = 2;
      cells[0].material_id = 0;

      cells[1].vertices[0] = 4;
      cells[1].vertices[1] = 2;
      cells[1].vertices[2] = 5;
      cells[1].vertices[3] = 3;
      cells[1].material_id = 0;

      cells[2].vertices[0] = 0;
      cells[2].vertices[1] = 4;
      cells[2].vertices[2] = 6;
      cells[2].vertices[3] = 5;
      cells[2].material_id = 0;

      Triangulation<dim>     triangulation;
      triangulation.create_triangulation (points_glob, cells, SubCellData());

      switch (i)
        {
        case 0:
          triangulation.begin_active()->set_refine_flag();
          break;
        case 1:
          (++(triangulation.begin_active()))->set_refine_flag();
          break;
        case 2:
          (++(++(triangulation.begin_active())))->set_refine_flag();
          break;
        default:
          Assert (false, ExcNotImplemented());
        }
      triangulation.execute_coarsening_and_refinement ();

      do_project (triangulation, fe, p, order_difference);
    }
}



int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1>();
  test<2>();
  test<3>();
}

