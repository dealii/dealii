// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// make sure that we treat FE_Q elements correctly if face_flip==true ||
// face_rotation==true. for p>=3, we need to revert the order of dofs somehow
// between the two sides of the face. The same applies for lines in non-standard
// orientation. this test is derived from deal.II/project_q_03

char logname[] = "output";


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
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_q.h>

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



DeclException1 (ExcFailedProjection,
                double,
                << "The projection was supposed to exactly represent the "
                << "original function, but the relative residual is "
                << arg1);


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
        Assert (error.l2_norm() <= 1e-10*projection.l2_norm(),
                ExcFailedProjection(error.l2_norm() / projection.l2_norm()));
    }
}



// test with a 3d grid that has cells with face_rotation==false ||
// face_flip==false and hanging nodes. maybe it triggers all sorts of
// assumptions that may be hidden in places
//
// the mesh we use is a 7 cell moebius-type in 3d, with the first and last
// cell refined in turn. that then makes 2 meshes with 14 active cells
// each. this also cycles through all possibilities of coarser or finer cell
// having face_rotation==false || face_flip==false
//
// the whole procedure is repeated for three combinations of face_rotation and
// face_flip (the standard case in which evereything is as usual is left out)
template <int dim>
void test_with_wrong_face_orientation (const FiniteElement<dim> &fe,
                                       const unsigned int        p,
                                       const unsigned int        order_difference = 0)
{
  if (dim != 3)
    return;
  for (unsigned int j=1; j<4; ++j)
    // j=1: face_rotation=true,  face_flip=true
    // j=2: face_rotation=false, face_flip=true
    // j=3: face_rotation=true,  face_flip=false
    for (unsigned int i=0; i<2; ++i)
      {
        Triangulation<dim>     triangulation;
        GridGenerator::moebius (triangulation,7,j,1.0,0.2);
        typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active();
        if (i==0)
          {
            std::advance (cell,6);
            deallog<<"face_rotation="<<cell->face_rotation(2)
                   <<", face_flip="<<cell->face_flip(2)<<std::endl;
            cell = triangulation.begin_active();
          }
        std::advance (cell, i*6);
        cell->set_refine_flag ();
        triangulation.execute_coarsening_and_refinement ();

        do_project (triangulation, fe, p, order_difference);
      }
}




int main ()
{
  std::ofstream logfile(logname);
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<3>();
}



template <int dim>
void test ()
{
  test_with_wrong_face_orientation (FE_Q<dim>(3), 3);
}
