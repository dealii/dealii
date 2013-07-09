//----------------------------  project_q_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  project_q_01.cc  ---------------------------


// check that VectorTools::project works for Q elements correctly

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

char logname[] = "project_q_01/output";


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
          {
        deallog << "Projection failed with relative error "
                << error.l2_norm() / projection.l2_norm()
                << std::endl;
        Assert (false, ExcInternalError());
          }
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



template <int dim>
void test ()
{
  for (unsigned int p=1; p<6-dim; ++p)
    test_no_hanging_nodes (FE_Q<dim>(p), p);
}


int main ()
{
  std::ofstream logfile(logname);
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();
}
