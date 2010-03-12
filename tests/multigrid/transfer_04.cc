//----------------------------------------------------------------------------
//    transfer.cc,v 1.13 2005/12/30 16:07:03 guido Exp
//    Version:
//
//    Copyright (C) 2000 - 2007, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

// like _01 but on adaptively refined grid

//TODO:[GK] Add checks for RT again!

#include "../tests.h"
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_out.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_tools.h>
#include <fe/fe_dgp.h>
#include <fe/fe_dgq.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/mg_level_object.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;


template <typename number>
class ScalingMatrix : public Subscriptor
{
  public:
				     /**
				      * Constructor setting the
				      * scaling factor. Default is
				      * constructing the identity
				      * matrix.
				      */
    ScalingMatrix(number scaling_factor = 1.);
    				     /**
				      * Apply preconditioner.
				      */
    template<class VECTOR>
    void vmult (VECTOR&, const VECTOR&) const;

				     /**
				      * Apply transpose
				      * preconditioner. Since this is
				      * the identity, this function is
				      * the same as
				      * vmult().
				      */
    template<class VECTOR>
    void Tvmult (VECTOR&, const VECTOR&) const;
				     /**
				      * Apply preconditioner, adding to the previous value.
				      */
    template<class VECTOR>
    void vmult_add (VECTOR&, const VECTOR&) const;

				     /**
				      * Apply transpose
				      * preconditioner, adding. Since this is
				      * the identity, this function is
				      * the same as
				      * vmult_add().
				      */
    template<class VECTOR>
    void Tvmult_add (VECTOR&, const VECTOR&) const;

  private:
    number factor;
};


//----------------------------------------------------------------------//

template<typename number>
ScalingMatrix<number>::ScalingMatrix(number factor)
		:
		factor(factor)
{}


template<typename number>
template<class VECTOR>
inline void
ScalingMatrix<number>::vmult (VECTOR &dst, const VECTOR &src) const
{
  dst.equ(factor, src);
}

template<typename number>
template<class VECTOR>
inline void
ScalingMatrix<number>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  dst.equ(factor, src);
}

template<typename number>
template<class VECTOR>
inline void
ScalingMatrix<number>::vmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(factor, src);
}



template<typename number>
template<class VECTOR>
inline void
ScalingMatrix<number>::Tvmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(factor, src);
}

template <int dim, typename number, int spacedim>
void
reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
	       MGLevelObject<dealii::Vector<number> > &v)
{
  for (unsigned int level=v.get_minlevel();
       level<=v.get_maxlevel();++level)
    {
      unsigned int n = mg_dof.n_dofs (level);
      v[level].reinit(n);
    }

}


template<int dim>
void refine_mesh (Triangulation<dim> &triangulation)
{
  bool cell_refined = false;
  for (typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active();
      cell != triangulation.end(); ++cell)
  {
   for (unsigned int vertex=0;
          vertex < GeometryInfo<dim>::vertices_per_cell;
          ++vertex)
      {
        const Point<dim> p = cell->vertex(vertex);
        const Point<dim> origin = (dim == 2 ?
                                    Point<dim>(0,0) :
                                    Point<dim>(0,0,0));
        const double dist = p.distance(origin);
        if(dist<0.25/M_PI)
        {
          cell->set_refine_flag ();
          cell_refined = true;
          break;
        }
      }
  }
  if(!cell_refined)//if no cell was selected for refinement, refine global
    for (typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active();
        cell != triangulation.end(); ++cell)
      cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement ();
}

template <int dim>
void initialize (const MGDoFHandler<dim> &dof,
    Vector<double> &u)
{
  unsigned int counter=0;
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices(dofs_per_cell);
  for (typename MGDoFHandler<dim>::active_cell_iterator
      cell = dof.begin_active();
      cell != dof.end(); ++cell)
  {
    cell->get_dof_indices(dof_indices);
    for(unsigned int i=0; i<dofs_per_cell; ++i)
      u(dof_indices[i]) = ++counter;
  }
}

template <int dim>
void initialize (const MGDoFHandler<dim> &dof,
    MGLevelObject<Vector<double> > &u)
{
  unsigned int counter=0;
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices(dofs_per_cell);
  typename MGDoFHandler<dim>::cell_iterator
      cell = dof.begin(0);
    cell->get_mg_dof_indices(dof_indices);
    for(unsigned int i=0; i<dofs_per_cell; ++i)
      u[0](dof_indices[i]) = ++counter;
}

template <int dim>
void print (const MGDoFHandler<dim> &dof,
    MGLevelObject<Vector<double> > &u)
{
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices(dofs_per_cell);
  for(unsigned int l=0; l<dof.get_tria().n_levels(); ++l)
  {
    deallog << std::endl;
    deallog << "Level " << l << std::endl;
    for (typename MGDoFHandler<dim>::cell_iterator
        cell = dof.begin(l);
        cell != dof.end(l); ++cell)
    {
      cell->get_mg_dof_indices(dof_indices);
      for(unsigned int i=0; i<dofs_per_cell; ++i)
        deallog << ' ' << (int)u[l](dof_indices[i]);
    }
  }
}

template <int dim>
void print_diff (const MGDoFHandler<dim> &dof_1, const MGDoFHandler<dim> &dof_2,
    const Vector<double> &u, const Vector<double> &v)
{
  Vector<double> diff;
  diff.reinit (u);
  const unsigned int dofs_per_cell = dof_1.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices_1(dofs_per_cell);
  std::vector<unsigned int> dof_indices_2(dofs_per_cell);
  for (typename MGDoFHandler<dim>::active_cell_iterator
      cell_1 = dof_1.begin_active(), cell_2 = dof_2.begin_active();
      cell_1 != dof_1.end(); ++cell_1, ++cell_2)
  {
    cell_1->get_dof_indices(dof_indices_1);
    cell_2->get_dof_indices(dof_indices_2);
    for(unsigned int i=0; i<dofs_per_cell; ++i)
      diff(dof_indices_1[i]) = u(dof_indices_1[i]) - v(dof_indices_2[i]);
  }
  deallog << std::endl;
  deallog << "diff " << diff.l2_norm() << std::endl;
}

template <int dim>
void check_smoother(const FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr,-1,1);
  refine_mesh(tr);
  refine_mesh(tr);

  //std::ostringstream out_filename;
  //out_filename << "gitter.eps";
  //std::ofstream grid_output (out_filename.str().c_str());
  //GridOut grid_out;
  //grid_out.write_eps (tr, grid_output);


  MGDoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);

  MGDoFHandler<dim> mgdof_renumbered(tr);
  mgdof_renumbered.distribute_dofs(fe);

  std::vector<unsigned int> block_component (2,0);
  block_component[1] = 1;

  DoFRenumbering::component_wise (mgdof_renumbered, block_component);
  for (unsigned int level=0;level<tr.n_levels();++level)
    DoFRenumbering::component_wise (mgdof_renumbered, level, block_component);

  ScalingMatrix<double> s1(-1);
  ScalingMatrix<double> s2(2.);
  ScalingMatrix<double> s8(8.);
  
  MGLevelObject<ScalingMatrix<double> > matrices(0, tr.n_levels()-1);
  matrices[0] = s1;
  matrices[1] = s2;
  matrices[2] = s8;

  typedef PreconditionJacobi<ScalingMatrix<double> > RELAX;
  MGSmootherPrecondition<ScalingMatrix, RELAX, Vector<double> > smoother(mem);
  RELAX::AdditionalData smoother_data;
  smoother.initialize(matrices, smoother_data);
  smoother.set_steps(1);

  smoother.smooth(level,u,rhs);

  Vector<double> u(mgdof.n_dofs());
  initialize(mgdof,u);

  MGLevelObject<Vector<double> > v(0, tr.n_levels()-1);
  reinit_vector(mgdof, v);

  transfer.copy_to_mg(mgdof, v, u);
  print(mgdof, v);

  u=0;
  initialize(mgdof, v);
  for(unsigned int l=0; l<tr.n_levels()-1; ++l)
    transfer.prolongate (l+1, v[l+1], v[l]);

  transfer.copy_from_mg(mgdof, u, v);
  Vector<double> diff = u;

  initialize(mgdof_renumbered,u);
  reinit_vector(mgdof_renumbered, v);
  transfer_renumbered.copy_to_mg(mgdof_renumbered, v, u);
  print(mgdof_renumbered, v);

  u=0;
  initialize(mgdof_renumbered, v);
  for(unsigned int l=0; l<tr.n_levels()-1; ++l)
    transfer_renumbered.prolongate (l+1, v[l+1], v[l]);

  transfer_renumbered.copy_from_mg(mgdof_renumbered, u, v);
  print_diff(mgdof_renumbered, mgdof, u, diff);
}


int main()
{
  std::ofstream logfile("transfer_04/output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check_simple (FESystem<2>(FE_Q<2>(1), 2));
}
