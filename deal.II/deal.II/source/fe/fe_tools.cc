//----------------------------  fe_tools.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools.cc  ---------------------------


#include <base/quadrature.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_iterator.h>
#include <fe/fe_tools.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>

// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif



template <int dim, typename number>
void FETools::get_interpolation_matrix(const FiniteElement<dim> &fe1,
				       const FiniteElement<dim> &fe2,
				       FullMatrix<number> &interpolation_matrix)
{
  Assert (fe1.n_components() == fe2.n_components(),
	  ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
  Assert(interpolation_matrix.m()==fe2.dofs_per_cell &&
	 interpolation_matrix.n()==fe1.dofs_per_cell,
	 ExcMatrixDimensionMismatch(interpolation_matrix.m(),
				    interpolation_matrix.n(),
				    fe2.dofs_per_cell,
				    fe1.dofs_per_cell));

				   // Initialize FEValues for fe1 at
				   // the unit support points of the
				   // fe2 element.
  const typename std::vector<Point<dim> > &
    fe2_support_points = fe2.get_unit_support_points ();

  Assert(fe2_support_points.size()==fe2.dofs_per_cell,
	 typename FiniteElementBase<dim>::ExcFEHasNoSupportPoints());

  for (unsigned int i=0; i<fe2.dofs_per_cell; ++i)	
    {
      const unsigned int i1 = fe2.system_to_component_index(i).first;
      for (unsigned int j=0; j<fe1.dofs_per_cell; ++j)
	{
	  const unsigned int j1 = fe1.system_to_component_index(j).first;
	  if (i1==j1)
	    interpolation_matrix(i,j) = fe1.shape_value (j,fe2_support_points[i]);
	  else
	    interpolation_matrix(i,j)=0.;
	}  
    }
}



template <int dim, typename number>
void FETools::get_back_interpolation_matrix(const FiniteElement<dim> &fe1,
					    const FiniteElement<dim> &fe2,
					    FullMatrix<number> &interpolation_matrix)
{
  Assert (fe1.n_components() == fe2.n_components(),
	  ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
  Assert(interpolation_matrix.m()==fe1.dofs_per_cell &&
	 interpolation_matrix.n()==fe1.dofs_per_cell, 
	 ExcMatrixDimensionMismatch(interpolation_matrix.m(),
				    interpolation_matrix.n(),
				    fe1.dofs_per_cell,
				    fe1.dofs_per_cell));
    
  FullMatrix<number> first_matrix (fe2.dofs_per_cell, fe1.dofs_per_cell);
  FullMatrix<number> second_matrix(fe1.dofs_per_cell, fe2.dofs_per_cell);
  
  get_interpolation_matrix(fe1, fe2, first_matrix);
  get_interpolation_matrix(fe2, fe1, second_matrix);

				   // int_matrix=second_matrix*first_matrix
  second_matrix.mmult(interpolation_matrix, first_matrix);
}



template <int dim, typename number>
void FETools::get_interpolation_difference_matrix(const FiniteElement<dim> &fe1,
						  const FiniteElement<dim> &fe2,
						  FullMatrix<number> &difference_matrix)
{
  Assert (fe1.n_components() == fe2.n_components(),
	  ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
  Assert(difference_matrix.m()==fe1.dofs_per_cell &&
	 difference_matrix.n()==fe1.dofs_per_cell, 
	 ExcMatrixDimensionMismatch(difference_matrix.m(),
				    difference_matrix.n(),
				    fe1.dofs_per_cell,
				    fe1.dofs_per_cell));
   
  FullMatrix<number> interpolation_matrix(fe1.dofs_per_cell);
  get_back_interpolation_matrix(fe1, fe2, interpolation_matrix);
  
  for (unsigned int i=0; i<fe1.dofs_per_cell; ++i)
    difference_matrix(i,i) = 1.;
  
				   // compute difference
  difference_matrix.add (-1, interpolation_matrix);
}



template <int dim, typename number>
void FETools::interpolate(const DoFHandler<dim> &dof1,
			  const Vector<number> &u1,
			  const DoFHandler<dim> &dof2,
			  Vector<number> &u2)
{
  Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());
  Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
  Assert(u2.size()==dof2.n_dofs(), ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

  const unsigned int dofs_per_cell1=dof1.get_fe().dofs_per_cell;
  const unsigned int dofs_per_cell2=dof2.get_fe().dofs_per_cell;
  
  Vector<number> u1_local(dofs_per_cell1);
  Vector<number> u2_local(dofs_per_cell2);

  FullMatrix<number> interpolation_matrix(dofs_per_cell2,
					  dofs_per_cell1);
  FETools::get_interpolation_matrix(dof1.get_fe(), dof2.get_fe(),
				    interpolation_matrix);
  
  typename DoFHandler<dim>::active_cell_iterator cell1 = dof1.begin_active(),
						 endc1 = dof1.end(),
						 cell2 = dof2.begin_active(),
						 endc2 = dof2.end();

  std::vector<unsigned int> index_multiplicity(dof2.n_dofs(),0);
  std::vector<unsigned int> dofs (dofs_per_cell2);
  u2.clear ();
  
  for (; cell1!=endc1, cell2!=endc2; ++cell1, ++cell2) 
    {
      cell1->get_dof_values(u1, u1_local);
      interpolation_matrix.vmult(u2_local, u1_local);
      cell2->get_dof_indices(dofs);
      for (unsigned int i=0; i<dofs_per_cell2; ++i)
	{
	  u2(dofs[i])+=u2_local(i);
	  ++index_multiplicity[dofs[i]];
	}
    }
  for (unsigned int i=0; i<dof2.n_dofs(); ++i)
    {
      Assert(index_multiplicity[i]!=0, ExcInternalError());
      u2(i) /= index_multiplicity[i];
    }
}



template <int dim, typename number>
void FETools::back_interpolate(const DoFHandler<dim> &dof1,
			       const Vector<number> &u1,
			       const FiniteElement<dim> &fe2,
			       Vector<number> &u1_interpolated)
{
  Assert(dof1.get_fe().n_components() == fe2.n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
  Assert(u1_interpolated.size()==dof1.n_dofs(),
	 ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));  

  const unsigned int dofs_per_cell1=dof1.get_fe().dofs_per_cell;

  Vector<number> u1_local(dofs_per_cell1);
  Vector<number> u1_int_local(dofs_per_cell1);
  
  FullMatrix<number> interpolation_matrix(dofs_per_cell1, dofs_per_cell1);
  FETools::get_back_interpolation_matrix(dof1.get_fe(), fe2,
					 interpolation_matrix);

  typename DoFHandler<dim>::active_cell_iterator cell1 = dof1.begin_active(),
						 endc1 = dof1.end();
  
  for (; cell1!=endc1; ++cell1) 
    {
      cell1->get_dof_values(u1, u1_local);
      interpolation_matrix.vmult(u1_int_local, u1_local);
      cell1->set_dof_values(u1_int_local, u1_interpolated);
    }
}



template <int dim, typename number>
void FETools::interpolation_difference(const DoFHandler<dim> &dof1,
				       const Vector<number> &u1,
				       const FiniteElement<dim> &fe2,
				       Vector<number> &u1_difference)
{
  Assert(dof1.get_fe().n_components() == fe2.n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
  Assert(u1_difference.size()==dof1.n_dofs(),
	 ExcDimensionMismatch(u1_difference.size(), dof1.n_dofs()));  

  const unsigned int dofs_per_cell1=dof1.get_fe().dofs_per_cell;

  Vector<number> u1_local(dofs_per_cell1);
  Vector<number> u1_diff_local(dofs_per_cell1);
  
  FullMatrix<number> difference_matrix(dofs_per_cell1, dofs_per_cell1);
  FETools::get_interpolation_difference_matrix(dof1.get_fe(), fe2,
					       difference_matrix);
  
  typename DoFHandler<dim>::active_cell_iterator cell1 = dof1.begin_active(),
						 endc1 = dof1.end();
  
  for (; cell1!=endc1; ++cell1) 
    {
      cell1->get_dof_values(u1, u1_local);
      difference_matrix.vmult(u1_diff_local, u1_local);
      cell1->set_dof_values(u1_diff_local, u1_difference);
    }
}



template <int dim, typename number>
void FETools::extrapolate(const DoFHandler<dim> &dof1,
			  const Vector<number> &u1,
			  const DoFHandler<dim> &dof2,
			  Vector<number> &u2)
{
  Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
  Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
  Assert(u2.size()==dof2.n_dofs(), ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

  Vector<number> u3(dof2.n_dofs());
  interpolate(dof1, u1, dof2, u3);

  const unsigned int dofs_per_cell  = dof2.get_fe().dofs_per_cell;
  const Triangulation<dim> &tria=dof1.get_tria();
  Vector<number> dof_values(dofs_per_cell);
  
  for (unsigned int level=0; level<tria.n_levels()-1; ++level)
    {
      typename DoFHandler<dim>::cell_iterator cell=dof2.begin(level),
					      endc=dof2.end(level);

      for (; cell!=endc; ++cell)
	if (!cell->active())
	  {
					     // check whether this
					     // cell has active
					     // children
	    bool active_children=false;
	    for (unsigned int child_n=0;
		 child_n<GeometryInfo<dim>::children_per_cell; ++child_n)
	      if (cell->child(child_n)->active())
		{
		  active_children=true;
		  break;
		}

					     // if there are active
					     // children, the we have
					     // to work on this
					     // cell. get the data
					     // from the one vector
					     // and set it on the
					     // other
	    if (active_children)
	      {
		cell->get_interpolated_dof_values(u3, dof_values);
		cell->set_dof_values_by_interpolation(dof_values, u2);
	      }
	  }
    }
}




/*-------------- Explicit Instantiations -------------------------------*/

template
void FETools::get_interpolation_matrix(const FiniteElement<deal_II_dimension> &,
				       const FiniteElement<deal_II_dimension> &,
				       FullMatrix<double> &);
template
void FETools::get_back_interpolation_matrix(const FiniteElement<deal_II_dimension> &,
					    const FiniteElement<deal_II_dimension> &,
					    FullMatrix<double> &);
template
void FETools::get_interpolation_difference_matrix(const FiniteElement<deal_II_dimension> &,
						  const FiniteElement<deal_II_dimension> &,
						  FullMatrix<double> &);
template
void FETools::interpolate(const DoFHandler<deal_II_dimension> &,
			  const Vector<double> &,
			  const DoFHandler<deal_II_dimension> &,
			  Vector<double> &);
template
void FETools::back_interpolate(const DoFHandler<deal_II_dimension> &,
			       const Vector<double> &,
			       const FiniteElement<deal_II_dimension> &,
			       Vector<double> &);
template
void FETools::interpolation_difference(const DoFHandler<deal_II_dimension> &,
				       const Vector<double> &,
				       const FiniteElement<deal_II_dimension> &,
				       Vector<double> &);
template
void FETools::extrapolate(const DoFHandler<deal_II_dimension> &,
			  const Vector<double> &,
			  const DoFHandler<deal_II_dimension> &,
			  Vector<double> &);


template
void FETools::get_interpolation_matrix(const FiniteElement<deal_II_dimension> &,
				       const FiniteElement<deal_II_dimension> &,
				       FullMatrix<float> &);
template
void FETools::get_back_interpolation_matrix(const FiniteElement<deal_II_dimension> &,
					    const FiniteElement<deal_II_dimension> &,
					    FullMatrix<float> &);
template
void FETools::get_interpolation_difference_matrix(const FiniteElement<deal_II_dimension> &,
						  const FiniteElement<deal_II_dimension> &,
						  FullMatrix<float> &);
template
void FETools::interpolate(const DoFHandler<deal_II_dimension> &,
			  const Vector<float> &,
			  const DoFHandler<deal_II_dimension> &,
			  Vector<float> &);
template
void FETools::back_interpolate(const DoFHandler<deal_II_dimension> &,
			       const Vector<float> &,
			       const FiniteElement<deal_II_dimension> &,
			       Vector<float> &);
template
void FETools::interpolation_difference(const DoFHandler<deal_II_dimension> &,
				       const Vector<float> &,
				       const FiniteElement<deal_II_dimension> &,
				       Vector<float> &);
template
void FETools::extrapolate(const DoFHandler<deal_II_dimension> &,
			  const Vector<float> &,
			  const DoFHandler<deal_II_dimension> &,
			  Vector<float> &);

/*----------------------------   fe_tools.cc     ---------------------------*/
