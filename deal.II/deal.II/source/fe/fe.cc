/* $Id$ */

#include <fe/fe.h>
#include <fe/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>



template <int dim>
FEValues<dim>::UpdateStruct::UpdateStruct () :
		update_q_points(false),
		update_gradients(false),
		update_jacobians(false),
		update_JxW_values(false),
		update_ansatz_points(false) {};
		




/*------------------------------- FEValues -------------------------------*/


template <int dim>
FEValues<dim>::FEValues (const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &quadrature,
			 const UpdateStruct       &update_flags) :
		n_quadrature_points(quadrature.n_quadrature_points),
		total_dofs(fe.total_dofs),
		shape_values(fe.total_dofs, quadrature.n_quadrature_points),
		shape_gradients(fe.total_dofs,
				vector<Point<dim> >(quadrature.n_quadrature_points)),
		unit_shape_gradients(fe.total_dofs,
				     vector<Point<dim> >(quadrature.n_quadrature_points)),
		weights(quadrature.n_quadrature_points, 0),
		JxW_values(quadrature.n_quadrature_points, 0),
		quadrature_points(quadrature.n_quadrature_points, Point<dim>()),
		unit_quadrature_points(quadrature.n_quadrature_points, Point<dim>()),
		jacobi_matrices (quadrature.n_quadrature_points,
				 dFMatrix(dim,dim)),
		update_flags (update_flags)
{
  for (unsigned int i=0; i<fe.total_dofs; ++i)
    for (unsigned int j=0; j<n_quadrature_points; ++j) 
      {
	shape_values(i,j) = fe.shape_value(i, quadrature.quad_point(j));
	unit_shape_gradients[i][j]
	  = fe.shape_grad(i, quadrature.quad_point(j));
      };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      weights[i] = quadrature.weight(i);
      unit_quadrature_points[i] = quadrature.quad_point(i);
    };
};



template <int dim>
double FEValues<dim>::shape_value (const unsigned int i,
				   const unsigned int j) const {
  Assert (i<(unsigned int)shape_values.m(), ExcInvalidIndex (i, shape_values.m()));
  Assert (j<(unsigned int)shape_values.n(), ExcInvalidIndex (j, shape_values.n()));

  return shape_values(i,j);
};



template <int dim>
const Point<dim> &
FEValues<dim>::shape_grad (const unsigned int i,
			   const unsigned int j) const {
  Assert (i<(unsigned int)shape_values.m(), ExcInvalidIndex (i, shape_values.m()));
  Assert (j<(unsigned int)shape_values.n(), ExcInvalidIndex (j, shape_values.n()));
  Assert (update_flags.update_gradients, ExcAccessToUninitializedField());

  return shape_gradients[i][j];
};



template <int dim>
const Point<dim> & FEValues<dim>::quadrature_point (const unsigned int i) const {
  Assert (i<n_quadrature_points, ExcInvalidIndex(i, n_quadrature_points));
  Assert (update_flags.update_q_points, ExcAccessToUninitializedField());
  
  return quadrature_points[i];
};



template <int dim>
const Point<dim> & FEValues<dim>::ansatz_point (const unsigned int i) const {
  Assert (i<ansatz_points.size(), ExcInvalidIndex(i, ansatz_points.size()));
  Assert (update_flags.update_ansatz_points, ExcAccessToUninitializedField());
  
  return ansatz_points[i];
};



template <int dim>
double FEValues<dim>::JxW (const unsigned int i) const {
  Assert (i<n_quadrature_points, ExcInvalidIndex(i, n_quadrature_points));
  Assert (update_flags.update_JxW_values, ExcAccessToUninitializedField());
  
  return JxW_values[i];
};



template <int dim>
void FEValues<dim>::reinit (const typename Triangulation<dim>::cell_iterator &cell,
			    const FiniteElement<dim>                         &fe) {
				   // fill jacobi matrices and real
				   // quadrature points
  if (update_flags.update_jacobians || update_flags.update_q_points)
    fe.fill_fe_values (cell,
		       unit_quadrature_points,
		       jacobi_matrices,
		       update_flags.update_jacobians,
		       ansatz_points,
		       update_flags.update_ansatz_points,
		       quadrature_points,
		       update_flags.update_q_points);

				   // compute gradients on real element if
				   // requested
  if (update_flags.update_gradients) 
    {
      Assert (update_flags.update_jacobians, ExcCannotInitializeField());
      
      for (unsigned int i=0; i<fe.total_dofs; ++i)
	for (unsigned int j=0; j<n_quadrature_points; ++j)
	  for (unsigned int s=0; s<dim; ++s)
	    {
	      shape_gradients[i][j](s) = 0;
	      
					       // (grad psi)_s =
					       // (grad_{\xi\eta})_b J_{bs}
					       // with J_{bs}=(d\xi_b)/(dx_s)
	      for (unsigned int b=0; b<dim; ++b)
		shape_gradients[i][j](s)
		  +=
		  unit_shape_gradients[i][j](b) * jacobi_matrices[j](b,s);
	    };
    };
  
  
				   // compute Jacobi determinants in
				   // quadrature points.
				   // refer to the general doc for
				   // why we take the inverse of the
				   // determinant
  if (update_flags.update_JxW_values) 
    {
      Assert (update_flags.update_jacobians,
	      ExcCannotInitializeField());
      for (unsigned int i=0; i<n_quadrature_points; ++i)
	JxW_values[i] = weights[i] / jacobi_matrices[i].determinant();
    };
};







/*------------------------------- FiniteElementBase ----------------------*/

template <int dim>
FiniteElementBase<dim>::FiniteElementBase (const FiniteElementBase<dim> &f) :
		total_dofs(f.total_dofs),
		interface_constraints (f.interface_constraints)
{
  for (unsigned int i=0; i<(1<<dim); ++i) 
    {
      restriction[i] = f.restriction[i];
      prolongation[i] = f.prolongation[i];
    };
}



template <int dim>
const dFMatrix &
FiniteElementBase<dim>::restrict (const unsigned int child) const {
  Assert (child<(1<<dim), ExcInvalidIndex(child));
  return restriction[child];
};



template <int dim>
const dFMatrix &
FiniteElementBase<dim>::prolongate (const unsigned int child) const {
  Assert (child<(1<<dim), ExcInvalidIndex(child));
  return prolongation[child];
};



template <int dim>
const dFMatrix &
FiniteElementBase<dim>::constraints () const {
  return interface_constraints;
};



template <int dim>
bool FiniteElementBase<dim>::operator == (const FiniteElementBase<dim> &f) const {
  return ((total_dofs == f.total_dofs) &&
	  (interface_constraints == f.interface_constraints));
};



template <int dim>
double FiniteElementBase<dim>::shape_value (const unsigned int,
					    const Point<dim> &) const {
  Assert (false, ExcPureFunctionCalled());
  return 0;
};



template <int dim>
Point<dim> FiniteElementBase<dim>::shape_grad (const unsigned int,
					       const Point<dim> &) const {
  Assert (false, ExcPureFunctionCalled());
  return Point<dim>();
};




template <int dim>
void FiniteElementBase<dim>::fill_fe_values (const typename Triangulation<dim>::cell_iterator &,
					     const vector<Point<dim> > &,
					     vector<dFMatrix> &,
					     const bool,
					     vector<Point<dim> > &,
					     const bool,
					     vector<Point<dim> > &,
					     const bool) const {
  Assert (false, ExcPureFunctionCalled());
};



/*------------------------------- FiniteElement ----------------------*/


bool FiniteElement<1>::operator == (const FiniteElement<1> &f) const {
  return ((dofs_per_vertex == f.dofs_per_vertex) &&
	  (dofs_per_line == f.dofs_per_line) &&
	  (FiniteElementBase<1>::operator == (f)));
};



void FiniteElement<1>::fill_fe_values (const Triangulation<1>::cell_iterator &cell,
				       const vector<Point<1> > &unit_points,
				       vector<dFMatrix>  &jacobians,
				       const bool         compute_jacobians,
				       vector<Point<1> > &ansatz_points,
				       const bool         compute_ansatz_points,
				       vector<Point<1> > &q_points,
				       const bool         compute_q_points) const {
				   // local mesh width
  const double h=(cell->vertex(1)(0) - cell->vertex(0)(0));

  for (unsigned int i=0; i<q_points.size(); ++i) 
    {
      if (compute_jacobians)
	jacobians[i](0,0) = 1./h;
      if (compute_q_points)
	q_points[i] = cell->vertex(0) + h*unit_points[i];
    };

				   // compute ansatz points. The first ones
				   // belong to vertex one, the second ones
				   // to vertex two, all following are
				   // equally spaced along the line
  if (compute_ansatz_points) 
    {
      ansatz_points.erase (ansatz_points.begin(), ansatz_points.end());
      ansatz_points.reserve (total_dofs);

				       // first the dofs in the vertices
      for (unsigned int vertex=0; vertex<2; vertex++) 
	for (unsigned int i=0; i<dofs_per_vertex; ++i)
	  ansatz_points.push_back (cell->vertex(vertex));
      
				       // now dofs on line
      for (unsigned int i=0; i<dofs_per_line; ++i) 
	ansatz_points.push_back (cell->vertex(0) +
				 Point<1>((i+1.0)/(total_dofs+1.0)*h));
    };
};



bool FiniteElement<2>::operator == (const FiniteElement<2> &f) const {
  return ((dofs_per_vertex == f.dofs_per_vertex) &&
	  (dofs_per_line == f.dofs_per_line) &&
	  (dofs_per_quad == f.dofs_per_quad) &&
	  (FiniteElementBase<2>::operator == (f)));
};



void FiniteElement<2>::fill_fe_values (const Triangulation<2>::cell_iterator &,
				       const vector<Point<2> > &,
				       vector<dFMatrix> &,
				       const bool,
				       vector<Point<2> > &,
				       const bool,
				       vector<Point<2> > &,
				       const bool) const {
  Assert (false, ExcPureFunctionCalled());
};





/*------------------------------- Explicit Instantiations -------------*/

template struct FEValues<1>::UpdateStruct;
template struct FEValues<2>::UpdateStruct;

template class FEValues<1>;
template class FEValues<2>;

template class FiniteElementBase<1>;
template class FiniteElementBase<2>;

template class FiniteElement<1>;
template class FiniteElement<2>;
