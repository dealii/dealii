/* $Id$ */

#include <fe/fe.h>
#include <fe/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary.h>





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
					     const bool,
					     const Boundary<dim> &) const {
  Assert (false, ExcPureFunctionCalled());
};



template <int dim>
void FiniteElementBase<dim>::face_ansatz_points (const typename Triangulation<dim>::face_iterator &,
						 const Boundary<dim> &,
						 vector<Point<dim> > &) const {
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
				       const bool         compute_q_points,
				       const Boundary<1> &) const {
  Assert (jacobians.size() == unit_points.size(),
	  ExcWrongFieldDimension(jacobians.size(), unit_points.size()));
  Assert (q_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(q_points.size(), unit_points.size()));
  Assert (ansatz_points.size() == total_dofs,
	  ExcWrongFieldDimension(ansatz_points.size(), total_dofs));


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
      unsigned int next = 0;
				       // first the dofs in the vertices
      for (unsigned int vertex=0; vertex<2; vertex++) 
	for (unsigned int i=0; i<dofs_per_vertex; ++i)
	  ansatz_points[next++] = cell->vertex(vertex);
      
				       // now dofs on line
      for (unsigned int i=0; i<dofs_per_line; ++i) 
	ansatz_points[next++] = cell->vertex(0) +
				Point<1>((i+1.0)/(total_dofs+1.0)*h);
    };
};



void FiniteElement<1>::face_ansatz_points (const typename Triangulation<1>::face_iterator &,
					   const Boundary<1> &,
					   vector<Point<1> > &) const {
				   // is this function useful in 1D?
  Assert (false, ExcPureFunctionCalled());
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
				       const bool,
				       const Boundary<2> &) const {
  Assert (false, ExcPureFunctionCalled());
};



void FiniteElement<2>::face_ansatz_points (const typename Triangulation<2>::face_iterator &,
					   const Boundary<2> &,
					   vector<Point<2> > &) const {
				   // is this function useful in 1D?
  Assert (false, ExcPureFunctionCalled());
};





/*------------------------------- Explicit Instantiations -------------*/

template class FiniteElementBase<1>;
template class FiniteElementBase<2>;

template class FiniteElement<1>;
template class FiniteElement<2>;
