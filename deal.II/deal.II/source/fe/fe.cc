/* $Id$ */

#include <fe/fe.h>
#include <fe/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/tria_boundary.h>





/*------------------------------- FiniteElementData ----------------------*/


bool FiniteElementData<1>::operator== (const FiniteElementData<1> &f) const {
  return ((dofs_per_vertex == f.dofs_per_vertex) &&
	  (dofs_per_line == f.dofs_per_line) &&
	  (total_dofs == f.total_dofs));
};



bool FiniteElementData<2>::operator== (const FiniteElementData<2> &f) const {
  return ((dofs_per_vertex == f.dofs_per_vertex) &&
	  (dofs_per_line == f.dofs_per_line) &&
	  (dofs_per_quad == f.dofs_per_quad) &&
	  (total_dofs == f.total_dofs));
};





/*------------------------------- FiniteElementBase ----------------------*/




FiniteElementBase<1>::FiniteElementBase (const unsigned int dofs_per_vertex,
					 const unsigned int dofs_per_line,
					 const unsigned int dofs_per_quad) :
		FiniteElementData<1> (dofs_per_vertex,
				      dofs_per_line)
{
  Assert (dofs_per_quad==0, ExcInternalError());

  const unsigned int dim=1;
  for (unsigned int i=0; i<(1<<dim); ++i) 
    {
      restriction[i].reinit (total_dofs, total_dofs);
      prolongation[i].reinit (total_dofs, total_dofs);
    };
  interface_constraints.reinit (1,1);
  interface_constraints(0,0)=1.;
};



FiniteElementBase<2>::FiniteElementBase (const unsigned int dofs_per_vertex,
					 const unsigned int dofs_per_line,
					 const unsigned int dofs_per_quad) :
		FiniteElementData<2> (dofs_per_vertex,
				      dofs_per_line,
				      dofs_per_quad)
{
  const unsigned int dim=2;
  for (unsigned int i=0; i<(1<<dim); ++i) 
    {
      restriction[i].reinit (total_dofs, total_dofs);
      prolongation[i].reinit (total_dofs, total_dofs);
    };
  interface_constraints.reinit (dofs_per_vertex+2*dofs_per_line,
				2*dofs_per_vertex+dofs_per_line);
};



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
  if (dim==1)
    Assert ((interface_constraints.m()==1) && (interface_constraints.n()==1),
	    ExcWrongInterfaceMatrixSize(interface_constraints.m(),
					interface_constraints.n()));
  
  return interface_constraints;
};



template <int dim>
bool FiniteElementBase<dim>::operator == (const FiniteElementBase<dim> &f) const {
  return ((static_cast<FiniteElementData<dim> >(*this) ==
	   static_cast<FiniteElementData<dim> >(f)) &&
	  (interface_constraints == f.interface_constraints));
};





/*------------------------------- FiniteElement ----------------------*/


void FiniteElement<1>::fill_fe_values (const DoFHandler<1>::cell_iterator &cell,
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


template <int dim>
void FiniteElement<dim>::fill_fe_values (const DoFHandler<dim>::cell_iterator &,
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



void FiniteElement<1>::fill_fe_face_values (const DoFHandler<1>::cell_iterator &,
					    const unsigned int       ,
					    const vector<Point<0> > &,
					    const vector<Point<1> > &,
					    vector<dFMatrix>        &,
					    const bool               ,
					    vector<Point<1> >       &,
					    const bool               ,
					    vector<Point<1> >       &,
					    const bool               ,
					    vector<double>          &,
					    const bool              ,
					    const Boundary<1>       &) const {
  Assert (false, ExcNotImplemented());
}



template <int dim>
void FiniteElement<dim>::fill_fe_face_values (const DoFHandler<dim>::cell_iterator &cell,
					      const unsigned int           face_no,
					      const vector<Point<dim-1> > &unit_points,
					      const vector<Point<dim> > &global_unit_points,
					      vector<dFMatrix>    &jacobians,
					      const bool           compute_jacobians,
					      vector<Point<dim> > &ansatz_points,
					      const bool           compute_ansatz_points,
					      vector<Point<dim> > &q_points,
					      const bool           compute_q_points,
					      vector<double>      &face_jacobi_determinants,
					      const bool           compute_face_jacobians,
					      const Boundary<dim> &boundary) const {
  Assert (jacobians.size() == unit_points.size(),
	  ExcWrongFieldDimension(jacobians.size(), unit_points.size()));
  Assert (q_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(q_points.size(), unit_points.size()));
  Assert (global_unit_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(global_unit_points.size(), unit_points.size()));
  Assert (ansatz_points.size() == dofs_per_face,
	  ExcWrongFieldDimension(ansatz_points.size(), dofs_per_face));
  
  static vector<Point<dim> > dummy(total_dofs);
  fill_fe_values (cell, global_unit_points,
		  jacobians, compute_jacobians,
		  dummy, false,
		  q_points, compute_q_points,
		  boundary);
  
  cout << "Global unit points:\n";
  for (unsigned int p=0; p<unit_points.size(); ++p)
    cout << "    " << global_unit_points[p] << endl;

  if (compute_ansatz_points)
    get_face_ansatz_points (cell->face(face_no), boundary, ansatz_points);

  if (compute_face_jacobians)
    get_face_jacobians (cell->face(face_no), boundary,
			unit_points, face_jacobi_determinants);
  
  cout << "Global ansatz points:\n";
  for (unsigned int p=0; p<unit_points.size(); ++p)
    cout << "    " << ansatz_points[p] << endl;
};



/*------------------------------- Explicit Instantiations -------------*/

template class FiniteElementData<1>;
template class FiniteElementData<2>;

template class FiniteElementBase<1>;
template class FiniteElementBase<2>;

template class FiniteElement<1>;
template class FiniteElement<2>;



