/* $Id$ */

#include <fe/fe.h>
#include <fe/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>

extern TriaIterator<1,CellAccessor<1> > __dummy2687; // for gcc2.8



/*------------------------------- FEValues -------------------------------*/


template <int dim>
FEValues<dim>::FEValues (const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &quadrature) :
		shape_values(fe.total_dofs, quadrature.n_quad_points),
		shape_gradients(fe.total_dofs,
				vector<Point<dim> >(quadrature.n_quad_points)),
		unit_shape_gradients(fe.total_dofs,
				     vector<Point<dim> >(quadrature.n_quad_points)),
		weights(quadrature.n_quad_points, 0),
		JxW_values(quadrature.n_quad_points, 0),
		quadrature_points(quadrature.n_quad_points, Point<dim>()),
		unit_quadrature_points(quadrature.n_quad_points, Point<dim>()),
		jacobi_matrices (quadrature.n_quad_points)
{
  for (unsigned int i=0; i<fe.total_dofs; ++i)
    for (unsigned int j=0; j<quadrature.n_quad_points; ++j) 
      {
	shape_values(i,j) = fe.shape_value(i, quadrature.quad_point(j));
	unit_shape_gradients[i][j]
	  = fe.shape_grad(i, quadrature.quad_point(j));
      };

  for (unsigned int i=0; i<weights.size(); ++i) 
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

  return shape_gradients[i][j];
};



void FEValues<1>::reinit (const Triangulation<1>::cell_iterator &cell,
			  const FiniteElement<1>                &fe) {
  const unsigned int dim=1;
				   // fill jacobi matrices and real
				   //quadrature points
  fe.fill_fe_values (cell,
		     unit_quadrature_points,
		     jacobi_matrices,
		     quadrature_points);
				   // compute gradients on real element
  for (unsigned int i=0; i<fe.total_dofs; ++i)
    for (unsigned int j=0; j<quadrature_points.size(); ++j)
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
  
				   // compute Jacobi determinants in
				   // quadrature points.
				   // refer to the general doc for
				   // why we take the inverse of the
				   // determinant
  for (unsigned int i=0; i<quadrature_points.size(); ++i)
    JxW_values[i] = weights[i] / jacobi_matrices[i].determinant();
};



void FEValues<2>::reinit (const Triangulation<2>::cell_iterator &cell,
			  const FiniteElement<2>                &fe) {
  const unsigned int dim=2;
				   // fill jacobi matrices and real
				   //quadrature points
  fe.fill_fe_values (cell,
		     unit_quadrature_points,
		     jacobi_matrices,
		     quadrature_points);
				   // compute gradients on real element
  for (unsigned int i=0; i<fe.total_dofs; ++i)
    for (unsigned int j=0; j<quadrature_points.size(); ++j)
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
  
				   // refer to the general doc for
				   // why we take the inverse of the
				   // determinant
  for (unsigned int i=0; i<quadrature_points.size(); ++i)
    JxW_values[i] = weights[i] / jacobi_matrices[i].determinant();
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



void FiniteElementBase<1>::fill_fe_values (const Triangulation<1>::cell_iterator &,
					   const vector<Point<1> > &,
					   vector<dFMatrix> &,
					   vector<Point<1> > &) const {
  Assert (false, ExcPureFunctionCalled());
};



void FiniteElementBase<2>::fill_fe_values (const Triangulation<2>::cell_iterator &,
					   const vector<Point<2> > &,
					   vector<dFMatrix> &,
					   vector<Point<2> > &) const {
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
				       vector<Point<1> > &points) const {
				   // local mesh width
  double h=(cell->vertex(1)(0) - cell->vertex(0)(0));

  unsigned int n_points = unit_points.size();
  for (unsigned int i=0; i<n_points; ++i) 
    {
      jacobians[i](0,0) = 1./h;
      points[i] = cell->vertex(0) + h*unit_points[i];
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
				       vector<Point<2> > &) const {
  Assert (false, ExcPureFunctionCalled());
};





/*------------------------------- Explicit Instantiations -------------*/

template class FEValues<1>;
template class FEValues<2>;

template class FiniteElementBase<1>;
template class FiniteElementBase<2>;

template class FiniteElement<1>;
template class FiniteElement<2>;
