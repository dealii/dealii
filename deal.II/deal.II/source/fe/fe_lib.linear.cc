/* $Id$ */

#include <fe/fe_lib.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>

extern TriaIterator<1,CellAccessor<1> > __dummy2687; // for gcc2.8
extern TriaIterator<2,CellAccessor<2> > __dummy2688; // for gcc2.8


FELinear<1>::FELinear () :
		FiniteElement<1> (1, 0)
{
  restriction[0].reinit (2,2);
  restriction[1].reinit (2,2);

				   // for restriction and prolongation matrices:
				   // note that we do not add up all the
				   // contributions since then we would get
				   // two summands per vertex in 1d (four
				   // in 2d, etc), but only one per line dof.
				   // We could accomplish for that by dividing
				   // the vertex dof values by 2 (4, etc), but
				   // would get into trouble at the boundary
				   // of the domain since there only one
				   // cell contributes to a vertex. Rather,
				   // we do not add up the contributions but
				   // set them right into the matrices!
  restriction[0](0,0) = 1.0;
  restriction[0](0,1) = 1./2.;
  restriction[0](1,1) = 1./2.;

  restriction[1](0,0) = 1./2.;
  restriction[1](1,0) = 1./2.;
  restriction[1](1,1) = 1.0;


  prolongation[0].reinit (2,2);
  prolongation[1].reinit (2,2);
  
  prolongation[0](0,0) = 1.0;
  prolongation[0](1,0) = 1./2.;
  prolongation[0](1,1) = 1./2.;

  prolongation[1](0,0) = 1./2.;
  prolongation[1](0,1) = 1./2.;
  prolongation[1](1,1) = 1.0;
};



double
FELinear<1>::shape_value(const unsigned int i,
			       const Point<1>& p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  switch (i)
    {
    case 0: return 1.-p(0);
    case 1: return p(0);
    }
  return 0.;
}



Point<1>
FELinear<1>::shape_grad(const unsigned int i,
			const Point<1>&) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  switch (i)
    {
    case 0: return Point<1>(-1.);
    case 1: return Point<1>(1.);
    }
  return Point<1>();
}



void FELinear<1>::fill_fe_values (const Triangulation<1>::cell_iterator &cell,
				  const vector<Point<1> >               &unit_points,
				  vector<dFMatrix>  &jacobians,
				  vector<Point<1> > &points) const {
				   // simply pass down
  FiniteElement<1>::fill_fe_values (cell, unit_points, jacobians, points);
};





FELinear<2>::FELinear () :
		FiniteElement<2> (1, 0, 0)
{
  interface_constraints.reinit(1,2);
  interface_constraints(0,0) = 1./2.;
  interface_constraints(0,1) = 1./2.;

  restriction[0].reinit(4,4);
  restriction[1].reinit(4,4);
  restriction[2].reinit(4,4);
  restriction[3].reinit(4,4);

  prolongation[0].reinit(4,4);
  prolongation[1].reinit(4,4);
  prolongation[2].reinit(4,4);
  prolongation[3].reinit(4,4);

  restriction[0](0,0) = 1.0;
  restriction[0](0,1) = 1./2.;
  restriction[0](1,1) = 1./2.;
  restriction[0](0,3) = 1./2.;
  restriction[0](3,3) = 1./2.;
  restriction[0](0,2) = 1./4.;
  restriction[0](1,2) = 1./4.;
  restriction[0](2,2) = 1./4.;
  restriction[0](3,2) = 1./4.;

  restriction[1](1,1) = 1.0;
  restriction[1](0,0) = 1./2.;
  restriction[1](1,0) = 1./2.;
  restriction[1](1,2) = 1./2.;
  restriction[1](2,2) = 1./2.;
  restriction[1](0,3) = 1./4.;
  restriction[1](1,3) = 1./4.;
  restriction[1](2,3) = 1./4.;
  restriction[1](3,3) = 1./4.;

  restriction[2](2,2) = 1.0;
  restriction[2](2,1) = 1./2.;
  restriction[2](1,1) = 1./2.;
  restriction[2](2,3) = 1./2.;
  restriction[2](3,3) = 1./2.;
  restriction[2](0,0) = 1./4.;
  restriction[2](1,0) = 1./4.;
  restriction[2](2,0) = 1./4.;
  restriction[2](3,0) = 1./4.;

  restriction[3](3,3) = 1.0;
  restriction[3](0,0) = 1./2.;
  restriction[3](3,0) = 1./2.;
  restriction[3](2,2) = 1./2.;
  restriction[3](3,2) = 1./2.;
  restriction[3](0,1) = 1./4.;
  restriction[3](1,1) = 1./4.;
  restriction[3](2,1) = 1./4.;
  restriction[3](3,1) = 1./4.;

  prolongation[0](0,0) = 1.0;
  prolongation[0](1,0) = 1./2.;
  prolongation[0](1,1) = 1./2.;
  prolongation[0](3,0) = 1./2.;
  prolongation[0](3,3) = 1./2.;
  prolongation[0](2,0) = 1./4.;
  prolongation[0](2,1) = 1./4.;
  prolongation[0](2,2) = 1./4.;
  prolongation[0](2,3) = 1./4.;

  prolongation[1](1,1) = 1.0;
  prolongation[1](0,0) = 1./2.;
  prolongation[1](0,1) = 1./2.;
  prolongation[1](2,1) = 1./2.;
  prolongation[1](2,2) = 1./2.;
  prolongation[1](3,0) = 1./4.;
  prolongation[1](3,1) = 1./4.;
  prolongation[1](3,2) = 1./4.;
  prolongation[1](3,3) = 1./4.;

  prolongation[2](2,2) = 1.0;
  prolongation[2](1,2) = 1./2.;
  prolongation[2](1,1) = 1./2.;
  prolongation[2](3,2) = 1./2.;
  prolongation[2](3,3) = 1./2.;
  prolongation[2](0,0) = 1./4.;
  prolongation[2](0,1) = 1./4.;
  prolongation[2](0,2) = 1./4.;
  prolongation[2](0,3) = 1./4.;

  prolongation[3](3,3) = 1.0;
  prolongation[3](0,0) = 1./2.;
  prolongation[3](0,3) = 1./2.;
  prolongation[3](2,2) = 1./2.;
  prolongation[3](2,3) = 1./2.;
  prolongation[3](1,0) = 1./4.;
  prolongation[3](1,1) = 1./4.;
  prolongation[3](1,2) = 1./4.;
  prolongation[3](1,3) = 1./4.;
};



double
FELinear<2>::shape_value (const unsigned int i,
			  const Point<2>& p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  switch (i)
    {
    case 0: return (1.-p(0)) * (1.-p(1));
    case 1: return p(0) * (1.-p(1));
    case 2: return p(0) * p(1);
    case 3: return (1.-p(0)) * p(1);
    }
  return 0.;
}



Point<2>
FELinear<2>::shape_grad (const unsigned int i,
			 const Point<2>& p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  switch (i)
    {
    case 0: return Point<2> (p(1)-1., p(0)-1.);
    case 1: return Point<2> (1.-p(1), -p(0));
    case 2: return Point<2> (p(1), p(0));
    case 3: return Point<2> (-p(1), 1.-p(0));
    }
  return Point<2> ();
}



// this function may be generalised to three or more dimensions with gcc2.8
// you will have to change th number of vertices
void FELinear<2>::fill_fe_values (const Triangulation<2>::cell_iterator &cell,
				  const vector<Point<2> >               &unit_points,
				  vector<dFMatrix>  &jacobians,
				  vector<Point<2> > &points) const {
  const unsigned int dim=2;
  const unsigned int n_vertices=4;
  
  unsigned int n_points=unit_points.size();

				   // initialize points to zero
  for (unsigned int i=0; i<n_points; ++i)
    points[i] = Point<dim> ();
  
				   // note: let x_l be the vector of the
				   // lth quadrature point in real space and
				   // xi_l that on the unit cell, let further
				   // p_j be the vector of the jth vertex
				   // of the cell in real space and
				   // N_j(xi_l) be the value of the associated
				   // basis function at xi_l, then
				   // x_l(xi_l) = sum_j p_j N_j(xi_l)
  for (unsigned int j=0; j<n_vertices; ++j) 
    for (unsigned int l=0; l<n_points; ++l) 
      points[l] += cell->vertex(j) * shape_value(j, unit_points[l]);

// computation of jacobian still missing
};






FEQuadratic<1>::FEQuadratic () :
		FiniteElement<1> (1, 1) {};



void FEQuadratic<1>::fill_fe_values (const Triangulation<1>::cell_iterator &cell,
				     const vector<Point<1> >               &unit_points,
				     vector<dFMatrix>  &jacobians,
				     vector<Point<1> > &points) const {
				   // simply pass down
  FiniteElement<1>::fill_fe_values (cell, unit_points, jacobians, points);
};



FEQuadratic<2>::FEQuadratic () :
		FiniteElement<2> (1, 1, 1)
{
  interface_constraints.reinit(3,3);
  interface_constraints(0,2) = 1.0;
  interface_constraints(1,0) = 3./8.;
  interface_constraints(1,1) = -1./8.;
  interface_constraints(1,2) = 3./4.;
  interface_constraints(2,0) = -1./8.;
  interface_constraints(2,1) = 3./8.;
  interface_constraints(2,2) = 3./4.;
};





FECubic<1>::FECubic () :
		FiniteElement<1> (1, 2) {};



void FECubic<1>::fill_fe_values (const Triangulation<1>::cell_iterator &cell,
				 const vector<Point<1> >               &unit_points,
				 vector<dFMatrix>  &jacobians,
				 vector<Point<1> > &points) const {
				   // simply pass down
  FiniteElement<1>::fill_fe_values (cell, unit_points, jacobians, points);
};




FECubic<2>::FECubic () :
		FiniteElement<2> (1, 2, 4) {};

