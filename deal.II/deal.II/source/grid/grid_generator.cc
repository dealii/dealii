//----------------------------  grid_generator.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_generator.cc  ---------------------------


#include <grid/grid_generator.h>
#include <grid/tria.h>
#include <cmath>


template <int dim>
void GridGenerator::hyper_rectangle (Triangulation<dim>& tria,
				     const Point<dim>& p_1,
				     const Point<dim>& p_2)
{
				   // First, normalize input such that
				   // p1 is lower in all coordinate directions.
  Point<dim> p1(p_1);
  Point<dim> p2(p_2);
  
  for (unsigned int i=0;i<dim;++i)
    if (p1(i) > p2(i))
      {
	double h = p1(i);
	p1(i) = p2(i);
	p2(i) = h;
      }
  
  std::vector<Point<dim> > vertices (GeometryInfo<dim>::vertices_per_cell);
  switch (dim)
    {
      case 1:
	vertices[0](0) = p1(0);
	vertices[1](0) = p2(0);
	break;
      case 2:
	vertices[0](0) = p1(0);
	vertices[1](0) = p2(0);
	vertices[2](0) = p2(0);
	vertices[3](0) = p1(0);
	vertices[0](1) = p1(1);
	vertices[1](1) = p1(1);
	vertices[2](1) = p2(1);
	vertices[3](1) = p2(1);
	break;
      case 3:
	vertices[0](0) = p1(0);
	vertices[1](0) = p2(0);
	vertices[2](0) = p2(0);
	vertices[3](0) = p1(0);
	vertices[4](0) = p1(0);
	vertices[5](0) = p2(0);
	vertices[6](0) = p2(0);
	vertices[7](0) = p1(0);
	vertices[0](1) = p1(0);
	vertices[1](1) = p1(0);
	vertices[2](1) = p1(0);
	vertices[3](1) = p1(0);
	vertices[4](1) = p2(0);
	vertices[5](1) = p2(0);
	vertices[6](1) = p2(0);
	vertices[7](1) = p2(0);
	vertices[0](2) = p1(0);
	vertices[1](2) = p1(0);
	vertices[2](2) = p2(0);
	vertices[3](2) = p2(0);
	vertices[4](2) = p1(0);
	vertices[5](2) = p1(0);
	vertices[6](2) = p2(0);
	vertices[7](2) = p2(0);
	break;
      default:
	Assert (false, ExcNotImplemented ());
    }

				   // Prepare cell data
  std::vector<CellData<dim> > cells (1);
  for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i)
    cells[0].vertices[i] = i;
  cells[0].material_id = 0;

  tria.create_triangulation (vertices, cells, SubCellData());
}



template <int dim>
void GridGenerator::hyper_cube(Triangulation<dim> &tria,
			       const double left,
			       const double right)
{
  Point<dim> p1;
  Point<dim> p2;
  for (unsigned int i=0;i<dim;++i)
    {
      p1(i) = left;
      p2(i) = right;
    }
  hyper_rectangle (tria, p1, p2);
}


#if deal_II_dimension == 1

template <>
void GridGenerator::hyper_cube_slit<> (Triangulation<1> &,
				       const double,
				       const double) {
  Assert (false, ExcInternalError());
};



template <>
void GridGenerator::hyper_L<> (Triangulation<1> &,
			       const double,
			       const double) {
  Assert (false, ExcInternalError());
};



template <>
void GridGenerator::hyper_ball<> (Triangulation<1> &,
				  const Point<1> &,
				  const double) {
  Assert (false, ExcInternalError());
};



template <>
void GridGenerator::hyper_shell<> (Triangulation<1> &,
				   const Point<1> &,
				   const double,
				   const double,
				   const unsigned int) {
  Assert (false, ExcInternalError());
};

#endif



#if deal_II_dimension == 2

template<>
void GridGenerator::enclosed_hyper_cube (Triangulation<2> &tria,
					 const double l,
					 const double r,
					 const double d,
					 bool colorize)
{
  std::vector<Point<2> > vertices(16);
  double coords[4];
  coords[0] = l-d;
  coords[1] = l;
  coords[2] = r;
  coords[3] = r+d;

  unsigned int k=0;
  for (unsigned int i0=0;i0<4;++i0)
    for (unsigned int i1=0;i1<4;++i1)
      vertices[k++] = Point<2>(coords[i1], coords[i0]);

  const unsigned char materials[9] = { 5, 4, 6,
				       1, 0, 2,
				       9, 8,10
  };
  
  std::vector<CellData<2> > cells(9);
  k = 0;
  for (unsigned int i0=0;i0<3;++i0)
    for (unsigned int i1=0;i1<3;++i1)
      {
	cells[k].vertices[0] = i1+4*i0;
	cells[k].vertices[1] = i1+4*i0+1;
	cells[k].vertices[2] = i1+4*i0+5;
	cells[k].vertices[3] = i1+4*i0+4;
	if (colorize)
	  cells[k].material_id = materials[k];
	++k;
      }
  tria.create_triangulation (vertices,
			     cells,
			     SubCellData());       // no boundary information
}



template <>
void GridGenerator::hyper_cube_slit<> (Triangulation<2> &tria,
				       const double left,
				       const double right) {
  const double rl2=(right+left)/2;
  const Point<2> vertices[10] = { Point<2>(left, left ),
				    Point<2>(rl2,  left ),
				    Point<2>(rl2,  rl2  ),
				    Point<2>(left, rl2  ),
				    Point<2>(right,left ),
				    Point<2>(right,rl2  ),
				    Point<2>(rl2,  right),
				    Point<2>(left, right),
				    Point<2>(right,right),
				    Point<2>(rl2,  left ) };
  const int cell_vertices[4][4] = { { 0,1,2,3 },
				    { 9,4,5,2 },
				    { 3,2,6,7 },
				    { 2,5,8,6 } };
  std::vector<CellData<2> > cells (4, CellData<2>());
  for (unsigned int i=0; i<4; ++i)
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  tria.create_triangulation (std::vector<Point<2> >(&vertices[0], &vertices[10]),
			     cells,
			     SubCellData());       // no boundary information
};



template <>
void GridGenerator::hyper_L<> (Triangulation<2> &tria,
			       const double a,
			       const double b) {
  const unsigned int dim=2;
  const Point<dim> vertices[8] = { Point<dim> (a,a),
				     Point<dim> ((a+b)/2,a),
				     Point<dim> (b,a),
				     Point<dim> (a,(a+b)/2),
				     Point<dim> ((a+b)/2,(a+b)/2),
				     Point<dim> (b,(a+b)/2),
				     Point<dim> (a,b),
				     Point<dim> ((a+b)/2,b)  };
  const int cell_vertices[3][4] = {{0, 1, 4, 3},
				   {1, 2, 5, 4},
				   {3, 4, 7, 6}};

  std::vector<CellData<2> > cells (3, CellData<2>());
  
  for (unsigned int i=0; i<3; ++i) 
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (std::vector<Point<dim> >(&vertices[0], &vertices[8]),
			     cells,
			     SubCellData());       // no boundary information
};



template <>
void GridGenerator::hyper_ball<> (Triangulation<2> &tria,
				  const Point<2>    &p,
				  const double      radius) {
  const double a = 1./(1+sqrt(2.0)); // equilibrate cell sizes at transition
				   // from the inner part to the radial
				   // cells
  const Point<2> vertices[8] = { p+Point<2>(-1,-1)*(radius/sqrt(2.0)),
				   p+Point<2>(+1,-1)*(radius/sqrt(2.0)),
				   p+Point<2>(-1,-1)*(radius/sqrt(2.0)*a),
				   p+Point<2>(+1,-1)*(radius/sqrt(2.0)*a),
				   p+Point<2>(-1,+1)*(radius/sqrt(2.0)*a),
				   p+Point<2>(+1,+1)*(radius/sqrt(2.0)*a),
				   p+Point<2>(-1,+1)*(radius/sqrt(2.0)),
				   p+Point<2>(+1,+1)*(radius/sqrt(2.0)) };
  
  const int cell_vertices[5][4] = {{0, 1, 3, 2},
				   {0, 2, 4, 6},
				   {2, 3, 5, 4},
				   {1, 7, 5, 3},
				   {6, 4, 5, 7}};

  std::vector<CellData<2> > cells (5, CellData<2>());
  
  for (unsigned int i=0; i<5; ++i) 
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (std::vector<Point<2> >(&vertices[0], &vertices[8]),
			     cells,
			     SubCellData());       // no boundary information
};



template <>
void GridGenerator::hyper_shell<> (Triangulation<2>   &tria,
				   const Point<2>     &center,
				   const double        inner_radius,
				   const double        outer_radius,
				   const unsigned int  n_cells)
{
  Assert ((inner_radius > 0) && (inner_radius < outer_radius),
	  ExcInvalidRadii ());
  
//TODO:[?] Unify the various places where PI is defined to a central instance  
  const double pi = 3.141592653589793238462;
  
				   // determine the number of cells
				   // for the grid. if not provided by
				   // the user determine it such that
				   // the length of each cell on the
				   // median (in the middle between
				   // the two circles) is equal to its
				   // radial extent (which is the
				   // difference between the two
				   // radii)
  const unsigned int N = (n_cells == 0 ?
			  static_cast<unsigned int>
			  (ceil((2*pi* (outer_radius + inner_radius)/2) /
				(outer_radius - inner_radius))) :
			  n_cells);

				   // set up N vertices on the
				   // outer and N vertices on
				   // the inner circle. the
				   // first N ones are on the
				   // outer one, and all are
				   // numbered counter-clockwise
  std::vector<Point<2> > vertices(2*N);
  for (unsigned int i=0; i<N; ++i)
    {
      vertices[i] = Point<2>( cos(2*pi * i/N),
			      sin(2*pi * i/N)) * outer_radius;
      vertices[i+N] = vertices[i] * (inner_radius/outer_radius);

      vertices[i]   += center;
      vertices[i+N] += center;
    };

  std::vector<CellData<2> > cells (N, CellData<2>());
	
  for (unsigned int i=0; i<N; ++i) 
    {
      cells[i].vertices[0] = i;
      cells[i].vertices[1] = (i+1)%N;
      cells[i].vertices[2] = N+((i+1)%N);
      cells[i].vertices[3] = N+i;
	    
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (vertices, cells, SubCellData());
};



template <>
void GridGenerator::half_hyper_shell<> (Triangulation<2>   &tria,
					const Point<2>     &center,
					const double        inner_radius,
					const double        outer_radius,
					const unsigned int  n_cells)
{
  Assert ((inner_radius > 0) && (inner_radius < outer_radius),
	  ExcInvalidRadii ());
  
  const double pi     = 3.14159265359;
				   // determine the number of cells
				   // for the grid. if not provided by
				   // the user determine it such that
				   // the length of each cell on the
				   // median (in the middle between
				   // the two circles) is equal to its
				   // radial extent (which is the
				   // difference between the two
				   // radii)
  const unsigned int N = (n_cells == 0 ?
			  static_cast<unsigned int>
			  (ceil((pi* (outer_radius + inner_radius)/2) /
				(outer_radius - inner_radius))) :
			  n_cells);

				   // set up N+1 vertices on the
				   // outer and N+1 vertices on
				   // the inner circle. the
				   // first N+1 ones are on the
				   // outer one, and all are
				   // numbered counter-clockwise
  std::vector<Point<2> > vertices(2*(N+1));
  for (unsigned int i=0; i<=N; ++i)
    {
				       // enforce that the x-coordinates
				       // of the first and last point of
				       // each half-circle are exactly
				       // zero (contrary to what we may
				       // compute using the imprecise
				       // value of pi)
      vertices[i] =  Point<2>( ( (i==0) || (i==N) ?
				 0 :
				 cos(pi * i/N - pi/2) ),
			       sin(pi * i/N - pi/2)) * outer_radius;
      vertices[i+N+1] = vertices[i] * (inner_radius/outer_radius);

      vertices[i]     += center;
      vertices[i+N+1] += center;
    };


  std::vector<CellData<2> > cells (N, CellData<2>());
	
  for (unsigned int i=0; i<N; ++i) 
    {
      cells[i].vertices[0] = i;
      cells[i].vertices[1] = (i+1)%(N+1);
      cells[i].vertices[2] = N+1+((i+1)%(N+1));
      cells[i].vertices[3] = N+1+i;
	    
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (vertices, cells, SubCellData());
};



#endif


#if deal_II_dimension == 3

template <>
void GridGenerator::hyper_cube_slit<> (Triangulation<3> &,
				       const double,
				       const double) {
  Assert (false, ExcNotImplemented());
};



template<>
void GridGenerator::enclosed_hyper_cube (Triangulation<3> &tria,
					 const double l,
					 const double r,
					 const double d,
					 bool colorize)
{
  std::vector<Point<3> > vertices(64);
  double coords[4];
  coords[0] = l-d;
  coords[1] = l;
  coords[2] = r;
  coords[3] = r+d;

  unsigned int k=0;
  for (unsigned int i0=0;i0<4;++i0)
    for (unsigned int i1=0;i1<4;++i1)
      for (unsigned int i2=0;i2<4;++i2)
	vertices[k++] = Point<3>(coords[i2], coords[i1], coords[i0]);

  const unsigned char materials[27] = {
	21,20,22,
	17,16,18,
	25,24,26,
	5 , 4, 6,
	1 , 0, 2,
	9 , 8,10,
	37,36,38,
	33,32,34,
	41,40,42
  };
  
  std::vector<CellData<3> > cells(27);
  k = 0;
  for (unsigned int i0=0;i0<3;++i0)
    for (unsigned int i1=0;i1<3;++i1)
      for (unsigned int i2=0;i2<3;++i2)
	{
	  cells[k].vertices[0] = i2+4*i1+16*i0;
	  cells[k].vertices[1] = i2+4*i1+16*i0+1;
	  cells[k].vertices[2] = i2+4*i1+16*i0+5;
	  cells[k].vertices[3] = i2+4*i1+16*i0+4;
	  cells[k].vertices[4] = i2+4*i1+16*i0+16;
	  cells[k].vertices[5] = i2+4*i1+16*i0+17;
	  cells[k].vertices[6] = i2+4*i1+16*i0+21;
	  cells[k].vertices[7] = i2+4*i1+16*i0+20;
	  if (colorize)
	    cells[k].material_id = materials[k];
	  ++k;
	}
  tria.create_triangulation (vertices,
			     cells,
			     SubCellData());       // no boundary information
}
  
template <>
void GridGenerator::hyper_L<> (Triangulation<3> &tria,
			       const double a,
			       const double b) {
  const unsigned int dim=3;
				   // we slice out the top back right
				   // part of the cube
  const Point<dim> vertices[26]
    = {
				     // front face of the big cube
    Point<dim> (a,      a,a),
      Point<dim> ((a+b)/2,a,a),
      Point<dim> (b,      a,a),
      Point<dim> (a,      a,(a+b)/2),
      Point<dim> ((a+b)/2,a,(a+b)/2),
      Point<dim> (b,      a,(a+b)/2),
      Point<dim> (a,      a,b),
      Point<dim> ((a+b)/2,a,b),
      Point<dim> (b,      a,b),
				       // middle face of the big cube
      Point<dim> (a,      (a+b)/2,a),
      Point<dim> ((a+b)/2,(a+b)/2,a),
      Point<dim> (b,      (a+b)/2,a),
      Point<dim> (a,      (a+b)/2,(a+b)/2),
      Point<dim> ((a+b)/2,(a+b)/2,(a+b)/2),
      Point<dim> (b,      (a+b)/2,(a+b)/2),
      Point<dim> (a,      (a+b)/2,b),
      Point<dim> ((a+b)/2,(a+b)/2,b),
      Point<dim> (b,      (a+b)/2,b),
				       // back face of the big cube
				       // last (top right) point is missing
      Point<dim> (a,      b,a),
      Point<dim> ((a+b)/2,b,a),
      Point<dim> (b,      b,a),
      Point<dim> (a,      b,(a+b)/2),
      Point<dim> ((a+b)/2,b,(a+b)/2),
      Point<dim> (b,      b,(a+b)/2),
      Point<dim> (a,      b,b),
      Point<dim> ((a+b)/2,b,b)
      };
  const int cell_vertices[7][8] = {{0, 1, 4, 3, 9, 10, 13, 12},
				   {1, 2, 5, 4, 10, 11, 14, 13},
				   {3, 4, 7, 6, 12, 13, 16, 15},
				   {4, 5, 8, 7, 13, 14, 17, 16},
				   {9, 10, 13, 12, 18, 19, 22, 21},
				   {10, 11, 14, 13, 19, 20, 23, 22},
				   {12, 13, 16, 15, 21, 22, 25, 24}};

  std::vector<CellData<3> > cells (7, CellData<3>());
  
  for (unsigned int i=0; i<7; ++i) 
    {
      for (unsigned int j=0; j<8; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (std::vector<Point<dim> >(&vertices[0], &vertices[26]),
			     cells,
			     SubCellData());       // no boundary information
};



template <>
void GridGenerator::hyper_ball<> (Triangulation<3> &tria,
				  const Point<3> &p,
				  const double radius)
{
				   // this function used to be
				   // implemented by the code below,
				   // but it turned out that it didn't
				   // work as expected: there were
				   // faces that we used in both
				   // orientations, leading the
				   // triangulation function to
				   // believe that the faces were
				   // external, even though they were
				   // in fact internal to the
				   // ball. this leads to strange
				   // results.
				   //
				   // since we haven't found a working
				   // enumeration of cell vertices and
				   // faces, this function is disabled
				   // altogether
  Assert(false, ExcNotImplemented());
  
  const double a = 1./(1+sqrt(3)); // equilibrate cell sizes at transition
				   // from the inner part to the radial
				   // cells
  const unsigned int n_vertices = 16;
  const Point<3> vertices[n_vertices]
    = {
				     // first the vertices of the inner
				     // cell
      p+Point<3>(-1,-1,-1)*(radius/sqrt(3)*a),
      p+Point<3>(+1,-1,-1)*(radius/sqrt(3)*a),
      p+Point<3>(+1,+1,-1)*(radius/sqrt(3)*a),
      p+Point<3>(-1,+1,-1)*(radius/sqrt(3)*a),
      p+Point<3>(-1,-1,+1)*(radius/sqrt(3)*a),
      p+Point<3>(+1,-1,+1)*(radius/sqrt(3)*a),
      p+Point<3>(+1,+1,+1)*(radius/sqrt(3)*a),
      p+Point<3>(-1,+1,+1)*(radius/sqrt(3)*a),
				       // now the eight vertices at
				       // the outer sphere
      p+Point<3>(-1,-1,-1)*(radius/sqrt(3)),
      p+Point<3>(+1,-1,-1)*(radius/sqrt(3)),
      p+Point<3>(+1,+1,-1)*(radius/sqrt(3)),
      p+Point<3>(-1,+1,-1)*(radius/sqrt(3)),
      p+Point<3>(-1,-1,+1)*(radius/sqrt(3)),
      p+Point<3>(+1,-1,+1)*(radius/sqrt(3)),
      p+Point<3>(+1,+1,+1)*(radius/sqrt(3)),
      p+Point<3>(-1,+1,+1)*(radius/sqrt(3))
      };

				   // one needs to draw the seven cubes to
				   // understand what's going on here
  const unsigned int n_cells = 7;
  const int cell_vertices[n_cells][8] = {{0, 1, 2, 3, 4, 5, 6, 7},
					 {8, 9, 10, 11, 0, 1, 2, 3},
					 {9, 13, 14, 10, 1, 5, 6, 2},
					 {12, 4, 7, 15, 13, 5, 6, 14},
					 {8, 0, 3, 11, 12, 4, 7, 15},
					 {11, 10,14, 15, 3, 2, 6, 7},
					 {8, 9, 13, 12, 0, 1, 5, 4}};
  
  std::vector<CellData<3> > cells (n_cells, CellData<3>());
  
  for (unsigned int i=0; i<n_cells; ++i) 
    {
      for (unsigned int j=0; j<GeometryInfo<3>::vertices_per_cell; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (std::vector<Point<3> >(&vertices[0], &vertices[n_vertices]),
			     cells,
			     SubCellData());       // no boundary information
};


template <>
void GridGenerator::hyper_shell<> (Triangulation<3>   &,
				   const Point<3>     &,
				   const double        ,
				   const double        ,
				   const unsigned int  )
{
  Assert (false, ExcNotImplemented());
};


#endif


// explicit instantiations
template void
GridGenerator::hyper_rectangle (Triangulation<deal_II_dimension> &,
				const Point<deal_II_dimension>&,
				const Point<deal_II_dimension>&);
template void
GridGenerator::hyper_cube (Triangulation<deal_II_dimension> &,
			   const double,
			   const double);

// template void GridGenerator::hyper_L (Triangulation<deal_II_dimension> &,
// 				      const double,
// 				      const double);
// template void GridGenerator::hyper_ball (Triangulation<deal_II_dimension> &,
// 					 const Point<deal_II_dimension> &,
// 					 const double);

