/* $Id$ */


#include <grid/grid_generator.h>
#include <grid/tria.h>
#include <cmath>





#if deal_II_dimension == 1

template <>
void GridGenerator::hyper_cube<> (Triangulation<1> &tria,
				  const double left,
				  const double right) {
  const Point<1> vertices[2] = { Point<1>(left), Point<1>(right) };
  vector<CellData<1> > cells (1, CellData<1>());
  cells[0].vertices[0] = 0;
  cells[0].vertices[1] = 1;
  cells[0].material_id = 0;
  
  tria.create_triangulation (vector<Point<1> >(&vertices[0], &vertices[2]),
			     cells,
			     SubCellData());       // no boundary information
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

#endif



#if deal_II_dimension == 2

template <>
void GridGenerator::hyper_cube<> (Triangulation<2> &tria,
				  const double left,
				  const double right) {
  const Point<2> vertices[4] = { Point<2>(left,left),
				 Point<2>(right,left),
				 Point<2>(right,right),
				 Point<2>(left,right)  };
  const int cell_vertices[1][4] = { { 0,1,2,3 } };
  vector<CellData<2> > cells (1, CellData<2>());
  for (unsigned int j=0; j<4; ++j)
    cells[0].vertices[j] = cell_vertices[0][j];
  cells[0].material_id = 0;
  
  tria.create_triangulation (vector<Point<2> >(&vertices[0], &vertices[4]),
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

  vector<CellData<2> > cells (3, CellData<2>());
  
  for (unsigned int i=0; i<3; ++i) 
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (vector<Point<dim> >(&vertices[0], &vertices[8]),
			     cells,
			     SubCellData());       // no boundary information
};



template <>
void GridGenerator::hyper_ball<> (Triangulation<2> &tria,
				  const Point<2>    &p,
				  const double      radius) {
  const double a = 1./(1+sqrt(2)); // equilibrate cell sizes at transition
				   // from the inner part to the radial
				   // cells
  const Point<2> vertices[8] = { p+Point<2>(-1,-1)*(radius/sqrt(2)),
				 p+Point<2>(+1,-1)*(radius/sqrt(2)),
				 p+Point<2>(-1,-1)*(radius/sqrt(2)*a),
				 p+Point<2>(+1,-1)*(radius/sqrt(2)*a),
				 p+Point<2>(-1,+1)*(radius/sqrt(2)*a),
				 p+Point<2>(+1,+1)*(radius/sqrt(2)*a),
				 p+Point<2>(-1,+1)*(radius/sqrt(2)),
				 p+Point<2>(+1,+1)*(radius/sqrt(2)) };
  
  const int cell_vertices[5][4] = {{0, 1, 3, 2},
				   {0, 2, 4, 6},
				   {2, 3, 5, 4},
				   {1, 7, 5, 3},
				   {6, 4, 5, 7}};

  vector<CellData<2> > cells (5, CellData<2>());
  
  for (unsigned int i=0; i<5; ++i) 
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (vector<Point<2> >(&vertices[0], &vertices[8]),
			     cells,
			     SubCellData());       // no boundary information
};

#endif


#if deal_II_dimension == 3

template <>
void GridGenerator::hyper_cube<> (Triangulation<3> &tria,
				  const double left,
				  const double right) {
  const Point<3> vertices[8] = { Point<3>(left,left,left),
				 Point<3>(right,left,left),
				 Point<3>(right,left,right),
				 Point<3>(left,left,right),

				 Point<3>(left,right,left),				 
				 Point<3>(right,right,left),
				 Point<3>(right,right,right),
				 Point<3>(left,right,right)};
  const int cell_vertices[1][8] = { { 0,1,2,3,4,5,6,7 } };
  vector<CellData<3> > cells (1, CellData<3>());
  for (unsigned int j=0; j<8; ++j)
    cells[0].vertices[j] = cell_vertices[0][j];
  cells[0].material_id = 0;
  
  tria.create_triangulation (vector<Point<3> >(&vertices[0], &vertices[8]),
			cells,
			SubCellData());       // no boundary information
};



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

  vector<CellData<3> > cells (7, CellData<3>());
  
  for (unsigned int i=0; i<7; ++i) 
    {
      for (unsigned int j=0; j<8; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (vector<Point<dim> >(&vertices[0], &vertices[26]),
			     cells,
			     SubCellData());       // no boundary information
};



template <>
void GridGenerator::hyper_ball<> (Triangulation<3> &tria,
				  const Point<3> &p,
				  const double radius) {
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

  vector<CellData<3> > cells (n_cells, CellData<3>());
  
  for (unsigned int i=0; i<n_cells; ++i) 
    {
      for (unsigned int j=0; j<GeometryInfo<3>::vertices_per_cell; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (vector<Point<3> >(&vertices[0], &vertices[n_vertices]),
			     cells,
			     SubCellData());       // no boundary information
};



#endif




// explicit instantiations
// template void GridGenerator::hyper_cube (Triangulation<deal_II_dimension> &,
// 					 const double,
// 					 const double);
// template void GridGenerator::hyper_L (Triangulation<deal_II_dimension> &,
// 				      const double,
// 				      const double);
// template void GridGenerator::hyper_ball (Triangulation<deal_II_dimension> &,
// 					 const Point<deal_II_dimension> &,
// 					 const double);

