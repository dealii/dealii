// $Id$
// (c) Guido Kanschat
//
// Compute support points

#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_lib.dg.h>
#include <fe/fe_values.h>
#include <vector>
#include <fstream>
#include <string>

template <int dim>
inline void
check_support (FiniteElement<dim>& finel, const char* name)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);
  DoFHandler<dim> dof (tr);
  dof.distribute_dofs (finel);

  vector<Point<dim> > cell_points (finel.dofs_per_cell);
  vector<Point<dim> > face_points (finel.dofs_per_face);
  
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  finel.get_support_points (cell, cell_points);

  cout << name << '<' << dim << '>' << " cell support points" << endl;
  
  for (unsigned int k=0;k<cell_points.size();++k)
    cout << cell_points[k] << endl;
  
  for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
    {
      cout << name << '<' << dim << '>' << " face " << i << " support points" << endl;
      DoFHandler<dim>::active_face_iterator face = cell->face(i);
      finel.get_face_support_points (face, face_points);
      
      for (unsigned int k=0;k<face_points.size();++k)
	cout << face_points[k] << endl;
    }
}

template <int dim>
inline void
check_matrices (FiniteElement<dim>& fe, const char* name)
{
  for (unsigned int i=0;i<GeometryInfo<dim>::children_per_cell;++i)
    {
      cout << name << '<' << dim << '>' << " restriction " << i << endl;
      fe.restrict(i).print_formatted (cout, 3, false, 6, "~");
      cout << name << '<' << dim << '>' << " embedding " << i << endl;
      fe.prolongate(i).print_formatted (cout, 3, false, 6, "~");
    }
}


#define CHECK_S(EL,dim) FE ## EL<dim> EL; check_support(EL, #EL);
#define CHECK_M(EL,dim) FE ## EL<dim> EL; check_matrices(EL, #EL);
#define CHECK_ALL(EL,dim) FE ## EL<dim> EL; check_support(EL, #EL); check_matrices(EL,#EL)

int
main()
{
  if (true)
    {
      CHECK_ALL(Q1,2);
      CHECK_ALL(Q2,2);
      CHECK_ALL(Q3,2);
      CHECK_ALL(Q4,2);
      CHECK_ALL(DG_Q0,2);
      CHECK_ALL(DG_Q1,2);
      CHECK_ALL(DG_Q2,2);
      CHECK_ALL(DG_Q3,2);
    }
  if (true)
    {
      CHECK_ALL(Q1,3);
      CHECK_ALL(Q2,3);
      CHECK_ALL(DG_Q0,3);
      CHECK_ALL(DG_Q1,3);
    }
  
  return 0;
}
