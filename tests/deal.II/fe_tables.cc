// $Id$

// deal_II_libraries.g=-ldeal_II_2d.g -ldeal_II_1d.g
// deal_II_libraries=-ldeal_II_2d

#include <base/exceptions.h>
#include <base/point.h>
#include <lac/vector.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_lib.dg.h>
#include <fe/fe_system.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/grid_generator.h>
#include <iomanip>
#include <fstream>

#include <base/logstream.h>

#define TEST_ELEMENT(e) { deallog.push(#e); e el;\
  print_fe_statistics(el); deallog.pop(); deallog << endl; }
#define TEST_MULTIPLE(e,n,d) { deallog.push(#e "x" #n); e eb; FESystem<d> el(eb,n); \
  print_fe_statistics(el); deallog.pop(); deallog << endl; }
#define TEST_MIXED2(e1,n1,e2,n2,d) { deallog.push(#e1 "x" #n1 "-" #e2 "x" #n2);\
  e1 eb1; e2 eb2; FESystem<d> el(eb1,n1,eb2,n2);\
  print_fe_statistics(el); deallog.pop(); deallog << endl; }


template<int dim>
inline void
print_fe_statistics(const FiniteElement<dim>& fe)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr,-1,1);
  DoFHandler<dim> dof(&tr);
  dof.distribute_dofs(fe);
  StraightBoundary<dim> boundary;
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  DoFHandler<dim>::active_face_iterator face = dof.begin_active_face();

  vector<Point<dim> > unit_points(fe.dofs_per_cell);
  vector<Point<dim> > support_points(fe.dofs_per_cell);
  vector<Point<dim> > face_support_points(fe.dofs_per_face);

  fe.get_unit_support_points(unit_points);
  fe.get_support_points(cell, support_points);
  fe.get_face_support_points(face, face_support_points);
  
  deallog << "dofs_per_cell" << " " << fe.dofs_per_cell;
  deallog << ": vertex" << " " << fe.dofs_per_vertex;
  deallog << "  line" << " " << fe.dofs_per_line;
  deallog << "  quad" <<  " " <<fe.dofs_per_quad << endl;
  deallog << "n_transform_fct " << fe.n_transform_functions() << endl;
  deallog << "n_components " << fe.n_components() << endl;
  deallog.push("components");
  for (unsigned i=0;i<fe.dofs_per_cell;++i)
    {
      pair<unsigned,unsigned> p = fe.system_to_component_index(i);
      deallog << "Index " << i << " ("
	      << p.first << "," << p.second << ") -> "
	      << fe.component_to_system_index(p.first, p.second)
	      << " support " << support_points[i] << " unit: " << unit_points[i]
	      << endl;
    }
  for (unsigned i=0;i<fe.dofs_per_face;++i)
    {
      pair<unsigned,unsigned> p = fe.face_system_to_component_index(i);
      deallog << "FaceIndex " << i << " ("
	      << p.first << "," << p.second << ") -> " 
	      << fe.face_component_to_system_index(p.first, p.second)
	      << " support " << face_support_points[i]
	      << endl;
    }
  deallog.pop();
}

int main()
{
  ofstream logfile("fe_tables.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  logfile.precision(4);
  
  deallog.push("GeometryInfo");

  deallog.push("1D");
  deallog << " vertices: " << GeometryInfo<1>::vertices_per_cell
	  << " lines: " << GeometryInfo<1>::lines_per_cell
	  << " quads: " << GeometryInfo<1>::quads_per_cell
	  << " hexes: " << GeometryInfo<1>::hexes_per_cell
	  << endl;
  deallog.pop();
  
  deallog.push("2D");
  deallog << " vertices: " << GeometryInfo<2>::vertices_per_cell
	  << " lines: " << GeometryInfo<2>::lines_per_cell
	  << " quads: " << GeometryInfo<2>::quads_per_cell
	  << " hexes: " << GeometryInfo<2>::hexes_per_cell
	  << endl;
  deallog.pop();
  
  deallog.push("3D");
  deallog << " vertices: " << GeometryInfo<3>::vertices_per_cell
	  << " lines: " << GeometryInfo<3>::lines_per_cell
	  << " quads: " << GeometryInfo<3>::quads_per_cell
	  << " hexes: " << GeometryInfo<3>::hexes_per_cell
	  << endl;
  deallog.pop();
  
  deallog.push("4D");
  deallog << " vertices: " << GeometryInfo<4>::vertices_per_cell
	  << " lines: " << GeometryInfo<4>::lines_per_cell
	  << " quads: " << GeometryInfo<4>::quads_per_cell
	  << " hexes: " << GeometryInfo<4>::hexes_per_cell
	  << endl;
  deallog.pop();
  
  deallog.pop();
  
  TEST_ELEMENT(FEDG_Q0<2>);
  TEST_ELEMENT(FEDG_Q1<2>);
  
  TEST_ELEMENT(FEQ1<2>);
  TEST_ELEMENT(FEQ2<2>);
  TEST_ELEMENT(FEQ3<2>);
  TEST_ELEMENT(FEQ4<2>);

  TEST_MULTIPLE(FEQ1<2>,3,2);
  TEST_MULTIPLE(FEQ2<2>,3,2);
  TEST_MULTIPLE(FEQ3<2>,3,2);
  
  TEST_MIXED2(FEQ1<2>,1,FEDG_Q0<2>,1,2);
  TEST_MIXED2(FEQ2<2>,3,FEQ1<2>,1,2);
  TEST_MIXED2(FEQ3<2>,3,FEQ2<2>,2,2);
}
