// 	$Id$

// test for correctness of gradients on a given cell. originally, this
// file did the testing itself, by comparing the computed results with
// builtin values; this would be unnecessary now, by simply comparing
// against a prebuilt output file, but it works this way also. if the
// computations fail at some time, the "ok" will be replaced by
// "WRONG!", which would let the comparison against the stored .expect
// file fail also.

// deal_II_libraries.g=-ldeal_II_2d.g
// deal_II_libraries=-ldeal_II_2d


#include <grid/tria.h>
#include <grid/tria_boundary.h>
#include <grid/dof.h>
#include <grid/grid_generator.h>
#include <fe/fe_values.h>
#include <fe/fe_lib.lagrange.h>
#include <base/quadrature_lib.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <lac/vector.h>




int main () {
  Triangulation<2> tria;
  GridGenerator::hyper_cube (tria,0,1);
  tria.begin_active()->vertex(2)(0) = 2;

  FEQ1<2> fe;
  DoFHandler<2> dof(&tria);
  dof.distribute_dofs(fe);

  QTrapez<2> q;
  FEValues<2> fevalues(fe,q,update_gradients);
  fevalues.reinit (dof.begin_active());
  
  
  Vector<double> val(4);

  cout << "------------------------------------------------------" << endl
       << "Testing transformation of gradients of shape function:" << endl;
  
				   // test for each of the four
				   // shape functions
  bool testcase_succeeded = true;
  for (unsigned int vertex=0; vertex<4; ++vertex)
    {
      val.clear ();
      val(vertex) = 1;

      vector<Tensor<1,2> > grads(4);
      fevalues.get_function_grads (val, grads);

      
      bool ok;
      switch (vertex) 
	{
	  case 0:
		ok = ((grads[0] == Point<2>(-1,-1)) &&
		      (grads[1] == Point<2>(0,-1)) &&
		      (grads[2] == Point<2>(-1,1)) &&
		      (grads[3] == Point<2>(0,0)));
		break;
	  case 1:
		ok = ((grads[0] == Point<2>(1,0)) &&
		      (grads[1] == Point<2>(0,0)) &&
		      (grads[2] == Point<2>(1,-2)) &&
		      (grads[3] == Point<2>(0,-1)));
		break;
	  case 2:
		ok = ((grads[0] == Point<2>(0,0)) &&
		      (grads[1] == Point<2>(0.5,0)) &&
		      (grads[2] == Point<2>(0,1)) &&
		      (grads[3] == Point<2>(0.5,0.5)));
		break;
	  case 3:
		ok = ((grads[0] == Point<2>(0,1)) &&
		      (grads[1] == Point<2>(-0.5,1)) &&
		      (grads[2] == Point<2>(0,0)) &&
		      (grads[3] == Point<2>(-0.5,0.5)));
		break;
	};

      cout << "  Shape function " << vertex
	   << ": "
	   << (ok ? "OK" : "WRONG!")
	   << endl;

      if (!ok)
	testcase_succeeded = false;
    };

  cout << "------------------------------------------------------" << endl;

  if (testcase_succeeded)
    return 0;
  else
    return 1;
};
