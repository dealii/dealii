#include <base/logstream.h>
#include <lac/sparse_matrix_ez.h>

#include <fstream>


int main ()
{
  std::ofstream logfile("sparse_matrix_ez_iterator.output");
  logfile.setf(std::ios::fixed);
  logfile.precision(2);
  deallog.attach(logfile);

  SparseMatrixEZ<double> m (4,4,4);
  m.set (0,0,1);
  m.set (1,0,2);
  m.set (2,0,3);
  m.set (2,2,0);  // should be ignored
  m.set (0,1,1);
  m.set (0,2,2);
  m.set (0,3,3);
  
  for (SparseMatrixEZ<double>::const_iterator i = m.begin(); i != m.end(); ++i)
    deallog << i->row() << ' ' << i->column() << ' ' << i->value()
	    << std::endl;
}

