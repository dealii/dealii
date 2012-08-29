/**********************************************************************
 * $Id: constraint_test.cc -1   $
 *
 * Copyright Guido Kanschat, 2010, 2012
 *
 **********************************************************************/

// reduced test case from constraint_c1.cc, causes a hang in close()


#include <deal.II/base/job_identifier.h>
#include <deal.II/lac/constraint_matrix.h>

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace dealii;



void
run()
{
  ConstraintMatrix     constraints;

  {
  constraints.add_line(7);
  std::vector<std::pair<unsigned int, double> > rhs;
  rhs.push_back(std::pair<unsigned int, double>(41,1.0));
  rhs.push_back(std::pair<unsigned int, double>(42,1.0));
  constraints.add_entries(7, rhs);
  }

  {
  constraints.add_line(41);
  std::vector<std::pair<unsigned int, double> > rhs;
  rhs.push_back(std::pair<unsigned int, double>(42,1.0));
  constraints.add_entries(41, rhs);
  }

  {
  constraints.add_line(42);
  std::vector<std::pair<unsigned int, double> > rhs;
  rhs.push_back(std::pair<unsigned int, double>(41,1.0));
  constraints.add_entries(42, rhs);
  }

  deallog << "Closing" << std::endl;
  constraints.close();
  deallog << "Closed" << std::endl;

}


int main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  
  run();
}
