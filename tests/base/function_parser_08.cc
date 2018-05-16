#include "../tests.h"
#include <map>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/function_parser.h>

// simple function parser test for erfc functionality.

void
test1()
{
  // set up problem:
  std::string variables = "x,y";
  std::string expression = "erfc(x)+erfc(y)";
  std::map<std::string,double> constants;

  // FunctionParser with 2 variables and 1 component:
  FunctionParser<2> fp(1);
  fp.initialize(variables,
                expression,
                constants);

  // Point at which we want to evaluate the function
  Point<2> point(0.1, -0.7);

  // evaluate the expression at 'point':
  double result = fp.value(point);

  deallog << "Function '" << expression << "'"
          << " @ " << point
          << " is " << result << std::endl;
}


void
test2()
{
  // Define some constants that will be used by the function parser
  std::map<std::string,double> constants;
  constants["pi"] = numbers::PI;

  // Define the variables that will be used inside the expressions
  std::string variables = "x,y,z";

  // Define the expressions of the individual components of a
  // vector valued function with two components:
  std::vector<std::string> expressions(2);
  expressions[0] = "erfc(pi*z)";
  expressions[1] = "sin(2*pi*y)*erfc(x^2)";

  // function parser with 3 variables and 2 components
  FunctionParser<3> vector_function(2);

  // And populate it with the newly created objects.
  vector_function.initialize(variables,
                             expressions,
                             constants);

  // Point at which we want to evaluate the function
  Point<3> point(0.0, 1.0, 1.0);

  // This Vector will store the result
  Vector<double> result(2);

  // Fill 'result' by evaluating the function
  vector_function.vector_value(point, result);

  // We can also only evaluate the 2nd component:
  double c = vector_function.value(point, 1);

  Assert(c == result[1], ExcInternalError());

  // Output the evaluated function
  deallog << "Function '" << expressions[0] << "," << expressions[1] << "'"
          << " @ " << point
          << " is " << result << std::endl;
}


int
main()
{
  initlog();

  test1();
  test2();
}
