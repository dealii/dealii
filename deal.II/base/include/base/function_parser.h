//----------------------- function_parser.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------- function_parser.h  ---------------------------


#ifndef __deal2__function_parser_h
#define __deal2__function_parser_h

#include <base/config.h>
#include <base/exceptions.h>
#include <base/function.h>
#include <vector>
#include <map>
#include <functionparser/fparser.h>

template <int dim> class Point;
template <int rank_, int dim> class Tensor;
template <int dim> class Tensor<1,dim>;
template <typename number> class Vector;

/**
 * Function Parser. Wrapper class for the fparser library (see
 * http://www.students.tut.fi/~warp/FunctionParser/). This class lets
 * you evaluate strings such as "sqrt(1-x^2+y^2)" with given values of
 * 'x' and 'y'. Some of the informations contained here are copied
 * verbatim from the fparser.txt file that comes with the fparser
 * library. Please refer also to that file both for clarifications on
 * how this wrapper works, as well as for any issue regarding the
 * licence that applies to this class.
 *
 * By using this class you indicate that you accept the terms of the
 * licence that comes with the fparser library.
 *
 * The following example shows how to use this class:
 @verbatim

 
 // Define some constants that will be used by the function parser
 std::map<std::string> constants;
 constants["pi"] = M_PI;
 
 // Define the variables that will be used inside the expressions
 std::string variables = "x,y,z";
 
 // Define the expressions of the vector_valued function.
 std::vector<std::string> expressions(2);
 expressions[0] = "sin(2*pi*x)+sinh(pi*z)";
 expressions[1] = "sin(2*pi*y)*exp(x^2)";
 
 // Generate an empty function with two components
 ParsedFunction<3> vector_function(2);
 
 // And populate it with the newly created objects.
 vector_function.initialize(variables,
			    expressions,
			    constants);
 
 @endverbatim
 * 
 * See http://www.students.tut.fi/~warp/FunctionParser/ for an
 * explanation on how the underlining library works. 
 *
 
 From the fparser.txt file:
 @verbatim

 The library is intended to be very fast. It byte-compiles the
 function string at parse time and interpretes this byte-code at
 evaluation time.  The evaluation is straightforward and no recursions
 are done (uses stack arithmetic).  Empirical tests show that it indeed
 is very fast (specially compared to libraries which evaluate functions
 by just interpreting the raw function string).

 The library is made in ISO C++ and requires a standard-conforming C++
 compiler.

 @endverbatim
 * This object overloads for you the virtual methods value() and
 * vector_value() of the Function base class with the byte compiled
 * versions of the expressions given to the initialize() methods.
 *
 * @author Luca Heltai, 2005
 */
template <int dim>
class FunctionParser : public Function<dim>
{
  public:
  /** Constructor for Parsed functions. Its arguments are the same of
      the base class Function. The only difference is that this object
      needs to be initialized with initialize() method before you can
      use it. If an attempt to use this function is made before the
      initialize() method has been called, then an exception is
      thrown. */
  FunctionParser (const unsigned int n_components = 1,
		  const double       initial_time = 0.0); 
  
  /** Destructor. Explicitly delete the FunctionParser objects (there
      is one for each component of the function). */
  ~FunctionParser();
  
  /** Type for the constant map. Used by the initialize() method. */
  typedef std::map<std::string, double> ConstMap;
  
  /** Iterator for the constants map. Used by the initialize() method. */
  typedef ConstMap::iterator ConstMapIterator;

  /** Initialize the function. From the fparser.txt file:
      @verbatim
      The function string understood by the class is very similar to the C-syntax.
      Arithmetic float expressions can be created from float literals, variables
      or functions using the following operators in this order of precedence:

      ()             expressions in parentheses first
      -A             unary minus
      A^B            exponentiation (A raised to the power B)
      A*B  A/B  A%B  multiplication, division and modulo
      A+B  A-B       addition and subtraction
      A=B  A<B  A>B  comparison between A and B (result is either 0 or 1)
      A&B            result is 1 if int(A) and int(B) differ from 0, else 0.
      A|B            result is 1 if int(A) or int(B) differ from 0, else 0.

      Since the unary minus has higher precedence than any other operator, for
      example the following expression is valid: x*-y
      Note that the '=' comparison can be inaccurate due to floating point
      precision problems (eg. "sqrt(100)=10" probably returns 0, not 1).
      
      The class supports these functions:
      
      abs(A)    : Absolute value of A. If A is negative, returns -A otherwise
		  returns A.
      acos(A)   : Arc-cosine of A. Returns the angle, measured in radians,
		  whose cosine is A.
      acosh(A)  : Same as acos() but for hyperbolic cosine.
      asin(A)   : Arc-sine of A. Returns the angle, measured in radians, whose
		  sine is A.
      asinh(A)  : Same as asin() but for hyperbolic sine.
      atan(A)   : Arc-tangent of (A). Returns the angle, measured in radians,
                  whose tangent is (A).
      atan2(A,B): Arc-tangent of A/B. The two main differences to atan() is
                  that it will return the right angle depending on the signs of
                  A and B (atan() can only return values betwen -pi/2 and pi/2),
                  and that the return value of pi/2 and -pi/2 are possible.
      atanh(A)  : Same as atan() but for hyperbolic tangent.
      ceil(A)   : Ceiling of A. Returns the smallest integer greater than A.
                  Rounds up to the next higher integer.
      cos(A)    : Cosine of A. Returns the cosine of the angle A, where A is
                  measured in radians.
      cosh(A)   : Same as cos() but for hyperbolic cosine.
      cot(A)    : Cotangent of A (equivalent to 1/tan(A)).
      csc(A)    : Cosecant of A (equivalent to 1/sin(A)).
      eval(...) : This a recursive call to the function to be evaluated. The
                  number of parameters must be the same as the number of parameters
		  taken by the function. Usually called inside if() to avoid
		  infinite recursion.
      exp(A)    : Exponential of A. Returns the value of e raised to the power
                  A where e is the base of the natural logarithm, i.e. the
                  non-repeating value approximately equal to 2.71828182846.
      floor(A)  : Floor of A. Returns the largest integer less than A. Rounds
                  down to the next lower integer.
      if(A,B,C) : If int(A) differs from 0, the return value of this function is B,
                  else C. Only the parameter which needs to be evaluated is
                  evaluated, the other parameter is skipped; this makes it safe to
                  use eval() in them.
      int(A)    : Rounds A to the closest integer. 0.5 is rounded to 1.
      log(A)    : Natural (base e) logarithm of A.
      log10(A)  : Base 10 logarithm of A.
      max(A,B)  : If A>B, the result is A, else B.
      min(A,B)  : If A<B, the result is A, else B.
      sec(A)    : Secant of A (equivalent to 1/cos(A)).
      sin(A)    : Sine of A. Returns the sine of the angle A, where A is
                  measured in radians.
      sinh(A)   : Same as sin() but for hyperbolic sine.
      sqrt(A)   : Square root of A. Returns the value whose square is A.
      tan(A)    : Tangent of A. Returns the tangent of the angle A, where A
                  is measured in radians.
      tanh(A)   : Same as tan() but for hyperbolic tangent.


    Examples of function string understood by the class:

    "1+2"
    "x-1"
    "-sin(sqrt(x^2+y^2))"
    "sqrt(XCoord*XCoord + YCoord*YCoord)"
    
    An example of a recursive function is the factorial function:
    
    "if(n>1, n*eval(n-1), 1)"

  Note that a recursive call has some overhead, which makes it a bit slower
  than any other operation. It may be a good idea to avoid recursive functions
  in very time-critical applications. Recursion also takes some memory, so
  extremely deep recursions should be avoided (eg. millions of nested recursive
  calls).

  Also note that the if() function is the only place where making a recursive
  call is safe. In any other place it will cause an infinite recursion (which
  will make the program eventually run out of memory).
  @endverbatim
  
  This methods accepts the following parameters: 
  
  - <b>vars</b>. It contains a string with the variables that will be
    used by the expressions to be evaluated. Note that the variables
    can have any name (of course different from the function names
    defined above!), but the order IS important. The first variable
    will correspond to the first component of the point in which the
    function is evaluated, the second variable to the second component
    and so forth. If this function is also time dependent, then it is
    necessary to specify it by setting the #time_dependent parameter
    to true.     
    An Exception is thrown if the number of variables specified here
    is different from dim (if this function is not time-dependent) or
    from dim+1 (if it is time-dependent).
    
  - <b>expressions</b>. This is a list of strings containing the
    expressions that will be byte compiled by the internal parser
    (FunctionParser). Note that the size of this vector must match
    exactly the number of components of the FunctionParser. If this is
    not the case, an exception is thrown.
    
  - <b>constants</b>. The map of constants is used to pass any
    necessary constant that we want to specify in our expressions (in
    the example above the number pi). An expression is valid if and
    only if it contains only defined variables and defined constants
    (other than the functions specified above). If a constant is given
    whose name is not valid (eg: <tt>constants["sin"] = 1.5;</tt>) an
    exception is thrown.
    
  - <b>time_dependent</b>. If this is a time dependent function, then
    the last variable is assumed to be the time variable, and
    this->get_time() is passed to the function. Naturally the number
    of variables parsed by the initialize() method in this case is
    dim+1. It defaults to false, i.e.. do not consider time.

  - <b>use_degrees</b>. Parameter to decide if the trigonometric
    functions work in radians or degrees. The default for this
    parameter is false, i.e. use radians and not degrees.
    
    
  */
  void initialize(const std::string vars,
		  const std::vector<std::string> expressions,
		  const ConstMap constants = 0,
		  bool time_dependent = false,
		  bool use_degrees = false);
  
  /** Initialize the function. Same as above, but for scalar
      functions. An exception is thrown if this method is called and
      the function has more than one component. 
      
      An example of time dependent scalar function is the following:
      @verbatim
     
      // Empty constants object
      std::map<std::string> constants;
      
      // Variables that will be used inside the expressions
      std::string variables = "x,y,t";
      
      // Define the expression of the scalar time dependent function.
      std::string expression = "exp(y*x)*exp(-t)";
      
      // Generate an empty scalar function
      FunctionParser<2> function;
      
      // And populate it with the newly created objects.
      function.initialize(variables,
			  expression,
			  constants,
			  true);	// This tells the parser that 
				 	// it is a time-dependent function
					// and there is another variable
					// to be taken into account (t).
 
     @endverbatim
      
  */
  void initialize(const std::string vars,
		  const std::string expression,
		  const ConstMap constants = 0,
		  bool time_dependent = false,
		  bool use_degrees = false);
  

  /**
   * Return the value of the function at the given point. Unless there
   * is only one component (i.e. the function is scalar), you should
   * state the component you want to have evaluated; it defaults to
   * zero, i.e. the first component. */
  virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const;

  /**
   * Return all components of a vector-valued function at a given
   * point.
   *
   * <tt>values</tt> shall have the
   * right size beforehand,
   * i.e. #n_components.
   */
  virtual void vector_value (const Point<dim>   &p,
			     Vector<double>     &values) const;

  /** @addtogroup Exceptions
   * @{ */    
  DeclException2 (ExcParseError,
		  int, char*, 
		  << "Parsing Error at Column " << arg1
		  << ". The parser said: " << arg2);
    
  DeclException2 (ExcInvalidExpressionSize,
		  int, int, 
		  << "The number of components (" << arg1
		  << ") is not equal to the number of expressions (" 
		  << arg2 << ").");
    
  //@}
 private: 
  /** A pointer to the actual function parsers. */
 fparser::FunctionParser * fp;
    
 /** State of usability. This variable is checked every time the
     function is called for evaluation. It's set to true in the
     initialize() methods.*/
 bool initialized;
    
 /** Number of variables. If this is also a function of time, then the
     number of variables is dim+1, otherwhise it is dim. In the case
     that this is a time dependent function, the time is supposed to
     be the last variable. If #n_vars is not identical to the number
     of the variables parsed by the initialize() method, then an
     exception is thrown.   
 */
 unsigned int n_vars;
};


#endif

  
