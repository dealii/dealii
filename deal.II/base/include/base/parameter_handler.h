/**********************   parameter-handler.h     ****************************/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __parameter_handler_H
#define __parameter_handler_H
/**********************   parameter-handler.h     ****************************/


// public classes; to be declared below
class ParameterHandler;
class MultipleParameterLoop;


#include <map>
#include <vector>
#include <string>
#include <base/exceptions.h>

class istream;
class ostream;



/**
 * List of possible output formats.
 */
enum OutputStyle {
      Text, LaTeX, HTML
};



/**
 * Declare some regexps which
 * may be used to define patterns.
 */ 
struct Patterns {
  public:
				     /**
				      * Base class to declare common
				      * interface.
				      */
    class PatternBase {
      public:
					 /**
					  * Make destructor of this and all
					  * derived classes virtual.
					  */
	virtual ~PatternBase () {};
	
					 /**
					  * Return true if the given string
					  * matches the pattern.
					  */
	virtual bool match (const string &test_string) const = 0;
	
					 /**
					  * Return a string describing the
					  * pattern.
					  */
	virtual string description () const = 0;
	
					 /**
					  * Return a pointer to an exact
					  * copy of the object. This is
					  * necessary since we want to store
					  * objects of this type in
					  * containers.
					  */
	virtual PatternBase * clone () const = 0;
    };
    
				     /**
				      * Test for the string being an
				      * integer.
				      */
    class Integer : public PatternBase {
      public:
	virtual bool match (const string &test_string) const;
	virtual string description () const;
	virtual PatternBase * clone () const;
    };
    
				     /**
				      * Test for the string being a double.
				      */
    class Double : public PatternBase {
      public:
	virtual bool match (const string &test_string) const;
	virtual string description () const;
	virtual PatternBase * clone () const;
    };
    
				     /**
				      * Test for the string being one of
				      * a sequence of values given like a
				      * regular expression. For example, if
				      * the string given to the constructor
				      * is "red|blue|black", then the #match#
				      * function returns #true# exactly if
				      * the string is either "red" or "blue"
				      * or "black". Spaces around the pipe
				      * signs do not matter and are
				      * eliminated.
				      */
    class Sequence : public PatternBase {
      public:
	Sequence (const string &seq);
	virtual bool match (const string &test_string) const;
	virtual string description () const;
	virtual PatternBase * clone () const;
      private:
	string sequence;
    };

				     /**
				      * Test for the string being either
				      * "true" or "false". This is mapped
				      * to the #Sequence# class.
				      */
    class Bool : public Sequence {
      public:
	Bool ();
	virtual PatternBase * clone () const;
    };
        
				     /**
				      * Always returns true when testing a
				      * string.
				      */
    class Anything : public PatternBase {
      public:
	virtual bool match (const string &test_string) const;
	virtual string description () const;
	virtual PatternBase * clone () const;
    };
};



      

/**
 *   The #ParameterHandler# class provides a standard interface to an input file
 *   which provides at run-time for program parameters such as time step sizes,
 *   geometries, right hand sides etc. The input for the program is given in files,
 *   streams or strings in memory using text like
 *   \begin{verbatim}
 *     set Time step size = 0.3
 *     set Geometry       = [0,1]x[0,3]
 *   \end{verbatim}
 *   Input may be sorted into subsection trees in order to give the input a logical
 *   structure.
 *
 *   
 *   \subsection{Declaration of entries}
 *   
 *   In order to use the facilities of a #ParameterHandler# object, one first has
 *   to make known the different entries the input file may or may not contain. This
 *   is done in the following way:
 *   \begin{verbatim}
 *     ...
 *     ParameterHandler prm;
 *     prm.declare_entry ("Time step size",
 *                       "0.2",
 *                       Patterns::Double());
 *     prm.declare_entry ("Geometry",
 *                       "[0,1]x[0,1]",
 *                       Patterns::Anything());
 *     ...
 *   \end{verbatim}
 *   Each entry is declared using the function #declare_entry#. The first parameter is
 *   the name of the entry (in short: the entry). The second is the default answer to
 *   be taken in case the entry is not specified in the input file. The third parameter
 *   is a regular expression which the input (and the default answer) has to match.
 *   Several such regular expressions are defined in #Patterns#.
 *
 *   Entries may be located in subsections which form a kind of input tree. For example
 *   input parameters for linear solver routines should be classified in a subsection
 *   named #Linear solver# or any other suitable name. This is accomplished in the
 *   following way:
 *   \begin{verbatim}
 *     ...
 *       LinEq eq;
 *       eq.declare_parameters (prm);
 *     ...
 *
 *     void LinEq::declare_parameters (ParameterHandler &prm) {
 *       prm.enter_subsection("Linear solver");
 *       prm.declare_entry ("Solver",
 *                          "CG",
 *		            Patterns::Sequence("CG|GMRES|GaussElim"));
 *       prm.declare_entry ("Maximum number of iterations",
 *                          "20",
 *		            ParameterHandler::RegularExpressions::Integer());
 *       ...
 *       prm.leave_subsection ();
 *     };
 *   \end{verbatim}
 *
 *   Subsections may be nested. For example a nonlinear solver may have a linear solver
 *   as member object. Then the function call tree would be something like (if the class
 *   #NonLinEq# has a member variables #eq# of type #LinEq#):
 *   \begin{verbatim}
 *     void NonLinEq::declare_parameters (ParameterHandler &prm) {
 *       prm.enter_subsection ("Nonlinear solver");
 *       prm.declare_entry ("Nonlinear method",
 *                          "Newton-Raphson",
 *		            ParameterHandler::RegularExpressions::Anything());
 *       eq.declare_parameters (prm);
 *       prm.leave_subsection ();
 *     };
 *   \end{verbatim}
 *
 *   For class member functions which declare the different entries we propose to use the
 *   common name #declare_parameters#. In normal cases this method can be #static# since the
 *   entries will not depend on any previous knowledge. Classes for which entries should
 *   logically be grouped into subsections should declare these subsections themselves. If
 *   a class has two or more member variables of the same type both of which should have
 *   their own parameters, this parent class' method #declare_parameters# is responsible to
 *   group them into different subsections:
 *   \begin{verbatim}
 *     void NonLinEq::declare_parameters (ParameterHandler &prm) {
 *       prm.enter_subsection ("Nonlinear solver");
 *       prm.enter_subsection ("Linear solver 1");
 *       eq1.declare_parameters (prm);
 *       prm.leave_subsection ();
 *
 *       prm.enter_subsection ("Linear solver 2");
 *       eq2.declare_parameters (prm);
 *       prm.leave_subsection ();
 *       prm.leave_subsection ();
 *     };	
 *   \end{verbatim}
 *
 *
 *   \subsection{Input files and special characters}
 *
 *   For the first example above the input file would look like the following:
 *   \begin{verbatim}
 *     ...
 *     subsection Nonlinear solver
 *       set Nonlinear method = Gradient
 *       subsection Linear solver
 *         set Solver                        = CG
 *         set Maxmimum number of iterations = 30
 *       end
 *     end
 *     ...                       # other stuff
 *   \end{verbatim}
 *   The words #subsection#, #set# and #end# may be either written in lowercase or uppercase
 *   letters. Leading and trailing whitespace is removed, multiple whitespace is condensed into
 *   only one. Since the latter applies also to the name of an entry, an entry name will not
 *   be recognised if in the declaration multiple whitespace is used.
 *
 *   In entry names and values the following characters are not allowed: \#, #{#, #}#, #|#.
 *   Their use is reserved for the \ref{MultipleParameterLoop} class.
 *   
 *   Comments starting with \# are skipped.
 *   
 *   We propose to use the following
 *   scheme to name entries: start the first word with a capital letter and use lowercase
 *   letters further on. The same applies to the possible entry values to the right of the
 *   #=# sign.
 *
 *   
 *   \subsection{Reading data from input sources}
 *   
 *   In order to read input you can use three possibilities: reading from an #istream# object,
 *   reading from a file of which the name is given and reading from a string in memory in
 *   which the lines are separated by #\n# characters. These possibilites are used as follows:
 *   \begin{verbatim}
 *     ParameterHandler prm;
 *     ...
 *     // declaration of entries
 *     ...
 *     prm.read_input (cin);         // read input from standard in,
 *     // or
 *     prm.read_input ("simulation.in");
 *     // or
 *     char *in = "set Time step size = 0.3 \n ...";
 *     prm.read_input (in);
 *     ...
 *   \end{verbatim}
 *   You can use several sources of input successively. Entries which are changed more than
 *   once will be overwritten everytime they are used. It is suggested to let the name of
 *   parameter input end in #.prm#.
 *
 *   You should not try to declare entries using #declare_entry# and #enter_subsection# with as
 *   yet unknown subsection names after using #read_input#. The results in this case are
 *   unspecified.
 *
 *   If an error occurs upon reading the input, error messages are written to #cerr#.
 *
 *   
 *   \subsection{Getting entry values out of a #ParameterHandler# object}
 *   
 *   Each class gets its data out of a #ParameterHandler# object by calling the #get (...)#
 *   member functions like this:
 *   \begin{verbatim}
 *      void NonLinEq::get_parameters (ParameterHandler &prm) {
 *       prm.enter_subsection ("Nonlinear solver");
 *       string method = prm.get ("Nonlinear method");
 *       eq.get_parameters (prm);
 *       prm.leave_subsection ();
 *     };
 *   \end{verbatim}
 *   #get()# returns the value of the given entry. If the entry was not specified in the input
 *   source(s), the default value is returned. You have to enter and leave subsections
 *   exactly as you did when declaring subsection. You may chose the order in which to
 *   transverse the subsection tree.
 *
 *   It is guaranteed that only entries matching the given regular expression are returned,
 *   i.e. an input entry value which does not match the regular expression is not stored.
 *
 *   You can use #get# to retrieve the parameter in text form, #get_integer# to get an integer
 *   or #get_double# to get a double. You can also use #get_bool#.
 *   It will cause an internal error if the string could not be 
 *   converted to an integer, double or a bool. This should, though, not
 *   happen if you correctly specified the regular expression for this entry; you should not
 *   try to get out an integer or a double from an entry for which no according regular
 *   expression was set. The internal error is raised through the #Assert()# macro family
 *   which only works in debug mode.
 *
 *   If you want to print out all user selectable features, use the
 *   #print_parameters# function. It is generally a good idea to print all parameters
 *   at the beginning of a log file, since this way input and output are together in
 *   one file which makes matching at a later time easier. Additionally, the function
 *   also print those entries which have not been modified in the input file und are
 *   thus set to default values; since default values may change in the process of
 *   program development, you cannot know the values of parameters not specified in the
 *   input file.
 *   
 *   
 *   \subsection{Style guide for data retrieval}
 *   
 *   We propose that every class which gets data out of a #ParameterHandler# object provides
 *   a function named #get_parameters#. This should be declared #virtual#. #get_parameters#
 *   functions in derived classes should call the #BaseClass::get_parameters# function.
 *
 *   
 *   \subsection{Possible future extensions}
 *   
 *   \begin{itemize}
 *   \item Allow long input lines to be broken by appending a backslash character
 *     (just like C macros and shell input).
 *   \item Provide an #input filename# command for the input file to enable users to put the
 *     most common parameters into separate files.
 *   \end{itemize}  
 *
 *
 *   
 *   \subsection{Worked Example}
 *
 *   This is the code:
 *   \begin{verbatim}
 *     #include <iostream.h>
 *     #include "../include/parameter_handler.h"
 *     
 *     
 *     class LinEq {
 *       public:
 *         static void declare_parameters (ParameterHandler &prm);
 *         void get_parameters (ParameterHandler &prm);
 *       private:
 *         string Method;
 *         int    MaxIterations;
 *     };
 *     
 *     
 *     class Problem {
 *       private:
 *         LinEq eq1, eq2;
 *         string Matrix1, Matrix2;
 *         string outfile;
 *       public:
 *         static void declare_parameters (ParameterHandler &prm);
 *         void get_parameters (ParameterHandler &prm);
 *     };
 *     
 *     
 *     
 *     void LinEq::declare_parameters (ParameterHandler &prm) {
 *                                        // declare parameters for the linear
 *                                        // solver in a subsection
 *       prm.enter_subsection ("Linear solver");
 *       prm.declare_entry ("Solver",
 *                          "CG",
 *                          Patterns::Sequence("CG|BiCGStab|GMRES"));
 *       prm.declare_entry ("Maximum number of iterations",
 *                          "20",
 *                          Patterns::Integer());
 *       prm.leave_subsection ();
 *     };
 *           
 *           
 *     void LinEq::get_parameters (ParameterHandler &prm) {
 *       prm.enter_subsection ("Linear solver");
 *       Method        = prm.get ("Solver");
 *       MaxIterations = prm.get_integer ("Maximum number of iterations");
 *       prm.leave_subsection ();
 *       cout << "  LinEq: Method=" << Method << ", MaxIterations=" << MaxIterations << endl;
 *     };
 *           
 *           
 *           
 *     void Problem::declare_parameters (ParameterHandler &prm) {
 *                                        // first some global parameter entries
 *       prm.declare_entry ("Output file",
 *                          "out",
 *                          Patterns::Anything());
 *       prm.declare_entry ("Equation 1",
 *                          "Laplace",
 *                          Patterns::Anything());
 *       prm.declare_entry ("Equation 2",
 *                          "Elasticity",
 *                          Patterns::Anything());
 *     
 *                                        // declare parameters for the
 *                                        // first equation
 *       prm.enter_subsection ("Equation 1");
 *       prm.declare_entry ("Matrix type",
 *                          "Sparse",
 *                          Patterns::Sequence("Full|Sparse|Diagonal"));
 *       LinEq::declare_parameters (prm);  // for eq1
 *       prm.leave_subsection ();
 *           
 *                                        // declare parameters for the
 *                                        // second equation
 *       prm.enter_subsection ("Equation 2");
 *       prm.declare_entry ("Matrix type",
 *                          "Sparse",
 *                          Patterns::Sequence("Full|Sparse|Diagonal"));
 *       LinEq::declare_parameters (prm);  // for eq2
 *       prm.leave_subsection ();
 *     };
 *           
 *           
 *     void Problem::get_parameters (ParameterHandler &prm) {
 *                                        // entries of the problem class
 *       outfile = prm.get ("Output file");
 *     
 *       string equation1 = prm.get ("Equation 1"),
 *              equation2 = prm.get ("Equation 2");
 *       
 *                                        // get parameters for the
 *                                        // first equation
 *       prm.enter_subsection ("Equation 1");
 *       Matrix1 = prm.get ("Matrix type");
 *       eq1.get_parameters (prm);         // for eq1
 *       prm.leave_subsection ();
 *           
 *                                        // get parameters for the
 *                                        // second equation
 *       prm.enter_subsection ("Equation 2");
 *       Matrix2 = prm.get ("Matrix type");
 *       eq2.get_parameters (prm);         // for eq2
 *       prm.leave_subsection ();
 *     
 *       cout << "  Problem: outfile=" << outfile << endl
 *            << "           eq1="     << equation1 << ", eq2=" << equation2 << endl
 *            << "           Matrix1=" << Matrix1 << ", Matrix2=" << Matrix2 << endl;
 *     };
 *          
 *           
 *           
 *           
 *     void main () {
 *       ParameterHandler prm;
 *       Problem p;
 *           
 *       p.declare_parameters (prm);
 *           
 *                                        // read input from "prmtest.prm"; giving
 *                                        // argv[1] would also be a good idea
 *       prm.read_input ("prmtest.prm");
 *           
 *                                        // print parameters to cout as ASCII text
 *       cout << endl << endl;
 *       prm.print_parameters (cout, Text);
 *           
 *                                        // get parameters into the program
 *       cout << endl << endl
 *            << "Getting parameters:" << endl;
 *       p.get_parameters (prm);
 *     };
 *   \end{verbatim}
 *
 *   
 *   This is the input file (named "prmtest.prm"):
 *   \begin{verbatim}
 *                                 # first declare the types of equations
 *     set Equation 1 = Poisson
 *     set Equation 2 = Navier-Stokes
 *
 *     subsection Equation 1
 *       set Matrix type = Sparse
 *       subsection Linear solver    # parameters for linear solver 1
 *         set Solver                       = Gauss-Seidel
 *         set Maximum number of iterations = 40
 *       end
 *     end
 *
 *     subsection Equation 2
 *       set Matrix type = Full
 *       subsection Linear solver
 *         set Solver                       = CG
 *         set Maximum number of iterations = 100
 *       end
 *     end
 *   \end{verbatim}
 *
 *   And here is the ouput of the program:
 *   \begin{verbatim}
 *     Line 8:
 *         The entry value
 *             Gauss-Seidel
 *         for the entry named
 *             Solver
 *         does not match the given regular expression
 *             CG|BiCGStab|GMRES
 *
 *
 *     Listing of Parameters
 *     ---------------------
 *       Equation 1  = Poisson  <Laplace>
 *       Equation 2  = Navier-Stokes  <Elasticity>
 *       Output file= out
 *       subsection Equation 1
 *         Matrix type = Sparse  <Sparse>
 *         subsection Linear solver
 *           Maximum number of iterations = 40  <20>
 *           Solver                      = CG
 *       subsection Equation 2
 *         Matrix type = Full  <Sparse>
 *         subsection Linear solver
 *           Maximum number of iterations = 100  <20>
 *           Solver                       = CG  <CG>
 *
 *
 *     Getting parameters:
 *       LinEq: Method=CG, MaxIterations=40
 *       LinEq: Method=CG, MaxIterations=100
 *       Problem: outfile=out
 *                eq1=Poisson, eq2=Navier-Stokes
 *                Matrix1=Sparse, Matrix2=Full
 *   \end{verbatim}
 *
 *   
 *   \subsection{References}
 *
 *   This class is inspired by the #MenuSystem# class of #DiffPack#.
 *
 *   @memo  This class provides a standard interface to an input file
 *   which provides at run-time for program parameters such as time step sizes,
 *   geometries, right hand sides etc.
 *   
 *   @author Wolfgang Bangerth, October 1997, revised February 1998
 *   @see MultipleParameterLoop
 */
class ParameterHandler {
  public:
				     /**
				      * Constructor.
				      */
    ParameterHandler ();

				     /**
				      * Destructor. Declare this only to have
				      * a virtual destructor, which is safer
				      * as we have virtual functions.
				      * It actually does nothing spectacular.
				      */
    virtual ~ParameterHandler ();
    
				     /**
				      * Read input from a stream until stream
				      * returns #eof# condition or error.
				      *
				      * Return whether the read was successful.
				      */
    virtual bool read_input (istream &input);
    
    				     /**
				      * Read input from a file the name of which
				      * is given.
				      *
				      * Return whether the read was successful.
				      */
    virtual bool read_input (const string &filename);
    
    				     /**
				      * Read input from a string in memory. The
				      * lines in memory have to be separated by
				      * #\n# characters.
				      *
				      * Return whether the read was successful.
				      */
    virtual bool read_input_from_string (const char *s);

				     /**
				      * Return status of this object:
				      * #true#=clean or #false#=error occured.
				      */
    bool ok() const   { return status; };

				     /**
				      * clear status bit and contents.
				      */
    void clear ();


				     /**
				      * Declare a new entry with name #entry#,
				      * default and for which
				      * any input has to match the #pattern#
				      * (default: any pattern).
				      * @return #false# if entry already exists or
				      * default value does not match the regular
				      * expression; #true# otherwise.
				      */
    bool declare_entry    (const string &entry,
			   const string &default_value,
			   const Patterns::PatternBase &pattern
			   = Patterns::Anything());
    
				     /**
				      * Enter a subsection; if not yet
				      * existent, declare it.
				      */
    void enter_subsection (const string &subsection);
    
				     /**
				      * Leave present subsection.
				      * @return #false# if there is no subsection
				      * to leave; true otherwise.
				      */
    bool leave_subsection ();

				     /**
				      * Return value of entry #entry_string#.
				      * If the entry was changed, then the changed
				      * value is returned, otherwise the default
				      * value. If the value of an undeclared entry
				      * is required, an empty string is returned and
				      * #assert# is used to check whether this entry
				      * was declared (therefore an exception may be
				      * thrown).
				      */
    const string & get (const string &entry_string) const;
    
				     /**
				      * Return value of entry #entry_string# as
				      * #long int#.
				      */
    long int       get_integer (const string &entry_string) const;
    
				     /**
				      * Return value of entry #entry_string# as
				      * #double#.
				      */
    double         get_double (const string &entry_string) const;

				     /**
				      * Return value of entry #entry_string# as
				      * #bool#.
				      */
    bool           get_bool (const string &entry_string) const;
				     /**
				      * Print all parameters with the given style
				      * to #out#. Presently only Text and LaTeX
				      * are implemented.
				      */
    ostream & print_parameters (ostream &out, const OutputStyle style);

				     /**
				      * Print out the parameters of the subsection
				      * given by the #subsection_path# member
				      * variable.
				      */
    void print_parameters_section (ostream &out,
				   const OutputStyle Style,
				   const unsigned int indent_level);
    
				   

				     /**
				      * Exception
				      */
    DeclException1 (ExcEntryAlreadyExists,
		    string,
		    << "The following entry already exists: " << arg1);
				     /**
				      * Exception
				      */
    DeclException2 (ExcDefaultDoesNotMatchPattern,
		    string, string,
		    << "The default string <" << arg1
		    << "> does not match the given pattern <" << arg2 << ">");  
				     /**
				      * Exception
				      */
    DeclException0 (ExcAlreadyAtTopLevel);
    				     /**
				      * Exception
				      */
    DeclException1 (ExcEntryUndeclared,
		    string,
		    << "You cant ask for entry <" << arg1 << "> you have not yet declared");  
    				     /**
				      * Exception
				      */
    DeclException1 (ExcConversionError,
		    string,
		    << "Error when trying to convert the following string: " << arg1);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);

    
  private:
				     /**
				      * ok bit
				      */
    bool status;

				     /**
				      * Whatever is in a section:
				      * map of entry names together with
				      * entry content and regexp, and
				      * list of subsections.
				      */
    struct Section {
	~Section ();

	typedef map<string, pair<string,Patterns::PatternBase*> > EntryType;
	
	EntryType             entries;
	map<string, Section*> subsections;
    };

				     /**
				      * Path of presently selected subsections;
				      * empty list means top level
				      */
    vector<string> subsection_path;

				     /**
				      * List of default values organized as a
				      * tree of subsections
				      */
    Section defaults;
    
				     /**
				      * Analogue list of changed entries. The
				      * tree of subsections is there even if there
				      * are no changed entry values in a
				      * subsection; therefore #enter_subsection#
				      * has to create the tree in both #Defaults#
				      * and #changed_entries.
				      */
    Section changed_entries;

				     /**
				      * Scan one line of input.
				      * #lineno# is the number of the line presently
				      * scanned (for the logs if there are messages).
				      * @return #false# if line contained stuff
				      * that could not be understood, the uppermost
				      * subsection was to be left by an #END# or
				      * #end# statement, a value for a non-declared
				      * entry was given or teh entry value did not
				      * match the regular expression. #true#
				      * otherwise
				      */
    bool scan_line (string line, const unsigned int lineno);

				     /**
				      * Get a pointer to the #Section# structure
				      * in the #Defaults# tree
				      * for the subsection we are presently in.
				      */
    Section*       get_present_defaults_subsection ();
    
				     /**
				      * Same, #const# version.
				      */ 
    const Section* get_present_defaults_subsection () const;

				     /**
				      * Get a pointer to the #Section# structure
				      * in the #changed_entries# tree
				      * for the subsection we are presently in.
				      */
    Section* get_present_changed_subsection ();
    
    				     /**
				      * Same, #const# version.
				      */
    const Section* get_present_changed_subsection () const;

    friend class MultipleParameterLoop;
};

	



/**
 *   The class #MultipleParameterLoop# offers an easy possibility to test several
 *   parameter sets during one run of the program. For this it uses the
 *   #ParameterHandler# class to read in data in a standardized form, searches for
 *   variant entry values and performs a loop over all combinations of parameters.
 *
 *   Variant entry values are given like this:
 *   \begin{verbatim}
 *     set Time step size = { 0.1 | 0.2 | 0.3 }
 *   \end{verbatim}
 *   The loop will then perform three runs of the program, one for each value
 *   of #Time step size#, while all other parameters are as specified or with their
 *   default value. If there are several variant entry values in the input a loop is
 *   performed for each combination of variant values:
 *   \begin{verbatim}
 *     set Time step size = { 0.1 | 0.2 }
 *     set Solver         = { CG  | GMRES }
 *   \end{verbatim}
 *   will result in four runs of the programs, with time step 0.1 and 0.2 for each
 *   of the two solvers.
 *
 *   Opposite to a variant entry, an array entry looks like this:
 *   \begin{verabtim}
 *     set Output file = ofile.{{ 1 | 2 | 3 | 4 }}
 *   \end{verbatim}
 *   This indicates that if there are variant entries producing a total of four
 *   different runs will write their results to the files #ofile.1#, #ofile.2#,
 *   #ofile.3# and #ofile.4#, respectively. Array entries do not generate multiple
 *   runs of the main loop themselves, but if there are variant entries, then in
 *   the #n#th run of the main loop, also the #n#th value of an array is returned.
 *
 *   Since the different variants are constructed in the order of declaration, not in
 *   the order in which the variat entries appear in the input file, it may be
 *   difficult to guess the mapping between the different variants and the appropriate
 *   entry in an array. You will have to check the order of declaration, or use
 *   only one variant entry.
 *   
 *   It is guaranteed that only selections which match the regular expression given
 *   upon declaration of an entry are given back to the program. If a variant value
 *   does not match the regular expression, the default value is stored and an error
 *   is issued. Before the first run of the loop, all possible values are checked
 *   for their conformance, so that the error is issued at the very beginning of the
 *   program.
 *
 *   
 *   \subsection{Usage}
 *   
 *   The usage of this class is similar to the #ParameterHandler# class. First the
 *   entries and subsections have to be declared, then a loop is performed in which
 *   the different parameter sets are set, a new instance of a user class is created
 *   which is then called. Taking the classes of the example for the
 *   #ParameterHandler# class, the extended program would look like this:
 *   \begin{verbatim}
 *     class HelperClass : public MultipleParameterLoop::UserClass {
 *       public:
 *         HelperClass ();
 *                     
 *         virtual void create_new (unsigned int runNo);
 *         virtual void declare_parameters (ParameterHandler &prm); 
 *         virtual void run (ParameterHandler &prm);
 *       private:
 *         Problem *p;
 *     };
 *           
 *           
 *     HelperClass::HelperClass () : p(0) {};
 *           
 *           
 *     void HelperClass::create_new (unsigned int runNo) {
 *       if (p) delete p;
 *       p = new Problem;
 *     };
 *                    
 *           
 *     void HelperClass::declare_parameters (ParameterHandler &prm) {
 *           				   // entries of the problem class
 *             			   // note: must be static member!
 *       Problem::declare_parameters (prm);
 *     };
 *           
 *           
 *     void HelperClass::run (ParameterHandler &prm) {
 *       p->get_parameters (prm);
 *     };
 *               
 *           
 *           
 *     void main () {
 *                                        // don't know why we have to write
 *     				   // "class" here, but it won't work
 *     				   // otherwise
 *       class MultipleParameterLoop prm;
 *       HelperClass h;
 *           
 *       h.declare_parameters (prm);
 *       prm.read_input ("prmtest.prm");
 *       prm.loop (h);
 *     };    
 *   \end{verbatim}
 *         
 *   As can be seen, first a new helper class has to be set up. This must contain
 *   a virtual constructor for a problem class. You can also derive your problem
 *   class from #MultipleParameterLoop::UserClass# and let #create_new# clear all
 *   member variables. If you have access to all inherited member variables in
 *   some way this is the recommended procedure. A third possibility is to use
 *   multiple inheritance and derive a helper class from both the
 *   #MultipleParameterLoop::UserClass# and the problem class. In any case,
 *   #create_new# has to provide a clean problem object which is the problem in
 *   the second and third possibility. However, if possible, the second way should
 *   be chosen.
 *
 *   The derived class also
 *   has to provide for member functions which declare the entries and which run
 *   the program. Running the program includes getting the parameters out of the
 *   #ParameterHandler# object.
 *
 *   After defining an object of this helper class and an object of the
 *   #MultipleParameterLoop# class, the entries have to be declared in the same way
 *   as for the #ParameterHandler# class. Then the input has to be read. Finally
 *   the loop is called. This executes the following steps:
 *   \begin{verbatim}
 *     for each combination
 *       {
 *         UserObject.create_new (runNo);
 *
 *  set parameters for this run
 *  
 *         UserObject.run (*this);
 *       };
 *   \end{verbatim}
 *   #UserObject# is the parameter to the #loop# function. #create_new# is given the number
 *   of the run (starting from one) to enable naming output files differently for each
 *   run.
 *
 *   
 *   \subsection{Syntax for variant and array entry values}
 *   
 *   Variant values are specified like #prefix{ v1 | v2 | v3 | ... }postfix#. Whitespace
 *   to the right of the opening brace #{# is ignored as well as to the left of the
 *   closing brace #}# while whitespace on the respectively other side is not ignored.
 *   Whitespace around the mid symbols #|# is also ignored. The empty selection
 *   #prefix{ v1 | }postfix# is also allowed and produces the strings #prefixv1postfix# and
 *   #prefixpostfix#.
 *
 *   The syntax for array values is equal, apart from the double braces:
 *   #prefix{{ v1 | v2 | v3 }}postfix#.
 *   
 *
 *   \subsection{Worked example}
 *   
 *   Given the above extensions to the example program for the #ParameterHandler# and the
 *   following input file
 *   \begin{verbatim}
 *     set Equation 1 = Poisson
 *     set Equation 2 = Navier-Stokes
 *     set Output file= results.{{ 1 | 2 | 3 | 4 | 5 | 6 }}
 *
 *     subsection Equation 1
 *       set Matrix type = Sparse
 *       subsection Linear solver 
 *         set Solver                       = CG
 *         set Maximum number of iterations = { 10 | 20 | 30 }
 *       end
 *     end
 *
 *     subsection Equation 2
 *       set Matrix type = Full
 *       subsection Linear solver
 *         set Solver                       = { BiCGStab | GMRES }
 *         set Maximum number of iterations = 100
 *       end
 *     end
 *   \end{verbatim}
 *   this is the output:
 *   \begin{verbatim}
 *     LinEq: Method=CG, MaxIterations=10
 *     LinEq: Method=BiCGStab, MaxIterations=100
 *     Problem: outfile=results.1
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *     LinEq: Method=CG, MaxIterations=20
 *     LinEq: Method=BiCGStab, MaxIterations=100
 *     Problem: outfile=results.2
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *     LinEq: Method=CG, MaxIterations=30
 *     LinEq: Method=BiCGStab, MaxIterations=100
 *     Problem: outfile=results.3
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *     LinEq: Method=CG, MaxIterations=10
 *     LinEq: Method=GMRES, MaxIterations=100
 *     Problem: outfile=results.4
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *     LinEq: Method=CG, MaxIterations=20
 *     LinEq: Method=GMRES, MaxIterations=100
 *     Problem: outfile=results.5
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *     LinEq: Method=CG, MaxIterations=30
 *     LinEq: Method=GMRES, MaxIterations=100
 *     Problem: outfile=results.6
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *   \end{verbatim}
 *   Since #create_new# gets the number of the run it would also be possible to output
 *   the number of the run.
 *   
 *   
 *   \subsection{References}
 *   This class is inspired by the #Multipleloop# class of #DiffPack#.
 *
 *   @memo  This class provides an interface to an input file which provides at
 *   run-time for multiple program parameters sets. The class performs a loop over
 *   all combinations of parameter sets.
 *   
 *   @author Wolfgang Bangerth, October 1997
 *   @version 1.0
 *   @see ParameterHandler
 */
class MultipleParameterLoop : public ParameterHandler {
  public:
				     /**
				      * This is the class the helper class or the
				      * problem class has to be derived of.
				      */
    class UserClass {
      public:
					 /**
					  * #create_new# must provide a clean
					  * object, either by creating a new one
					  * or by cleaning an old one.
					  */
	virtual void create_new (const unsigned int runNo) = 0;
	
					 /**
					  * This should declare parameters and call
					  * the #declare_parameters# function of the
					  * problem class.
					  */
	virtual void declare_parameters (ParameterHandler &prm) = 0;
	
					 /**
					  * Get the parameters and run any
					  * necessary action.
					  */
	virtual void run (ParameterHandler &prm) = 0;
    };

				     /**
				      * Constructor
				      */
    MultipleParameterLoop ();

				     /**
				      * Destructor. Declare this only to have
				      * a virtual destructor, which is safer
				      * as we have virtual functions.
				      * It actually does nothing spectacular.
				      */
    virtual ~MultipleParameterLoop ();

    				     /**
				      * Read input from a stream until stream
				      * returns #eof# condition or error.
				      */
    virtual bool read_input (istream &Input);
    
    				     /**
				      * Read input from a file the name of which
				      * is given.
				      */
    virtual bool read_input (const string &FileName);
    
    				     /**
				      * Read input from a string in memory. The
				      *  lines in memory have to be separated by
				      *  #\n# characters.
				      */
    virtual bool read_input_from_string (const char *s);

				     /**
				      * run the central loop.
				      */
    void loop (UserClass &uc);

  private:
				     /**
				      *	Declare what a multiple entry is: a variant
				      *	entry (in curly braces #{}#) or an
				      * array (in double curly braces #{{}}#).
				      */
    enum MultipleEntryType {
	  variant, array
    };

				     /**
				      * An object in the list of entries with
				      * multiple values.
				      */
    class Entry {
      public:
					 /**
					  * Constructor
					  */
	Entry () {};
	
					 /**
					  * Construct an object with given subsection
					  * path, name and value. The splitting up
					  * into the different variants is done
					  * later by #split_different_values.
					  */
	Entry (const vector<string> &Path, const string &Name, const string &Value);

					 /**
					  * Split the entry value into the different
					  * branches.
					  */
	void split_different_values ();

					 /**
					  * Path to variant entry.
					  */
	vector<string> subsection_path;

					 /**
					  * Name of entry.
					  */
	string         entry_name;

					 /**
					  * Original variant value.
					  */
	string         entry_value;
	
					 /**
					  * List of entry values constructed out of
					  * what was given in the input file (that
					  * is stored in #EntryValue#.
					  */
	vector<string> different_values;

					 /**
					  * Store whether this entry is a variant
					  * entry or an array.
					  */
	MultipleEntryType      type;
    };

				     /**
				      * List of variant entry values.
				      */
    vector<Entry> multiple_choices;
    
				     /**
				      * Number of branches constructed from the
				      * different combinations of the variants.
				      * This obviously equals the number of runs
				      * to be performed.
				      */
    int n_branches;

				     /**
				      * Initialize the different branches, i.e.
				      * construct the combinations.
				      */
    void init_branches ();
    
				     /**
				      * Initialize the branches in the given
				      * section.
				      */
    void init_branches_section (const ParameterHandler::Section &sec);
    
				     /**
				      * Transfer the entry values for one run
				      * to the entry tree.
				      */
    void fill_entry_values (const unsigned int run_no);
};




/**********************   parameter-handler.h     ****************************/
/* end of #ifndef __parameter_handler_H */
#endif
/**********************   parameter-handler.h     ****************************/
