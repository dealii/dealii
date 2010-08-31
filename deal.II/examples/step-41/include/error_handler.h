#ifndef LH_ERROR_HANDLER_H
#define LH_ERROR_HANDLER_H

#include <fstream>

#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <lac/vector.h>

#include <grid/tria.h>

// #include <numerics/error_estimator.h>
#include <base/function.h>
#include <numerics/solution_transfer.h>
#include <numerics/data_out.h>

#include <base/convergence_table.h>
#include <base/logstream.h>
#include <base/config.h>

#include <base/parameter_handler.h>

#include <map>

enum NormFlags { 
    None = 0x00, 
    Linfty = 0x01,
    L2 = 0x02, 
    W1infty = 0x04,
    H1 = 0x08,
    AddUp = 0x10
};

using namespace dealii;

template <int dim, typename VEC=Vector<double> >
class ErrorHandler : public Subscriptor
{
public:
    /** The constructor takes the mpi initialization stuff. */
    ErrorHandler ();
    
    /** Initialize the given values for the paramter file. */
    static void declare_parameters(ParameterHandler &prm, 
                                   unsigned int ntables=1);
    
    /** Parse the given parameter handler. */
    void parse_parameters(ParameterHandler &prm);
  
    /** Calculate the error of the numeric solution in variuous norms. Store
        the result in the given table. */
    void error_from_exact(const DoFHandler<dim> & vspace, 
                          const VEC &solution, 
                          const Function<dim> &exact,
                          unsigned int table_no = 0,
                          double dt=0.); 
  
    /** Difference between two solutions in two different vector spaces. */
    void difference(const DoFHandler<dim> &, const VEC &,
                    const DoFHandler<dim> &, const VEC &, 
                    unsigned int table_no = 0, double dt=0.); 
  
    /** Difference between two solutions in the same vector space. */
    void difference(const DoFHandler<dim> &, const VEC &, 
                    const VEC &, unsigned int table_no = 0, double dt=0.);
    
    /** By default output first table. */
    void output_table(const unsigned int table_no=0);
    
private:  
    /** Error results.*/
    std::vector<ConvergenceTable>  tables;  
    
    /** Headers for tables and output. Contains the name of the solution
        components. */
    std::vector<std::string> headers;
    
    /** Headers for latex tables. Contains the name of the solution
        components. */
    std::vector<std::string> latex_headers;
    
    /** Captions for latex. */
    std::vector<std::string> latex_captions;
    
    /** Names of the tables. */
    std::vector<std::string> names;
    
    /** Type of error to compute per components. */
    std::vector<std::vector<NormFlags> > types;
    
    /** The parameters have been read. */
    bool initialized;
    
    /** Compute the error. If this is false, all functions regarding
        errors are disabled and don't do anything.*/
    bool compute_error;
    
    /** Add convergence rates. */
    std::vector<bool> add_rates;
    
    /** Write the error files. */
    bool write_error;
    
    /** Output the error file also on screen. */
    bool output_error;
    
    /** The error file format. */
    std::string error_file_format;
    
    /** The extra column to add to the tables. */    
    std::vector<std::map<std::string, bool> > extras;

    /** Wether or not to calculate the rates according to the given keys. */
    std::vector<std::string> rate_keys;
};

/**
 * Global operator which returns an object in which all bits are set
 * which are either set in the first or the second argument. This
 * operator exists since if it did not then the result of the bit-or
 * <tt>operator |</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * NormFlags.
 */
inline
NormFlags
operator | (NormFlags f1, NormFlags f2)
{
  return static_cast<NormFlags> (
				 static_cast<unsigned int> (f1) |
				 static_cast<unsigned int> (f2));
}

/**
 * Global operator which sets the bits from the second argument also
 * in the first one.
 */
inline
NormFlags &
operator |= (NormFlags &f1, NormFlags f2)
{
  f1 = f1 | f2;
  return f1;
}


/**
 * Global operator which returns an object in which all bits are set
 * which are set in the first as well as the second argument. This
 * operator exists since if it did not then the result of the bit-and
 * <tt>operator &</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * NormFlags.
 */
inline
NormFlags
operator & (NormFlags f1, NormFlags f2)
{
  return static_cast<NormFlags> (
				 static_cast<unsigned int> (f1) &
				 static_cast<unsigned int> (f2));
}


/**
 * Global operator which clears all the bits in the first argument if
 * they are not also set in the second argument.
 */
inline
NormFlags &
operator &= (NormFlags &f1, NormFlags f2)
{
  f1 = f1 & f2;
  return f1;
}

#endif
