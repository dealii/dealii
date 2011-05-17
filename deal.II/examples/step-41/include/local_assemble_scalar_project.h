#ifndef LOCAL_ASSEMBLE_SCALAR_PROJECT
#define LOCAL_ASSEMBLE_SCALAR_PROJECT

#include "local_assemble.h"
#include <deal.II/base/parsed_function.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe.h>

#include <fstream>
#include <iostream>
#include <deal.II/base/parameter_handler.h>


template <int dim>
class LocalAssembleScalarProject : public LocalAssembleBase<dim>
{
 public:
  
  LocalAssembleScalarProject(); 
  
  ~LocalAssembleScalarProject(); 
  
  virtual void assemble_rhs_term
    (const typename MGDoFHandler<dim>::active_cell_iterator&,
     Vector<double> &);


  void reinit(FiniteElement<dim>&,
	      Table<2, double > &,
	      Table<2, double > &);
			
 private:
  /** A pointer to the finite element.*/
  SmartPointer<FiniteElement<dim> > fe;

  /** A pointer to the plastic strain */
  Table<dim,double >   hard_table;
  Table<dim,double >   iter_table;

  /** A pointer to fe_values objects. */
  SmartPointer<FEValues<dim> > fe_v;

  
};

#endif
