#ifndef LOCAL_ASSEMBLE_ELASTIC_MATRIX
#define LOCAL_ASSEMBLE_ELASTIC_MATRIX

#include "local_assemble.h"
#include "parsed_symmetric_tensor_function.h"
#include <base/parsed_function.h>

#include <base/logstream.h>
#include <base/smartpointer.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <fe/fe_values.h>
#include <fe/fe.h>

#include <fstream>
#include <iostream>
#include <base/parameter_handler.h>


template <int dim>
class LocalAssembleElasticMatrix : public LocalAssembleBase<dim>
{
 public:
  
  LocalAssembleElasticMatrix(); 
  
  ~LocalAssembleElasticMatrix(); 
  
  virtual void  assemble_cell_term
    (const typename MGDoFHandler<dim>::active_cell_iterator& cell, 
     FullMatrix<double> &);


  void reinit(FiniteElement<dim>&);

  void parameters(ParameterHandler &prm);
			
 private:
  /** A pointer to the finite element.*/
  SmartPointer<FiniteElement<dim> > fe;

  /** A pointer to the rhs function.*/
  ParsedSymmetricTensorFunction<4, dim>  C;

  /** A pointer to fe_values objects. */
  SmartPointer<FEValues<dim> > fe_v;
  
};

#endif
