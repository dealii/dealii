#ifndef LOCAL_ASSEMBLE_STANDARD
#define LOCAL_ASSEMBLE_STANDARD

#include "local_assemble_base.h"
#include <base/logstream.h>
#include <base/smartpointer.h>
#include <fe/fe_values.h>
#include <fe/fe.h>

#include <fstream>
#include <iostream>
#include <base/parameter_handler.h>


template <int dim, typename DH=MGDoFHandler<dim> >
  class LocalAssemble : public LocalAssembleBase<dim, DH>
{
public:
#if deal_II_dimension != 1
  LocalAssemble();
  
  /** Given two face finite elements, assemble the two matrices.*/
  virtual void  assemble_face_terms (FEFaceValuesBase<dim> &fe_v,
				     FEFaceValuesBase<dim> &fe_n_v,
				     FullMatrix<double> &,
				     FullMatrix<double> &);

  virtual void  assemble_face_term
  (const typename DH::active_cell_iterator& , 
   const unsigned int,
   const typename DH::active_cell_iterator& , 
   const unsigned int,
   FullMatrix<double> &,
   FullMatrix<double> &);

  virtual void  assemble_face_term
  (const typename DH::active_cell_iterator& , 
   const unsigned int,
   const unsigned int,
   const typename DH::active_cell_iterator& , 
   const unsigned int,
   FullMatrix<double> &,
   FullMatrix<double> &);

  virtual void  assemble_face_term
  (const typename DH::active_cell_iterator& , 
   const unsigned int,
   const typename DH::active_cell_iterator& , 
   const unsigned int,
   const unsigned int,
   FullMatrix<double> &,
   FullMatrix<double> &);

  /** A pointer to fe_face_values objects. */
  SmartPointer<FEFaceValues<dim> > fe_face_v;

  /** A pointer to fe_face_values objects. */
  SmartPointer<FEFaceValues<dim> > fe_face_n_v;

  /** A pointer to fe_face_values objects. */
  SmartPointer<FESubfaceValues<dim> > fe_sub_face_v;
#endif

};

#endif
