#ifndef __deal2__local_assemble_base_h
#define __deal2__local_assemble_base_h

#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <fe/fe_values.h>
#include "assemble_flags.h"

#include <fstream>
#include <iostream>

using namespace dealii;

/** Base class for local assemblers. This is the object you should
  derive your own classes from, if you want to use the VirtualMatrix
  class or the MyTools::assemble utility. It provides a common
  interface to the procedure for assembling locally the matrices. All
  functions are virtual and would throw an exception if called. You
  should really implement your own class deriving it from this one. An
  example is provided in the LocalAssembleLaplace class, which
  assemble the Laplace equations locally for continuous finite
  elements.

  All methods work basically in the same way. They take a FEValues
  object and a FullMatrix one as input, and fill the FullMatrix with
  appropriate values. Exceptions are the face and boundary terms,
  which take FEFaceValues as input.

  This object is supposed to be used in a Virtual Matrix via the 
  VirtualMatrix::enter method.

  Here the only thing which is actually done is to create the internal 
  AssembleFlags object upon construction, which is used to determine which
  of the methods will be called by the virtual matrix.
 */
template <int dim, typename DH=MGDoFHandler<dim> >
class LocalAssembleBase : public Subscriptor
{
    public:
	virtual ~LocalAssembleBase() {};
	/** This object will be called for each cell of a triangulation. */
	// virtual void assemble_cell_term(const FEValues<dim>& fe_v,
	// 	FullMatrix<double> &u_v_matrix) const;
	virtual void assemble_cell_term
	    (const typename DH::active_cell_iterator&, 
	     FullMatrix<double> &cell_m);

	/** This object will be called for each boundary face of a triangulation.*/
	virtual void assemble_boundary_term
	    (const typename DH::active_cell_iterator&, const unsigned int, 
	     FullMatrix<double> &cell_m);

	/** This object will be called for each face of a
	    triangulation. This one is called when the face is shared
	    by neighbors on the same level of refinement. */
	virtual void assemble_face_term
	  (const typename DH::active_cell_iterator &cell, 
	   const unsigned int face_no, 
	   const typename DH::active_cell_iterator &neighbor, 
	   const unsigned int neighbor_face_no, 
	   FullMatrix<double> &cell_m,
	   FullMatrix<double> &neighbor_cell_m);
	
	/** This object will be called for each face of a
	    triangulation. This one is called when the current face is
	    coarser. */
	virtual void assemble_face_term
	  (const typename DH::active_cell_iterator &cell, 
	   const unsigned int face_no, 
	   const unsigned int sub_face_no, 
	   const typename DH::active_cell_iterator &neighbor, 
	   const unsigned int neighbor_face_no, 
	   FullMatrix<double> &cell_m,
	   FullMatrix<double> &neighbor_cell_m);
	
	/** This object will be called for each face of a
	    triangulation. This one is called when the current face is
	    finer. */
	virtual void assemble_face_term
	  (const typename DH::active_cell_iterator &cell, 
	   const unsigned int face_no, 
	   const typename DH::active_cell_iterator &neighbor, 
	   const unsigned int neighbor_face_no, 
	   const unsigned int neighbor_subface_no,
	   FullMatrix<double> &cell_m,
	   FullMatrix<double> &neighbor_cell_m);
	
	/** Assemble rhs. This will be called for each cell.*/
	virtual void assemble_rhs_term
	  (const typename DH::active_cell_iterator&, 
	   Vector<double> &);

	/** This object will be called for each boundary face of a triangulation.*/
	virtual void assemble_rhs_boundary_term
	  (const typename DH::active_cell_iterator&, const unsigned int, 
	   Vector<double> &);
	
	AssembleFlags flags;
};

#endif
