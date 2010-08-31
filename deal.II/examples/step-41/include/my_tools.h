//---------------------------------------------------------------------------
//    $Id: fe.h,v 1.124 2005/09/17 09:19:19 guido Exp $
//    Version: $Name:  $
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__my_tools_h
#define __deal2__my_tools_h

#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <multigrid/mg_base.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <fe/fe.h>
#include <base/quadrature_lib.h>

#include <dofs/dof_constraints.h>
#include "assemble_flags.h"
#include "local_assemble.h"
#include "vector_space.h"

/**
 * Some assembly routines. These are an extension of the already
 * existing deal.II assembly routines that make use of the support
 * class LocalAssembleBase. These routines call functions of the local
 * assemblers passing local matrices or vectors and the cell. It is
 * responsability of the local assemblers to fill in the informations
 * correctly. Unlike step-12 of the deal.II library, the cell iterator
 * is passed, not the fe_values. This means that the local assembler
 * should build its own fe_values object. It is its responsability to
 * select all the outer details on how to actually compute these
 * matrices...
 * 
 *
 * @author Luca Heltai,
 * 2005, 2008
 */

namespace MyTools {
    using namespace std;

    /** DG-Like assembly routine for matrices. Hanging node
     * constraints are not taken into account, and one should do this
     * outside. However face terms are correctly handled, and this
     * method is capable of assemblying fully DG objects. */
    template <int dim, typename DH, typename MATRIX>
    void assemble(DH &dof_handler, MATRIX& system_matrix, 
		  LocalAssembleBase<dim, DH> &local, 
		  const unsigned int this_mpi_process = 0) 
    {
	const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
	std::vector<unsigned int> dofs (dofs_per_cell);
	std::vector<unsigned int> dofs_neighbor (dofs_per_cell);
    
	// Now we create the cell matrices and vectors. Here we need two
	// cell matrices, both for face terms that include test functions
	// <code>vi</code> (internal shape functions, i.e. shape functions
	// of the current cell). 
	// 
	// i stands for internal, e stands for external
	FullMatrix<double> ui_vi_matrix (dofs_per_cell, dofs_per_cell);
	FullMatrix<double> ui_ve_matrix (dofs_per_cell, dofs_per_cell);

	// Furthermore we need some cell iterators.
	typename DH::active_cell_iterator
	    cell = dof_handler.begin_active(),
	    endc = dof_handler.end();

	// Now we start the loop over all
	// active cells.
	for (;cell!=endc; ++cell) 
	    if( cell->subdomain_id() == this_mpi_process )
	    {
		// On each cell we need to reset the <code>ui_vi_matrix</code>
		// and <code>cell_vector</code> to zero, before assembling the
		// cell terms.
		ui_vi_matrix = 0;
      
		// and call the function that assembles the cell terms. The
		// first argument is the <code>FEValues</code> that was
		// previously reinit'ed on the current cell.
		if( local.flags & assemble_cell ) 
		    local.assemble_cell_term(cell, ui_vi_matrix);
	  
		// As in previous examples the vector `dofs' includes the
		// dof_indices.
		cell->get_dof_indices (dofs);
	  
		// This is the start of the nested loop over all faces. 
		if( local.flags & (assemble_boundary|assemble_face) ) {
		    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
		    {
			// First we set the face iterator
			typename DH::face_iterator face=cell->face(face_no);
	  
			// and clear the <code>ui_ve_matrix</code> on each face.
			ui_ve_matrix = 0;
		  
			// Now we distinguish the four different cases in the
			// ordering mentioned above. We start with faces belonging
			// to the boundary of the domain.
			if (face->at_boundary())
			{
			    // and assemble the corresponding face terms.
			    if(local.flags & assemble_boundary) 
				local.assemble_boundary_term(cell, face_no,
							     ui_vi_matrix);
			}
			else if(local.flags & assemble_face) 
			{
			    // Now we are not on the boundary of the domain,
			    // therefore there must exist a neighboring cell.
			    typename DH::cell_iterator neighbor=
				cell->neighbor(face_no);

			    // We proceed with the second and most complicated case:
			    // the neighboring cell is more refined than the current
			    // cell. As in deal.II neighboring cells are restricted
			    // to have a level difference of not more than one, the
			    // neighboring cell is known to be at most ONCE more
			    // refined than the current cell. Furthermore also the
			    // face is more refined, i.e. it has children. Here we
			    // note that the following part of code will not work
			    // for <code>dim==1</code>.
			    if (face->has_children())
			    {
				// First we store which number the current cell has
				// in the list of neighbors of the neighboring
				// cell. Hence,
				// neighbor-@>neighbor(neighbor_face_no) equals the
				// current cell <code>cell</code>.
				const unsigned int neighbor_face_no=
				    cell->neighbor_of_neighbor(face_no);
		  
		  
				// We loop over subfaces
				for (unsigned int subface_no=0;
				     subface_no<face->n_children(); ++subface_no)
				{
				    // and set the cell iterator
				    // <code>neighbor_child</code> to the cell
				    // placed `behind' the current subface.
				    typename DH::active_cell_iterator
					neighbor_child
					= cell->neighbor_child_on_subface (face_no, subface_no);
		      
				    // As these are quite complicated indirections
				    // which one does not usually get right at first
				    // attempt we check for the internal
				    // consistency.
				    Assert (neighbor_child->face(neighbor_face_no) == face->child(subface_no),
					    ExcInternalError());
				    Assert (!neighbor_child->has_children(), ExcInternalError());

				    // We need to reset the
				    // <code>ui_ve_matrix</code> on each subface
				    // because on each subface the <code>un</code>
				    // belong to different neighboring cells.
				    ui_ve_matrix = 0;
		      
				    // As already mentioned above for the current
				    // case (case 2) we employ the
				    // <code>FESubfaceValues</code> of the current
				    // cell (here reinited for the current cell,
				    // face and subface) and we employ the
				    // FEFaceValues of the neighboring child cell.
				    local.assemble_face_term(cell, face_no, subface_no,
							     neighbor_child, neighbor_face_no,
							     ui_vi_matrix,
							     ui_ve_matrix);
		      
				    // Then we get the dof indices of the
				    // neighbor_child cell
				    neighbor_child->get_dof_indices (dofs_neighbor);
		      						
				    // and distribute <code>ui_ve_matrix</code> to
				    // the system_matrix
				    for (unsigned int i=0; i<dofs_per_cell; ++i)
					for (unsigned int k=0; k<dofs_per_cell; ++k)
					    system_matrix.add(dofs_neighbor[i], dofs[k],
							      ui_ve_matrix(i,k));
				}
				// End of <code>if (face-@>has_children())</code>
			    }
			    else
			    {
				// We proceed with case 3, i.e. neighboring cell is
				// of the same refinement level as the current cell.
				if (neighbor->level() == cell->level()) 
				{
				    // Like before we store which number the current
				    // cell has in the list of neighbors of the
				    // neighboring cell.
				    const unsigned int neighbor_face_no=cell->neighbor_of_neighbor(face_no);

				    local.assemble_face_term(cell, face_no,
							     neighbor, neighbor_face_no,
							     ui_vi_matrix,
							     ui_ve_matrix);
				    // End of <code>if (neighbor-@>level() ==
				    // cell-@>level())</code>
				}
				else
				{
				    // Finally we consider case 4. When the
				    // neighboring cell is not finer and not of the
				    // same refinement level as the current cell it
				    // must be coarser.
				    Assert(neighbor->level() < cell->level(), ExcInternalError());

				    // Find out the how many'th face_no and
				    // subface_no the current face is w.r.t. the
				    // neighboring cell.
				    const std::pair<unsigned int, unsigned int> faceno_subfaceno=
					cell->neighbor_of_coarser_neighbor(face_no);
				    const unsigned int neighbor_face_no=faceno_subfaceno.first,
					neighbor_subface_no=faceno_subfaceno.second;

				    Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
										 neighbor_subface_no)
					    == cell,
					    ExcInternalError());

				    local.assemble_face_term(cell, face_no,
							     neighbor, 
							     neighbor_face_no,
							     neighbor_subface_no,
							     ui_vi_matrix,
							     ui_ve_matrix);
				}

				// Now we get the dof indices of the
				// <code>neighbor</code> cell,
				neighbor->get_dof_indices (dofs_neighbor);
		 						
				// and distribute the
				// <code>ui_ve_matrix</code>.
				for (unsigned int i=0; i<dofs_per_cell; ++i)
				    for (unsigned int k=0; k<dofs_per_cell; ++k)
					system_matrix.add(dofs_neighbor[i], dofs[k],
							  ui_ve_matrix(i,k));
			    }
			    // End of <code>face not at boundary</code>:
			}
			// End of loop over all faces:
		    }
		}
	    
		// Finally we distribute the
		// <code>ui_vi_matrix</code>
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		    for (unsigned int j=0; j<dofs_per_cell; ++j)
			system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i,j));
	    }
    }


    /** MG assembly routine for matrices. Hanging node constraints are
     * taken into account. However face terms are NOT correctly
     * handled, that is, this method is currently not capable of
     * handling fluxex. 
     
     Different levels can be assembled. If the level is -1, then the
     active level is used.
     */
    template <int dim, typename MATRIX>
    void assemble(MGDoFHandler<dim> &mg_dof_handler, ConstraintMatrix &hang,
		  int level, MATRIX& m, LocalAssembleBase<dim, MGDoFHandler<dim> > &local, 
		  const unsigned int this_mpi_process = 0) 
    {
	AssertThrow( !(local.flags & assemble_face), ExcNotImplemented());
	// Cell iterators.
	typename MGDoFHandler<dim>::active_cell_iterator cell, endc;
	if(level < 0) { 
	    cell = mg_dof_handler.begin_active();
	    endc = mg_dof_handler.end();
	} else {
	    cell = mg_dof_handler.begin_active(level);	
	    endc = mg_dof_handler.end(level);
	}

	unsigned int dofs_per_cell = 
	    mg_dof_handler.get_fe().dofs_per_cell;

	// dofs and neighbor dofs.
	vector<unsigned int> dofs(dofs_per_cell);
	vector<unsigned int> neighbor_dofs(dofs_per_cell);

	FullMatrix<double> cell_m(dofs_per_cell, dofs_per_cell);
	FullMatrix<double> cell_neighbor_m(dofs_per_cell, dofs_per_cell);

	for (; cell!=endc; ++cell)
	    if( cell->subdomain_id() == this_mpi_process )
	    {
		// Reset local matrices.
		cell_m = 0;
		cell_neighbor_m = 0;
	  
		// Update indices
		if(level<0)
		    cell->get_dof_indices(dofs);
		else
		    cell->get_mg_dof_indices (dofs);

		// Assemble cell term if needed
		if( local.flags & assemble_cell ) {

		  //
		  local.assemble_cell_term(cell, cell_m);
		}
		for(unsigned int n=0; n<dofs.size(); ++n)
		hang.distribute_local_to_global(cell_m, dofs, m);
#if (DEAL_II_DIM > 1)
		// Assemble boundary or face terms if needed
		if( local.flags & (assemble_boundary|assemble_face) ) {
	    
		    for (unsigned int face_no=0; 
			 face_no<GeometryInfo<dim>::faces_per_cell; 
			 ++face_no)
		    {
			typename MGDoFHandler<dim>::face_iterator face=cell->face(face_no);
		
			if (face->at_boundary())
			{
			    if(local.flags & assemble_boundary)   {
				cell_m = 0;
				local.assemble_boundary_term(cell, face_no, cell_m);
				hang.distribute_local_to_global(cell_m, dofs, m);
			    }
			} // End of Boundary terms
			else if(local.flags & assemble_face) {
			    typename MGDoFHandler<dim>::cell_iterator neighbor=
				cell->neighbor(face_no);;
		  
			    if (face->has_children()) // The neighbor is more refined
			    {
				const unsigned int n_face_no=
				    cell->neighbor_of_neighbor(face_no);
		      
		      
				for (unsigned int subface_no=0;
				     subface_no<face->n_children(); ++subface_no)
				{
				    typename MGDoFHandler<dim>::active_cell_iterator
					neighbor_child
					= cell->neighbor_child_on_subface (face_no, subface_no);
			  
				    Assert (neighbor_child->face(n_face_no) == face->child(subface_no),
					    ExcInternalError());
				    Assert (!neighbor_child->has_children(), ExcInternalError());
			  
				    cell_m = 0;
				    cell_neighbor_m = 0;
				    local.assemble_face_term(cell, face_no, subface_no,
							     neighbor_child, n_face_no,
							     cell_m, cell_neighbor_m);
				    if(level<0)
					neighbor_child->get_dof_indices (neighbor_dofs);
				    else 
					neighbor_child->get_mg_dof_indices (neighbor_dofs);
			  
				    hang.distribute_local_to_global(cell_m, dofs, m);
				    // TBA Distribute neighbors!!!
				    // hang.distribute_local_to_global(cell_m, neighbor_dofs, dofs, m);
				    //


				}
			    } // face has children
			    else
			    {
				if (neighbor->level() == cell->level())
				{
				    const unsigned int n_face_no=cell->neighbor_of_neighbor(face_no);
			  
				    cell_neighbor_m = 0;
				    cell_m = 0;
				    local.assemble_face_term(cell, face_no, 
							     neighbor, n_face_no, 
							     cell_m,
							     cell_neighbor_m);
			  
				    if(level<0)
					neighbor->get_dof_indices (neighbor_dofs);
				    else 
					neighbor->get_mg_dof_indices (neighbor_dofs);
			  
				    hang.distribute_local_to_global(cell_m, dofs, m);
				    // TBA: distribute neighbor dofs....
				    // hang.distribute_local_to_global(cell_m, neighbor_dofs, dofs, m);
				    // 
			  
				} // Same level on the two neighbors
				else
				{
				    Assert(neighbor->level() < cell->level(), ExcInternalError());
			  
				    const std::pair<unsigned int, unsigned int> faceno_subfaceno=
					cell->neighbor_of_coarser_neighbor(face_no);
				    const unsigned int neighbor_face_no=faceno_subfaceno.first,
					neighbor_subface_no=faceno_subfaceno.second;
			  
				    Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
										 neighbor_subface_no)
					    == cell,
					    ExcInternalError());
			  
				    cell_neighbor_m = 0;
				    cell_m = 0;
				    local.assemble_face_term(cell, face_no,
							     neighbor, neighbor_face_no, 
							     neighbor_subface_no,
							     cell_m, cell_neighbor_m);
			  
				    if(level<0)
					neighbor->get_dof_indices (neighbor_dofs);
				    else 
					neighbor->get_mg_dof_indices (neighbor_dofs);
			  
				    hang.distribute_local_to_global(cell_m, dofs, m);
				    // TBA Distribute neighbor dofs
				    // hang.distribute_local_to_global(cell_m, neighbor_dofs, dofs, m);
				    //

				} // The neighbor is coarser
			    } // else : if face has children
			} // skipped if no flags for face terms
		    } // faces loop
		} // Skipped if no flags on either faces or boundary 
#else
    //1D face terms
#endif 
	    } // cells loop
    }
  
    /** MG assembly routine for matrices and vectors. Hanging node
     * constraints are taken into account. However face terms are NOT
     * correctly handled, that is, this method is currently not
     * capable of handling fluxex.
     
     Different levels can be assembled. If the level is -1, then the
     active level is used.
     */
    template <int dim, typename MATRIX, typename VEC>
    void assemble(MGDoFHandler<dim> &mg_dof_handler, ConstraintMatrix &hang,
		  int level, MATRIX& m, VEC& rhs, LocalAssembleBase<dim, MGDoFHandler<dim> > &local, 
		  const unsigned int this_mpi_process = 0) 
    {
	// Cell iterators.
	typename MGDoFHandler<dim>::active_cell_iterator cell, endc;
	if(level < 0) { 
	    cell = mg_dof_handler.begin_active();
	    endc = mg_dof_handler.end();
	} else {
	    cell = mg_dof_handler.begin_active(level);	
	    endc = mg_dof_handler.end(level);
	}

	unsigned int dofs_per_cell = 
	    mg_dof_handler.get_fe().dofs_per_cell;

	// dofs and neighbor dofs.
	vector<unsigned int> dofs(dofs_per_cell);
	vector<unsigned int> neighbor_dofs(dofs_per_cell);

	FullMatrix<double> cell_m(dofs_per_cell, dofs_per_cell);
	FullMatrix<double> cell_neighbor_m(dofs_per_cell, dofs_per_cell);
	Vector<double>	   local_rhs(dofs_per_cell);

	for (; cell!=endc; ++cell)
	    if( cell->subdomain_id() == this_mpi_process )
	    {
		// Reset local matrices.
		cell_m = 0;
		cell_neighbor_m = 0;
	  
		// Update indices
		if(level<0)
		    cell->get_dof_indices(dofs);
		else
		    cell->get_mg_dof_indices (dofs);

		// Assemble cell term if needed
		if( local.flags & assemble_cell ) 
		    local.assemble_cell_term(cell, cell_m);
		if( local.flags & assemble_rhs_cell)
		    local.assemble_rhs_term(cell, local_rhs);
		
		hang.distribute_local_to_global(cell_m, dofs, m);
		hang.distribute_local_to_global(local_rhs, dofs, rhs);

#if (DEAL_II_DIM > 1)
		// Assemble boundary or face terms if needed
		if( local.flags & (assemble_boundary|assemble_face) ) {
	    
		    for (unsigned int face_no=0; 
			 face_no<GeometryInfo<dim>::faces_per_cell; 
			 ++face_no)
		    {
			typename MGDoFHandler<dim>::face_iterator face=cell->face(face_no);
		
			if (face->at_boundary())
			{
			    if(local.flags & assemble_boundary)   {
				cell_m = 0;
				local.assemble_boundary_term(cell, face_no, cell_m);
				hang.distribute_local_to_global(cell_m, dofs, m);
			    }
			    if(local.flags & assemble_rhs_boundary) {
				local_rhs =0;
				local.assemble_rhs_boundary_term(cell, face_no, local_rhs);
				hang.distribute_local_to_global(local_rhs, dofs, rhs);
			    }
			} // End of Boundary terms
			else if(local.flags & assemble_face) {
			    typename MGDoFHandler<dim>::cell_iterator neighbor=
				cell->neighbor(face_no);;
		  
			    if (face->has_children()) // The neighbor is more refined
			    {
				const unsigned int n_face_no=
				    cell->neighbor_of_neighbor(face_no);
		      
		      
				for (unsigned int subface_no=0;
				     subface_no<face->n_children(); ++subface_no)
				{
				    typename MGDoFHandler<dim>::active_cell_iterator
					neighbor_child
					= cell->neighbor_child_on_subface (face_no, subface_no);
			  
				    Assert (neighbor_child->face(n_face_no) == face->child(subface_no),
					    ExcInternalError());
				    Assert (!neighbor_child->has_children(), ExcInternalError());
			  
				    cell_m = 0;
				    cell_neighbor_m = 0;
				    local.assemble_face_term(cell, face_no, subface_no,
							     neighbor_child, n_face_no,
							     cell_m, cell_neighbor_m);
				    if(level<0)
					neighbor_child->get_dof_indices (neighbor_dofs);
				    else 
					neighbor_child->get_mg_dof_indices (neighbor_dofs);
			  
				    hang.distribute_local_to_global(cell_m, dofs, m);
				    // TBA Distribute neighbors!!!
				    // hang.distribute_local_to_global(cell_m, neighbor_dofs, dofs, m);
				    //


				}
			    } // face has children
			    else
			    {
				if (neighbor->level() == cell->level())
				{
				    const unsigned int n_face_no=cell->neighbor_of_neighbor(face_no);
			  
				    cell_neighbor_m = 0;
				    cell_m = 0;
				    local.assemble_face_term(cell, face_no, 
							     neighbor, n_face_no, 
							     cell_m,
							     cell_neighbor_m);
			  
				    if(level<0)
					neighbor->get_dof_indices (neighbor_dofs);
				    else 
					neighbor->get_mg_dof_indices (neighbor_dofs);
			  
				    hang.distribute_local_to_global(cell_m, dofs, m);
				    // TBA: distribute neighbor dofs....
				    // hang.distribute_local_to_global(cell_m, neighbor_dofs, dofs, m);
				    // 
			  
				} // Same level on the two neighbors
				else
				{
				    Assert(neighbor->level() < cell->level(), ExcInternalError());
			  
				    const std::pair<unsigned int, unsigned int> faceno_subfaceno=
					cell->neighbor_of_coarser_neighbor(face_no);
				    const unsigned int neighbor_face_no=faceno_subfaceno.first,
					neighbor_subface_no=faceno_subfaceno.second;
			  
				    Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
										 neighbor_subface_no)
					    == cell,
					    ExcInternalError());
			  
				    cell_neighbor_m = 0;
				    cell_m = 0;
				    local.assemble_face_term(cell, face_no,
							     neighbor, neighbor_face_no, 
							     neighbor_subface_no,
							     cell_m, cell_neighbor_m);
			  
				    if(level<0)
					neighbor->get_dof_indices (neighbor_dofs);
				    else 
					neighbor->get_mg_dof_indices (neighbor_dofs);
			  
				    hang.distribute_local_to_global(cell_m, dofs, m);
				    // TBA Distribute neighbor dofs
				    // hang.distribute_local_to_global(cell_m, neighbor_dofs, dofs, m);
				    //

				} // The neighbor is coarser
			    } // else : if face has children
			} // skipped if no flags for face terms
		    } // faces loop

		} // Skipped if no flags on either faces or boundary 

#else
    //1D face terms
#endif 
	    } // cells loop
    }


      /** MG assembly routine for rhs vectors. Hanging node
     * constraints are taken into account. However face terms are NOT
     * correctly handled, that is, this method is currently not
     * capable of handling fluxex.
     
     Different levels can be assembled. If the level is -1, then the
     active level is used.
     */
    template <int dim, typename VECTOR>
    void assemble_rhs(MGDoFHandler<dim> &mg_dof_handler, ConstraintMatrix &hang,
		      int level, VECTOR& rhs, LocalAssembleBase<dim, MGDoFHandler<dim> > &local,
		      const unsigned int this_mpi_process = 0) 
    {
	// Cell iterators.
	typename MGDoFHandler<dim>::active_cell_iterator cell, endc;
	if(level < 0) { 
	    cell = mg_dof_handler.begin_active();
	    endc = mg_dof_handler.end();
	} else {
	    cell = mg_dof_handler.begin_active(level);	
	    endc = mg_dof_handler.end(level);
	}
      
	unsigned int dofs_per_cell = 
	    mg_dof_handler.get_fe().dofs_per_cell;
      
	// dofs and neighbor dofs.
	vector<unsigned int> dofs(dofs_per_cell);
      
	Vector<double> cell_rhs(dofs_per_cell);
      
	for (; cell!=endc; ++cell)
	    if( cell->subdomain_id() == this_mpi_process )
	    {
		// Update indices
		if(level<0)
		    cell->get_dof_indices(dofs);
		else
		    cell->get_mg_dof_indices (dofs);
	  
		// Assemble cell term if needed
		if( local.flags & assemble_rhs_cell ) {
		    // Reset the cell term
		    cell_rhs = 0;
		    local.assemble_rhs_term(cell, cell_rhs);
		    // Distribute the local rhs to the global matrix.
		    hang.distribute_local_to_global(cell_rhs, dofs, rhs);
		}

		// Assemble boundary or face terms if needed
		if( local.flags & (assemble_rhs_boundary/* |assemble_rhs_face */) )
		    for (unsigned int face_no=0; 
			 face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
		    {
			typename MGDoFHandler<dim>::face_iterator face=cell->face(face_no);
			if (face->at_boundary())
			{ 
			    cell_rhs = 0;
			    local.assemble_rhs_boundary_term(cell, face_no, cell_rhs);
			    // Distribute the local rhs to the global matrix.
			    hang.distribute_local_to_global(cell_rhs, dofs, rhs);
			}
		    } // face loop	  
	    } // Cell loop
    }
    
    /** DG-Like assembly routine for matrices and vectors. Hanging
     * node constraints are not taken into account, and one should do
     * this outside. However face terms are correctly handled, and
     * this method is capable of assemblying fully DG objects. */
    template <int dim, typename DH, typename MATRIX, typename VEC>
      void assemble(DH &dof_handler, MATRIX& system_matrix, VEC& system_rhs,
		  LocalAssembleBase<dim, DH> &local, 
		  const unsigned int this_mpi_process = 0) 
    {
	const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
	std::vector<unsigned int> dofs (dofs_per_cell);
	std::vector<unsigned int> dofs_neighbor (dofs_per_cell);
    
	// Now we create the cell matrices and vectors. Here we need two
	// cell matrices, both for face terms that include test functions
	// <code>vi</code> (internal shape functions, i.e. shape functions
	// of the current cell). 
	// 
	// i stands for internal, e stands for external
	FullMatrix<double> ui_vi_matrix (dofs_per_cell, dofs_per_cell);
	FullMatrix<double> ui_ve_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double>     local_rhs(dofs_per_cell);

	// Furthermore we need some cell iterators.
	typename DH::active_cell_iterator
	    cell = dof_handler.begin_active(),
	    endc = dof_handler.end();

	// Now we start the loop over all
	// active cells.
	for (;cell!=endc; ++cell) 
	    if( cell->subdomain_id() == this_mpi_process )
	    {
		// On each cell we need to reset the <code>ui_vi_matrix</code>
		// and <code>cell_vector</code> to zero, before assembling the
		// cell terms.
		ui_vi_matrix = 0;
		local_rhs = 0;
      
		// and call the function that assembles the cell terms. The
		// first argument is the <code>FEValues</code> that was
		// previously reinit'ed on the current cell.
		if( local.flags & assemble_cell ) 
		    local.assemble_cell_term(cell, ui_vi_matrix);
		
		if( local.flags & assemble_rhs_cell )
		    local.assemble_rhs_term(cell, local_rhs);
	  
		// As in previous examples the vector `dofs' includes the
		// dof_indices.
		cell->get_dof_indices (dofs);
	  
		// This is the start of the nested loop over all faces. 
		if( local.flags & (assemble_boundary|assemble_face) ) {
		    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
		    {
			// First we set the face iterator
			typename DH::face_iterator face=cell->face(face_no);
	  
			// and clear the <code>ui_ve_matrix</code> on each face.
			ui_ve_matrix = 0;
		  
			// Now we distinguish the four different cases in the
			// ordering mentioned above. We start with faces belonging
			// to the boundary of the domain.
			if (face->at_boundary())
			{
			    // and assemble the corresponding face terms.
			    if(local.flags & assemble_boundary) 
				local.assemble_boundary_term(cell, face_no,
							     ui_vi_matrix);
			    if(local.flags & assemble_rhs_boundary)
				local.assemble_rhs_boundary_term(cell, face_no,
								 local_rhs);
			}
			else if(local.flags & assemble_face) 
			{
			    // Now we are not on the boundary of the domain,
			    // therefore there must exist a neighboring cell.
			    typename DH::cell_iterator neighbor=
				cell->neighbor(face_no);

			    // We proceed with the second and most complicated case:
			    // the neighboring cell is more refined than the current
			    // cell. As in deal.II neighboring cells are restricted
			    // to have a level difference of not more than one, the
			    // neighboring cell is known to be at most ONCE more
			    // refined than the current cell. Furthermore also the
			    // face is more refined, i.e. it has children. Here we
			    // note that the following part of code will not work
			    // for <code>dim==1</code>.
			    if (face->has_children())
			    {
				// First we store which number the current cell has
				// in the list of neighbors of the neighboring
				// cell. Hence,
				// neighbor-@>neighbor(neighbor_face_no) equals the
				// current cell <code>cell</code>.
				const unsigned int neighbor_face_no=
				    cell->neighbor_of_neighbor(face_no);
		  
		  
				// We loop over subfaces
				for (unsigned int subface_no=0;
				     subface_no<face->n_children(); ++subface_no)
				{
				    // and set the cell iterator
				    // <code>neighbor_child</code> to the cell
				    // placed `behind' the current subface.
				    typename DH::active_cell_iterator
					neighbor_child
					= cell->neighbor_child_on_subface (face_no, subface_no);
		      
				    // As these are quite complicated indirections
				    // which one does not usually get right at first
				    // attempt we check for the internal
				    // consistency.
				    Assert (neighbor_child->face(neighbor_face_no) == face->child(subface_no),
					    ExcInternalError());
				    Assert (!neighbor_child->has_children(), ExcInternalError());

				    // We need to reset the
				    // <code>ui_ve_matrix</code> on each subface
				    // because on each subface the <code>un</code>
				    // belong to different neighboring cells.
				    ui_ve_matrix = 0;
		      
				    // As already mentioned above for the current
				    // case (case 2) we employ the
				    // <code>FESubfaceValues</code> of the current
				    // cell (here reinited for the current cell,
				    // face and subface) and we employ the
				    // FEFaceValues of the neighboring child cell.
				    local.assemble_face_term(cell, face_no, subface_no,
							     neighbor_child, neighbor_face_no,
							     ui_vi_matrix,
							     ui_ve_matrix);
		      
				    // Then we get the dof indices of the
				    // neighbor_child cell
				    neighbor_child->get_dof_indices (dofs_neighbor);
		      						
				    // and distribute <code>ui_ve_matrix</code> to
				    // the system_matrix
				    for (unsigned int i=0; i<dofs_per_cell; ++i)
					for (unsigned int k=0; k<dofs_per_cell; ++k)
					    system_matrix.add(dofs_neighbor[i], dofs[k],
							      ui_ve_matrix(i,k));
				}
				// End of <code>if (face-@>has_children())</code>
			    }
			    else
			    {
				// We proceed with case 3, i.e. neighboring cell is
				// of the same refinement level as the current cell.
				if (neighbor->level() == cell->level()) 
				{
				    // Like before we store which number the current
				    // cell has in the list of neighbors of the
				    // neighboring cell.
				    const unsigned int neighbor_face_no=cell->neighbor_of_neighbor(face_no);

				    local.assemble_face_term(cell, face_no,
							     neighbor, neighbor_face_no,
							     ui_vi_matrix,
							     ui_ve_matrix);
				    // End of <code>if (neighbor-@>level() ==
				    // cell-@>level())</code>
				}
				else
				{
				    // Finally we consider case 4. When the
				    // neighboring cell is not finer and not of the
				    // same refinement level as the current cell it
				    // must be coarser.
				    Assert(neighbor->level() < cell->level(), ExcInternalError());

				    // Find out the how many'th face_no and
				    // subface_no the current face is w.r.t. the
				    // neighboring cell.
				    const std::pair<unsigned int, unsigned int> faceno_subfaceno=
					cell->neighbor_of_coarser_neighbor(face_no);
				    const unsigned int neighbor_face_no=faceno_subfaceno.first,
					neighbor_subface_no=faceno_subfaceno.second;

				    Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
										 neighbor_subface_no)
					    == cell,
					    ExcInternalError());

				    local.assemble_face_term(cell, face_no,
							     neighbor, 
							     neighbor_face_no,
							     neighbor_subface_no,
							     ui_vi_matrix,
							     ui_ve_matrix);
				}

				// Now we get the dof indices of the
				// <code>neighbor</code> cell,
				neighbor->get_dof_indices (dofs_neighbor);
		 						
				// and distribute the
				// <code>ui_ve_matrix</code>.
				for (unsigned int i=0; i<dofs_per_cell; ++i)
				    for (unsigned int k=0; k<dofs_per_cell; ++k)
					system_matrix.add(dofs_neighbor[i], dofs[k],
							  ui_ve_matrix(i,k));
			    }
			    // End of <code>face not at boundary</code>:
			}
			// End of loop over all faces:
		    }
		}
	    
		// Finally we distribute the
		// <code>ui_vi_matrix</code> 
		for (unsigned int i=0; i<dofs_per_cell; ++i) {
		    for (unsigned int j=0; j<dofs_per_cell; ++j)
			system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i,j));
		    system_rhs(dofs[i]) += local_rhs(i);
		}
	    }
    }
    
    /** Assembly routine for rhs. Hanging node constraints are not
     * taken into account, and one should do this outside. */
    template <int dim, typename DH, typename VEC>
      void assemble_rhs(DH &dof_handler, VEC& system_rhs,
		      LocalAssembleBase<dim, DH> &local, 
		      const unsigned int this_mpi_process = 0) 
    {
	const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
	std::vector<unsigned int> dofs (dofs_per_cell);
	
	Vector<double>     local_rhs(dofs_per_cell);

	// Furthermore we need some cell iterators.
	typename DH::active_cell_iterator
	    cell = dof_handler.begin_active(),
	    endc = dof_handler.end();

	// Now we start the loop over all
	// active cells.
	for (;cell!=endc; ++cell) 
	    if( cell->subdomain_id() == this_mpi_process )
	    {
		// On each cell we need to reset the <code>ui_vi_matrix</code>
		// and <code>cell_vector</code> to zero, before assembling the
		// cell terms.
		local_rhs = 0;
      
		// and call the function that assembles the cell terms. The
		// first argument is the <code>FEValues</code> that was
		// previously reinit'ed on the current cell.
		if( local.flags & assemble_rhs_cell )
		    local.assemble_rhs_term(cell, local_rhs);
	  
		// As in previous examples the vector `dofs' includes the
		// dof_indices.
		cell->get_dof_indices (dofs);
	  
		// This is the start of the nested loop over all faces. 
		if( local.flags & assemble_rhs_boundary ) {
		    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
		    {
			// First we set the face iterator
			typename DH::face_iterator face=cell->face(face_no);
		  
			// Now we distinguish the four different cases in the
			// ordering mentioned above. We start with faces belonging
			// to the boundary of the domain.
			if (face->at_boundary())
			{
			    if(local.flags)
				local.assemble_rhs_boundary_term(cell, face_no,
								 local_rhs);
			}
		    }
		}
	    
		// Finally we distribute the
		// <code>rhs</code> 
		for (unsigned int i=0; i<dofs_per_cell; ++i) {
		    system_rhs(dofs[i]) += local_rhs(i);
		}
	    }
    }
}

#endif
