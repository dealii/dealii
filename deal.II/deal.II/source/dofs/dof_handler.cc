//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/memory_consumption.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_levels.h>
#include <dofs/dof_faces.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_levels.h>
#include <grid/tria.h>
#include <base/geometry_info.h>
#include <fe/fe.h>

#include <set>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN


//TODO[WB]: do not use a plain pointer for DoFHandler::faces, but rather an
//auto_ptr or some such thing. alternatively, why not use the DoFFaces object
//right away?

template<int dim, int spacedim>
const unsigned int DoFHandler<dim,spacedim>::dimension;

template<int dim, int spacedim>
const unsigned int DoFHandler<dim,spacedim>::space_dimension;

template <int dim, int spacedim>
const unsigned int DoFHandler<dim,spacedim>::invalid_dof_index;

template <int dim, int spacedim>
const unsigned int DoFHandler<dim,spacedim>::default_fe_index;


// reference the invalid_dof_index variable explicitely to work around
// a bug in the icc8 compiler
namespace internal
{
  template <int dim, int spacedim>
  const unsigned int * dummy ()
  {
    return &dealii::DoFHandler<dim,spacedim>::invalid_dof_index;
  }

  template const unsigned int * dummy<deal_II_dimension,deal_II_dimension> ();
}



namespace internal
{
  namespace DoFHandler
  {
				     // access class
				     // dealii::DoFHandler instead of
				     // namespace internal::DoFHandler
    using dealii::DoFHandler;

/**
 * A class with the same purpose as the similarly named class of the
 * Triangulation class. See there for more information.
 */
    struct Implementation
    {
					 /**
					  * Distribute dofs on the given cell,
					  * with new dofs starting with index
					  * @p next_free_dof. Return the next
					  * unused index number. The finite
					  * element used is the one given to
					  * @p distribute_dofs, which is copied
					  * to @p selected_fe.
					  *
					  * This function is excluded from the
					  * @p distribute_dofs function since
					  * it can not be implemented dimension
					  * independent.
					  */
	template <int spacedim>
	static
	unsigned int
	distribute_dofs_on_cell (const DoFHandler<1,spacedim> &dof_handler,
				 typename DoFHandler<1,spacedim>::active_cell_iterator &cell,
				 unsigned int          next_free_dof)
	  {

					     // distribute dofs of vertices
	    for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
	      {
		typename DoFHandler<1,spacedim>::cell_iterator
		  neighbor = cell->neighbor(v);

		if (neighbor.state() == IteratorState::valid)
		  {
						     // find true neighbor; may be its
						     // a child of @p{neighbor}
		    while (neighbor->has_children())
		      neighbor = neighbor->child(v==0 ? 1 : 0);

						     // has neighbor already been processed?
		    if (neighbor->user_flag_set())
						       // copy dofs
		      {
			if (v==0)
			  for (unsigned int d=0;
			       d<dof_handler.selected_fe->dofs_per_vertex; ++d)
			    cell->set_vertex_dof_index (0, d,
							neighbor->vertex_dof_index (1, d));
			else
			  for (unsigned int d=0;
			       d<dof_handler.selected_fe->dofs_per_vertex; ++d)
			    cell->set_vertex_dof_index (1, d,
							neighbor->vertex_dof_index (0, d));

							 // next neighbor
			continue;
		      }
		  }

						 // otherwise: create dofs newly
		for (unsigned int d=0;
		     d<dof_handler.selected_fe->dofs_per_vertex; ++d)
		  cell->set_vertex_dof_index (v, d, next_free_dof++);
	      }

					     // dofs of line
	    for (unsigned int d=0;
		 d<dof_handler.selected_fe->dofs_per_line; ++d)
	      cell->set_dof_index (d, next_free_dof++);

					     // note that this cell has been
					     // processed
	    cell->set_user_flag ();

	    return next_free_dof;
	  }



	template <int spacedim>
	static
	unsigned int
	distribute_dofs_on_cell (const DoFHandler<2,spacedim> &dof_handler,
				 typename DoFHandler<2,spacedim>::active_cell_iterator &cell,
				 unsigned int          next_free_dof)
	  {
	    if (dof_handler.selected_fe->dofs_per_vertex > 0)
					       // number dofs on vertices
	      for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
						 // check whether dofs for this
						 // vertex have been distributed
						 // (only check the first dof)
		if (cell->vertex_dof_index(vertex, 0) == DoFHandler<2,spacedim>::invalid_dof_index)
		  for (unsigned int d=0; d<dof_handler.selected_fe->dofs_per_vertex; ++d)
		    cell->set_vertex_dof_index (vertex, d, next_free_dof++);

					     // for the four sides
	    if (dof_handler.selected_fe->dofs_per_line > 0)
	      for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
		{
		  typename DoFHandler<2,spacedim>::line_iterator
		    line = cell->line(side);

						   // distribute dofs if necessary:
						   // check whether line dof is already
						   // numbered (check only first dof)
		  if (line->dof_index(0) == DoFHandler<2,spacedim>::invalid_dof_index)
						     // if not: distribute dofs
		    for (unsigned int d=0; d<dof_handler.selected_fe->dofs_per_line; ++d)
		      line->set_dof_index (d, next_free_dof++);
		}


					     // dofs of quad
	    if (dof_handler.selected_fe->dofs_per_quad > 0)
	      for (unsigned int d=0; d<dof_handler.selected_fe->dofs_per_quad; ++d)
		cell->set_dof_index (d, next_free_dof++);


					     // note that this cell has been processed
	    cell->set_user_flag ();

	    return next_free_dof;
	  }


	template <int spacedim>
	static
	unsigned int
	distribute_dofs_on_cell (const DoFHandler<3,spacedim> &dof_handler,
				 typename DoFHandler<3,spacedim>::active_cell_iterator &cell,
				 unsigned int          next_free_dof)
	  {
	    if (dof_handler.selected_fe->dofs_per_vertex > 0)
					       // number dofs on vertices
	      for (unsigned int vertex=0; vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
						 // check whether dofs for this
						 // vertex have been distributed
						 // (only check the first dof)
		if (cell->vertex_dof_index(vertex, 0) == DoFHandler<3,spacedim>::invalid_dof_index)
		  for (unsigned int d=0; d<dof_handler.selected_fe->dofs_per_vertex; ++d)
		    cell->set_vertex_dof_index (vertex, d, next_free_dof++);

					     // for the lines
	    if (dof_handler.selected_fe->dofs_per_line > 0)
	      for (unsigned int l=0; l<GeometryInfo<3>::lines_per_cell; ++l)
		{
		  typename DoFHandler<3,spacedim>::line_iterator
		    line = cell->line(l);

						   // distribute dofs if necessary:
						   // check whether line dof is already
						   // numbered (check only first dof)
		  if (line->dof_index(0) == DoFHandler<3,spacedim>::invalid_dof_index)
						     // if not: distribute dofs
		    for (unsigned int d=0; d<dof_handler.selected_fe->dofs_per_line; ++d)
		      line->set_dof_index (d, next_free_dof++);
		}

					     // for the quads
	    if (dof_handler.selected_fe->dofs_per_quad > 0)
	      for (unsigned int q=0; q<GeometryInfo<3>::quads_per_cell; ++q)
		{
		  typename DoFHandler<3,spacedim>::quad_iterator
		    quad = cell->quad(q);

						   // distribute dofs if necessary:
						   // check whether quad dof is already
						   // numbered (check only first dof)
		  if (quad->dof_index(0) == DoFHandler<3,spacedim>::invalid_dof_index)
						     // if not: distribute dofs
		    for (unsigned int d=0; d<dof_handler.selected_fe->dofs_per_quad; ++d)
		      quad->set_dof_index (d, next_free_dof++);
		}


					     // dofs of hex
	    if (dof_handler.selected_fe->dofs_per_hex > 0)
	      for (unsigned int d=0; d<dof_handler.selected_fe->dofs_per_hex; ++d)
		cell->set_dof_index (d, next_free_dof++);


					     // note that this cell has been
					     // processed
	    cell->set_user_flag ();

	    return next_free_dof;
	  }


					 /**
					  * Implementation of the general template
					  * of same name.
					  */
	template <int spacedim>
	static
	void renumber_dofs (const std::vector<unsigned int> &new_numbers,
			    DoFHandler<1,spacedim>          &dof_handler)
	  {
					     // note that we can not use cell
					     // iterators in this function since
					     // then we would renumber the dofs on
					     // the interface of two cells more
					     // than once. Anyway, this way it's
					     // not only more correct but also
					     // faster; note, however, that dof
					     // numbers may be invalid_dof_index,
					     // namely when the appropriate
					     // vertex/line/etc is unused
	    for (std::vector<unsigned int>::iterator
		   i=dof_handler.vertex_dofs.begin();
		 i!=dof_handler.vertex_dofs.end(); ++i)
	      if (*i != DoFHandler<1,spacedim>::invalid_dof_index)
		*i = new_numbers[*i];
	      else
						 // if index is
						 // invalid_dof_index:
						 // check if this one
						 // really is unused
		Assert (dof_handler.tria
			->vertex_used((i-dof_handler.vertex_dofs.begin()) /
				      dof_handler.selected_fe->dofs_per_vertex)
			== false,
			ExcInternalError ());

	    for (unsigned int level=0; level<dof_handler.levels.size(); ++level)
	      for (std::vector<unsigned int>::iterator
		     i=dof_handler.levels[level]->lines.dofs.begin();
		   i!=dof_handler.levels[level]->lines.dofs.end(); ++i)
		if (*i != DoFHandler<1,spacedim>::invalid_dof_index)
		  *i = new_numbers[*i];
	  }



	template <int spacedim>
	static
	void renumber_dofs (const std::vector<unsigned int> &new_numbers,
			    DoFHandler<2,spacedim>          &dof_handler)
	  {
					     // note that we can not use cell
					     // iterators in this function since
					     // then we would renumber the dofs on
					     // the interface of two cells more
					     // than once. Anyway, this way it's
					     // not only more correct but also
					     // faster; note, however, that dof
					     // numbers may be invalid_dof_index,
					     // namely when the appropriate
					     // vertex/line/etc is unused
	    for (std::vector<unsigned int>::iterator
		   i=dof_handler.vertex_dofs.begin();
		 i!=dof_handler.vertex_dofs.end(); ++i)
	      if (*i != DoFHandler<2,spacedim>::invalid_dof_index)
		*i = new_numbers[*i];
	      else
						 // if index is invalid_dof_index:
						 // check if this one really is
						 // unused
		Assert (dof_handler.tria
			->vertex_used((i-dof_handler.vertex_dofs.begin()) /
				      dof_handler.selected_fe->dofs_per_vertex)
			== false,
			ExcInternalError ());

	    for (std::vector<unsigned int>::iterator
		   i=dof_handler.faces->lines.dofs.begin();
		 i!=dof_handler.faces->lines.dofs.end(); ++i)
	      if (*i != DoFHandler<2,spacedim>::invalid_dof_index)
		*i = new_numbers[*i];

	    for (unsigned int level=0; level<dof_handler.levels.size(); ++level)
	      {
		for (std::vector<unsigned int>::iterator
		       i=dof_handler.levels[level]->quads.dofs.begin();
		     i!=dof_handler.levels[level]->quads.dofs.end(); ++i)
		  if (*i != DoFHandler<2,spacedim>::invalid_dof_index)
		    *i = new_numbers[*i];
	      }
	  }


	template <int spacedim>
	static
	void renumber_dofs (const std::vector<unsigned int> &new_numbers,
			    DoFHandler<3,spacedim>          &dof_handler)
	  {
					     // note that we can not use cell
					     // iterators in this function since
					     // then we would renumber the dofs on
					     // the interface of two cells more
					     // than once. Anyway, this way it's
					     // not only more correct but also
					     // faster; note, however, that dof
					     // numbers may be invalid_dof_index,
					     // namely when the appropriate
					     // vertex/line/etc is unused
	    for (std::vector<unsigned int>::iterator
		   i=dof_handler.vertex_dofs.begin();
		 i!=dof_handler.vertex_dofs.end(); ++i)
	      if (*i != DoFHandler<3,spacedim>::invalid_dof_index)
		*i = new_numbers[*i];
	      else
						 // if index is invalid_dof_index:
						 // check if this one really is
						 // unused
		Assert (dof_handler.tria
			->vertex_used((i-dof_handler.vertex_dofs.begin()) /
				      dof_handler.selected_fe->dofs_per_vertex)
			== false,
			ExcInternalError ());

	    for (std::vector<unsigned int>::iterator
		   i=dof_handler.faces->lines.dofs.begin();
		 i!=dof_handler.faces->lines.dofs.end(); ++i)
	      if (*i != DoFHandler<3,spacedim>::invalid_dof_index)
		*i = new_numbers[*i];
	    for (std::vector<unsigned int>::iterator
		   i=dof_handler.faces->quads.dofs.begin();
		 i!=dof_handler.faces->quads.dofs.end(); ++i)
	      if (*i != DoFHandler<3,spacedim>::invalid_dof_index)
		*i = new_numbers[*i];

	    for (unsigned int level=0; level<dof_handler.levels.size(); ++level)
	      {
		for (std::vector<unsigned int>::iterator
		       i=dof_handler.levels[level]->hexes.dofs.begin();
		     i!=dof_handler.levels[level]->hexes.dofs.end(); ++i)
		  if (*i != DoFHandler<3,spacedim>::invalid_dof_index)
		    *i = new_numbers[*i];
	      }
	  }



					 /**
					  * Implement the function of same name in
					  * the mother class.
					  */
	template <int spacedim>
	static
	unsigned int
	max_couplings_between_dofs (const DoFHandler<1,spacedim> &dof_handler)
	  {
	    return std::min(3*dof_handler.selected_fe->dofs_per_vertex +
			    2*dof_handler.selected_fe->dofs_per_line,
			    dof_handler.n_dofs());
	  }



	template <int spacedim>
	static
	unsigned int
	max_couplings_between_dofs (const DoFHandler<2,spacedim> &dof_handler)
	  {

					     // get these numbers by drawing pictures
					     // and counting...
					     // example:
					     //   |     |     |
					     // --x-----x--x--X--
					     //   |     |  |  |
					     //   |     x--x--x
					     //   |     |  |  |
					     // --x--x--*--x--x--
					     //   |  |  |     |
					     //   x--x--x     |
					     //   |  |  |     |
					     // --X--x--x-----x--
					     //   |     |     |
					     // x = vertices connected with center vertex *;
					     //   = total of 19
					     // (the X vertices are connected with * if
					     // the vertices adjacent to X are hanging
					     // nodes)
					     // count lines -> 28 (don't forget to count
					     // mother and children separately!)
	    unsigned int max_couplings;
	    switch (dof_handler.tria->max_adjacent_cells())
	      {
		case 4:
		      max_couplings=19*dof_handler.selected_fe->dofs_per_vertex +
				    28*dof_handler.selected_fe->dofs_per_line +
				    8*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 5:
		      max_couplings=21*dof_handler.selected_fe->dofs_per_vertex +
				    31*dof_handler.selected_fe->dofs_per_line +
				    9*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 6:
		      max_couplings=28*dof_handler.selected_fe->dofs_per_vertex +
				    42*dof_handler.selected_fe->dofs_per_line +
				    12*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 7:
		      max_couplings=30*dof_handler.selected_fe->dofs_per_vertex +
				    45*dof_handler.selected_fe->dofs_per_line +
				    13*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 8:
		      max_couplings=37*dof_handler.selected_fe->dofs_per_vertex +
				    56*dof_handler.selected_fe->dofs_per_line +
				    16*dof_handler.selected_fe->dofs_per_quad;
		      break;

						       // the following
						       // numbers are not
						       // based on actual
						       // counting but by
						       // extrapolating the
						       // number sequences
						       // from the previous
						       // ones (for example,
						       // for dofs_per_vertex,
						       // the sequence above
						       // is 19, 21, 28, 30,
						       // 37, and is continued
						       // as follows):
		case 9:
		      max_couplings=39*dof_handler.selected_fe->dofs_per_vertex +
				    59*dof_handler.selected_fe->dofs_per_line +
				    17*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 10:
		      max_couplings=46*dof_handler.selected_fe->dofs_per_vertex +
				    70*dof_handler.selected_fe->dofs_per_line +
				    20*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 11:
		      max_couplings=48*dof_handler.selected_fe->dofs_per_vertex +
				    73*dof_handler.selected_fe->dofs_per_line +
				    21*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 12:
		      max_couplings=55*dof_handler.selected_fe->dofs_per_vertex +
				    84*dof_handler.selected_fe->dofs_per_line +
				    24*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 13:
		      max_couplings=57*dof_handler.selected_fe->dofs_per_vertex +
				    87*dof_handler.selected_fe->dofs_per_line +
				    25*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 14:
		      max_couplings=63*dof_handler.selected_fe->dofs_per_vertex +
				    98*dof_handler.selected_fe->dofs_per_line +
				    28*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 15:
		      max_couplings=65*dof_handler.selected_fe->dofs_per_vertex +
				    103*dof_handler.selected_fe->dofs_per_line +
				    29*dof_handler.selected_fe->dofs_per_quad;
		      break;
		case 16:
		      max_couplings=72*dof_handler.selected_fe->dofs_per_vertex +
				    114*dof_handler.selected_fe->dofs_per_line +
				    32*dof_handler.selected_fe->dofs_per_quad;
		      break;

		default:
		      Assert (false, ExcNotImplemented());
		      max_couplings=0;
	      }
	    return std::min(max_couplings,dof_handler.n_dofs());
	  }


	template <int spacedim>
	static
	unsigned int
	max_couplings_between_dofs (const DoFHandler<3,spacedim> &dof_handler)
	  {
//TODO:[?] Invent significantly better estimates than the ones in this function

					     // doing the same thing here is a
					     // rather complicated thing, compared
					     // to the 2d case, since it is hard
					     // to draw pictures with several
					     // refined hexahedra :-) so I
					     // presently only give a coarse
					     // estimate for the case that at most
					     // 8 hexes meet at each vertex
					     //
					     // can anyone give better estimate
					     // here?
	    const unsigned int max_adjacent_cells
	      = dof_handler.tria->max_adjacent_cells();

	    unsigned int max_couplings;
	    if (max_adjacent_cells <= 8)
	      max_couplings=7*7*7*dof_handler.selected_fe->dofs_per_vertex +
			    7*6*7*3*dof_handler.selected_fe->dofs_per_line +
			    9*4*7*3*dof_handler.selected_fe->dofs_per_quad +
			    27*dof_handler.selected_fe->dofs_per_hex;
	    else
	      {
		Assert (false, ExcNotImplemented());
		max_couplings=0;
	      }

	    return std::min(max_couplings,dof_handler.n_dofs());
	  }


					 /**
					  * Reserve enough space in the
					  * <tt>levels[]</tt> objects to store the
					  * numbers of the degrees of freedom
					  * needed for the given element. The
					  * given element is that one which
					  * was selected when calling
					  * @p distribute_dofs the last time.
					  */
	template <int spacedim>
	static
	void reserve_space (DoFHandler<1,spacedim> &dof_handler)
	  {
	    dof_handler.vertex_dofs
	      .resize(dof_handler.tria->n_vertices() *
		      dof_handler.selected_fe->dofs_per_vertex,
		      DoFHandler<1,spacedim>::invalid_dof_index);

	    for (unsigned int i=0; i<dof_handler.tria->n_levels(); ++i)
	      {
		dof_handler.levels
		  .push_back (new internal::DoFHandler::DoFLevel<1>);

		dof_handler.levels.back()->lines.dofs
		  .resize (dof_handler.tria->n_raw_lines(i) *
			   dof_handler.selected_fe->dofs_per_line,
			   DoFHandler<1,spacedim>::invalid_dof_index);

		dof_handler.levels.back()->cell_dof_indices_cache
		  .resize (dof_handler.tria->n_raw_lines(i) *
			   dof_handler.selected_fe->dofs_per_cell,
			   DoFHandler<1,spacedim>::invalid_dof_index);
	      }
	  }


	template <int spacedim>
	static
	void reserve_space (DoFHandler<2,spacedim> &dof_handler)
	  {
	    dof_handler.vertex_dofs
	      .resize(dof_handler.tria->n_vertices() *
		      dof_handler.selected_fe->dofs_per_vertex,
		      DoFHandler<2,spacedim>::invalid_dof_index);

	    for (unsigned int i=0; i<dof_handler.tria->n_levels(); ++i)
	      {
		dof_handler.levels.push_back (new internal::DoFHandler::DoFLevel<2>);

		dof_handler.levels.back()->quads.dofs
		  .resize (dof_handler.tria->n_raw_quads(i) *
			   dof_handler.selected_fe->dofs_per_quad,
			   DoFHandler<2,spacedim>::invalid_dof_index);

		dof_handler.levels.back()->cell_dof_indices_cache
		  .resize (dof_handler.tria->n_raw_quads(i) *
			   dof_handler.selected_fe->dofs_per_cell,
			   DoFHandler<2,spacedim>::invalid_dof_index);
	      }

	    dof_handler.faces = new internal::DoFHandler::DoFFaces<2>;
	    dof_handler.faces->lines.dofs
	      .resize (dof_handler.tria->n_raw_lines() *
		       dof_handler.selected_fe->dofs_per_line,
		       DoFHandler<2,spacedim>::invalid_dof_index);
	  }


	template <int spacedim>
	static
	void reserve_space (DoFHandler<3,spacedim> &dof_handler)
	  {
	    dof_handler.vertex_dofs
	      .resize(dof_handler.tria->n_vertices() *
		      dof_handler.selected_fe->dofs_per_vertex,
		      DoFHandler<3,spacedim>::invalid_dof_index);

	    for (unsigned int i=0; i<dof_handler.tria->n_levels(); ++i)
	      {
		dof_handler.levels.push_back (new internal::DoFHandler::DoFLevel<3>);

		dof_handler.levels.back()->hexes.dofs
		  .resize (dof_handler.tria->n_raw_hexs(i) *
			   dof_handler.selected_fe->dofs_per_hex,
			   DoFHandler<3,spacedim>::invalid_dof_index);

		dof_handler.levels.back()->cell_dof_indices_cache
		  .resize (dof_handler.tria->n_raw_hexs(i) *
			   dof_handler.selected_fe->dofs_per_cell,
			   DoFHandler<3,spacedim>::invalid_dof_index);
	      }
	    dof_handler.faces = new internal::DoFHandler::DoFFaces<3>;

	    dof_handler.faces->lines.dofs
	      .resize (dof_handler.tria->n_raw_lines() *
		       dof_handler.selected_fe->dofs_per_line,
		       DoFHandler<3,spacedim>::invalid_dof_index);
	    dof_handler.faces->quads.dofs
	      .resize (dof_handler.tria->n_raw_quads() *
		       dof_handler.selected_fe->dofs_per_quad,
		       DoFHandler<3,spacedim>::invalid_dof_index);
	  }
    };
  }
}



template<int dim, int spacedim>
DoFHandler<dim,spacedim>::DoFHandler (const Triangulation<dim,spacedim> &tria)
		:
		tria(&tria, typeid(*this).name()),
		selected_fe(0, typeid(*this).name()),
		faces(NULL),
		used_dofs (0)
{}


template <int dim, int spacedim>
DoFHandler<dim,spacedim>::~DoFHandler ()
{
				   // release allocated memory
  clear ();
}


/*------------------------ Cell iterator functions ------------------------*/


template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_cell_iterator
DoFHandler<dim,spacedim>::begin_raw (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return begin_raw_line (level);
      case 2:
	    return begin_raw_quad (level);
      case 3:
	    return begin_raw_hex (level);
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_cell_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::cell_iterator
DoFHandler<dim,spacedim>::begin (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return begin_line (level);
      case 2:
	    return begin_quad (level);
      case 3:
	    return begin_hex (level);
      default:
	    Assert (false, ExcImpossibleInDim(dim));
	    return cell_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::active_cell_iterator
DoFHandler<dim,spacedim>::begin_active (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return begin_active_line (level);
      case 2:
	    return begin_active_quad (level);
      case 3:
	    return begin_active_hex (level);
      default:
	    Assert (false, ExcNotImplemented());
	    return active_cell_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_cell_iterator
DoFHandler<dim,spacedim>::last_raw () const
{
  switch (dim)
    {
      case 1:
	    return last_raw_line ();
      case 2:
	    return last_raw_quad ();
      case 3:
	    return last_raw_hex ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_cell_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_cell_iterator
DoFHandler<dim,spacedim>::last_raw (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return last_raw_line (level);
      case 2:
	    return last_raw_quad (level);
      case 3:
	    return last_raw_hex (level);
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_cell_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::cell_iterator
DoFHandler<dim,spacedim>::last () const
{
  switch (dim)
    {
      case 1:
	    return last_line ();
      case 2:
	    return last_quad ();
      case 3:
	    return last_hex ();
      default:
	    Assert (false, ExcNotImplemented());
	    return cell_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::cell_iterator
DoFHandler<dim,spacedim>::last (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return last_line (level);
      case 2:
	    return last_quad (level);
      case 3:
	    return last_hex (level);
      default:
	    Assert (false, ExcNotImplemented());
	    return cell_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::active_cell_iterator
DoFHandler<dim,spacedim>::last_active () const
{
  switch (dim)
    {
      case 1:
	    return last_active_line ();
      case 2:
	    return last_active_quad ();
      case 3:
	    return last_active_hex ();
      default:
	    Assert (false, ExcNotImplemented());
	    return active_cell_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::active_cell_iterator
DoFHandler<dim,spacedim>::last_active (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return last_active_line (level);
      case 2:
	    return last_active_quad (level);
      case 3:
	    return last_active_hex (level);
      default:
	    Assert (false, ExcNotImplemented());
	    return active_cell_iterator();
    }
}


template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_cell_iterator
DoFHandler<dim,spacedim>::end () const
{
  switch (dim)
    {
      case 1:
	    return end_line();
      case 2:
	    return end_quad();
      case 3:
	    return end_hex();
      default:
	    Assert (false, ExcImpossibleInDim(dim));
	    return raw_cell_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_cell_iterator
DoFHandler<dim, spacedim>::end_raw (const unsigned int level) const
{
  return (level == tria->n_levels()-1 ?
	  end() :
	  begin_raw (level+1));
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::cell_iterator
DoFHandler<dim, spacedim>::end (const unsigned int level) const
{
  return (level == tria->n_levels()-1 ?
	  cell_iterator(end()) :
	  begin (level+1));
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_cell_iterator
DoFHandler<dim, spacedim>::end_active (const unsigned int level) const
{
  return (level == tria->n_levels()-1 ?
	  active_cell_iterator(end()) :
	  begin_active (level+1));
}


/*------------------------ Face iterator functions ------------------------*/


template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_face_iterator
DoFHandler<dim,spacedim>::begin_raw_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return begin_raw_line ();
      case 3:
	    return begin_raw_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_face_iterator ();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::face_iterator
DoFHandler<dim,spacedim>::begin_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return begin_line ();
      case 3:
	    return begin_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return face_iterator ();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::active_face_iterator
DoFHandler<dim,spacedim>::begin_active_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return begin_active_line ();
      case 3:
	    return begin_active_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return active_face_iterator ();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_face_iterator
DoFHandler<dim, spacedim>::end_raw_face () const
{
  return end_face();
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_face_iterator
DoFHandler<dim,spacedim>::end_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return end_line ();
      case 3:
	    return end_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_face_iterator ();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_face_iterator
DoFHandler<dim, spacedim>::end_active_face () const
{
  return active_face_iterator(end_face());
}





template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_face_iterator
DoFHandler<dim,spacedim>::last_raw_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return last_raw_line ();
      case 3:
	    return last_raw_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_face_iterator ();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::face_iterator
DoFHandler<dim,spacedim>::last_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return last_line ();
      case 3:
	    return last_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_face_iterator ();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::active_face_iterator
DoFHandler<dim,spacedim>::last_active_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return last_active_line ();
      case 3:
	    return last_active_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_face_iterator ();
    }
}


/*------------------------ Line iterator functions ------------------------*/



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_line_iterator
DoFHandler<dim, spacedim>::begin_raw_line (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    Assert (level<tria->n_levels(), ExcInvalidLevel(level));

	    if (tria->n_raw_lines(level) == 0)
	      return end_line ();

	    return raw_line_iterator (tria,
				      level,
				      0,
				      this);

      default:
	    Assert (level == 0, ExcFacesHaveNoLevel());
	    return raw_line_iterator (tria,
				      0,
				      0,
				      this);
    }
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::line_iterator
DoFHandler<dim, spacedim>::begin_line (const unsigned int level) const
{
  				   // level is checked in begin_raw
  raw_line_iterator ri = begin_raw_line (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_line_iterator
DoFHandler<dim, spacedim>::begin_active_line (const unsigned int level) const
{
  				   // level is checked in begin_raw
  line_iterator i = begin_line (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_line_iterator
DoFHandler<dim, spacedim>::end_line () const
{
  return raw_line_iterator (tria,
			    -1,
			    -1,
			    this);
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_line_iterator
DoFHandler<dim, spacedim>::last_raw_line (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    Assert (level<tria->n_levels(), ExcInvalidLevel(level));
	    Assert (tria->n_raw_lines(level) != 0,
		    ExcEmptyLevel (level));

	    return raw_line_iterator (tria,
				      level,
				      tria->n_raw_lines(level)-1,
				      this);

      default:
	    Assert (level == 0, ExcFacesHaveNoLevel());
	    return raw_line_iterator (tria,
				      0,
				      tria->n_raw_lines()-1,
				      this);
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_line_iterator
DoFHandler<dim, spacedim>::last_raw_line () const
{
  if (dim == 1)
    return last_raw_line (tria->n_levels()-1);
  else
    return last_raw_line (0);
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::line_iterator
DoFHandler<dim, spacedim>::last_line (const unsigned int level) const
{
  				   // level is checked in last_raw
  raw_line_iterator ri = last_raw_line(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == IteratorState::valid)
    if (ri->used()==true)
      return ri;
  return ri;
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::line_iterator
DoFHandler<dim, spacedim>::last_line () const
{
  if (dim == 1)
    return last_line (tria->n_levels()-1);
  else
    return last_line (0);
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_line_iterator
DoFHandler<dim, spacedim>::last_active_line (const unsigned int level) const
{
				   // level is checked in last_raw
  line_iterator i = last_line(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == IteratorState::valid)
    if (i->has_children()==false)
      return i;
  return i;
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_line_iterator
DoFHandler<dim, spacedim>::last_active_line () const
{
  if (dim == 1)
    return last_active_line (tria->n_levels()-1);
  else
    return last_active_line (0);
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_line_iterator
DoFHandler<dim, spacedim>::end_raw_line (const unsigned int level) const
{
  Assert (dim == 1 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 1)
    return (level == tria->n_levels()-1 ?
	    end_line() :
	    begin_raw_line (level+1));
  else
    return end_line();
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::line_iterator
DoFHandler<dim, spacedim>::end_line (const unsigned int level) const
{
  Assert (dim == 1 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 1)
    return (level == tria->n_levels()-1 ?
	    line_iterator(end_line()) :
	    begin_line (level+1));
  else
    return line_iterator(end_line());
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_line_iterator
DoFHandler<dim, spacedim>::end_active_line (const unsigned int level) const
{
  Assert (dim == 1 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 1)
    return (level == tria->n_levels()-1 ?
	    active_line_iterator(end_line()) :
	    begin_active_line (level+1));
  else
    return active_line_iterator(end_line());
}



/*------------------------ Quad iterator functions ------------------------*/


template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_quad_iterator
DoFHandler<dim,spacedim>::begin_raw_quad (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_hex_iterator();
      case 2:
      {
	Assert (level<tria->n_levels(), ExcInvalidLevel(level));

	if (tria->n_raw_quads(level) == 0)
	  return end_quad();

	return raw_quad_iterator (tria,
				  level,
				  0,
				  this);
      }

      case 3:
      {
	Assert (level == 0, ExcFacesHaveNoLevel());

	return raw_quad_iterator (tria,
				  0,
				  0,
				  this);
      }


      default:
	    Assert (false, ExcNotImplemented());
	    return raw_hex_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::quad_iterator
DoFHandler<dim,spacedim>::begin_quad (const unsigned int level) const
{
				   // level is checked in begin_raw
  raw_quad_iterator ri = begin_raw_quad (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::active_quad_iterator
DoFHandler<dim,spacedim>::begin_active_quad (const unsigned int level) const
{
				   // level is checked in begin_raw
  quad_iterator i = begin_quad (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_quad_iterator
DoFHandler<dim, spacedim>::end_raw_quad (const unsigned int level) const
{
  Assert (dim == 2 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 2)
    return (level == tria->n_levels()-1 ?
	    end_quad() :
	    begin_raw_quad (level+1));
  else
    return end_quad();
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::quad_iterator
DoFHandler<dim, spacedim>::end_quad (const unsigned int level) const
{
  Assert (dim == 2 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 2)
    return (level == tria->n_levels()-1 ?
	    quad_iterator(end_quad()) :
	    begin_quad (level+1));
  else
    return quad_iterator(end_quad());
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_quad_iterator
DoFHandler<dim, spacedim>::end_active_quad (const unsigned int level) const
{
  Assert(dim == 2 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 2)
    return (level == tria->n_levels()-1 ?
	    active_quad_iterator(end_quad()) :
	    begin_active_quad (level+1));
  else
    return active_quad_iterator(end_quad());
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_quad_iterator
DoFHandler<dim,spacedim>::end_quad () const
{
  return raw_quad_iterator (tria,
			    -1,
			    -1,
			    this);
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_quad_iterator
DoFHandler<dim,spacedim>::last_raw_quad (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_quad_iterator();
      case 2:
	    Assert (level<tria->n_levels(),
		    ExcInvalidLevel(level));
	    Assert (tria->n_raw_quads(level) != 0,
		    ExcEmptyLevel (level));
	    return raw_quad_iterator (tria,
				      level,
				      tria->n_raw_quads(level)-1,
				      this);
      case 3:
	    Assert (level == 0, ExcFacesHaveNoLevel());
	    return raw_quad_iterator (tria,
				      0,
				      tria->n_raw_quads()-1,
				      this);
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_quad_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_quad_iterator
DoFHandler<dim,spacedim>::last_raw_quad () const
{
  if (dim == 2)
    return last_raw_quad (tria->n_levels()-1);
  else
    return last_raw_quad (0);
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::quad_iterator
DoFHandler<dim,spacedim>::last_quad (const unsigned int level) const
{
				   // level is checked in last_raw
  raw_quad_iterator ri = last_raw_quad(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == IteratorState::valid)
    if (ri->used()==true)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::quad_iterator
DoFHandler<dim,spacedim>::last_quad () const
{
  if (dim == 2)
    return last_quad (tria->n_levels()-1);
  else
    return last_quad (0);
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::active_quad_iterator
DoFHandler<dim,spacedim>::last_active_quad (const unsigned int level) const
{
				   // level is checked in last_raw
  quad_iterator i = last_quad(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == IteratorState::valid)
    if (i->has_children()==false)
      return i;
  return i;
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::active_quad_iterator
DoFHandler<dim,spacedim>::last_active_quad () const
{
  if (dim == 2)
    return last_active_quad (tria->n_levels()-1);
  else
    return last_active_quad (0);
}


/*------------------------ Hex iterator functions ------------------------*/


template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::raw_hex_iterator
DoFHandler<dim,spacedim>::begin_raw_hex (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
      case 2:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_hex_iterator();
      case 3:
      {
	Assert (level<tria->n_levels(), ExcInvalidLevel(level));

	if (tria->n_raw_hexs(level) == 0)
	  return end_hex();

	return raw_hex_iterator (tria,
				 level,
				 0,
				 this);
      }

      default:
	    Assert (false, ExcNotImplemented());
	    return raw_hex_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim,spacedim>::hex_iterator
DoFHandler<dim,spacedim>::begin_hex (const unsigned int level) const
{
				   // level is checked in begin_raw
  raw_hex_iterator ri = begin_raw_hex (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_hex_iterator
DoFHandler<dim, spacedim>::begin_active_hex (const unsigned int level) const
{
				   // level is checked in begin_raw
  hex_iterator i = begin_hex (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_hex_iterator
DoFHandler<dim, spacedim>::end_raw_hex (const unsigned int level) const
{
  return (level == tria->n_levels()-1 ?
	  end_hex() :
	  begin_raw_hex (level+1));
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::hex_iterator
DoFHandler<dim, spacedim>::end_hex (const unsigned int level) const
{
  return (level == tria->n_levels()-1 ?
	  hex_iterator(end_hex()) :
	  begin_hex (level+1));
}


template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_hex_iterator
DoFHandler<dim, spacedim>::end_active_hex (const unsigned int level) const
{
  return (level == tria->n_levels()-1 ?
	  active_hex_iterator(end_hex()) :
	  begin_active_hex (level+1));
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_hex_iterator
DoFHandler<dim, spacedim>::end_hex () const
{
  return raw_hex_iterator (tria,
			   -1,
			   -1,
			   this);
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_hex_iterator
DoFHandler<dim, spacedim>::last_raw_hex (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
      case 2:
	    Assert (false, ExcImpossibleInDim(dim));
	    return raw_hex_iterator();

      case 3:
	    Assert (level<tria->n_levels(),
		    ExcInvalidLevel(level));
	    Assert (tria->n_raw_hexs(level) != 0,
		    ExcEmptyLevel (level));

	    return raw_hex_iterator (tria,
				     level,
				     tria->n_raw_hexs(level)-1,
				     this);
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_hex_iterator();
    }
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::raw_hex_iterator
DoFHandler<dim, spacedim>::last_raw_hex () const
{
  return last_raw_hex (tria->n_levels()-1);
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::hex_iterator
DoFHandler<dim, spacedim>::last_hex (const unsigned int level) const
{
				   // level is checked in last_raw
  raw_hex_iterator ri = last_raw_hex(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == IteratorState::valid)
    if (ri->used()==true)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::hex_iterator
DoFHandler<dim, spacedim>::last_hex () const
{
  return last_hex (tria->n_levels()-1);
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_hex_iterator
DoFHandler<dim, spacedim>::last_active_hex (const unsigned int level) const
{
				   // level is checked in last_raw
  hex_iterator i = last_hex(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == IteratorState::valid)
    if (i->has_children()==false)
      return i;
  return i;
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_hex_iterator
DoFHandler<dim, spacedim>::last_active_hex () const
{
  return last_active_hex (tria->n_levels()-1);
}




//---------------------------------------------------------------------------



#if deal_II_dimension == 1

template <>
unsigned int DoFHandler<1>::n_boundary_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  return 2*selected_fe->dofs_per_vertex;
}



template <>
unsigned int DoFHandler<1>::n_boundary_dofs (const FunctionMap &boundary_indicators) const
{
  Assert (selected_fe != 0, ExcNoFESelected());

				   // check that only boundary
				   // indicators 0 and 1 are allowed
				   // in 1d
  for (FunctionMap::const_iterator i=boundary_indicators.begin();
       i!=boundary_indicators.end(); ++i)
    Assert ((i->first == 0) || (i->first == 1),
	    ExcInvalidBoundaryIndicator());

  return boundary_indicators.size()*selected_fe->dofs_per_vertex;
}



template <>
unsigned int DoFHandler<1>::n_boundary_dofs (const std::set<unsigned char> &boundary_indicators) const
{
  Assert (selected_fe != 0, ExcNoFESelected());

				   // check that only boundary
				   // indicators 0 and 1 are allowed
				   // in 1d
  for (std::set<unsigned char>::const_iterator i=boundary_indicators.begin();
       i!=boundary_indicators.end(); ++i)
    Assert ((*i == 0) || (*i == 1),
	    ExcInvalidBoundaryIndicator());

  return boundary_indicators.size()*selected_fe->dofs_per_vertex;
}


template <>
unsigned int DoFHandler<1,2>::n_boundary_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  return 2*selected_fe->dofs_per_vertex;
}



template <>
unsigned int DoFHandler<1,2>::n_boundary_dofs (const FunctionMap &boundary_indicators) const
{
  Assert (selected_fe != 0, ExcNoFESelected());

				   // check that only boundary
				   // indicators 0 and 1 are allowed
				   // in 1d
  for (FunctionMap::const_iterator i=boundary_indicators.begin();
       i!=boundary_indicators.end(); ++i)
    Assert ((i->first == 0) || (i->first == 1),
	    ExcInvalidBoundaryIndicator());

  return boundary_indicators.size()*selected_fe->dofs_per_vertex;
}



template <>
unsigned int DoFHandler<1,2>::n_boundary_dofs (const std::set<unsigned char> &boundary_indicators) const
{
  Assert (selected_fe != 0, ExcNoFESelected());

				   // check that only boundary
				   // indicators 0 and 1 are allowed
				   // in 1d
  for (std::set<unsigned char>::const_iterator i=boundary_indicators.begin();
       i!=boundary_indicators.end(); ++i)
    Assert ((*i == 0) || (*i == 1),
	    ExcInvalidBoundaryIndicator());

  return boundary_indicators.size()*selected_fe->dofs_per_vertex;
}

#endif


template<int dim, int spacedim>
unsigned int DoFHandler<dim,spacedim>::n_boundary_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());

  std::set<int> boundary_dofs;

  const unsigned int dofs_per_face = selected_fe->dofs_per_face;
  std::vector<unsigned int> dofs_on_face(dofs_per_face);

				   // loop over all faces to check
				   // whether they are at a
				   // boundary. note that we need not
				   // take special care of single
				   // lines (using
				   // @p{cell->has_boundary_lines}),
				   // since we do not support
				   // boundaries of dimension dim-2,
				   // and so every boundary line is
				   // also part of a boundary face.
  active_face_iterator face = begin_active_face (),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (face->at_boundary())
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  boundary_dofs.insert(dofs_on_face[i]);
      }
  return boundary_dofs.size();
}



template<int dim, int spacedim>
unsigned int
DoFHandler<dim,spacedim>::n_boundary_dofs (const FunctionMap &boundary_indicators) const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (boundary_indicators.find(255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());

  std::set<int> boundary_dofs;

  const unsigned int dofs_per_face = selected_fe->dofs_per_face;
  std::vector<unsigned int> dofs_on_face(dofs_per_face);
  active_face_iterator face = begin_active_face (),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (boundary_indicators.find(face->boundary_indicator()) !=
	boundary_indicators.end())
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  boundary_dofs.insert(dofs_on_face[i]);
      }
  return boundary_dofs.size();
}



template<int dim, int spacedim>
unsigned int
DoFHandler<dim,spacedim>::n_boundary_dofs (const std::set<unsigned char> &boundary_indicators) const
{
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (boundary_indicators.find (255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());

  std::set<int> boundary_dofs;

  const unsigned int dofs_per_face = selected_fe->dofs_per_face;
  std::vector<unsigned int> dofs_on_face(dofs_per_face);
  active_face_iterator face = begin_active_face (),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (std::find (boundary_indicators.begin(),
		   boundary_indicators.end(),
		   face->boundary_indicator()) !=
	boundary_indicators.end())
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  boundary_dofs.insert(dofs_on_face[i]);
      }
  return boundary_dofs.size();
}



template<int dim, int spacedim>
unsigned int
DoFHandler<dim,spacedim>::memory_consumption () const
{
  unsigned int mem = (MemoryConsumption::memory_consumption (tria) +
		      MemoryConsumption::memory_consumption (selected_fe) +
		      MemoryConsumption::memory_consumption (block_info_object) +
		      MemoryConsumption::memory_consumption (levels) +
		      MemoryConsumption::memory_consumption (*faces) +
		      MemoryConsumption::memory_consumption (faces) +
		      MemoryConsumption::memory_consumption (used_dofs) +
		      MemoryConsumption::memory_consumption (vertex_dofs));
  for (unsigned int i=0; i<levels.size(); ++i)
    mem += MemoryConsumption::memory_consumption (*levels[i]);

  return mem;
}



template<int dim, int spacedim>
void DoFHandler<dim,spacedim>::
  distribute_dofs (const FiniteElement<dim,spacedim> &ff,
		   const unsigned int                 offset)
{
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());

  selected_fe = &ff;

                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  clear_space ();
  internal::DoFHandler::Implementation::reserve_space (*this);

				   // Clear user flags because we will
				   // need them. But first we save
				   // them and make sure that we
				   // restore them later such that at
				   // the end of this function the
				   // Triangulation will be in the
				   // same state as it was at the
				   // beginning of this function.
  std::vector<bool> user_flags;
  tria->save_user_flags(user_flags);
  const_cast<Triangulation<dim,spacedim> &>(*tria).clear_user_flags ();

  unsigned int next_free_dof = offset;
  active_cell_iterator cell = begin_active(),
		       endc = end();

  for (; cell != endc; ++cell)
    next_free_dof
      = internal::DoFHandler::Implementation::
        distribute_dofs_on_cell (*this, cell, next_free_dof);

  used_dofs = next_free_dof;

				   // update the cache used
				   // for cell dof indices
  for (cell_iterator cell = begin(); cell != end(); ++cell)
    cell->update_cell_dof_indices_cache ();

				   // finally restore the user flags
  const_cast<Triangulation<dim,spacedim> &>(*tria).load_user_flags(user_flags);

  block_info_object.initialize(*this);
}


template<int dim, int spacedim>
void DoFHandler<dim,spacedim>::initialize_local_block_info ()
{
  block_info_object.initialize_local(*this);
}


template<int dim, int spacedim>
void DoFHandler<dim,spacedim>::clear ()
{
				   // release lock to old fe
  selected_fe = 0;

				   // release memory
  clear_space ();
}



template <int dim, int spacedim>
void
DoFHandler<dim,spacedim>::
renumber_dofs (const std::vector<unsigned int> &new_numbers)
{
  Assert (new_numbers.size() == n_dofs(), ExcRenumberingIncomplete());

#ifdef DEBUG
				   // assert that the new indices are
				   // consecutively numbered
  if (true)
    {
      std::vector<unsigned int> tmp(new_numbers);
      std::sort (tmp.begin(), tmp.end());
      std::vector<unsigned int>::const_iterator p = tmp.begin();
      unsigned int                         i = 0;
      for (; p!=tmp.end(); ++p, ++i)
	Assert (*p == i, ExcNewNumbersNotConsecutive(i));
    }
#endif


  internal::DoFHandler::Implementation::renumber_dofs (new_numbers, *this);

				   // update the cache used for cell dof
				   // indices
  for (cell_iterator cell = begin(); cell != end(); ++cell)
    cell->update_cell_dof_indices_cache ();

}



template <int dim, int spacedim>
unsigned int
DoFHandler<dim,spacedim>::max_couplings_between_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());

  return internal::DoFHandler::Implementation::max_couplings_between_dofs (*this);
}



template <int dim, int spacedim>
unsigned int
DoFHandler<dim,spacedim>::max_couplings_between_boundary_dofs () const
{
  Assert (selected_fe != 0, ExcNoFESelected());

  switch (dim)
    {
      case 1:
	    return selected_fe->dofs_per_vertex;
      case 2:
	    return (3*selected_fe->dofs_per_vertex +
		    2*selected_fe->dofs_per_line);
      case 3:
					     // we need to take refinement of
					     // one boundary face into
					     // consideration here; in fact,
					     // this function returns what
					     // #max_coupling_between_dofs<2>
					     // returns
					     //
					     // we assume here, that only four
					     // faces meet at the boundary;
					     // this assumption is not
					     // justified and needs to be
					     // fixed some time. fortunately,
					     // ommitting it for now does no
					     // harm since the matrix will cry
					     // foul if its requirements are
					     // not satisfied
	    return (19*selected_fe->dofs_per_vertex +
		    28*selected_fe->dofs_per_line +
		    8*selected_fe->dofs_per_quad);
      default:
	    Assert (false, ExcNotImplemented());
	    return numbers::invalid_unsigned_int;
    }
}



template<int dim, int spacedim>
void DoFHandler<dim,spacedim>::clear_space ()
{
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.resize (0);

  delete faces;
  faces = 0;

  std::vector<unsigned int> tmp;
  std::swap (vertex_dofs, tmp);
}


/*-------------- Explicit Instantiations -------------------------------*/
template class DoFHandler<deal_II_dimension>;

#if deal_II_dimension==1 || deal_II_dimension==2
template class DoFHandler<deal_II_dimension,deal_II_dimension+1>;
#endif


DEAL_II_NAMESPACE_CLOSE
