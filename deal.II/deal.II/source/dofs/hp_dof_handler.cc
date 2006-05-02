//----------------------------  hp_dof_handler.cc  ------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_dof_handler.cc  ------------------------


#include <base/memory_consumption.h>
#include <dofs/hp_dof_handler.h>
#include <dofs/hp_dof_levels.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_levels.h>
#include <grid/tria.h>
#include <grid/geometry_info.h>
#include <fe/fe.h>

#include <set>
#include <algorithm>
#include <functional>


namespace hp
{
  
  template <int dim>
  const unsigned int DoFHandler<dim>::dimension;

  template <int dim>
  const unsigned int DoFHandler<dim>::invalid_dof_index;

  template <int dim>
  const unsigned int DoFHandler<dim>::default_fe_index;
  


  template <int dim>
  DoFHandler<dim>::DoFHandler (const Triangulation<dim> &tria)
                  :
                  tria(&tria),
                  used_dofs (0)
  {
    create_active_fe_table ();
    tria.add_refinement_listener (*this);
  }


  template <int dim>
  DoFHandler<dim>::~DoFHandler ()
  {
                                     // unsubscribe as a listener to refinement
                                     // of the underlying triangulation
    tria->remove_refinement_listener (*this);
  
                                     // ...and release allocated memory
    clear ();
  }


#if deal_II_dimension == 1

  template <>
  DoFHandler<1>::raw_cell_iterator
  DoFHandler<1>::begin_raw (const unsigned int level) const
  {
    return begin_raw_line (level);
  }


  template <>
  DoFHandler<1>::cell_iterator
  DoFHandler<1>::begin (const unsigned int level) const
  {
    return begin_line (level);
  }


  template <>
  DoFHandler<1>::active_cell_iterator
  DoFHandler<1>::begin_active (const unsigned int level) const
  {
    return begin_active_line (level);
  }


  template <>
  DoFHandler<1>::raw_cell_iterator
  DoFHandler<1>::end () const
  {
    return end_line ();
  }


  template <>
  DoFHandler<1>::raw_cell_iterator
  DoFHandler<1>::last_raw () const
  {
    return last_raw_line ();
  }


  template <>
  DoFHandler<1>::raw_cell_iterator
  DoFHandler<1>::last_raw (const unsigned int level) const
  {
    return last_raw_line (level);
  }


  template <>
  DoFHandler<1>::cell_iterator
  DoFHandler<1>::last () const
  {
    return last_line ();
  }


  template <>
  DoFHandler<1>::cell_iterator
  DoFHandler<1>::last (const unsigned int level) const
  {
    return last_line (level);
  }


  template <>
  DoFHandler<1>::active_cell_iterator
  DoFHandler<1>::last_active () const
  {
    return last_active_line ();
  }


  template <>
  DoFHandler<1>::active_cell_iterator
  DoFHandler<1>::last_active (const unsigned int level) const
  {
    return last_active_line (level);
  }


  template <>
  DoFHandler<1>::raw_face_iterator
  DoFHandler<1>::begin_raw_face (const unsigned int) const
  {
    Assert (false, ExcFunctionNotUseful());
    return 0;
  }


  template <>
  DoFHandler<1>::face_iterator
  DoFHandler<1>::begin_face (const unsigned int) const
  {
    Assert (false, ExcFunctionNotUseful());
    return 0;
  }


  template <>
  DoFHandler<1>::active_face_iterator
  DoFHandler<1>::begin_active_face (const unsigned int) const
  {
    Assert (false, ExcFunctionNotUseful());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_face_iterator
  DoFHandler<1>::end_face () const
  {
    Assert (false, ExcFunctionNotUseful());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_face_iterator
  DoFHandler<1>::last_raw_face () const
  {
    Assert (false, ExcFunctionNotUseful());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_face_iterator
  DoFHandler<1>::last_raw_face (const unsigned int) const
  {
    Assert (false, ExcFunctionNotUseful());
    return 0;
  }


  template <>
  DoFHandler<1>::face_iterator
  DoFHandler<1>::last_face () const
  {
    Assert (false, ExcFunctionNotUseful());
    return 0;
  }


  template <>
  DoFHandler<1>::face_iterator
  DoFHandler<1>::last_face (const unsigned int) const
  {
    Assert (false, ExcFunctionNotUseful());
    return 0;
  }


  template <>
  DoFHandler<1>::active_face_iterator
  DoFHandler<1>::last_active_face () const
  {
    Assert (false, ExcFunctionNotUseful());
    return 0;
  }


  template <>
  DoFHandler<1>::active_face_iterator
  DoFHandler<1>::last_active_face (const unsigned int) const
  {
    Assert (false, ExcFunctionNotUseful());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_quad_iterator
  DoFHandler<1>::begin_raw_quad (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::quad_iterator
  DoFHandler<1>::begin_quad (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::active_quad_iterator
  DoFHandler<1>::begin_active_quad (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_quad_iterator
  DoFHandler<1>::end_quad () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_quad_iterator
  DoFHandler<1>::last_raw_quad (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::quad_iterator
  DoFHandler<1>::last_quad (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::active_quad_iterator
  DoFHandler<1>::last_active_quad (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_quad_iterator
  DoFHandler<1>::last_raw_quad () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::quad_iterator
  DoFHandler<1>::last_quad () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::active_quad_iterator
  DoFHandler<1>::last_active_quad () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_hex_iterator
  DoFHandler<1>::begin_raw_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::hex_iterator
  DoFHandler<1>::begin_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::active_hex_iterator
  DoFHandler<1>::begin_active_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_hex_iterator
  DoFHandler<1>::end_hex () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_hex_iterator
  DoFHandler<1>::last_raw_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::raw_hex_iterator
  DoFHandler<1>::last_raw_hex () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::hex_iterator
  DoFHandler<1>::last_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::hex_iterator
  DoFHandler<1>::last_hex () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::active_hex_iterator
  DoFHandler<1>::last_active_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<1>::active_hex_iterator
  DoFHandler<1>::last_active_hex () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }

#endif


#if deal_II_dimension == 2

  template <>
  DoFHandler<2>::raw_cell_iterator
  DoFHandler<2>::begin_raw (const unsigned int level) const
  {
    return begin_raw_quad (level);
  }


  template <>
  DoFHandler<2>::cell_iterator
  DoFHandler<2>::begin (const unsigned int level) const
  {
    return begin_quad (level);
  }


  template <>
  DoFHandler<2>::active_cell_iterator
  DoFHandler<2>::begin_active (const unsigned int level) const
  {
    return begin_active_quad (level);
  }


  template <>
  DoFHandler<2>::raw_cell_iterator
  DoFHandler<2>::end () const
  {
    return end_quad ();
  }


  template <>
  DoFHandler<2>::raw_cell_iterator
  DoFHandler<2>::last_raw () const
  {
    return last_raw_quad ();
  }


  template <>
  DoFHandler<2>::raw_cell_iterator
  DoFHandler<2>::last_raw (const unsigned int level) const
  {
    return last_raw_quad (level);
  }


  template <>
  DoFHandler<2>::cell_iterator
  DoFHandler<2>::last () const
  {
    return last_quad ();
  }


  template <>
  DoFHandler<2>::cell_iterator
  DoFHandler<2>::last (const unsigned int level) const
  {
    return last_quad (level);
  }


  template <>
  DoFHandler<2>::active_cell_iterator
  DoFHandler<2>::last_active () const
  {
    return last_active_quad ();
  }


  template <>
  DoFHandler<2>::active_cell_iterator
  DoFHandler<2>::last_active (const unsigned int level) const
  {
    return last_active_quad (level);
  }


  template <>
  DoFHandler<2>::raw_face_iterator
  DoFHandler<2>::begin_raw_face (const unsigned int level) const
  {
    return begin_raw_line (level);
  }


  template <>
  DoFHandler<2>::face_iterator
  DoFHandler<2>::begin_face (const unsigned int level) const
  {
    return begin_line (level);
  }


  template <>
  DoFHandler<2>::active_face_iterator
  DoFHandler<2>::begin_active_face (const unsigned int level) const
  {
    return begin_active_line (level);
  }


  template <>
  DoFHandler<2>::raw_face_iterator
  DoFHandler<2>::end_face () const
  {
    return end_line ();
  }


  template <>
  DoFHandler<2>::raw_face_iterator
  DoFHandler<2>::last_raw_face () const
  {
    return last_raw_line ();
  }


  template <>
  DoFHandler<2>::raw_face_iterator
  DoFHandler<2>::last_raw_face (const unsigned int level) const
  {
    return last_raw_line (level);
  }


  template <>
  DoFHandler<2>::face_iterator
  DoFHandler<2>::last_face () const
  {
    return last_line ();
  }


  template <>
  DoFHandler<2>::face_iterator
  DoFHandler<2>::last_face (const unsigned int level) const
  {
    return last_line (level);
  }


  template <>
  DoFHandler<2>::active_face_iterator
  DoFHandler<2>::last_active_face () const
  {
    return last_active_line ();
  }


  template <>
  DoFHandler<2>::active_face_iterator
  DoFHandler<2>::last_active_face (const unsigned int level) const
  {
    return last_active_line (level);
  }


  template <>
  DoFHandler<2>::raw_hex_iterator
  DoFHandler<2>::begin_raw_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<2>::hex_iterator
  DoFHandler<2>::begin_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<2>::active_hex_iterator
  DoFHandler<2>::begin_active_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<2>::raw_hex_iterator
  DoFHandler<2>::end_hex () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<2>::raw_hex_iterator
  DoFHandler<2>::last_raw_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<2>::raw_hex_iterator
  DoFHandler<2>::last_raw_hex () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<2>::hex_iterator
  DoFHandler<2>::last_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<2>::hex_iterator
  DoFHandler<2>::last_hex () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<2>::active_hex_iterator
  DoFHandler<2>::last_active_hex (const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


  template <>
  DoFHandler<2>::active_hex_iterator
  DoFHandler<2>::last_active_hex () const
  {
    Assert (false, ExcNotImplemented());
    return 0;
  }


#endif


#if deal_II_dimension == 3

  template <>
  DoFHandler<3>::raw_cell_iterator
  DoFHandler<3>::begin_raw (const unsigned int level) const
  {
    return begin_raw_hex (level);
  }


  template <>
  DoFHandler<3>::cell_iterator
  DoFHandler<3>::begin (const unsigned int level) const
  {
    return begin_hex (level);
  }


  template <>
  DoFHandler<3>::active_cell_iterator
  DoFHandler<3>::begin_active (const unsigned int level) const
  {
    return begin_active_hex (level);
  }


  template <>
  DoFHandler<3>::raw_cell_iterator
  DoFHandler<3>::end () const
  {
    return end_hex ();
  }


  template <>
  DoFHandler<3>::raw_cell_iterator
  DoFHandler<3>::last_raw () const
  {
    return last_raw_hex ();
  }


  template <>
  DoFHandler<3>::raw_cell_iterator
  DoFHandler<3>::last_raw (const unsigned int level) const
  {
    return last_raw_hex (level);
  }


  template <>
  DoFHandler<3>::cell_iterator
  DoFHandler<3>::last () const
  {
    return last_hex ();
  }


  template <>
  DoFHandler<3>::cell_iterator
  DoFHandler<3>::last (const unsigned int level) const
  {
    return last_hex (level);
  }


  template <>
  DoFHandler<3>::active_cell_iterator
  DoFHandler<3>::last_active () const
  {
    return last_active_hex ();
  }


  template <>
  DoFHandler<3>::active_cell_iterator
  DoFHandler<3>::last_active (const unsigned int level) const
  {
    return last_active_hex (level);
  }


  template <>
  DoFHandler<3>::raw_face_iterator
  DoFHandler<3>::begin_raw_face (const unsigned int level) const
  {
    return begin_raw_quad (level);
  }


  template <>
  DoFHandler<3>::face_iterator
  DoFHandler<3>::begin_face (const unsigned int level) const
  {
    return begin_quad (level);
  }


  template <>
  DoFHandler<3>::active_face_iterator
  DoFHandler<3>::begin_active_face (const unsigned int level) const
  {
    return begin_active_quad (level);
  }


  template <>
  DoFHandler<3>::raw_face_iterator
  DoFHandler<3>::end_face () const
  {
    return end_quad ();
  }


  template <>
  DoFHandler<3>::raw_face_iterator
  DoFHandler<3>::last_raw_face () const
  {
    return last_raw_quad ();
  }


  template <>
  DoFHandler<3>::raw_face_iterator
  DoFHandler<3>::last_raw_face (const unsigned int level) const
  {
    return last_raw_quad (level);
  }


  template <>
  DoFHandler<3>::face_iterator
  DoFHandler<3>::last_face () const
  {
    return last_quad ();
  }


  template <>
  DoFHandler<3>::face_iterator
  DoFHandler<3>::last_face (const unsigned int level) const
  {
    return last_quad (level);
  }


  template <>
  DoFHandler<3>::active_face_iterator
  DoFHandler<3>::last_active_face () const
  {
    return last_active_quad ();
  }


  template <>
  DoFHandler<3>::active_face_iterator
  DoFHandler<3>::last_active_face (const unsigned int level) const
  {
    return last_active_quad (level);
  }


#endif


  template <int dim>
  typename DoFHandler<dim>::raw_line_iterator
  DoFHandler<dim>::begin_raw_line (const unsigned int level) const
  {
    return raw_line_iterator (tria,
                              tria->begin_raw_line(level)->level(),
                              tria->begin_raw_line(level)->index(),
                              this);
  }


  template <int dim>
  typename DoFHandler<dim>::line_iterator
  DoFHandler<dim>::begin_line (const unsigned int level) const
  {
    return line_iterator (tria,
                          tria->begin_line(level)->level(),
                          tria->begin_line(level)->index(),
                          this);
  }


  template <int dim>
  typename DoFHandler<dim>::active_line_iterator
  DoFHandler<dim>::begin_active_line (const unsigned int level) const
  {
    return active_line_iterator (tria,
                                 tria->begin_active_line(level)->level(),
                                 tria->begin_active_line(level)->index(),
                                 this);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_quad_iterator
  DoFHandler<dim>::begin_raw_quad (const unsigned int level) const
  {
    return raw_quad_iterator (tria,
                              tria->begin_raw_quad(level)->level(),
                              tria->begin_raw_quad(level)->index(),
                              this);
  }


  template <int dim>
  typename DoFHandler<dim>::quad_iterator
  DoFHandler<dim>::begin_quad (const unsigned int level) const
  {
    return quad_iterator (tria,
                          tria->begin_quad(level)->level(),
                          tria->begin_quad(level)->index(),
                          this);
  }


  template <int dim>
  typename DoFHandler<dim>::active_quad_iterator
  DoFHandler<dim>::begin_active_quad (const unsigned int level) const
  {
    return active_quad_iterator (tria,
                                 tria->begin_active_quad(level)->level(),
                                 tria->begin_active_quad(level)->index(),
                                 this);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_hex_iterator
  DoFHandler<dim>::begin_raw_hex (const unsigned int level) const
  {
    return raw_hex_iterator (tria,
                             tria->begin_raw_hex(level)->level(),
                             tria->begin_raw_hex(level)->index(),
                             this);
  }


  template <int dim>
  typename DoFHandler<dim>::hex_iterator
  DoFHandler<dim>::begin_hex (const unsigned int level) const
  {
    return hex_iterator (tria,
                         tria->begin_hex(level)->level(),
                         tria->begin_hex(level)->index(),
                         this);
  }


  template <int dim>
  typename DoFHandler<dim>::active_hex_iterator
  DoFHandler<dim>::begin_active_hex (const unsigned int level) const
  {
    return active_hex_iterator (tria,
                                tria->begin_active_hex(level)->level(),
                                tria->begin_active_hex(level)->index(),
                                this);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_line_iterator
  DoFHandler<dim>::end_line () const
  {
    return raw_line_iterator (tria, -1, -1, this);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_quad_iterator
  DoFHandler<dim>::end_quad () const
  {
    return raw_quad_iterator (tria, -1, -1, this);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_hex_iterator
  DoFHandler<dim>::end_hex () const
  {
    return raw_hex_iterator (tria, -1, -1, this);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_cell_iterator
  DoFHandler<dim>::end_raw (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            end() :
            begin_raw (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::cell_iterator
  DoFHandler<dim>::end (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            cell_iterator(end()) :
            begin (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::active_cell_iterator
  DoFHandler<dim>::end_active (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            active_cell_iterator(end()) :
            begin_active (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::raw_face_iterator
  DoFHandler<dim>::end_raw_face (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            end_face() :
            begin_raw_face (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::face_iterator
  DoFHandler<dim>::end_face (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            face_iterator(end_face()) :
            begin_face (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::active_face_iterator
  DoFHandler<dim>::end_active_face (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            active_face_iterator(end_face()) :
            begin_active_face (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::raw_line_iterator
  DoFHandler<dim>::end_raw_line (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            end_line() :
            begin_raw_line (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::line_iterator
  DoFHandler<dim>::end_line (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            line_iterator(end_line()) :
            begin_line (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::active_line_iterator
  DoFHandler<dim>::end_active_line (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            active_line_iterator(end_line()) :
            begin_active_line (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::raw_quad_iterator
  DoFHandler<dim>::end_raw_quad (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            end_quad() :
            begin_raw_quad (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::quad_iterator
  DoFHandler<dim>::end_quad (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            quad_iterator(end_quad()) :
            begin_quad (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::active_quad_iterator
  DoFHandler<dim>::end_active_quad (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            active_quad_iterator(end_quad()) :
            begin_active_quad (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::raw_hex_iterator
  DoFHandler<dim>::end_raw_hex (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            end_hex() :
            begin_raw_hex (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::hex_iterator
  DoFHandler<dim>::end_hex (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            hex_iterator(end_hex()) :
            begin_hex (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::active_hex_iterator
  DoFHandler<dim>::end_active_hex (const unsigned int level) const
  {
    return (level == levels.size()-1 ?
            active_hex_iterator(end_hex()) :
            begin_active_hex (level+1));
  }


  template <int dim>
  typename DoFHandler<dim>::raw_line_iterator
  DoFHandler<dim>::last_raw_line (const unsigned int level) const
  {
    return raw_line_iterator (tria,
                              tria->last_raw_line(level)->level(),
                              tria->last_raw_line(level)->index(),
                              this);
  }


  template <int dim>
  typename DoFHandler<dim>::line_iterator
  DoFHandler<dim>::last_line (const unsigned int level) const
  {
    return line_iterator (tria,
                          tria->last_line(level)->level(),
                          tria->last_line(level)->index(),
                          this);
  }


  template <int dim>
  typename DoFHandler<dim>::active_line_iterator
  DoFHandler<dim>::last_active_line (const unsigned int level) const
  {
    return active_line_iterator (tria,
                                 tria->last_active_line(level)->level(),
                                 tria->last_active_line(level)->index(),
                                 this);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_quad_iterator
  DoFHandler<dim>::last_raw_quad (const unsigned int level) const
  {
    return raw_quad_iterator (tria,
                              tria->last_raw_quad(level)->level(),
                              tria->last_raw_quad(level)->index(),
                              this);
  }


  template <int dim>
  typename DoFHandler<dim>::quad_iterator
  DoFHandler<dim>::last_quad (const unsigned int level) const
  {
    return quad_iterator (tria,
                          tria->last_quad(level)->level(),
                          tria->last_quad(level)->index(),
                          this);
  }


  template <int dim>
  typename DoFHandler<dim>::active_quad_iterator
  DoFHandler<dim>::last_active_quad (const unsigned int level) const
  {
    return active_quad_iterator (tria,
                                 tria->last_active_quad(level)->level(),
                                 tria->last_active_quad(level)->index(),
                                 this);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_hex_iterator
  DoFHandler<dim>::last_raw_hex (const unsigned int level) const
  {
    return raw_hex_iterator (tria,
                             tria->last_raw_hex(level)->level(),
                             tria->last_raw_hex(level)->index(),
                             this);
  }


  template <int dim>
  typename DoFHandler<dim>::hex_iterator
  DoFHandler<dim>::last_hex (const unsigned int level) const
  {
    return hex_iterator (tria,
                         tria->last_hex(level)->level(),
                         tria->last_hex(level)->index(),
                         this);
  }


  template <int dim>
  typename DoFHandler<dim>::active_hex_iterator
  DoFHandler<dim>::last_active_hex (const unsigned int level) const
  {
    return active_hex_iterator (tria,
                                tria->last_active_hex(level)->level(),
                                tria->last_active_hex(level)->index(),
                                this);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_line_iterator
  DoFHandler<dim>::last_raw_line () const
  {
    return last_raw_line (levels.size()-1);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_quad_iterator
  DoFHandler<dim>::last_raw_quad () const
  {
    return last_raw_quad (levels.size()-1);
  }


  template <int dim>
  typename DoFHandler<dim>::raw_hex_iterator
  DoFHandler<dim>::last_raw_hex () const
  {
    return last_raw_hex (levels.size()-1);
  }


  template <int dim>
  typename DoFHandler<dim>::line_iterator
  DoFHandler<dim>::last_line () const
  {
    return last_line (levels.size()-1);
  }


  template <int dim>
  typename DoFHandler<dim>::quad_iterator
  DoFHandler<dim>::last_quad () const
  {
    return last_quad (levels.size()-1);
  }


  template <int dim>
  typename DoFHandler<dim>::hex_iterator
  DoFHandler<dim>::last_hex () const
  {
    return last_hex (levels.size()-1);
  }


  template <int dim>
  typename DoFHandler<dim>::active_line_iterator
  DoFHandler<dim>::last_active_line () const
  {
    return last_active_line (levels.size()-1);
  }


  template <int dim>
  typename DoFHandler<dim>::active_quad_iterator
  DoFHandler<dim>::last_active_quad () const
  {
    return last_active_quad (levels.size()-1);
  }


  template <int dim>
  typename DoFHandler<dim>::active_hex_iterator
  DoFHandler<dim>::last_active_hex () const
  {
    return last_active_hex (levels.size()-1);
  }


//------------------------------------------------------------------



#if deal_II_dimension == 1

  template <>
  unsigned int DoFHandler<1>::n_boundary_dofs () const
  {
    Assert (finite_elements != 0, ExcNoFESelected());

    DoFHandler<1>::cell_iterator cell;
    unsigned int n = 0;
  
                                     // search left-most cell
    cell = this->begin_active();
    while (!cell->at_boundary(0))
      cell = cell->neighbor(0);
    n += cell->get_fe().dofs_per_vertex;
  
                                     // same with right-most cell
    cell = this->begin_active();
    while (!cell->at_boundary(1))
      cell = cell->neighbor(1);
    n += cell->get_fe().dofs_per_vertex;

    return n;
  }



  template <>
  unsigned int DoFHandler<1>::n_boundary_dofs (const FunctionMap &boundary_indicators) const
  {
    Assert (finite_elements != 0, ExcNoFESelected());

                                     // check that only boundary
                                     // indicators 0 and 1 are allowed
                                     // in 1d
    for (FunctionMap::const_iterator i=boundary_indicators.begin();
         i!=boundary_indicators.end(); ++i)
      Assert ((i->first == 0) || (i->first == 1),
              ExcInvalidBoundaryIndicator());
  
    DoFHandler<1>::active_cell_iterator cell;
    unsigned int n = 0;
  
                                     // search left-most cell
    if (boundary_indicators.find (0) != boundary_indicators.end())
      {
        cell = this->begin_active();
        while (!cell->at_boundary(0))
          cell = cell->neighbor(0);
        n += cell->get_fe().dofs_per_vertex;
      }
  
                                     // same with right-most cell
    if (boundary_indicators.find (1) != boundary_indicators.end())
      {
        cell = this->begin_active();
        while (!cell->at_boundary(1))
          cell = cell->neighbor(1);
        n += cell->get_fe().dofs_per_vertex;
      }
  
    return n;
  }



  template <>
  unsigned int DoFHandler<1>::n_boundary_dofs (const std::set<unsigned char> &boundary_indicators) const
  {
    Assert (finite_elements != 0, ExcNoFESelected());

                                     // check that only boundary
                                     // indicators 0 and 1 are allowed
                                     // in 1d
    for (std::set<unsigned char>::const_iterator i=boundary_indicators.begin();
         i!=boundary_indicators.end(); ++i)
      Assert ((*i == 0) || (*i == 1),
              ExcInvalidBoundaryIndicator());
  
    DoFHandler<1>::active_cell_iterator cell;
    unsigned int n = 0;
  
                                     // search left-most cell
    if (boundary_indicators.find (0) != boundary_indicators.end())
      {
        cell = this->begin_active();
        while (!cell->at_boundary(0))
          cell = cell->neighbor(0);
        n += cell->get_fe().dofs_per_vertex;
      }
  
                                     // same with right-most cell
    if (boundary_indicators.find (1) != boundary_indicators.end())
      {
        cell = this->begin_active();
        while (!cell->at_boundary(1))
          cell = cell->neighbor(1);
        n += cell->get_fe().dofs_per_vertex;
      }
  
    return n;
  }

#endif


  template <int dim>
  unsigned int DoFHandler<dim>::n_boundary_dofs () const
  {
    Assert (finite_elements != 0, ExcNoFESelected());

    std::set<int> boundary_dofs;
    std::vector<unsigned int> dofs_on_face;
    dofs_on_face.reserve (this->get_fe ().max_dofs_per_face());

                                     // loop over all faces to check
                                     // whether they are at a
                                     // boundary. note that we need not
                                     // take special care of single
                                     // lines in 3d (using
                                     // @p{cell->has_boundary_lines}),
                                     // since we do not support
                                     // boundaries of dimension dim-2,
                                     // and so every boundary line is
                                     // also part of a boundary face.
    typename DoFHandler<dim>::active_cell_iterator cell = this->begin_active (),
                                                   endc = this->end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f))
          {
            const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
            dofs_on_face.resize (dofs_per_face);
          
            cell->face(f)->get_dof_indices (dofs_on_face);
            for (unsigned int i=0; i<dofs_per_face; ++i)
              boundary_dofs.insert(dofs_on_face[i]);
          };
    return boundary_dofs.size();  
  }



  template <int dim>
  unsigned int
  DoFHandler<dim>::n_boundary_dofs (const FunctionMap &boundary_indicators) const
  {
    Assert (finite_elements != 0, ExcNoFESelected());
    Assert (boundary_indicators.find(255) == boundary_indicators.end(),
            ExcInvalidBoundaryIndicator());
  
                                     // same as above, but with
                                     // additional checks for set of
                                     // boundary indicators
    std::set<int> boundary_dofs;
    std::vector<unsigned int> dofs_on_face;
    dofs_on_face.reserve (this->get_fe ().max_dofs_per_face());

    typename DoFHandler<dim>::active_cell_iterator cell = this->begin_active (),
                                                   endc = this->end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f) &&
            (boundary_indicators.find(cell->face(f)->boundary_indicator()) !=
             boundary_indicators.end()))
          {
            const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
            dofs_on_face.resize (dofs_per_face);
          
            cell->face(f)->get_dof_indices (dofs_on_face);
            for (unsigned int i=0; i<dofs_per_face; ++i)
              boundary_dofs.insert(dofs_on_face[i]);
          };
    return boundary_dofs.size();
  }



  template <int dim>
  unsigned int
  DoFHandler<dim>::n_boundary_dofs (const std::set<unsigned char> &boundary_indicators) const
  {
    Assert (finite_elements != 0, ExcNoFESelected());
    Assert (boundary_indicators.find (255) == boundary_indicators.end(),
            ExcInvalidBoundaryIndicator());

                                     // same as above, but with
                                     // additional checks for set of
                                     // boundary indicators
    std::set<int> boundary_dofs;
    std::vector<unsigned int> dofs_on_face;
    dofs_on_face.reserve (this->get_fe ().max_dofs_per_face());

    typename DoFHandler<dim>::active_cell_iterator cell = this->begin_active (),
                                                   endc = this->end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f) &&
            (boundary_indicators.find(cell->face(f)->boundary_indicator()) !=
             boundary_indicators.end()))
          {
            const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
            dofs_on_face.resize (dofs_per_face);
          
            cell->face(f)->get_dof_indices (dofs_on_face);
            for (unsigned int i=0; i<dofs_per_face; ++i)
              boundary_dofs.insert(dofs_on_face[i]);
          };
    return boundary_dofs.size();
  }


  template <int dim>
  unsigned int
  DoFHandler<dim>::memory_consumption () const
  {
    unsigned int mem = (MemoryConsumption::memory_consumption (tria) +
                        MemoryConsumption::memory_consumption (finite_elements) +
                        MemoryConsumption::memory_consumption (tria) +
                        MemoryConsumption::memory_consumption (levels) +
                        MemoryConsumption::memory_consumption (used_dofs) +
                        MemoryConsumption::memory_consumption (vertex_dofs));
    for (unsigned int i=0; i<levels.size(); ++i)
      mem += MemoryConsumption::memory_consumption (*levels[i]);
  
    return mem;
  }



  template <int dim>
  void DoFHandler<dim>::distribute_dofs (const hp::FECollection<dim> &ff)
  {
    Assert (tria->n_levels() > 0, ExcInvalidTriangulation());

    finite_elements = &ff;
  
                                     // This call ensures that the
                                     // active_fe_indices vectors are
                                     // initialized correctly.
    create_active_fe_table ();

    reserve_space ();

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
    const_cast<Triangulation<dim> &>(*tria).clear_user_flags ();
  
    unsigned int next_free_dof = 0;
    active_cell_iterator cell = begin_active(),
                         endc = end();

    for (; cell != endc; ++cell) 
      next_free_dof = distribute_dofs_on_cell (cell, next_free_dof);

    used_dofs = next_free_dof;

                                     // finally restore the user flags
    const_cast<Triangulation<dim> &>(*tria).load_user_flags(user_flags);
  }


  template <int dim>
  void DoFHandler<dim>::clear ()
  {
                                     // release lock to old fe
    finite_elements = 0;

                                     // release memory
    clear_space ();
  }


#if deal_II_dimension == 1

  template <>
  unsigned int
  DoFHandler<1>::distribute_dofs_on_cell (active_cell_iterator &cell,
					  unsigned int          next_free_dof)
  {
//TODO[WB]: this may be completely bogus -- look through it    
                                     // get the fe for this cell
    const FiniteElement<1> &fe = cell->get_fe();

                                     // distribute dofs of vertices
    for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
      {
        cell_iterator neighbor = cell->neighbor(v);

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
                  for (unsigned int d=0; d<fe.dofs_per_vertex; ++d)
                    cell->set_vertex_dof_index (0, d,
                                                neighbor->vertex_dof_index (1, d),
						cell->active_fe_index());
                else
                  for (unsigned int d=0; d<fe.dofs_per_vertex; ++d)
                    cell->set_vertex_dof_index (1, d,
                                                neighbor->vertex_dof_index (0, d),
						cell->active_fe_index());

                                                 // next neighbor
                continue;
              };
          };
            
                                         // otherwise: create dofs newly
        for (unsigned int d=0; d<fe.dofs_per_vertex; ++d)
          cell->set_vertex_dof_index (v, d, next_free_dof++,
				      cell->active_fe_index());
      };
  
                                     // dofs of line
    for (unsigned int d=0; d<fe.dofs_per_line; ++d)
      cell->set_dof_index (d, next_free_dof++,
			   cell->active_fe_index());

                                     // note that this cell has been processed
    cell->set_user_flag ();
  
    return next_free_dof;
  }

#endif


#if deal_II_dimension == 2

  template <>
  unsigned int
  DoFHandler<2>::distribute_dofs_on_cell (active_cell_iterator &cell,
					  unsigned int          next_free_dof)
  {
//TODO[WB]: Check for correctness    
                                     // get the fe for this cell
    const FiniteElement<2> &fe = cell->get_fe();

    if (fe.dofs_per_vertex > 0)
                                       // number dofs on vertices
      for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
                                         // check whether dofs for this
                                         // vertex have been distributed
                                         // (only check the first dof)
        if (cell->vertex_dof_index(vertex, 0, cell->active_fe_index())
	    == invalid_dof_index)
          for (unsigned int d=0; d<fe.dofs_per_vertex; ++d)
            cell->set_vertex_dof_index (vertex, d, next_free_dof++,
					cell->active_fe_index());
    
                                     // for the four sides
//TODO[?] Does not work for continuous FEs
    if (fe.dofs_per_line > 0)
      for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
        {
          line_iterator line = cell->line(side);

                                           // distribute dofs if necessary:
                                           // check whether line dof is already
                                           // numbered (check only first dof)
          if (line->dof_index(0,cell->active_fe_index()) == invalid_dof_index)
                                             // if not: distribute dofs
            for (unsigned int d=0; d<fe.dofs_per_line; ++d)
              line->set_dof_index (d, next_free_dof++,
				   cell->active_fe_index());
        };


                                     // dofs of quad
    if (fe.dofs_per_quad > 0)
      for (unsigned int d=0; d<fe.dofs_per_quad; ++d)
        cell->set_dof_index (d, next_free_dof++,
			     cell->active_fe_index());


                                     // note that this cell has been processed
    cell->set_user_flag ();
  
    return next_free_dof;
  }

#endif


#if deal_II_dimension == 3

  template <>
  unsigned int
  DoFHandler<3>::distribute_dofs_on_cell (active_cell_iterator &cell,
					  unsigned int          next_free_dof)
  {
//TODO[WB]: check for correctness    
                                     // get the fe for this cell
    const FiniteElement<3> &fe = cell->get_fe();

    if (fe.dofs_per_vertex > 0)
                                       // number dofs on vertices
      for (unsigned int vertex=0; vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
                                         // check whether dofs for this
                                         // vertex have been distributed
                                         // (only check the first dof)
        if (cell->vertex_dof_index(vertex, 0,cell->active_fe_index())
	    == invalid_dof_index)
          for (unsigned int d=0; d<fe.dofs_per_vertex; ++d)
            cell->set_vertex_dof_index (vertex, d, next_free_dof++,
					cell->active_fe_index());
    
                                     // for the lines
    if (fe.dofs_per_line > 0)
      for (unsigned int l=0; l<GeometryInfo<3>::lines_per_cell; ++l)
        {
          line_iterator line = cell->line(l);
	
                                           // distribute dofs if necessary:
                                           // check whether line dof is already
                                           // numbered (check only first dof)
          if (line->dof_index(0, cell->active_fe_index()) == invalid_dof_index)
                                             // if not: distribute dofs
            for (unsigned int d=0; d<fe.dofs_per_line; ++d)
              line->set_dof_index (d, next_free_dof++,
				   cell->active_fe_index());
        };

                                     // for the quads
    if (fe.dofs_per_quad > 0)
      for (unsigned int q=0; q<GeometryInfo<3>::quads_per_cell; ++q)
        {
          quad_iterator quad = cell->quad(q);
	
                                           // distribute dofs if necessary:
                                           // check whether quad dof is already
                                           // numbered (check only first dof)
          if (quad->dof_index(0, cell->active_fe_index()) == invalid_dof_index)
                                             // if not: distribute dofs
            for (unsigned int d=0; d<fe.dofs_per_quad; ++d)
              quad->set_dof_index (d, next_free_dof++,
				   cell->active_fe_index());
        };


                                     // dofs of hex
    if (fe.dofs_per_hex > 0)
      for (unsigned int d=0; d<fe.dofs_per_hex; ++d)
        cell->set_dof_index (d, next_free_dof++,
			     cell->active_fe_index());


                                     // note that this cell has been processed
    cell->set_user_flag ();
  
    return next_free_dof;
  }

#endif


#if deal_II_dimension == 1

  template <>
  void DoFHandler<1>::renumber_dofs (const std::vector<unsigned int> &new_numbers)
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
      };
#endif

                                     // note that we can not use cell iterators
                                     // in this function since then we would
                                     // renumber the dofs on the interface of
                                     // two cells more than once. Anyway, this
                                     // way it's not only more correct but also
                                     // faster; note, however, that dof numbers
                                     // may be invalid_dof_index, namely when the appropriate
                                     // vertex/line/etc is unused
    for (std::vector<unsigned int>::iterator i=vertex_dofs.begin(); i!=vertex_dofs.end(); ++i)
      if (*i != invalid_dof_index)
        *i = new_numbers[*i];
  
    for (unsigned int level=0; level<levels.size(); ++level) 
      for (std::vector<unsigned int>::iterator i=levels[level]->line_dofs.begin();
           i!=levels[level]->line_dofs.end(); ++i)
        if (*i != invalid_dof_index)
          *i = new_numbers[*i];
  }

#endif


#if deal_II_dimension == 2

  template <>
  void DoFHandler<2>::renumber_dofs (const std::vector<unsigned int> &new_numbers)
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
      };
#endif

                                     // note that we can not use cell iterators
                                     // in this function since then we would
                                     // renumber the dofs on the interface of
                                     // two cells more than once. Anyway, this
                                     // way it's not only more correct but also
                                     // faster; note, however, that dof numbers
                                     // may be invalid_dof_index, namely when the appropriate
                                     // vertex/line/etc is unused
    for (std::vector<unsigned int>::iterator i=vertex_dofs.begin(); i!=vertex_dofs.end(); ++i)
      if (*i != invalid_dof_index)
        *i = new_numbers[*i];
  
    for (unsigned int level=0; level<levels.size(); ++level) 
      {
        for (std::vector<unsigned int>::iterator i=levels[level]->line_dofs.begin();
             i!=levels[level]->line_dofs.end(); ++i)
          if (*i != invalid_dof_index)
            *i = new_numbers[*i];
        for (std::vector<unsigned int>::iterator i=levels[level]->quad_dofs.begin();
             i!=levels[level]->quad_dofs.end(); ++i)
          if (*i != invalid_dof_index)
            *i = new_numbers[*i];
      };
  }

#endif


#if deal_II_dimension == 3

  template <>
  void DoFHandler<3>::renumber_dofs (const std::vector<unsigned int> &new_numbers)
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
        unsigned int                              i = 0;
        for (; p!=tmp.end(); ++p, ++i)
          Assert (*p == i, ExcNewNumbersNotConsecutive(i));
      };
#endif

                                     // note that we can not use cell iterators
                                     // in this function since then we would
                                     // renumber the dofs on the interface of
                                     // two cells more than once. Anyway, this
                                     // way it's not only more correct but also
                                     // faster; note, however, that dof numbers
                                     // may be invalid_dof_index, namely when the appropriate
                                     // vertex/line/etc is unused
    for (std::vector<unsigned int>::iterator i=vertex_dofs.begin(); i!=vertex_dofs.end(); ++i)
      if (*i != invalid_dof_index)
        *i = new_numbers[*i];
  
    for (unsigned int level=0; level<levels.size(); ++level) 
      {
        for (std::vector<unsigned int>::iterator i=levels[level]->line_dofs.begin();
             i!=levels[level]->line_dofs.end(); ++i)
          if (*i != invalid_dof_index)
            *i = new_numbers[*i];
        for (std::vector<unsigned int>::iterator i=levels[level]->quad_dofs.begin();
             i!=levels[level]->quad_dofs.end(); ++i)
          if (*i != invalid_dof_index)
            *i = new_numbers[*i];
        for (std::vector<unsigned int>::iterator i=levels[level]->hex_dofs.begin();
             i!=levels[level]->hex_dofs.end(); ++i)
          if (*i != invalid_dof_index)
            *i = new_numbers[*i];
      };
  }

#endif


#if deal_II_dimension == 1

  template <>
  unsigned int
  DoFHandler<1>::max_couplings_between_dofs () const
  {
    Assert (finite_elements != 0, ExcNoFESelected());
    return std::min(3*finite_elements->max_dofs_per_vertex() +
                    2*finite_elements->max_dofs_per_line(), n_dofs());
  }



  template <>
  unsigned int
  DoFHandler<1>::max_couplings_between_boundary_dofs () const
  {
    Assert (finite_elements != 0, ExcNoFESelected());
    return finite_elements->max_dofs_per_vertex();
  }

#endif


#if deal_II_dimension == 2

  template <>
  unsigned int
  DoFHandler<2>::max_couplings_between_dofs () const
  {
    Assert (finite_elements != 0, ExcNoFESelected());

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
    switch (tria->max_adjacent_cells())
      {
        case 4:
              max_couplings=19*finite_elements->max_dofs_per_vertex() +
                            28*finite_elements->max_dofs_per_line() +
                            8*finite_elements->max_dofs_per_quad();
              break;
        case 5:
              max_couplings=21*finite_elements->max_dofs_per_vertex() +
                            31*finite_elements->max_dofs_per_line() +
                            9*finite_elements->max_dofs_per_quad();
              break;
        case 6:
              max_couplings=28*finite_elements->max_dofs_per_vertex() +
                            42*finite_elements->max_dofs_per_line() +
                            12*finite_elements->max_dofs_per_quad();
              break;
        case 7:
              max_couplings=30*finite_elements->max_dofs_per_vertex() +
                            45*finite_elements->max_dofs_per_line() +
                            13*finite_elements->max_dofs_per_quad();
              break;
        case 8:
              max_couplings=37*finite_elements->max_dofs_per_vertex() +
                            56*finite_elements->max_dofs_per_line() +
                            16*finite_elements->max_dofs_per_quad();
              break;
        default:
              Assert (false, ExcNotImplemented());
              max_couplings=0;
      };
    return std::min(max_couplings,n_dofs());
  }



  template <>
  unsigned int
  DoFHandler<2>::max_couplings_between_boundary_dofs () const
  {
    Assert (finite_elements != 0, ExcNoFESelected());
    return (3*finite_elements->max_dofs_per_vertex()
            +
            2*finite_elements->max_dofs_per_line());
  }

#endif


#if deal_II_dimension == 3

  template <>
  unsigned int
  DoFHandler<3>::max_couplings_between_dofs () const
  {
//TODO:[?] Invent significantly better estimates than the ones in this function  
    Assert (finite_elements != 0, ExcNoFESelected());

                                     // doing the same thing here is a rather
                                     // complicated thing, compared to the 2d
                                     // case, since it is hard to draw pictures
                                     // with several refined hexahedra :-) so I
                                     // presently only give a coarse estimate
                                     // for the case that at most 8 hexes meet
                                     // at each vertex
                                     //
                                     // can anyone give better estimate here?
    const unsigned int max_adjacent_cells = tria->max_adjacent_cells();

    unsigned int max_couplings;
    if (max_adjacent_cells <= 8)
      max_couplings=7*7*7*finite_elements->max_dofs_per_vertex() +
                    7*6*7*3*finite_elements->max_dofs_per_line() +
                    9*4*7*3*finite_elements->max_dofs_per_quad() +
                    27*finite_elements->max_dofs_per_hex();
    else
      {
        Assert (false, ExcNotImplemented());
        max_couplings=0;
      }
  
    return std::min(max_couplings,n_dofs());
  }


  template <>
  unsigned int
  DoFHandler<3>::max_couplings_between_boundary_dofs () const
  {
    Assert (finite_elements != 0, ExcNoFESelected());
                                     // we need to take refinement of
                                     // one boundary face into consideration
                                     // here; in fact, this function returns
                                     // what #max_coupling_between_dofs<2>
                                     // returns
                                     //
                                     // we assume here, that only four faces
                                     // meet at the boundary; this assumption
                                     // is not justified and needs to be
                                     // fixed some time. fortunately, ommitting
                                     // it for now does no harm since the
                                     // matrix will cry foul if its requirements
                                     // are not satisfied
    return (19*finite_elements->max_dofs_per_vertex() +
            28*finite_elements->max_dofs_per_line() +
            8*finite_elements->max_dofs_per_quad());
  }


#endif



#if deal_II_dimension == 1

  template <>
  void DoFHandler<1>::reserve_space ()
  {
    Assert (finite_elements != 0, ExcNoFESelected());
    Assert (finite_elements->size() > 0, ExcNoFESelected());
    Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
    Assert (tria->n_levels() == levels.size (), ExcInternalError ());

    Assert (finite_elements->max_dofs_per_vertex() == 0,
            ExcMessage ("hp finite elements are presently only supported "
                        "for discontinuous elements"));

                                     // Release all space except the
                                     // active_fe_indices field which
                                     // we have to backup before
    {
      std::vector<std::vector<unsigned int> > active_fe_backup(levels.size ());
      for (unsigned int level = 0; level<levels.size (); ++level)
        std::swap (levels[level]->active_fe_indices, active_fe_backup[level]);

                                       // delete all levels and set them up
                                       // newly, since vectors are
                                       // troublesome if you want to change
                                       // their size
      clear_space ();
  
      for (unsigned int level=0; level<tria->n_levels(); ++level) 
        {
          levels.push_back (new internal::hp::DoFLevel<1>);
          std::swap (active_fe_backup[level], levels[level]->active_fe_indices);
        }
    }

                                     // count how much space we need
                                     // on each level and set the
                                     // dof_*_offsets
                                     // data. initially set the latter
                                     // to an invalid index, and only
                                     // later set it to something
                                     // reasonable for active cells
    for (unsigned int level=0; level<tria->n_levels(); ++level) 
      {
        levels[level]->line_dof_offsets
          = std::vector<unsigned int> (tria->n_raw_lines(level),
                                       invalid_dof_index);

        unsigned int next_free_line_dof = 0;
        for (active_cell_iterator cell=begin_active(level);
             cell!=end_active(level); ++cell)
          if (!cell->has_children())
          {
            levels[level]->line_dof_offsets[cell->index()] = next_free_line_dof;
            next_free_line_dof +=
              (*finite_elements)[cell->active_fe_index()].dofs_per_line;
          }

        levels[level]->line_dofs = std::vector<unsigned int> (next_free_line_dof,
                                                              invalid_dof_index);
      }


                                     // safety check: make sure that
                                     // the number of DoFs we
                                     // allocated is actually correct
                                     // (above we have also set the
                                     // dof_*_offsets field, so
                                     // we couldn't use this simpler
                                     // algorithm)
#ifdef DEBUG
    for (unsigned int level=0; level<tria->n_levels(); ++level)
      {
        unsigned int n_line_dofs = 0;
        for (cell_iterator cell=begin_active(level);
             cell!=end_active(level); ++cell)
          if (!cell->has_children())
            n_line_dofs +=
              (*finite_elements)[cell->active_fe_index()].dofs_per_line;
        
        Assert (levels[level]->line_dofs.size() == n_line_dofs, ExcInternalError());
        Assert (static_cast<unsigned int>
                (std::count (levels[level]->line_dof_offsets.begin(),
                             levels[level]->line_dof_offsets.end(),
                             invalid_dof_index))
                ==
                tria->n_raw_lines(level) - tria->n_active_lines(level),
                ExcInternalError());
      }
#endif
  }


#endif


#if deal_II_dimension == 2

  template <>
  void DoFHandler<2>::reserve_space ()
  {
    const unsigned int dim = 2;
    
    Assert (finite_elements != 0, ExcNoFESelected());
    Assert (finite_elements->size() > 0, ExcNoFESelected());
    Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
    Assert (tria->n_levels() == levels.size (), ExcInternalError ());

                                     // Release all space except the
                                     // active_fe_indices field which
                                     // we have to backup before
    {
      std::vector<std::vector<unsigned int> > active_fe_backup(levels.size ());
      for (unsigned int level = 0; level<levels.size (); ++level)
        std::swap (levels[level]->active_fe_indices, active_fe_backup[level]);

                                       // delete all levels and set them up
                                       // newly, since vectors are
                                       // troublesome if you want to change
                                       // their size
      clear_space ();
  
      for (unsigned int level=0; level<tria->n_levels(); ++level) 
        {
          levels.push_back (new internal::hp::DoFLevel<2>);
          std::swap (active_fe_backup[level],
                     levels[level]->active_fe_indices);
        }
    }

                                     // QUAD (CELL) DOFs
                                     //
                                     // count how much space we need
                                     // on each level for the quad
                                     // dofs and set the
                                     // dof_*_offsets
                                     // data. initially set the latter
                                     // to an invalid index, and only
                                     // later set it to something
                                     // reasonable for active cells
    for (unsigned int level=0; level<tria->n_levels(); ++level) 
      {
        levels[level]->quad_dof_offsets
          = std::vector<unsigned int> (tria->n_raw_quads(level),
                                       invalid_dof_index);

        unsigned int next_free_quad_dof = 0;
        for (active_cell_iterator cell=begin_active(level);
             cell!=end_active(level); ++cell)
          if (!cell->has_children())
          {
            levels[level]->quad_dof_offsets[cell->index()]
              = next_free_quad_dof;
            next_free_quad_dof
              += (*finite_elements)[cell->active_fe_index()].dofs_per_quad;
          }

        levels[level]->quad_dofs
          = std::vector<unsigned int> (next_free_quad_dof,
                                       invalid_dof_index);
      }

                                     // safety check: make sure that
                                     // the number of DoFs we
                                     // allocated is actually correct
                                     // (above we have also set the
                                     // dof_*_offsets field, so
                                     // we couldn't use this simpler
                                     // algorithm)
#ifdef DEBUG
    for (unsigned int level=0; level<tria->n_levels(); ++level)
      {
        unsigned int n_quad_dofs = 0;
        for (cell_iterator cell=begin_active(level);
             cell!=end_active(level); ++cell)
          if (!cell->has_children())
            n_quad_dofs +=
              (*finite_elements)[cell->active_fe_index()].dofs_per_quad;
        
        Assert (levels[level]->quad_dofs.size() == n_quad_dofs,
                ExcInternalError());
        Assert (static_cast<unsigned int>
                (std::count (levels[level]->quad_dof_offsets.begin(),
                             levels[level]->quad_dof_offsets.end(),
                             invalid_dof_index))
                ==
                tria->n_raw_quads(level) - tria->n_active_quads(level),
                ExcInternalError());
      }
#endif

    
                                     // LINE DOFS
                                     //
                                     // same here: count line dofs,
                                     // then allocate as much space as
                                     // we need and prime the linked
                                     // list for lines (see the
                                     // description in hp::DoFLevels)
                                     // with the indices we will
                                     // need. note that our task is
                                     // more complicated since two
                                     // adjacent cells may have
                                     // different active_fe_indices,
                                     // in which case we need to
                                     // allocate *two* sets of line
                                     // dofs for the same line
                                     //
                                     // the way we do things is that
                                     // we loop over all active cells
                                     // (these are the ones that have
                                     // DoFs only anyway) and all
                                     // their faces. We note in the
                                     // user flags whether we have
                                     // previously visited a face and
                                     // if so skip it (consequently,
                                     // we have to save and later
                                     // restore the line flags)
    {
      std::vector<bool> saved_line_user_flags;
      const_cast<Triangulation<dim>&>(*tria)
        .save_user_flags_line (saved_line_user_flags);
      const_cast<Triangulation<dim>&>(*tria).clear_user_flags_line ();

                                       // an array to hold how many
                                       // slots (see the hp::DoFLevel
                                       // class) we will have to store
                                       // on each level
      std::vector<unsigned int> n_line_slots (tria->n_levels(), 0);
    
      for (active_cell_iterator cell=begin_active(); cell!=end(); ++cell)
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (! cell->face(face)->user_flag_set())
            {
                                               // ok, face has not been
                                               // visited. so we need to
                                               // allocate space for
                                               // it. let's see how much
                                               // we need: we need one
                                               // set if a) there is no
                                               // neighbor behind this
                                               // face, or b) the
                                               // neighbor is not on the
                                               // same level or further
                                               // refined, or c) the
                                               // neighbor is on the
                                               // same level, but
                                               // happens to have the
                                               // same active_fe_index:
              if (cell->at_boundary(face)
                  ||
                  (cell->neighbor(face)->level() < cell->level())
                  ||
                  cell->neighbor(face)->has_children()
                  ||
                  ((cell->neighbor(face)->level() == cell->level())
                   &&
                   !cell->neighbor(face)->has_children()
                   &&
                   (cell->active_fe_index() == cell->neighbor(face)->active_fe_index())))
                                                 // ok, one set of
                                                 // dofs. that makes
                                                 // one index, 1 times
                                                 // dofs_per_line
                                                 // dofs, and one stop
                                                 // index
                n_line_slots[cell->level()]
                  += (*finite_elements)[cell->active_fe_index()].dofs_per_line + 2;

                                               // otherwise we do
                                               // indeed need two
                                               // sets, i.e. two
                                               // indices, two sets of
                                               // dofs, and one stop
                                               // index:
              else
                n_line_slots[cell->level()]
                  += ((*finite_elements)[cell->active_fe_index()].dofs_per_line
                      +
                      (*finite_elements)[cell->neighbor(face)->active_fe_index()]
                      .dofs_per_line
                      +
                      3);
            
                                               // mark this face as
                                               // visited
              cell->face(face)->set_user_flag ();
            }

                                       // now that we know how many
                                       // line dofs we will have to
                                       // have on each level, allocate
                                       // the memory. note that we
                                       // allocate offsets for all
                                       // lines, though only the
                                       // active ones will have a
                                       // non-invalid value later on
      for (unsigned int level=0; level<tria->n_levels(); ++level) 
        {
          levels[level]->line_dof_offsets
            = std::vector<unsigned int> (tria->n_raw_lines(level),
                                         invalid_dof_index);
          levels[level]->line_dofs
            = std::vector<unsigned int> (n_line_slots[level],
                                         invalid_dof_index);        
        }

                                       // with the memory now
                                       // allocated, loop over the
                                       // cells again and prime the
                                       // _offset values as well as
                                       // the fe_index fields
      const_cast<Triangulation<dim>&>(*tria).clear_user_flags_line ();

      std::vector<unsigned int> next_free_line_slot (tria->n_levels(), 0);
      
      for (active_cell_iterator cell=begin_active(); cell!=end(); ++cell)
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (! cell->face(face)->user_flag_set())
            {
                                               // same decision tree
                                               // as before
              if (cell->at_boundary(face)
                  ||
                  (cell->neighbor(face)->level() < cell->level())
                  ||
                  cell->neighbor(face)->has_children()
                  ||
                  ((cell->neighbor(face)->level() == cell->level())
                   &&
                   !cell->neighbor(face)->has_children()
                   &&
                   (cell->active_fe_index() == cell->neighbor(face)->active_fe_index())))
                {
                  levels[cell->level()]
                    ->line_dof_offsets[cell->face(face)->index()]
                    = next_free_line_slot[cell->level()];

                                                   // set first slot
                                                   // for this line to
                                                   // active_fe_index
                                                   // of this face
                  levels[cell->level()]
                    ->line_dofs[next_free_line_slot[cell->level()]]
                    = cell->active_fe_index();

                                                   // the next
                                                   // dofs_per_line
                                                   // indices remain
                                                   // unset for the
                                                   // moment (i.e. at
                                                   // invalid_dof_index).
                                                   // following this
                                                   // comes the stop
                                                   // index, which
                                                   // also is
                                                   // invalid_dof_index
                                                   // and therefore
                                                   // does not have to
                                                   // be explicitly
                                                   // set
                  
                                                   // finally, mark
                                                   // those slots as
                                                   // used
                  next_free_line_slot[cell->level()]
                    += (*finite_elements)[cell->active_fe_index()].dofs_per_line + 2;
                }
              else
                {
                  levels[cell->level()]
                    ->line_dof_offsets[cell->face(face)->index()]
                    = next_free_line_slot[cell->level()];

                                                   // set first slot
                                                   // for this line to
                                                   // active_fe_index
                                                   // of this face
                  levels[cell->level()]
                    ->line_dofs[next_free_line_slot[cell->level()]]
                    = cell->active_fe_index();

                                                   // the next
                                                   // dofs_per_line
                                                   // indices remain
                                                   // unset for the
                                                   // moment (i.e. at
                                                   // invalid_dof_index).
                                                   // 
                                                   // then comes the
                                                   // fe_index for the
                                                   // neighboring
                                                   // cell:
                  levels[cell->level()]
                    ->line_dofs[next_free_line_slot[cell->level()]
                                +
                                (*finite_elements)[cell->active_fe_index()].dofs_per_line
                                +
                                1]
                    = cell->neighbor(face)->active_fe_index();
                                                   // then again a set
                                                   // of dofs that we
                                                   // need not set
                                                   // right now
                                                   // 
                                                   // following this
                                                   // comes the stop
                                                   // index, which
                                                   // also is
                                                   // invalid_dof_index
                                                   // and therefore
                                                   // does not have to
                                                   // be explicitly
                                                   // set
                  
                                                   // finally, mark
                                                   // those slots as
                                                   // used
                  next_free_line_slot[cell->level()]
                    += ((*finite_elements)[cell->active_fe_index()].dofs_per_line
                        +
                        (*finite_elements)[cell->neighbor(face)->active_fe_index()]
                        .dofs_per_line
                        +
                        3);
                }
              
                                               // mark this face as
                                               // visited
              cell->face(face)->set_user_flag ();
            }

                                       // we should have moved the
                                       // cursor for each level to the
                                       // total number of dofs on that
                                       // level. check that
      for (unsigned int level=0; level<tria->n_levels(); ++level)
        Assert (next_free_line_slot[level] == n_line_slots[level],
                ExcInternalError());
      
                                       // at the end, restore the user
                                       // flags for the lines
      const_cast<Triangulation<dim>&>(*tria)
        .load_user_flags_line (saved_line_user_flags);
    }
    

				     // VERTEX DOFS
    
				     // The final step is to allocate
				     // vertex dof information. since
				     // vertices are sequentially
				     // numbered, what we do first is
				     // to set up an array in which we
				     // record whether a vertex is
				     // associated with any of the
				     // given fe's, by setting a
				     // bit. in a later step, we then
				     // actually allocate memory for
				     // the required dofs
    {
      std::vector<std::vector<bool> >
	vertex_fe_association (finite_elements->size(),
			       std::vector<bool> (tria->n_vertices(), false));
      
      for (active_cell_iterator cell=begin_active(); cell!=end(); ++cell)
	for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	  vertex_fe_association[cell->active_fe_index()][cell->vertex_index(v)]
	    = true;
      
				       // in debug mode, make sure
				       // that each vertex is
				       // associated with at least one
				       // fe (note that except for
				       // unused vertices, all
				       // vertices are actually
				       // active)
#ifdef DEBUG
      for (unsigned int v=0; v<tria->n_vertices(); ++v)
	if (tria->vertex_used(v) == true)
	  {
	    unsigned int fe=0;
	    for (; fe<finite_elements->size(); ++fe)
	      if (vertex_fe_association[fe][v] == true)
		break;
	    Assert (fe != finite_elements->size(), ExcInternalError());
	  }
#endif

				       // next count how much memory
				       // we actually need. for each
				       // vertex, we need one slot per
				       // fe to store the fe_index,
				       // plus dofs_per_vertex for
				       // this fe. in addition, we
				       // need one slot as the end
				       // marker for the
				       // fe_indices. at the same time
				       // already fill the
				       // vertex_dofs_offsets field
      vertex_dofs_offsets.resize (tria->n_vertices(),
				  deal_II_numbers::invalid_unsigned_int);
      
      unsigned int vertex_slots_needed = 0;
      for (unsigned int v=0; v<tria->n_vertices(); ++v)
	if (tria->vertex_used(v) == true)
	  {
	    vertex_dofs_offsets[v] = vertex_slots_needed;
	    
	    for (unsigned int fe=0; fe<finite_elements->size(); ++fe)
	      if (vertex_fe_association[fe][v] == true)
		vertex_slots_needed += (*finite_elements)[fe].dofs_per_vertex + 1;
	    ++vertex_slots_needed;
	  }

				       // now allocate the space we
				       // have determined we need, and
				       // set up the linked lists for
				       // each of the vertices
      vertex_dofs.resize (vertex_slots_needed, invalid_dof_index);
      for (unsigned int v=0; v<tria->n_vertices(); ++v)
	if (tria->vertex_used(v) == true)
	  {
	    unsigned int pointer = vertex_dofs_offsets[v];
	    for (unsigned int fe=0; fe<finite_elements->size(); ++fe)
	      if (vertex_fe_association[fe][v] == true)
		{
						   // if this vertex
						   // uses this fe,
						   // then set the
						   // fe_index and
						   // move the pointer
						   // ahead
		  vertex_dofs[pointer] = fe;
		  pointer += (*finite_elements)[fe].dofs_per_vertex + 1;
		}
					     // finally place the end
					     // marker
	    vertex_dofs[pointer] = deal_II_numbers::invalid_unsigned_int;
	    ++pointer;
	  }
    }
    
  }

#endif


#if deal_II_dimension == 3

  template <>
  void DoFHandler<3>::reserve_space ()
  {
    Assert (finite_elements != 0, ExcNoFESelected());
    Assert (finite_elements->size() > 0, ExcNoFESelected());
    Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
    Assert (tria->n_levels() == levels.size (), ExcInternalError ());

    Assert (finite_elements->max_dofs_per_face() == 0,
            ExcMessage ("hp finite elements are presently only supported "
                        "for discontinuous elements"));

                                     // Release all space except the
                                     // active_fe_indices field which
                                     // we have to backup before
    {
      std::vector<std::vector<unsigned int> > active_fe_backup(levels.size ());
      for (unsigned int level = 0; level<levels.size (); ++level)
        std::swap (levels[level]->active_fe_indices, active_fe_backup[level]);

                                       // delete all levels and set them up
                                       // newly, since vectors are
                                       // troublesome if you want to change
                                       // their size
      clear_space ();
  
      for (unsigned int level=0; level<tria->n_levels(); ++level) 
        {
          levels.push_back (new internal::hp::DoFLevel<3>);
          std::swap (active_fe_backup[level], levels[level]->active_fe_indices);
        }
    }

                                     // count how much space we need
                                     // on each level and set the
                                     // dof_*_offsets
                                     // data. initially set the latter
                                     // to an invalid index, and only
                                     // later set it to something
                                     // reasonable for active cells
    for (unsigned int level=0; level<tria->n_levels(); ++level) 
      {
        levels[level]->hex_dof_offsets
          = std::vector<unsigned int> (tria->n_raw_hexs(level),
                                       invalid_dof_index);

        unsigned int next_free_hex_dof = 0;
        for (active_cell_iterator cell=begin_active(level);
             cell!=end_active(level); ++cell)
          if (!cell->has_children())
          {
            levels[level]->hex_dof_offsets[cell->index()] = next_free_hex_dof;
            next_free_hex_dof +=
              (*finite_elements)[cell->active_fe_index()].dofs_per_hex;
          }

        levels[level]->hex_dofs = std::vector<unsigned int> (next_free_hex_dof,
                                                              invalid_dof_index);
      }


                                     // safety check: make sure that
                                     // the number of DoFs we
                                     // allocated is actually correct
                                     // (above we have also set the
                                     // dof_*_offsets field, so
                                     // we couldn't use this simpler
                                     // algorithm)
#ifdef DEBUG
    for (unsigned int level=0; level<tria->n_levels(); ++level)
      {
        unsigned int n_hex_dofs = 0;
        for (cell_iterator cell=begin_active(level);
             cell!=end_active(level); ++cell)
          if (!cell->has_children())
            n_hex_dofs +=
              (*finite_elements)[cell->active_fe_index()].dofs_per_hex;
        
        Assert (levels[level]->hex_dofs.size() == n_hex_dofs, ExcInternalError());
        Assert (static_cast<unsigned int>
                (std::count (levels[level]->hex_dof_offsets.begin(),
                             levels[level]->hex_dof_offsets.end(),
                             invalid_dof_index))
                ==
                tria->n_raw_hexs(level) - tria->n_active_hexs(level),
                ExcInternalError());
      }
#endif
  }

#endif



  template <int dim>
  void DoFHandler<dim>::create_active_fe_table ()
  {
                                     // Create sufficiently many
                                     // hp::DoFLevels.
    while (levels.size () < tria->n_levels ())
      levels.push_back (new internal::hp::DoFLevel<dim>);

                                     // then make sure that on each
                                     // level we have the appropriate
                                     // size of active_fe_indices;
                                     // preset them to zero, i.e. the
                                     // default FE
    for (unsigned int level=0; level<levels.size(); ++level)
      {
        if (levels[level]->active_fe_indices.size () == 0)
          levels[level]->active_fe_indices.resize (tria->n_raw_cells(level),
                                                   0);
        else
          {
                                             // Either the
                                             // active_fe_indices have
                                             // size zero because they
                                             // were just created, or
                                             // the correct
                                             // size. Other sizes
                                             // indicate that
                                             // something went wrong.
            Assert (levels[level]->active_fe_indices.size () ==
                    tria->n_raw_cells(level),
                    ExcInternalError ());
          }
      }
  }


  template <int dim>
  void DoFHandler<dim>::pre_refinement_notification (const Triangulation<dim> &)
  {
    create_active_fe_table ();
  }


#if deal_II_dimension == 1
  template <>
  void DoFHandler<1>::pre_refinement_notification (const Triangulation<1> &tria)
  {
    create_active_fe_table ();

                                     // Remember if the cells have already
                                     // children. That will make the transfer
                                     // of the active_fe_index to the finer
                                     // levels easier.
    Assert (has_children.size () == 0, ExcInternalError ());
    for (unsigned int i=0; i<levels.size(); ++i)
      {
	const unsigned int lines_on_level = tria.n_raw_lines(i);
	std::vector<bool> *has_children_level =
          new std::vector<bool> (lines_on_level);

                                         // Check for each cell, if it has children.
	transform (tria.levels[i]->lines.children.begin (),
		   tria.levels[i]->lines.children.end (),
		   has_children_level->begin (),
		   std::bind2nd (std::not_equal_to<int>(), -1));

	has_children.push_back (has_children_level);
      }
  }
#endif

  
#if deal_II_dimension == 2
  template <>
  void DoFHandler<2>::pre_refinement_notification (const Triangulation<2> &tria)
  {
    create_active_fe_table ();

                                     // Remember if the cells have already
                                     // children. That will make the transfer
                                     // of the active_fe_index to the finer
                                     // levels easier.
    Assert (has_children.size () == 0, ExcInternalError ());
    for (unsigned int i=0; i<levels.size(); ++i)
      {
	const unsigned int quads_on_level = tria.n_raw_quads (i);
	std::vector<bool> *has_children_level =
          new std::vector<bool> (quads_on_level);

                                         // Check for each cell, if it has children.
	std::transform (tria.levels[i]->quads.children.begin (),
                        tria.levels[i]->quads.children.end (),
                        has_children_level->begin (),
                        std::bind2nd (std::not_equal_to<int>(), -1));

	has_children.push_back (has_children_level);
      }
  }
#endif

#if deal_II_dimension == 3
  template <>
  void DoFHandler<3>::pre_refinement_notification (const Triangulation<3> &tria)
  {
    create_active_fe_table ();

                                     // Remember if the cells have already
                                     // children. That will make the transfer
                                     // of the active_fe_index to the finer
                                     // levels easier.
    Assert (has_children.size () == 0, ExcInternalError ());
    for (unsigned int i=0; i<levels.size(); ++i)
      {
	const unsigned int hexes_on_level = tria.n_raw_hexs(i);
	std::vector<bool> *has_children_level =
          new std::vector<bool> (hexes_on_level);

                                         // Check for each cell, if it
                                         // has children.
	transform (tria.levels[i]->hexes.children.begin (),
		   tria.levels[i]->hexes.children.end (),
		   has_children_level->begin (),
		   std::bind2nd (std::not_equal_to<int>(), -1));

	has_children.push_back (has_children_level);
      }
  }
#endif

  template <int dim>
  void DoFHandler<dim>::post_refinement_notification (const Triangulation<dim> &tria)
  {
    Assert (has_children.size () == levels.size (), ExcInternalError ());

                                     // In each refinement at most one new
                                     // new level may appear. If that happened
                                     // it is appended to the DoFHandler
                                     // levels.
    if (levels.size () < tria.n_levels ())
      levels.push_back (new internal::hp::DoFLevel<dim>);

                                     // Coarsening can lead to the loss
                                     // of levels. Hence remove them.
    while (levels.size () > tria.n_levels ())
      { 
	delete levels[levels.size ()-1];
	levels.pop_back ();
      }

                                     // Resize active_fe_indices
                                     // vectors. use zero indicator to
                                     // extend
    for (unsigned int i=0; i<levels.size(); ++i)
      levels[i]->active_fe_indices.resize (tria.n_raw_cells(i), 0);

                                     // if a finite element collection
                                     // has already been set, then
                                     // actually try to set
                                     // active_fe_indices for child
                                     // cells of refined cells to the
                                     // active_fe_index of the mother
                                     // cell. if no finite element
                                     // collection has been assigned
                                     // yet, then all indicators are
                                     // zero anyway, and there is no
                                     // point trying to set anything
                                     // (besides, we would trip over
                                     // an assertion in
                                     // set_active_fe_index)
    if (finite_elements != 0)
      {
        cell_iterator cell = begin(),
                      endc = end ();
        for (; cell != endc; ++cell)
          {
                                             // Look if the cell got children during
                                             // refinement
                                             // Note: Although one level is added to
                                             // the DoFHandler levels, when the
                                             // triangulation got one, for the buffer
                                             // has_children this new level is not
                                             // required, because the cells on the
                                             // finest level never have children. Hence
                                             // cell->has_children () will always return
                                             // false on that level, which would cause
                                             // shortcut evaluation of the following
                                             // expression. Thus an index error in
                                             // has_children should never occur.
            if (cell->has_children () &&
                !(*has_children [cell->level ()])[cell->index ()])
              {
                                                 // Set active_fe_index in children to the
                                                 // same value as in the parent cell.
                for (unsigned int i = 0; i < cell->n_children(); ++i)
                  cell->child (i)->set_active_fe_index (cell->active_fe_index ());
              }
          }
      }
    
                                     // Free buffer objects
    std::vector<std::vector<bool> *>::iterator children_level;
    for (children_level = has_children.begin ();
	 children_level != has_children.end ();
	 ++children_level)
      delete (*children_level);
                                     /*
                                       for_each (has_children.begin (),
                                       has_children.end (),
                                       delete());
                                     */
    has_children.clear ();
  }


  template <int dim>
  void DoFHandler<dim>::clear_space ()
  {  
    for (unsigned int i=0; i<levels.size(); ++i)
      delete levels[i];
    levels.resize (0);

    {
      std::vector<unsigned int> tmp;
      std::swap (vertex_dofs, tmp);
    }
    
    {
      std::vector<unsigned int> tmp;
      std::swap (vertex_dofs_offsets, tmp);
    }
  }


/*-------------- Explicit Instantiations -------------------------------*/
  template class DoFHandler<deal_II_dimension>;
}
