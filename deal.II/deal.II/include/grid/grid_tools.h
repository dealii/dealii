//----------------------------  grid_tools.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_tools.h  ---------------------------
#ifndef __deal2__grid_tools_H
#define __deal2__grid_tools_H


template <int dim> class Triangulation;



/**
 * This class is a collection of algorithms working on
 * triangulations. See the descriptions of the individual functions
 * for more information.
 *
 * @author Wolfgang Bangerth, 2001
 */
class GridTools
{
  public:
				     /**
				      * Return the diameter of a
				      * triangulation. The diameter is
				      * computed using only the
				      * vertices, i.e. if the diameter
				      * should be larger than the
				      * maximal distance between
				      * boundary vertices due to a
				      * higher order mapping, then
				      * this function will not catch
				      * this.
				      */
    template <int dim>
    static
    double diameter (const Triangulation<dim> &tria);

				     /**
				      * Same function, but for 1d.
				      */
    static
    double diameter (const Triangulation<1> &tria);
};


/*----------------------------   grid_tools.h     ---------------------------*/
/* end of #ifndef __deal2__grid_tools_H */
#endif
/*----------------------------   grid_tools.h     ---------------------------*/
