//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__tria_object_h
#define __deal2__tria_object_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/geometry_info.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace Triangulation
  {
    
/**
 * Class template for the <tt>structdim</tt>-dimensional cells constituting
 * a Triangulation of dimension <tt>structdim</tt> or lower dimensional
 * objects of higher dimensions.  They are characterized by the
 * (global) indices of their faces, which are cells of dimension
 * <tt>structdim-1</tt> or vertices if <tt>structdim=1</tt>.
 *
 * @author Guido Kanschat, 2007
 */
    template <int structdim>
    class TriaObject
    {
      public:
	static const unsigned int dimension = structdim;
	
					 /**
					  * Default constructor,
					  * setting all face indices
					  * to invalid values.
					  */				  
	TriaObject ();
	
					 /**
					  * Constructor for a line
					  * object with the numbers of
					  * its two end points.
					  *
					  * Throws an exception if
					  * dimension is not one.
					  */
	TriaObject (const int i0, const int i1);

					 /**
					  * Constructor for a quadrilateral
					  * object with the numbers of
					  * its four lines.
					  *
					  * Throws an exception if
					  * dimension is not two.
					  */
	TriaObject (const int i0, const int i1,
		    const int i2, const int i3);

					 /**
					  * Constructor for a hexahedron
					  * object with the numbers of
					  * its six quadrilaterals.
					  *
					  * Throws an exception if
					  * dimension is not two.
					  */
	TriaObject (const int i0, const int i1,
                    const int i2, const int i3,
                    const int i4, const int i5);
    

					 /**
					  * Return the index of the
					  * ith face object.
					  */
	int face (const unsigned int i) const;

					 /**
					  * Set the index of the ith
					  * face object.
					  */
	void set_face (const unsigned int i, const int index);

					 /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        static unsigned int memory_consumption ();
	
      protected:
                                         /**
                                          *  Global indices of the
                                          *  face iterators bounding
                                          *  this cell if dim@>1, and
                                          *  the two vertex indices in
                                          *  1d.
                                          */
        int faces[GeometryInfo<structdim>::faces_per_cell];
    };

//----------------------------------------------------------------------//

    template <int structdim>
    inline
    TriaObject<structdim>::TriaObject ()
    {
      for (unsigned int i=0; i<GeometryInfo<structdim>::faces_per_cell; ++i)
	faces[i] = -1;
    }


    template <int structdim>
    inline
    TriaObject<structdim>::TriaObject (const int i0,
				       const int i1)
    {
      Assert (structdim==1, ExcImpossibleInDim(structdim));
      faces[0] = i0;
      faces[1] = i1;
    }


    template <int structdim>
    inline
    TriaObject<structdim>::TriaObject (const int i0,
				       const int i1,
				       const int i2,
				       const int i3)
    {
      Assert (structdim==2, ExcImpossibleInDim(structdim));
      faces[0] = i0;
      faces[1] = i1;
      faces[2] = i2;
      faces[3] = i3;
    }


    template <int structdim>
    inline
    TriaObject<structdim>::TriaObject (const int i0,
				       const int i1,
				       const int i2,
				       const int i3,
				       const int i4,
				       const int i5)
    {
      Assert (structdim==3, ExcImpossibleInDim(structdim));
      faces[0] = i0;
      faces[1] = i1;
      faces[2] = i2;
      faces[3] = i3;
      faces[4] = i4;
      faces[5] = i5;
    }


    template <int structdim>
    inline
    int TriaObject<structdim>::face (const unsigned int i) const
    {
      Assert (i<GeometryInfo<structdim>::faces_per_cell,
              ExcIndexRange(i,0,GeometryInfo<structdim>::faces_per_cell));
      return faces[i];
    }
    
    
    
    template <int structdim>
    inline
    void TriaObject<structdim>::set_face (const unsigned int i, const int index)
    {
      Assert (i<GeometryInfo<structdim>::faces_per_cell,
              ExcIndexRange(i,0,GeometryInfo<structdim>::faces_per_cell));
      faces[i] = index;
    }
    
    
    
    template <int structdim>
    inline
    unsigned int
    TriaObject<structdim>::memory_consumption ()
    {
      return sizeof(TriaObject<structdim>);
    }
    
    
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
