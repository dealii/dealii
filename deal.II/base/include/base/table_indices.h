//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__table_indices_h
#define __deal2__table_indices_h


#include <base/config.h>
#include <base/exceptions.h>



/**
 * Base class for an array of indices of fixed size used for the TableBase and
 * SymmetricTensor classes. Actually, this class serves a dual purpose, as it
 * not only stores indices into the TableBase class, but also the sizes of the
 * table in its various coordinates.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <int N>
class TableIndicesBase
{
  public:
                                     /**
                                      * Access the value of the
                                      * <tt>i</tt>th index.
                                      */
    unsigned int operator[] (const unsigned int i) const;

                                     /**
                                      * Compare two index fields for
                                      * equality.
                                      */
    bool operator == (const TableIndicesBase<N> &other) const;

                                     /**
                                      * Compare two index fields for
                                      * inequality.
                                      */
    bool operator != (const TableIndicesBase<N> &other) const;

                                     /**
                                      * Sort the indices in ascending
                                      * order. While this operation is not
                                      * very useful for Table objects, it is
                                      * used for the SymmetricTensor class.
                                      */
    void sort ();
      
  protected:
                                     /**
                                      * Store the indices in an array.
                                      */
    unsigned indices[N];
};


/**
 * Array of indices of fixed size used for the TableBase
 * class.
 *
 * This is the general template, and has no implementation. There are
 * a number of specializations that are actually implemented (one for
 * each used value of <tt>N</tt>), which only differ in the way they
 * implement their constructors (they take <tt>N</tt> arguments, something
 * that cannot be represented by a general template). Actual storage
 * of and access to data is done by the TableIndicesBase base
 * class of a specializations.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <int N>
class TableIndices
{
};



/**
 * Array of indices of fixed size used for the TableBase
 * class.
 *
 * This is the specialization for a one-dimensional table, i.e. a
 * vector. This class only differs in the non-default constructors
 * from the other specializations. Actual storage of and access to
 * data is done by the TableIndicesBase base class of a
 * specializations.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <>
class TableIndices<1> : public TableIndicesBase<1>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1);
};



/**
 * Array of indices of fixed size used for the TableBase
 * class.
 *
 * This is the specialization for a two-dimensional table. This class
 * only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the TableIndicesBase base class of a specializations.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <>
class TableIndices<2> : public TableIndicesBase<2>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2);
};



/**
 * Array of indices of fixed size used for the TableBase
 * class.
 *
 * This is the specialization for a three-dimensional table. This class
 * only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the TableIndicesBase base class of a specializations.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <>
class TableIndices<3> : public TableIndicesBase<3>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2,
                  const unsigned int index3);
};


/**
 * Array of indices of fixed size used for the TableBase
 * class.
 *
 * This is the specialization for a four-dimensional table. This class
 * only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the TableIndicesBase base class of a specializations.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <>
class TableIndices<4> : public TableIndicesBase<4>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2,
                  const unsigned int index3,
                  const unsigned int index4);
};


/**
 * Array of indices of fixed size used for the TableBase
 * class.
 *
 * This is the specialization for a five-dimensional table. This class
 * only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the TableIndicesBase base class of a specializations.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <>
class TableIndices<5> : public TableIndicesBase<5>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2,
                  const unsigned int index3,
                  const unsigned int index4,
                  const unsigned int index5);
};


/**
 * Array of indices of fixed size used for the TableBase
 * class.
 *
 * This is the specialization for a six-dimensional table. This class
 * only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the TableIndicesBase base class of a specializations.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <>
class TableIndices<6> : public TableIndicesBase<6>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2,
                  const unsigned int index3,
                  const unsigned int index4,
                  const unsigned int index5,
                  const unsigned int index6);
};


/**
 * Array of indices of fixed size used for the TableBase
 * class.
 *
 * This is the specialization for a seven-dimensional table. This
 * class only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the TableIndicesBase base class of a specializations.
 *
 * @author Wolfgang Bangerth, 2002, Ralf Hartmann 2004
 */
template <>
class TableIndices<7> : public TableIndicesBase<7>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2,
                  const unsigned int index3,
                  const unsigned int index4,
                  const unsigned int index5,
                  const unsigned int index6,
                  const unsigned int index7);
};


/* --------------------- Template and inline functions ---------------- */


template <int N>
inline
unsigned int
TableIndicesBase<N>::operator [] (const unsigned int i) const 
{
  Assert (i < N, ExcIndexRange (i, 0, N));
  return indices[i];
}



template <int N>
inline
bool
TableIndicesBase<N>::operator == (const TableIndicesBase<N> &other) const 
{
  for (unsigned int i=0; i<N; ++i)
    if (indices[i] != other.indices[i])
      return false;
  return true;
}



template <int N>
inline
bool
TableIndicesBase<N>::operator != (const TableIndicesBase<N> &other) const 
{
  return !(*this == other);
}



template <>
inline
void
TableIndicesBase<1>::sort ()
{}



template <>
inline
void
TableIndicesBase<2>::sort ()
{
  if (indices[1] < indices[0])
    std::swap (indices[0], indices[1]);
}



template <>
inline
void
TableIndicesBase<3>::sort ()
{
                                   // bubble sort for 3 elements
  if (indices[1] < indices[0])
    std::swap (indices[0], indices[1]);
  if (indices[2] < indices[1])
    std::swap (indices[1], indices[2]);
  if (indices[1] < indices[0])
    std::swap (indices[0], indices[1]);
}




inline
TableIndices<1>::TableIndices () 
{
  this->indices[0] = 0;
}



inline
TableIndices<1>::TableIndices (const unsigned int index1)
{
  this->indices[0] = index1;
}



inline
TableIndices<2>::TableIndices () 
{
  this->indices[0] = this->indices[1] = 0;
}



inline
TableIndices<2>::TableIndices (const unsigned int index1,
                               const unsigned int index2)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
}



inline
TableIndices<3>::TableIndices () 
{
  this->indices[0] = this->indices[1] = this->indices[2] = 0;
}



inline
TableIndices<3>::TableIndices (const unsigned int index1,
                               const unsigned int index2,
                               const unsigned int index3)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
  this->indices[2] = index3;
}


inline
TableIndices<4>::TableIndices ()
{
  this->indices[0] = this->indices[1] = this->indices[2] = this->indices[3] = 0;
}



inline
TableIndices<4>::TableIndices (const unsigned int index1,
                               const unsigned int index2,
                               const unsigned int index3,
			       const unsigned int index4)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
  this->indices[2] = index3;
  this->indices[3] = index4;
}



inline
TableIndices<5>::TableIndices ()
{
  this->indices[0] = this->indices[1]
                   = this->indices[2]
                   = this->indices[3]
                   = this->indices[4] = 0;
}



inline
TableIndices<5>::TableIndices (const unsigned int index1,
                               const unsigned int index2,
                               const unsigned int index3,
			       const unsigned int index4,
			       const unsigned int index5)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
  this->indices[2] = index3;
  this->indices[3] = index4;
  this->indices[4] = index5;
}



inline
TableIndices<6>::TableIndices ()
{
  this->indices[0] = this->indices[1] = this->indices[2]
		   = this->indices[3] = this->indices[4]
		   = this->indices[5] = 0;
}



inline
TableIndices<6>::TableIndices (const unsigned int index1,
                               const unsigned int index2,
                               const unsigned int index3,
			       const unsigned int index4,
			       const unsigned int index5,
			       const unsigned int index6)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
  this->indices[2] = index3;
  this->indices[3] = index4;
  this->indices[4] = index5;
  this->indices[5] = index6;
}



inline
TableIndices<7>::TableIndices ()
{
  this->indices[0] = this->indices[1] = this->indices[2]
		   = this->indices[3] = this->indices[4]
		   = this->indices[5] = this->indices[6] = 0;
}



inline
TableIndices<7>::TableIndices (const unsigned int index1,
                               const unsigned int index2,
                               const unsigned int index3,
			       const unsigned int index4,
			       const unsigned int index5,
			       const unsigned int index6,
			       const unsigned int index7)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
  this->indices[2] = index3;
  this->indices[3] = index4;
  this->indices[4] = index5;
  this->indices[5] = index6;
  this->indices[6] = index7;
}




#endif
