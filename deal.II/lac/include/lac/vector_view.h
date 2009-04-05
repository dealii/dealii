//---------------------------------------------------------------------------
//    $Id: vector.h 16829 2008-09-15 18:11:22Z heltai $
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __dealii__vector_view_h
#define __dealii__vector_view_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <lac/vector.h>

#include <cstdio>

DEAL_II_NAMESPACE_OPEN


/*! @addtogroup Vectors
 *@{
 */

/**
 * View of a numerical vector of data. This class provides an
 * interface compatible with the Vector<double> class (from which it
 * is inherited), that allows fast access to locations of memory
 * already allocated with arrays of type Number.
 *
 * This is in the same style of the vector view in the Trilinos
 * library.
 *
 * Notice that NO CHECKS are performed on the actual memory, and if
 * you try to access illegal areas of memory, your computer will
 * suffer from it. Use this class ONLY if you know exactly what you
 * are doing. Two constructors are provided. One for read-write
 * access, and one for read only access. 
 *
 * You are allowed to use this class on objects of type const
 * Vector<Number>, however you should be aware of the fact that the
 * constness of pointed to array is casted away, which means that you
 * should only use the const constructor when the actual object you
 * are constructing is itself a constant object.
 *
 * You WILL be allowed to change the vector afterwards, so make sure
 * you know what you are doing, before you change data without
 * realizing.
 *
 * Since this class does not own the memory that you are accessing,
 * you have to make sure that the life of the section of memory you
 * are viewing is longer than this object. No attempt is made to
 * ensure that this is the case.
 *
 * An example usage of this class is the following:
 *
 <code>
 // Create an array of length 5;
 double * array = new double[5];
 // Now create a view of the above array that is compatible with the
 // Vector<double> class
 VectorView<double> view(5, array);

 view(1) = 4;
 
 // The following line should output 4.
 cout << array[1] << endl;

 // If debug mode is on, then the following triggers an execption:
 view(6) = 4;
 
 // But notice that no checks are performed, so this is legal but WILL
 // NOT work
 VectorView<double> wrong_view(10, array);
 
 // Now no assert will be thrown if you type wrong_view(6), but most
 // likely a seg fault will occur.
 view(6) = 4;

 // Notice that this construction is legal. It will create a copy of
 // the array.
 const Vector<double> const_copy(view);
 
 // Now this is the correct way to instantiate a constant view of the
 // above vector:
 const VectorView<double> correct_const_copy_view(5, const_copy.begin());
 
 // While this will compile, BUT WILL NOT COMPLAIN if you try to write
 // on it!
 VectorView<double> wrong_const_copy_view(5, const_copy.begin());

 // Now writing to elements of wrong_const_copy_view is allowed, and
 // will change the same memory as the const_copy object.
 wrong_const_copy_view(1) = 5;
 
 if(copy_view(1) == wrong_const_copy_view(1)) cout << "Tautology";
 
 </code>
 *
 *
 * @note Instantiations for this template are provided for
 * <tt>@<float@>, @<double@>, @<long double@>,
 * @<std::complex@<float@>@>, @<std::complex@<double@>@>,
 * @<std::complex@<long double@>@></tt>; others can be generated in
 * application programs (see the section on @ref Instantiations in the
 * manual).
 * 
 * @author Luca Heltai
 */
template<typename Number>
class VectorView : public Vector<Number> {
public:    
    /** Read write constructor. Takes the size of the vector, just
     * like the standard one, but the data is picked starting from the
     * location of the pointer @p ptr.
     */
    VectorView(const unsigned int new_size, Number *ptr);
    
    /** The constant constructor is the same as above, however you
     * will not be able to access the data for write access.
     *
     * You should only use this class by constructing it as a const
     * VectorView<double>(size, ptr) object.
     *
     * Undefined behavior will occur if you construct it as a non
     * const object or attempt to write on it.
     */
    VectorView(const unsigned int new_size, const Number *ptr);
    
    /** This desctructor will only reset the internal sizes and the
     * interanl pointers,  but it will NOT clear the memory. */
    ~VectorView();
};



/*@}*/
/*----------------------- Inline functions ----------------------------------*/



template<typename Number>
inline
VectorView<Number>::VectorView(const unsigned int new_size, Number * ptr)
{
    this->vec_size	= new_size;
    this->max_vec_size	= new_size;
    this->val		= ptr;
}



template<typename Number>
inline
VectorView<Number>::VectorView(const unsigned int new_size, const Number * ptr) 
{
    this->vec_size	= new_size;
    this->max_vec_size	= new_size;
    this->val		= const_cast<Number*>(ptr);
}



template<typename Number>
inline
VectorView<Number>::~VectorView() {
    this->vec_size = 0;
    this->max_vec_size = 0;
    this->val = 0;
}

DEAL_II_NAMESPACE_CLOSE

#endif
