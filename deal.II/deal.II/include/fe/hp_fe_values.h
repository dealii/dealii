//----------------------------  hp_fe_values.h  ---------------------------
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
//----------------------------  hp_fe_values.h  ---------------------------
#ifndef __deal2__hp_fe_values_h
#define __deal2__hp_fe_values_h

#include <base/config.h>
#include <fe/fe.h>
#include <fe/fe_collection.h>
#include <fe/q_collection.h>
#include <fe/mapping_collection.h>
#include <fe/fe_values.h>

#include <map>
#include <boost/shared_ptr.hpp>

template <int dim> class MappingQ1;
template <int dim> class FiniteElement;


namespace internal
{
  namespace hp
  {
/**
 * Base class for the <tt>hp::FE*Values</tt> classes, storing the data
 * that is common to them. The main task of this class is to provide a
 * table where for every combination of finite element, mapping, and
 * quadrature object from their corresponding collection objects there
 * is a matching ::FEValues, ::FEFaceValues, or ::FESubfaceValues
 * object. To make things more efficient, however, these FE*Values
 * objects are only created once requested (lazy allocation).
 *
 * The first template parameter denotes the space dimension we are in,
 * the second the dimensionality of the object that we integrate on,
 * i.e. for usual @p hp::FEValues it is equal to the first one, while
 * for face integration it is one less. The third template parameter
 * indicates the type of underlying non-hp FE*Values base type,
 * i.e. it could either be ::FEValues, ::FEFaceValues, or
 * ::FESubfaceValues.
 * 
 * @ingroup hp
 *
 * @author Wolfgang Bangerth, 2003
 */
    template <int dim, int q_dim, class FEValues>
    class FEValuesBase 
    {
      public:        
                                         /**
                                          * Constructor. Set the fields
                                          * of this class to the values
                                          * indicated by the parameters
                                          * to the constructor.
                                          */
        FEValuesBase (const ::hp::MappingCollection<dim> &mapping_collection,
                      const ::hp::FECollection<dim>      &fe_collection,
                      const ::hp::QCollection<q_dim>     &q_collection,
                      const UpdateFlags             update_flags);
                                         /**
                                          * Constructor. Set the fields
                                          * of this class to the values
                                          * indicated by the parameters
                                          * to the constructor, and
                                          * choose a @p MappingQ1
                                          * object for the mapping
                                          * object.
                                          */
        FEValuesBase (const ::hp::FECollection<dim> &fe_collection,
                      const ::hp::QCollection<q_dim> &q_collection,
                      const UpdateFlags         update_flags);
        
                                         /**
                                          * Return a reference to the
                                          * @p FEValues object
                                          * selected by the last call
                                          * to
                                          * select_fe_values(). select_fe_values()
                                          * in turn is called when you
                                          * called the @p reinit
                                          * function of the
                                          * <tt>hp::FE*Values</tt>
                                          * class the last time.
                                          */
        const FEValues & get_present_fe_values () const;

      protected:

                                         /**
                                          * Select a FEValues object
                                          * suitable for the given FE,
                                          * quadrature, and mapping
                                          * indices. If such an object
                                          * doesn't yet exist, create
                                          * one.
                                          *
                                          * The function returns a
                                          * writable reference so that
                                          * derived classes can also
                                          * reinit() the selected
                                          * FEValues object.
                                          */
        FEValues &
        select_fe_values (const unsigned int fe_index,
                          const unsigned int mapping_index,
                          const unsigned int q_index);

      protected:
                                         /**
                                          * A pointer to the
                                          * collection of finite
                                          * elements to be used.
                                          */
        const SmartPointer<const ::hp::FECollection<dim> > fe_collection;
        
                                         /**
                                          * A pointer to the
                                          * collection of mappings to
                                          * be used.
                                          */
        const SmartPointer<const ::hp::MappingCollection<dim> > mapping_collection;
    
                                         /**
                                          * Copy of the quadrature
                                          * collection object
                                          * provided to the
                                          * constructor.
                                          */
        const ::hp::QCollection<q_dim> q_collection;

      private:
                                         /**
                                          * A table in which we store
                                          * pointers to fe_values
                                          * objects for different
                                          * finite element, mapping,
                                          * and quadrature objects
                                          * from our collection. The
                                          * first index indicates the
                                          * index of the finite
                                          * element within the
                                          * fe_collection, the second
                                          * the index of the mapping
                                          * within the mapping
                                          * collection, and the last
                                          * one the index of the
                                          * quadrature formula within
                                          * the q_collection.
                                          *
                                          * Initially, all entries
                                          * have zero pointers, and we
                                          * will allocate them lazily
                                          * as needed in
                                          * select_fe_values().
                                          */
        Table<3,boost::shared_ptr<FEValues> > fe_values_table;

                                         /**
                                          * Set of indices pointing at
                                          * the fe_values object
                                          * selected last time the
                                          * select_fe_value() function
                                          * was called.
                                          */
        TableIndices<3> present_fe_values_index;
        
                                         /**
                                          * Values of the update flags as
                                          * given to the constructor.
                                          */
        const UpdateFlags update_flags;
    };

  }
  
}


namespace hp
{
  
/**
 * 
 * @ingroup hp
 */  
  template <int dim>
  class FEValues : public internal::hp::FEValuesBase<dim,dim,::FEValues<dim> >
  {
    public:
                                       /**
                                        * Constructor. Initialize this
                                        * object with the given
                                        * parameters.
                                        *
                                        * The finite element
                                        * collection parameter is
                                        * actually ignored, but is in
                                        * the signature of this
                                        * function to make it
                                        * compatible with the
                                        * signature of the respective
                                        * constructor of the usual
                                        * FEValues object, with
                                        * the respective parameter in
                                        * that function also being the
                                        * return value of the
                                        * <tt>DoFHandler::get_fe()</tt>
                                        * function.
                                        */
      FEValues (const ::hp::MappingCollection<dim> &mapping_collection,
                const ::hp::FECollection<dim>  &fe_collection,
                const ::hp::QCollection<dim>       &q_collection,
                const UpdateFlags             update_flags);


                                       /**
                                        * Constructor. Initialize this
                                        * object with the given
                                        * parameters, and choose a
                                        * @p MappingQ1 object for the
                                        * mapping object.
                                        *
                                        * The finite element
                                        * collection parameter is
                                        * actually ignored, but is in
                                        * the signature of this
                                        * function to make it
                                        * compatible with the
                                        * signature of the respective
                                        * constructor of the usual
                                        * FEValues object, with
                                        * the respective parameter in
                                        * that function also being the
                                        * return value of the
                                        * <tt>DoFHandler::get_fe()</tt>
                                        * function.
                                        */
      FEValues (const hp::FECollection<dim> &fe_collection,
                const hp::QCollection<dim>      &q_collection,
                const UpdateFlags            update_flags);


                                       /**
                                        * Reinitialize the object for
                                        * the given cell.
                                        *
                                        * After the call, you can get
                                        * an FEValues object using the
                                        * get_present_fe_values()
                                        * function that corresponds to
                                        * the present cell. For this
                                        * FEValues object, we use the
                                        * additional arguments
                                        * described below to determine
                                        * which finite element,
                                        * mapping, and quadrature
                                        * formula to use. They are
                                        * order in such a way that the
                                        * arguments one may want to
                                        * change most frequently come
                                        * first. The rules for these
                                        * arguments are as follows:
                                        *
                                        * If the @p fe_index argument
                                        * to this function is left at
                                        * its default value, then we
                                        * use that finite element
                                        * within the hp::FECollection
                                        * passed to the constructor of
                                        * this class with index given
                                        * by
                                        * <code>cell-@>active_fe_index()</code>. Consequently,
                                        * the hp::FECollection
                                        * argument given to this
                                        * object should really be the
                                        * same as that used in the
                                        * construction of the
                                        * hp::DofHandler associated
                                        * with the present cell. On
                                        * the other hand, if a value
                                        * is given for this argument,
                                        * it overrides the choice of
                                        * <code>cell-@>active_fe_index()</code>.
                                        *
                                        * If the @p q_index argument
                                        * is left at its default
                                        * value, then we use that
                                        * quadrature formula within
                                        * the hp::QCollection passed
                                        * to the constructor of this
                                        * class with index given by
                                        * <code>cell-@>active_fe_index()</code>,
                                        * i.e. the same index as that
                                        * of the finite element. In
                                        * this case, there should be a
                                        * corresponding quadrature
                                        * formula for each finite
                                        * element in the
                                        * hp::FECollection. As a
                                        * special case, if the
                                        * quadrature collection
                                        * contains only a single
                                        * element (a frequent case if
                                        * one wants to use the same
                                        * quadrature object for all
                                        * finite elements in an hp
                                        * discretization, even if that
                                        * may not be the most
                                        * efficient), then this single
                                        * quadrature is used unless a
                                        * different value for this
                                        * argument is specified. On
                                        * the other hand, if a value
                                        * is given for this argument,
                                        * it overrides the choice of
                                        * <code>cell-@>active_fe_index()</code>
                                        * or the choice for the single
                                        * quadrature.
                                        *
                                        * If the @p mapping_index
                                        * argument is left at its
                                        * default value, then we use
                                        * that mapping object within
                                        * the hp::MappingCollection
                                        * passed to the constructor of
                                        * this class with index given
                                        * by
                                        * <code>cell-@>active_fe_index()</code>,
                                        * i.e. the same index as that
                                        * of the finite
                                        * element. As above, if the
                                        * mapping collection contains
                                        * only a single element (a
                                        * frequent case if one wants
                                        * to use a MappingQ1 object
                                        * for all finite elements in
                                        * an hp discretization), then
                                        * this single mapping is used
                                        * unless a different value for
                                        * this argument is specified.
                                        */
      void
      reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell,
              const unsigned int q_index = deal_II_numbers::invalid_unsigned_int,
              const unsigned int mapping_index = deal_II_numbers::invalid_unsigned_int,
              const unsigned int fe_index = deal_II_numbers::invalid_unsigned_int);
  };



/**
 * 
 * @ingroup hp
 */  
  template <int dim>
  class FEFaceValues : public internal::hp::FEValuesBase<dim,dim-1,::FEFaceValues<dim> >
  {
    public:
                                       /**
                                        * Constructor. Initialize this
                                        * object with the given
                                        * parameters.
                                        *
                                        * The finite element
                                        * collection parameter is
                                        * actually ignored, but is in
                                        * the signature of this
                                        * function to make it
                                        * compatible with the
                                        * signature of the respective
                                        * constructor of the usual
                                        * FEValues object, with
                                        * the respective parameter in
                                        * that function also being the
                                        * return value of the
                                        * <tt>DoFHandler::get_fe()</tt>
                                        * function.
                                        */
      FEFaceValues (const hp::MappingCollection<dim> &mapping_collection,
                    const hp::FECollection<dim>  &fe_collection,
                    const hp::QCollection<dim-1>     &q_collection,
                    const UpdateFlags             update_flags);


                                       /**
                                        * Constructor. Initialize this
                                        * object with the given
                                        * parameters, and choose a
                                        * @p MappingQ1 object for the
                                        * mapping object.
                                        *
                                        * The finite element
                                        * collection parameter is
                                        * actually ignored, but is in
                                        * the signature of this
                                        * function to make it
                                        * compatible with the
                                        * signature of the respective
                                        * constructor of the usual
                                        * FEValues object, with
                                        * the respective parameter in
                                        * that function also being the
                                        * return value of the
                                        * <tt>DoFHandler::get_fe()</tt>
                                        * function.
                                        */
      FEFaceValues (const hp::FECollection<dim>  &fe_collection,
                    const hp::QCollection<dim-1> &q_collection,
                    const UpdateFlags             update_flags);

                                       /**
                                        * Reinitialize the object for
                                        * the given cell and face.
                                        *
                                        * After the call, you can get
                                        * an FEFaceValues object using the
                                        * get_present_fe_values()
                                        * function that corresponds to
                                        * the present cell. For this
                                        * FEFaceValues object, we use the
                                        * additional arguments
                                        * described below to determine
                                        * which finite element,
                                        * mapping, and quadrature
                                        * formula to use. They are
                                        * order in such a way that the
                                        * arguments one may want to
                                        * change most frequently come
                                        * first. The rules for these
                                        * arguments are as follows:
                                        *
                                        * If the @p fe_index argument
                                        * to this function is left at
                                        * its default value, then we
                                        * use that finite element
                                        * within the hp::FECollection
                                        * passed to the constructor of
                                        * this class with index given
                                        * by
                                        * <code>cell-@>active_fe_index()</code>. Consequently,
                                        * the hp::FECollection
                                        * argument given to this
                                        * object should really be the
                                        * same as that used in the
                                        * construction of the
                                        * hp::DofHandler associated
                                        * with the present cell. On
                                        * the other hand, if a value
                                        * is given for this argument,
                                        * it overrides the choice of
                                        * <code>cell-@>active_fe_index()</code>.
                                        *
                                        * If the @p q_index argument
                                        * is left at its default
                                        * value, then we use that
                                        * quadrature formula within
                                        * the hp::QCollection passed
                                        * to the constructor of this
                                        * class with index given by
                                        * <code>cell-@>active_fe_index()</code>,
                                        * i.e. the same index as that
                                        * of the finite element. In
                                        * this case, there should be a
                                        * corresponding quadrature
                                        * formula for each finite
                                        * element in the
                                        * hp::FECollection. As a
                                        * special case, if the
                                        * quadrature collection
                                        * contains only a single
                                        * element (a frequent case if
                                        * one wants to use the same
                                        * quadrature object for all
                                        * finite elements in an hp
                                        * discretization, even if that
                                        * may not be the most
                                        * efficient), then this single
                                        * quadrature is used unless a
                                        * different value for this
                                        * argument is specified. On
                                        * the other hand, if a value
                                        * is given for this argument,
                                        * it overrides the choice of
                                        * <code>cell-@>active_fe_index()</code>
                                        * or the choice for the single
                                        * quadrature.
                                        *
                                        * If the @p mapping_index
                                        * argument is left at its
                                        * default value, then we use
                                        * that mapping object within
                                        * the hp::MappingCollection
                                        * passed to the constructor of
                                        * this class with index given
                                        * by
                                        * <code>cell-@>active_fe_index()</code>,
                                        * i.e. the same index as that
                                        * of the finite
                                        * element. As above, if the
                                        * mapping collection contains
                                        * only a single element (a
                                        * frequent case if one wants
                                        * to use a MappingQ1 object
                                        * for all finite elements in
                                        * an hp discretization), then
                                        * this single mapping is used
                                        * unless a different value for
                                        * this argument is specified.
                                        */
      void
      reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell,
              const unsigned int face_no,
              const unsigned int q_index = deal_II_numbers::invalid_unsigned_int,
              const unsigned int mapping_index = deal_II_numbers::invalid_unsigned_int,
              const unsigned int fe_index = deal_II_numbers::invalid_unsigned_int);
  };



/**
 * 
 * @ingroup hp
 */  
  template <int dim>
  class FESubfaceValues : public internal::hp::FEValuesBase<dim,dim-1,::FESubfaceValues<dim> >
  {
    public:
                                       /**
                                        * Constructor. Initialize this
                                        * object with the given
                                        * parameters.
                                        *
                                        * The finite element
                                        * collection parameter is
                                        * actually ignored, but is in
                                        * the signature of this
                                        * function to make it
                                        * compatible with the
                                        * signature of the respective
                                        * constructor of the usual
                                        * FEValues object, with
                                        * the respective parameter in
                                        * that function also being the
                                        * return value of the
                                        * <tt>DoFHandler::get_fe()</tt>
                                        * function.
                                        */
      FESubfaceValues (const hp::MappingCollection<dim> &mapping_collection,
                       const hp::FECollection<dim>  &fe_collection,
                       const hp::QCollection<dim-1>     &q_collection,
                       const UpdateFlags             update_flags);


                                       /**
                                        * Constructor. Initialize this
                                        * object with the given
                                        * parameters, and choose a
                                        * @p MappingQ1 object for the
                                        * mapping object.
                                        *
                                        * The finite element
                                        * collection parameter is
                                        * actually ignored, but is in
                                        * the signature of this
                                        * function to make it
                                        * compatible with the
                                        * signature of the respective
                                        * constructor of the usual
                                        * FEValues object, with
                                        * the respective parameter in
                                        * that function also being the
                                        * return value of the
                                        * <tt>DoFHandler::get_fe()</tt>
                                        * function.
                                        */
      FESubfaceValues (const hp::FECollection<dim> &fe_collection,
                       const hp::QCollection<dim-1>    &q_collection,
                       const UpdateFlags            update_flags);

                                       /**
                                        * Reinitialize the object for
                                        * the given cell, face, and subface.
                                        *
                                        * After the call, you can get
                                        * an FESubfaceValues object using the
                                        * get_present_fe_values()
                                        * function that corresponds to
                                        * the present cell. For this
                                        * FESubfaceValues object, we use the
                                        * additional arguments
                                        * described below to determine
                                        * which finite element,
                                        * mapping, and quadrature
                                        * formula to use. They are
                                        * order in such a way that the
                                        * arguments one may want to
                                        * change most frequently come
                                        * first. The rules for these
                                        * arguments are as follows:
                                        *
                                        * If the @p q_index argument
                                        * is left at its default
                                        * value, then we use that
                                        * quadrature formula within
                                        * the hp::QCollection passed
                                        * to the constructor of this
                                        * class with index given by
                                        * <code>cell-@>active_fe_index()</code>,
                                        * i.e. the same index as that
                                        * of the finite element. In
                                        * this case, there should be a
                                        * corresponding quadrature
                                        * formula for each finite
                                        * element in the
                                        * hp::FECollection. As a
                                        * special case, if the
                                        * quadrature collection
                                        * contains only a single
                                        * element (a frequent case if
                                        * one wants to use the same
                                        * quadrature object for all
                                        * finite elements in an hp
                                        * discretization, even if that
                                        * may not be the most
                                        * efficient), then this single
                                        * quadrature is used unless a
                                        * different value for this
                                        * argument is specified. On
                                        * the other hand, if a value
                                        * is given for this argument,
                                        * it overrides the choice of
                                        * <code>cell-@>active_fe_index()</code>
                                        * or the choice for the single
                                        * quadrature.
                                        *
                                        * If the @p mapping_index
                                        * argument is left at its
                                        * default value, then we use
                                        * that mapping object within
                                        * the hp::MappingCollection
                                        * passed to the constructor of
                                        * this class with index given
                                        * by
                                        * <code>cell-@>active_fe_index()</code>,
                                        * i.e. the same index as that
                                        * of the finite
                                        * element. As above, if the
                                        * mapping collection contains
                                        * only a single element (a
                                        * frequent case if one wants
                                        * to use a MappingQ1 object
                                        * for all finite elements in
                                        * an hp discretization), then
                                        * this single mapping is used
                                        * unless a different value for
                                        * this argument is specified.
                                        */
      void
      reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell,
              const unsigned int face_no,
              const unsigned int subface_no,
              const unsigned int q_index = deal_II_numbers::invalid_unsigned_int,
              const unsigned int mapping_index = deal_II_numbers::invalid_unsigned_int,
              const unsigned int fe_index = deal_II_numbers::invalid_unsigned_int);
  };
  
}


// -------------- inline and template functions --------------

namespace internal 
{
  namespace hp
  {
    template <int dim, int q_dim, class FEValues>
    inline
    const FEValues &
    FEValuesBase<dim,q_dim,FEValues>::get_present_fe_values () const
    {
      return *fe_values_table(present_fe_values_index);
    }
  }
  
}


#endif
