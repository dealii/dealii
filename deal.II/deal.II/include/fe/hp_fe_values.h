//----------------------------  hp_fe_values.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004 by the deal.II authors
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
#include <fe/fe_values.h>

#include <map>
#include <boost/shared_ptr.hpp>

template <int dim> class MappingQ1;
template <int dim> class FiniteElement;


namespace internal 
{
/**
 * Map between finite element objects and @p FEValues (the second
 * template parameter, which might be <tt>FEValues<dim></tt>,
 * <tt>FEFaceValues<dim></tt>, ...) objects. The <tt>hpFE*Values</tt> classes
 * use this to hold an <tt>FE*Values</tt> object for each finite element
 * that is used in the triangulation that it integrates on.
 */
  template <int dim, class FEValues>
  class FEValuesMap
  {
    public:
                                       /**
                                        * Make destructor virtual,
                                        * since we have virtual
                                        * functions in this class.
                                        */
      virtual ~FEValuesMap ();
      
                                       /**
                                        * Return a reference to the
                                        * @p FEValues object selected
                                        * by the last call to
                                        * @p select_fe_values. The
                                        * returned value is a constant
                                        * reference since the only
                                        * non-const function in the
                                        * underlying class would be
                                        * @p reinit, which you must
                                        * not call directly any way;
                                        * rather, use the @p reinit
                                        * function of the
                                        * <tt>hpFE*Values</tt> class.
                                        */
      const FEValues & get_present_fe_values () const;

    protected:
                                       /**
                                        * Create an object of type
                                        * @p FEValues for this
                                        * particular finite
                                        * element. Since the type of
                                        * @p FEValues is not known
                                        * here, we have to leave this
                                        * to derived classes.
                                        */
      virtual
      FEValues * create_fe_values (const FiniteElement<dim> &fe) const = 0;
      
                                       /**
                                        * Select the @p FEValues
                                        * object corresponding to the
                                        * finite element object given
                                        * as argument. If there is no
                                        * such @p FEValues object
                                        * yet, then one is created by
                                        * calling @p create_fe_values
                                        * with this finite element
                                        * object.
                                        *
                                        * A non-const reference to
                                        * this object is returned.
                                        */
      FEValues &
      select_fe_values (const FiniteElement<dim> &fe);
      
    private:
                                       /**
                                        * Have a map from pointers to
                                        * finite elements to objects
                                        * of type @p FEValues.
                                        *
                                        * Note that the shared pointer
                                        * will make sure that the
                                        * memory previously allocated
                                        * is released when the
                                        * destructor of this class is
                                        * run.
                                        */
      std::map<SmartPointer<const FiniteElement<dim> >,
               boost::shared_ptr<FEValues> > fe_to_fe_values_map;

                                       /**
                                        * Pointer to the @p FEValues
                                        * object that is used for the
                                        * present cell. This always
                                        * points to one of the objects
                                        * in the map above, and to
                                        * which one is determined by
                                        * the last call to
                                        * @p select_fe_values.
                                        */
      boost::shared_ptr<FEValues> present_fe_values;
  };


/**
 * Base class for the <tt>hpFE*Values</tt> classes, storing the data that
 * is common to them. The first template parameter denotes the space
 * dimension we are in, the second the dimensionality of the object
 * that we integrate on, i.e. for usual @p hpFEValues it is equal to
 * the first one, while for face integration it is one less.
 *
 * @author Wolfgang Bangerth, 2003
 */
  template <int dim, int q_dim>
  class hpFEValuesBase 
  {
    public:
                                       /**
                                        * Constructor. Set the fields
                                        * of this class to the values
                                        * indicated by the parameters
                                        * to the constructor.
                                        */
      hpFEValuesBase (const Mapping<dim>      &mapping,
                      const Quadrature<q_dim> &quadrature,
                      const UpdateFlags        update_flags);

                                       /**
                                        * Constructor. Set the fields
                                        * of this class to the values
                                        * indicated by the parameters
                                        * to the constructor, and
                                        * choose a @p MappingQ1
                                        * object for the mapping
                                        * object.
                                        */
      hpFEValuesBase (const Quadrature<q_dim> &quadrature,
                      const UpdateFlags        update_flags);
      
    protected:
                                       /**
                                        * Pointer to a mapping object as
                                        * provided to the constructor,
                                        * or a default @p MappingQ1
                                        * object if none was given.
                                        */
      const SmartPointer<const Mapping<dim> > mapping;
    
                                       /**
                                        * Pointer to a quadrature object
                                        * as provided to the
                                        * constructor.
                                        */
      const SmartPointer<const Quadrature<q_dim> > quadrature;

                                       /**
                                        * Values of the update flags as
                                        * given to the constructor.
                                        */
      const UpdateFlags update_flags;

    private:
                                       /**
                                        * Default mapping, in case
                                        * none was provided through
                                        * the call to the constructor.
                                        */
      static const MappingQ1<dim> default_mapping;
  };
}



template <int dim>
class hpFEValues : public internal::FEValuesMap<dim,FEValues<dim> >,
                   protected internal::hpFEValuesBase<dim,dim>
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
                                      * DoFHandler<tt>::get_fe()</tt>
                                      * function.
                                      */
    hpFEValues (const Mapping<dim>      &mapping,
                const FECollection<dim> &fe_collection,
                const Quadrature<dim>   &quadrature,
                const UpdateFlags        update_flags);

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
                                      * DoFHandler<tt>::get_fe()</tt>
                                      * function.
                                      */
    hpFEValues (const FECollection<dim> &fe_collection,
                const Quadrature<dim>   &quadrature,
                const UpdateFlags        update_flags);

                                     /**
                                      * Reinitialize the object for
                                      * the given cell. This selects
                                      * the right @p FEValues object
                                      * for the finite element in use
                                      * by the cell given, and calling
                                      * the @p reinit function on
                                      * this object.
                                      */
    void
    reinit (const typename hpDoFHandler<dim>::cell_iterator &cell);

  protected:
                                     /**
                                      * Create an object of type
                                      * @p FEValues for this
                                      * particular finite element.
                                      */
    virtual
    FEValues<dim> *
    create_fe_values (const FiniteElement<dim> &fe) const;
};



template <int dim>
class hpFEFaceValues : public internal::FEValuesMap<dim,FEFaceValues<dim> >,
                       protected internal::hpFEValuesBase<dim,dim-1>
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
                                      * DoFHandler<tt>::get_fe()</tt>
                                      * function.
                                      */
    hpFEFaceValues (const Mapping<dim>      &mapping,
                    const FECollection<dim> &fe_collection,
                    const Quadrature<dim-1> &quadrature,
                    const UpdateFlags        update_flags);

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
                                      * DoFHandler<tt>::get_fe()</tt>
                                      * function.
                                      */
    hpFEFaceValues (const FECollection<dim> &fe_collection,
                    const Quadrature<dim-1> &quadrature,
                    const UpdateFlags        update_flags);

                                     /**
                                      * Reinitialize the object for
                                      * the given cell. This selects
                                      * the right @p FEValues object
                                      * for the finite element in use
                                      * by the cell given, and calling
                                      * the @p reinit function on
                                      * this object.
                                      */
    void
    reinit (const typename hpDoFHandler<dim>::cell_iterator &cell,
            const unsigned int face_no);

  protected:
                                     /**
                                      * Create an object of type
                                      * @p FEValues for this
                                      * particular finite element.
                                      */
    virtual
    FEFaceValues<dim> *
    create_fe_values (const FiniteElement<dim> &fe) const;
};



template <int dim>
class hpFESubfaceValues : public internal::FEValuesMap<dim,FESubfaceValues<dim> >,
                          protected internal::hpFEValuesBase<dim,dim-1>
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
                                      * DoFHandler<tt>::get_fe()</tt>
                                      * function.
                                      */
    hpFESubfaceValues (const Mapping<dim>      &mapping,
                       const FECollection<dim> &fe_collection,
                       const Quadrature<dim-1> &quadrature,
                       const UpdateFlags        update_flags);

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
                                      * DoFHandler<tt>::get_fe()</tt>
                                      * function.
                                      */
    hpFESubfaceValues (const FECollection<dim> &fe_collection,
                       const Quadrature<dim-1> &quadrature,
                       const UpdateFlags        update_flags);

                                     /**
                                      * Reinitialize the object for
                                      * the given cell. This selects
                                      * the right @p FEValues object
                                      * for the finite element in use
                                      * by the cell given, and calling
                                      * the @p reinit function on
                                      * this object.
                                      */
    void
    reinit (const typename hpDoFHandler<dim>::cell_iterator &cell,
            const unsigned int face_no,
            const unsigned int subface_no);

  protected:
                                     /**
                                      * Create an object of type
                                      * @p FEValues for this
                                      * particular finite element.
                                      */
    virtual
    FESubfaceValues<dim> *
    create_fe_values (const FiniteElement<dim> &fe) const;
};


// -------------- inline and template functions --------------

namespace internal 
{
  template <int dim, class FEValues>
  inline
  const FEValues &
  FEValuesMap<dim,FEValues>::get_present_fe_values () const
  {
    return *present_fe_values;
  }
}




#endif
