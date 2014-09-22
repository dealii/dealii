// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__hp_fe_values_h
#define __deal2__hp_fe_values_h

#include <deal.II/base/config.h>
#include <deal.II/fe/fe.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/fe/fe_values.h>

#include <map>
#include <deal.II/base/std_cxx11/shared_ptr.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class MappingQ1;
template <int dim, int spacedim> class FiniteElement;



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
     * i.e. it could either be dealii::FEValues, dealii::FEFaceValues, or
     * dealii::FESubfaceValues.
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
       * Constructor. Set the fields of this class to the values indicated by
       * the parameters to the constructor.
       */
      FEValuesBase (const dealii::hp::MappingCollection<dim,FEValues::space_dimension> &mapping_collection,
                    const dealii::hp::FECollection<dim,FEValues::space_dimension>      &fe_collection,
                    const dealii::hp::QCollection<q_dim>     &q_collection,
                    const UpdateFlags             update_flags);
      /**
       * Constructor. Set the fields of this class to the values indicated by
       * the parameters to the constructor, and choose a @p MappingQ1 object
       * for the mapping object.
       */
      FEValuesBase (const dealii::hp::FECollection<dim,FEValues::space_dimension> &fe_collection,
                    const dealii::hp::QCollection<q_dim> &q_collection,
                    const UpdateFlags         update_flags);

      /**
       * Get a reference to the collection of finite element objects used
       * here.
       */
      const dealii::hp::FECollection<dim,FEValues::space_dimension> &
      get_fe_collection () const;

      /**
       * Get a reference to the collection of mapping objects used
       * here.
       */
      const dealii::hp::MappingCollection<dim,FEValues::space_dimension> &
      get_mapping_collection () const;

      /**
       * Get a reference to the collection of quadrature objects used
       * here.
       */
      const dealii::hp::QCollection<q_dim> &
      get_quadrature_collection () const;

      /**
       * Get the underlying update flags.
       */
      UpdateFlags get_update_flags() const;

      /**
       * Return a reference to the @p FEValues object selected by the last
       * call to select_fe_values(). select_fe_values() in turn is called when
       * you called the @p reinit function of the <tt>hp::FE*Values</tt> class
       * the last time.
       */
      const FEValues &get_present_fe_values () const;

    protected:

      /**
       * Select a FEValues object suitable for the given FE, quadrature, and
       * mapping indices. If such an object doesn't yet exist, create one.
       *
       * The function returns a writable reference so that derived classes can
       * also reinit() the selected FEValues object.
       */
      FEValues &
      select_fe_values (const unsigned int fe_index,
                        const unsigned int mapping_index,
                        const unsigned int q_index);

    protected:
      /**
       * A pointer to the collection of finite elements to be used.
       */
      const SmartPointer<const dealii::hp::FECollection<dim,FEValues::space_dimension>,
            FEValuesBase<dim,q_dim,FEValues> > fe_collection;

      /**
       * A pointer to the collection of mappings to be used.
       */
      const SmartPointer<const dealii::hp::MappingCollection<dim, FEValues::space_dimension>,
            FEValuesBase<dim,q_dim,FEValues> > mapping_collection;

      /**
       * Copy of the quadrature collection object provided to the constructor.
       */
      const dealii::hp::QCollection<q_dim> q_collection;

    private:
      /**
       * A table in which we store pointers to fe_values objects for different
       * finite element, mapping, and quadrature objects from our
       * collection. The first index indicates the index of the finite element
       * within the fe_collection, the second the index of the mapping within
       * the mapping collection, and the last one the index of the quadrature
       * formula within the q_collection.
       *
       * Initially, all entries have zero pointers, and we will allocate them
       * lazily as needed in select_fe_values().
       */
      dealii::Table<3,std_cxx11::shared_ptr<FEValues> > fe_values_table;

      /**
       * Set of indices pointing at the fe_values object selected last time
       * the select_fe_value() function was called.
       */
      TableIndices<3> present_fe_values_index;

      /**
       * Values of the update flags as given to the constructor.
       */
      const UpdateFlags update_flags;
    };

  }

}


namespace hp
{

  /**
   * An hp equivalent of the ::FEValues class. See the step-27
   * tutorial program for examples of use.
   *
   * The idea of this class is as follows: when one assembled matrices in the hp
   * finite element method, there may be different finite elements on different
   * cells, and consequently one may also want to use different quadrature
   * formulas for different cells. On the other hand, the ::FEValues efficiently
   * handles pre-evaluating whatever information is necessary for a single
   * finite element and quadrature object. This class brings these concepts
   * together: it provides a "collection" of ::FEValues objects.
   *
   * Upon construction, one passes not one finite element and quadrature object
   * (and possible a mapping), but a whole collection of type hp::FECollection
   * and hp::QCollection. Later on, when one sits on a concrete cell, one would
   * call the reinit() function for this particular cell, just as one does for a
   * regular ::FEValues object. The difference is that this time, the reinit()
   * function looks up the active_fe_index of that cell, if necessary creates a
   * ::FEValues object that matches the finite element and quadrature formulas
   * with that particular index in their collections, and then re-initializes it
   * for the current cell. The ::FEValues object that then fits the finite
   * element and quadrature formula for the current cell can then be accessed
   * using the get_present_fe_values() function, and one would work with it just
   * like with any ::FEValues object for non-hp DoF handler objects.
   *
   * The reinit() functions have additional arguments with default values. If
   * not specified, the function takes the index into the hp::FECollection,
   * hp::QCollection, and hp::MappingCollection objects from the active_fe_index
   * of the cell, as explained above. However, one can also select different
   * indices for a current cell. For example, by specifying a different index
   * into the hp::QCollection class, one does not need to sort the quadrature
   * objects in the quadrature collection so that they match one-to-one the
   * order of finite element objects in the FE collection (even though choosing
   * such an order is certainly convenient).
   *
   * Note that ::FEValues objects are created on the fly, i.e. only as they are
   * needed. This ensures that we do not create objects for every combination of
   * finite element, quadrature formula and mapping, but only those that will
   * actually be needed.
   *
   * This class has not yet been implemented for the use in the codimension
   * one case (<tt>spacedim != dim </tt>).
   *
   * @ingroup hp hpcollection
   * @author Wolfgang Bangerth, 2003
   */
  template <int dim, int spacedim=dim>
  class FEValues : public dealii::internal::hp::FEValuesBase<dim,dim,dealii::FEValues<dim,spacedim> >
  {
  public:

    static const unsigned int dimension = dim;

    static const unsigned int space_dimension = spacedim;

    /**
     * Constructor. Initialize this object with the given parameters.
     *
     * The finite element collection parameter is actually ignored, but is in
     * the signature of this function to make it compatible with the signature
     * of the respective constructor of the usual FEValues object, with the
     * respective parameter in that function also being the return value of
     * the <tt>DoFHandler::get_fe()</tt> function.
     */
    FEValues (const dealii::hp::MappingCollection<dim,spacedim> &mapping_collection,
              const dealii::hp::FECollection<dim,spacedim>  &fe_collection,
              const dealii::hp::QCollection<dim>       &q_collection,
              const UpdateFlags             update_flags);


    /**
     * Constructor. Initialize this object with the given parameters, and
     * choose a @p MappingQ1 object for the mapping object.
     *
     * The finite element collection parameter is actually ignored, but is in
     * the signature of this function to make it compatible with the signature
     * of the respective constructor of the usual FEValues object, with the
     * respective parameter in that function also being the return value of
     * the <tt>DoFHandler::get_fe()</tt> function.
     */
    FEValues (const hp::FECollection<dim,spacedim> &fe_collection,
              const hp::QCollection<dim>      &q_collection,
              const UpdateFlags            update_flags);


    /**
     * Reinitialize the object for the given cell.
     *
     * After the call, you can get an FEValues object using the
     * get_present_fe_values() function that corresponds to the present
     * cell. For this FEValues object, we use the additional arguments
     * described below to determine which finite element, mapping, and
     * quadrature formula to use. They are order in such a way that the
     * arguments one may want to change most frequently come first. The rules
     * for these arguments are as follows:
     *
     * If the @p fe_index argument to this function is left at its default
     * value, then we use that finite element within the hp::FECollection
     * passed to the constructor of this class with index given by
     * <code>cell-@>active_fe_index()</code>. Consequently, the
     * hp::FECollection argument given to this object should really be the
     * same as that used in the construction of the hp::DofHandler associated
     * with the present cell. On the other hand, if a value is given for this
     * argument, it overrides the choice of
     * <code>cell-@>active_fe_index()</code>.
     *
     * If the @p q_index argument is left at its default value, then we use
     * that quadrature formula within the hp::QCollection passed to the
     * constructor of this class with index given by
     * <code>cell-@>active_fe_index()</code>, i.e. the same index as that of
     * the finite element. In this case, there should be a corresponding
     * quadrature formula for each finite element in the hp::FECollection. As
     * a special case, if the quadrature collection contains only a single
     * element (a frequent case if one wants to use the same quadrature object
     * for all finite elements in an hp discretization, even if that may not
     * be the most efficient), then this single quadrature is used unless a
     * different value for this argument is specified. On the other hand, if a
     * value is given for this argument, it overrides the choice of
     * <code>cell-@>active_fe_index()</code> or the choice for the single
     * quadrature.
     *
     * If the @p mapping_index argument is left at its default value, then we
     * use that mapping object within the hp::MappingCollection passed to the
     * constructor of this class with index given by
     * <code>cell-@>active_fe_index()</code>, i.e. the same index as that of
     * the finite element. As above, if the mapping collection contains only a
     * single element (a frequent case if one wants to use a MappingQ1 object
     * for all finite elements in an hp discretization), then this single
     * mapping is used unless a different value for this argument is
     * specified.
     */
    template <class DH, bool lda>
    void
    reinit (const TriaIterator<DoFCellAccessor<DH,lda> > cell,
            const unsigned int q_index = numbers::invalid_unsigned_int,
            const unsigned int mapping_index = numbers::invalid_unsigned_int,
            const unsigned int fe_index = numbers::invalid_unsigned_int);

    /**
     * Like the previous function, but for non-hp iterators. The reason this
     * (and the other non-hp iterator) function exists is so that one can use
     * hp::FEValues not only for hp::DoFhandler objects, but for all sorts of
     * DoFHandler objects, and triangulations not associated with DoFHandlers
     * in general.
     *
     * Since <code>cell-@>active_fe_index()</code> doesn't make sense for
     * triangulation iterators, this function chooses the zero-th finite
     * element, mapping, and quadrature object from the relevant constructions
     * passed to the constructor of this object. The only exception is if you
     * specify a value different from the default value for any of these last
     * three arguments.
     */
    void
    reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
            const unsigned int q_index = numbers::invalid_unsigned_int,
            const unsigned int mapping_index = numbers::invalid_unsigned_int,
            const unsigned int fe_index = numbers::invalid_unsigned_int);


  };



  /**
   * This is the equivalent of the hp::FEValues class but for face integrations,
   * i.e. it is to hp::FEValues what ::FEFaceValues is to ::FEValues.
   *
   * The same comments apply as in the documentation of the hp::FEValues
   * class. However, it is important to note that it is here more common that
   * one would want to explicitly specify an index to a particular quadrature
   * formula in the reinit() functions. This is because the default index
   * corresponds to the finite element index on the current function. On the
   * other hand, integration on faces will typically have to happen with a
   * quadrature formula that is adjusted to the finite elements used on both
   * sides of a face. If one sorts the elements of the hp::FECollection with
   * ascending polynomial degree, and matches these finite elements with
   * corresponding quadrature formulas in the hp::QCollection passed to the
   * constructor, then the quadrature index passed to the reinit() function
   * should typically be something like <code>std::max
   * (cell-@>active_fe_index(), neighbor-@>active_fe_index()</code> to ensure
   * that a quadrature formula is chosen that is sufficiently accurate for
   * <em>both</em> finite elements.
   *
   * @ingroup hp hpcollection
   * @author Wolfgang Bangerth, 2003
   */
  template <int dim, int spacedim=dim>
  class FEFaceValues : public dealii::internal::hp::FEValuesBase<dim,dim-1,dealii::FEFaceValues<dim,spacedim> >
  {
  public:
    /**
     * Constructor. Initialize this object with the given parameters.
     *
     * The finite element collection parameter is actually ignored, but is in
     * the signature of this function to make it compatible with the signature
     * of the respective constructor of the usual FEValues object, with the
     * respective parameter in that function also being the return value of
     * the <tt>DoFHandler::get_fe()</tt> function.
     */
    FEFaceValues (const hp::MappingCollection<dim,spacedim> &mapping_collection,
                  const hp::FECollection<dim,spacedim>  &fe_collection,
                  const hp::QCollection<dim-1>     &q_collection,
                  const UpdateFlags             update_flags);


    /**
     * Constructor. Initialize this object with the given parameters, and
     * choose a @p MappingQ1 object for the mapping object.
     *
     * The finite element collection parameter is actually ignored, but is in
     * the signature of this function to make it compatible with the signature
     * of the respective constructor of the usual FEValues object, with the
     * respective parameter in that function also being the return value of
     * the <tt>DoFHandler::get_fe()</tt> function.
     */
    FEFaceValues (const hp::FECollection<dim,spacedim>  &fe_collection,
                  const hp::QCollection<dim-1> &q_collection,
                  const UpdateFlags             update_flags);

    /**
     * Reinitialize the object for the given cell and face.
     *
     * After the call, you can get an FEFaceValues object using the
     * get_present_fe_values() function that corresponds to the present
     * cell. For this FEFaceValues object, we use the additional arguments
     * described below to determine which finite element, mapping, and
     * quadrature formula to use. They are order in such a way that the
     * arguments one may want to change most frequently come first. The rules
     * for these arguments are as follows:
     *
     * If the @p fe_index argument to this function is left at its default
     * value, then we use that finite element within the hp::FECollection
     * passed to the constructor of this class with index given by
     * <code>cell-@>active_fe_index()</code>. Consequently, the
     * hp::FECollection argument given to this object should really be the
     * same as that used in the construction of the hp::DofHandler associated
     * with the present cell. On the other hand, if a value is given for this
     * argument, it overrides the choice of
     * <code>cell-@>active_fe_index()</code>.
     *
     * If the @p q_index argument is left at its default value, then we use
     * that quadrature formula within the hp::QCollection passed to the
     * constructor of this class with index given by
     * <code>cell-@>active_fe_index()</code>, i.e. the same index as that of
     * the finite element. In this case, there should be a corresponding
     * quadrature formula for each finite element in the hp::FECollection. As
     * a special case, if the quadrature collection contains only a single
     * element (a frequent case if one wants to use the same quadrature object
     * for all finite elements in an hp discretization, even if that may not
     * be the most efficient), then this single quadrature is used unless a
     * different value for this argument is specified. On the other hand, if a
     * value is given for this argument, it overrides the choice of
     * <code>cell-@>active_fe_index()</code> or the choice for the single
     * quadrature.
     *
     * If the @p mapping_index argument is left at its default value, then we
     * use that mapping object within the hp::MappingCollection passed to the
     * constructor of this class with index given by
     * <code>cell-@>active_fe_index()</code>, i.e. the same index as that of
     * the finite element. As above, if the mapping collection contains only a
     * single element (a frequent case if one wants to use a MappingQ1 object
     * for all finite elements in an hp discretization), then this single
     * mapping is used unless a different value for this argument is
     * specified.
     */
    template <class DH, bool lda>
    void
    reinit (const TriaIterator<DoFCellAccessor<DH,lda> > cell,
            const unsigned int face_no,
            const unsigned int q_index = numbers::invalid_unsigned_int,
            const unsigned int mapping_index = numbers::invalid_unsigned_int,
            const unsigned int fe_index = numbers::invalid_unsigned_int);

    /**
     * Like the previous function, but for non-hp iterators. The reason this
     * (and the other non-hp iterator) function exists is so that one can use
     * hp::FEValues not only for hp::DoFhandler objects, but for all sorts of
     * DoFHandler objects, and triangulations not associated with DoFHandlers
     * in general.
     *
     * Since <code>cell-@>active_fe_index()</code> doesn't make sense for
     * triangulation iterators, this function chooses the zero-th finite
     * element, mapping, and quadrature object from the relevant constructions
     * passed to the constructor of this object. The only exception is if you
     * specify a value different from the default value for any of these last
     * three arguments.
     */
    void
    reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
            const unsigned int face_no,
            const unsigned int q_index = numbers::invalid_unsigned_int,
            const unsigned int mapping_index = numbers::invalid_unsigned_int,
            const unsigned int fe_index = numbers::invalid_unsigned_int);
  };



  /**
   * This class implements for subfaces what hp::FEFaceValues does for
   * faces. See there for further documentation.
   *
   * @ingroup hp hpcollection
   * @author Wolfgang Bangerth, 2003
   */
  template <int dim, int spacedim=dim>
  class FESubfaceValues : public dealii::internal::hp::FEValuesBase<dim,dim-1,dealii::FESubfaceValues<dim,spacedim> >
  {
  public:
    /**
     * Constructor. Initialize this object with the given parameters.
     *
     * The finite element collection parameter is actually ignored, but is in
     * the signature of this function to make it compatible with the signature
     * of the respective constructor of the usual FEValues object, with the
     * respective parameter in that function also being the return value of
     * the <tt>DoFHandler::get_fe()</tt> function.
     */
    FESubfaceValues (const hp::MappingCollection<dim,spacedim> &mapping_collection,
                     const hp::FECollection<dim,spacedim>  &fe_collection,
                     const hp::QCollection<dim-1>     &q_collection,
                     const UpdateFlags             update_flags);


    /**
     * Constructor. Initialize this object with the given parameters, and
     * choose a @p MappingQ1 object for the mapping object.
     *
     * The finite element collection parameter is actually ignored, but is in
     * the signature of this function to make it compatible with the signature
     * of the respective constructor of the usual FEValues object, with the
     * respective parameter in that function also being the return value of
     * the <tt>DoFHandler::get_fe()</tt> function.
     */
    FESubfaceValues (const hp::FECollection<dim,spacedim> &fe_collection,
                     const hp::QCollection<dim-1>    &q_collection,
                     const UpdateFlags            update_flags);

    /**
     * Reinitialize the object for the given cell, face, and subface.
     *
     * After the call, you can get an FESubfaceValues object using the
     * get_present_fe_values() function that corresponds to the present
     * cell. For this FESubfaceValues object, we use the additional arguments
     * described below to determine which finite element, mapping, and
     * quadrature formula to use. They are order in such a way that the
     * arguments one may want to change most frequently come first. The rules
     * for these arguments are as follows:
     *
     * If the @p q_index argument is left at its default value, then we use
     * that quadrature formula within the hp::QCollection passed to the
     * constructor of this class with index given by
     * <code>cell-@>active_fe_index()</code>, i.e. the same index as that of
     * the finite element. In this case, there should be a corresponding
     * quadrature formula for each finite element in the hp::FECollection. As
     * a special case, if the quadrature collection contains only a single
     * element (a frequent case if one wants to use the same quadrature object
     * for all finite elements in an hp discretization, even if that may not
     * be the most efficient), then this single quadrature is used unless a
     * different value for this argument is specified. On the other hand, if a
     * value is given for this argument, it overrides the choice of
     * <code>cell-@>active_fe_index()</code> or the choice for the single
     * quadrature.
     *
     * If the @p mapping_index argument is left at its default value, then we
     * use that mapping object within the hp::MappingCollection passed to the
     * constructor of this class with index given by
     * <code>cell-@>active_fe_index()</code>, i.e. the same index as that of
     * the finite element. As above, if the mapping collection contains only a
     * single element (a frequent case if one wants to use a MappingQ1 object
     * for all finite elements in an hp discretization), then this single
     * mapping is used unless a different value for this argument is
     * specified.
     */
    template <class DH, bool lda>
    void
    reinit (const TriaIterator<DoFCellAccessor<DH,lda> > cell,
            const unsigned int face_no,
            const unsigned int subface_no,
            const unsigned int q_index = numbers::invalid_unsigned_int,
            const unsigned int mapping_index = numbers::invalid_unsigned_int,
            const unsigned int fe_index = numbers::invalid_unsigned_int);

    /**
     * Like the previous function, but for non-hp iterators. The reason this
     * (and the other non-hp iterator) function exists is so that one can use
     * hp::FEValues not only for hp::DoFhandler objects, but for all sorts of
     * DoFHandler objects, and triangulations not associated with DoFHandlers
     * in general.
     *
     * Since <code>cell-@>active_fe_index()</code> doesn't make sense for
     * triangulation iterators, this function chooses the zero-th finite
     * element, mapping, and quadrature object from the relevant constructions
     * passed to the constructor of this object. The only exception is if you
     * specify a value different from the default value for any of these last
     * three arguments.
     */
    void
    reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
            const unsigned int face_no,
            const unsigned int subface_no,
            const unsigned int q_index = numbers::invalid_unsigned_int,
            const unsigned int mapping_index = numbers::invalid_unsigned_int,
            const unsigned int fe_index = numbers::invalid_unsigned_int);
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



    template <int dim, int q_dim, class FEValues>
    inline
    const dealii::hp::FECollection<dim,FEValues::space_dimension> &
    FEValuesBase<dim,q_dim,FEValues>::get_fe_collection () const
    {
      return *fe_collection;
    }



    template <int dim, int q_dim, class FEValues>
    inline
    const dealii::hp::MappingCollection<dim,FEValues::space_dimension> &
    FEValuesBase<dim,q_dim,FEValues>::get_mapping_collection () const
    {
      return *mapping_collection;
    }



    template <int dim, int q_dim, class FEValues>
    inline
    const dealii::hp::QCollection<q_dim> &
    FEValuesBase<dim,q_dim,FEValues>::get_quadrature_collection () const
    {
      return q_collection;
    }



    template <int dim, int q_dim, class FEValues>
    inline
    dealii::UpdateFlags
    FEValuesBase<dim,q_dim,FEValues>::get_update_flags () const
    {
      return update_flags;
    }
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif
