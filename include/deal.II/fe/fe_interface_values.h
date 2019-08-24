#ifndef dealii_fe_interface_h
#define dealii_fe_interface_h

#include <deal.II/meshworker/mesh_loop.h>

/*
 * TODO:
 * - test hp
 * - tests
 * - use fefacet for boundaries?
 * - test GMG
 * - example with error estimator
 * - choose_gradient()
 * - add missing get_function_values(...)
 */

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
class FEInterfaceValues;

/**
 * Namespace for views you get from accessing FEInterfaceValues using an
 * Extractor.
 */
namespace FEInterfaceViews
{
  /**
   * Base class for the views.
   */
  template <int dim, int spacedim = dim>
  class Base
  {
  public:
    typedef double value_type;

    Base(const FEInterfaceValues<dim, spacedim> &fe_interface)
      : fe_interface(&fe_interface)
    {}

  protected:
    const FEInterfaceValues<dim, spacedim> *fe_interface;
  };

  /**
   * The view of a scalar variable.
   */
  template <int dim, int spacedim = dim>
  class Scalar : public Base<dim, spacedim>
  {
  public:
    /**
     * This is the type returned for values.
     */
    typedef double value_type;

    /**
     * This is the type returned for gradients, for example from
     * gradient_average
     */
    typedef dealii::Tensor<1, spacedim> gradient_type;

    /**
     * Constructor for an object that represents a single scalar component
     */
    Scalar(const FEInterfaceValues<dim, spacedim> &fe_interface,
           const unsigned int                      component)
      : Base<dim, spacedim>(fe_interface)
      , component(component)
      , extractor(component)
    {}

    /**
     * Return the jump $[u]=u_1 - u_2$ on the interface for the shape function
     * @p interface_dof_index in the quadrature point @p q_point.
     */
    value_type
    jump(const unsigned int interface_dof_index,
         const unsigned int q_point) const
    {
      const unsigned int shape_fct =
        this->fe_interface->interface_dof_index_to_fe_dof_index(
          interface_dof_index);
      const unsigned int fe_idx =
        this->fe_interface->interface_dof_index_to_fe_index(
          interface_dof_index);

      if (fe_idx == 0)
        return this->fe_interface->get_fe_values().shape_value_component(
          shape_fct, q_point, component);
      else
        {
          if (this->fe_interface->at_boundary())
            return 0.0;
          else
            return -this->fe_interface->get_fe_values_neighbor()
                      .shape_value_component(shape_fct, q_point, component);
        }
    }

    /**
     * Return the average value $\{u\}=\frac{1}{2}(u_1 + u_2)$ on the interface
     * for the shape
     * function @p interface_dof_index in the quadrature point @p q_point.
     */
    value_type
    average(const unsigned int interface_dof_index,
            const unsigned int q_point) const
    {
      const unsigned int shape_fct =
        this->fe_interface->interface_dof_index_to_fe_dof_index(
          interface_dof_index);
      const unsigned int fe_idx =
        this->fe_interface->interface_dof_index_to_fe_index(
          interface_dof_index);

      return 0.5 * this->fe_interface->get_fe_values(fe_idx)
                     .shape_value_component(shape_fct, q_point, component);
    }

    /**
     * Return the average of the gradient $\{\nabla u\}$ on the interface for
     * the shape
     * function @p interface_dof_index in the quadrature point @p q_point.
     */
    gradient_type
    gradient_average(const unsigned int interface_dof_index,
                     const unsigned int q_point) const
    {
      const unsigned int shape_fct =
        this->fe_interface->interface_dof_index_to_fe_dof_index(
          interface_dof_index);
      const unsigned int fe_idx =
        this->fe_interface->interface_dof_index_to_fe_index(
          interface_dof_index);

      return 0.5 * this->fe_interface->get_fe_values(fe_idx)[extractor]
                     .gradient(shape_fct, q_point);
    }

    /**
     * Return the current (if @p current_cell is true) or neighboring (if @p
     * current_cell is false) value on the interface for the shape function @p
     * interface_dof_index in the quadrature point @p q_point of component @p
     * component.
     */
    value_type
    choose(const bool         current_cell,
           const unsigned int interface_dof_index,
           const unsigned int q_point) const
    {
      const unsigned int shape_fct =
        this->fe_interface->interface_dof_index_to_fe_dof_index(
          interface_dof_index);
      const unsigned int fe_idx =
        this->fe_interface->interface_dof_index_to_fe_index(
          interface_dof_index);

      if (current_cell && fe_idx == 0)
        return this->fe_interface->get_fe_values().shape_value_component(
          shape_fct, q_point, component);
      if (!current_cell && fe_idx == 1)
        return this->fe_interface->get_fe_values_neighbor()
          .shape_value_component(shape_fct, q_point, component);

      return 0.0;
    }

  private:
    /**
     * component index of this view.
     */
    const unsigned int component;
    /**
     * The extractor for this view.
     */
    FEValuesExtractors::Scalar extractor;
  };

  template <int dim, int spacedim = dim>
  class Vector : public Base<dim, spacedim>
  {
  public:
    typedef dealii::Tensor<1, spacedim> value_type;

    typedef dealii::Tensor<2, spacedim> gradient_type;

    /**
     * Constructor for an object that represents a vector component
     */
    Vector(const FEInterfaceValues<dim, spacedim> &fe_interface,
           const unsigned int                      first_vector_component)
      : Base<dim, spacedim>(fe_interface)
      , first_vector_component(first_vector_component)
      , extractor(first_vector_component)
    {}

    value_type
    jump(const unsigned int idx, const unsigned int q_point) const
    {
      const unsigned int shape_fct =
        this->fe_interface->interface_dof_index_to_fe_dof_index(idx);
      if (shape_fct == idx)
        return this->fe_interface->get_fe_values()[extractor].value(shape_fct,
                                                                    q_point);
      else
        {
          if (this->fe_interface->at_boundary())
            {
              // we should not get here
              Assert(false, ExcInternalError());
              return value_type(); // return 0 tensor
            }
          else
            return -this->fe_interface->get_fe_values_neighbor()[extractor]
                      .value(shape_fct, q_point);
        }
    }

    value_type
    average(const unsigned int idx, const unsigned int q_point) const
    {
      const unsigned int shape_fct =
        this->fe_interface->interface_dof_index_to_fe_dof_index(idx);
      const unsigned int fe_idx =
        this->fe_interface->interface_dof_index_to_fe_index(idx);

      return (fe_idx == 0 ? 0.5 : -0.5) *
             this->fe_interface->get_fe_values(fe_idx)[extractor].value(
               shape_fct, q_point);
    }

    gradient_type
    gradient_average(const unsigned int idx, const unsigned int q_point) const
    {
      const unsigned int shape_fct =
        this->fe_interface->interface_dof_index_to_fe_dof_index(idx);
      const unsigned int fe_idx =
        this->fe_interface->interface_dof_index_to_fe_index(idx);
      return 0.5 * this->fe_interface->get_fe_values(fe_idx)[extractor]
                     .gradient(shape_fct, q_point);
    }

    value_type
    gradient_dot_n_average(const unsigned int idx,
                           const unsigned int q_point) const
    {
      const unsigned int shape_fct =
        this->fe_interface->interface_dof_index_to_fe_dof_index(idx);
      const unsigned int fe_idx =
        this->fe_interface->interface_dof_index_to_fe_index(idx);
      return 0.5 *
             this->fe_interface->get_fe_values(fe_idx)[extractor].gradient(
               shape_fct, q_point) *
             this->fe_interface->get_fe_values().normal_vector(q_point);
    }

    value_type
    choose(const bool         left,
           const unsigned int idx,
           const unsigned int q_point) const;

  private:
    const unsigned int         first_vector_component;
    FEValuesExtractors::Vector extractor;
  };

} // namespace FEInterfaceViews


/**
 * FEInterfaceValues is a data structure to access and assemble Finite Element
 * data on facets or interfaces of a mesh.
 *
 * Basically, it provides a way to access average and jump terms used in
 * Discontinuous Galerkin methods on faces between two neighboring cells.
 *
 * Internally, it provides an abstraction for two FEFaceValues.
 * The class introduces a new "interface dof index" that walks over the union of
 * the dof indices of the two FEFaceValues objects.
 *
 * The class is made to be used inside MeshWorker::mesh_loop. It is intended to
 * be a low level replacement for MeshWorker and LocalIntegrators.
 */
template <int dim, int spacedim = dim>
class FEInterfaceValues
{
public:
  /**
   * Construct an FEInterfaceValues from existing FEFaceValues.
   */
  FEInterfaceValues(const FEFaceValues<dim> &   fe,
                    const FESubfaceValues<dim> &fe_sub,
                    const FEFaceValues<dim> &   fe_neighbor,
                    const FESubfaceValues<dim> &fe_sub_neighbor)
    : internal_fe_face_values(fe)
    , internal_fe_subface_values(fe_sub)
    , internal_fe_face_values_neighbor(fe_neighbor)
    , internal_fe_subface_values_neighbor(fe_sub_neighbor)
    , fe_face_values(nullptr)
    , fe_face_values_neighbor(nullptr)
  {
    update_view_cache();
    n_dofs_fe          = fe.get_fe().n_dofs_per_cell();
    n_dofs_fe_neighbor = fe_neighbor.get_fe().n_dofs_per_cell();
  }

  /**
   * Construct the FeInterfaceValues with a single FiniteElement (same on both
   * sides of the facet). The FeFaceValues objects will be initialized with
   * the given @p mapping, @p quadrature, and @p update_flags.
   *
   */
  FEInterfaceValues(const Mapping<dim, spacedim> &      mapping,
                    const FiniteElement<dim, spacedim> &fe,
                    const Quadrature<dim - 1> &         quadrature,
                    const UpdateFlags                   update_flags)
    : internal_fe_face_values(mapping, fe, quadrature, update_flags)
    , internal_fe_subface_values(mapping, fe, quadrature, update_flags)
    , internal_fe_face_values_neighbor(mapping, fe, quadrature, update_flags)
    , internal_fe_subface_values_neighbor(mapping, fe, quadrature, update_flags)
    , fe_face_values(nullptr)
    , fe_face_values_neighbor(nullptr)
  {
    update_view_cache();
    n_dofs_fe          = fe.n_dofs_per_cell();
    n_dofs_fe_neighbor = fe.n_dofs_per_cell();
  }


  /**
   * Re-initialize this object to be used on a new interface given by two faces
   * of two neighboring cells.
   *
   * The arguments including their order is identical to the @p face_worker arguments
   * in MeshWorker::mesh_loop.
   */
  template <class Iterator>
  void
  reinit(const Iterator &    cell,
         const unsigned int &face_no,
         const unsigned int &sub_face_no,
         const Iterator &    cell_neighbor,
         const unsigned int &face_no_neighbor,
         const unsigned int &sub_face_no_neighbor)
  {
    if (sub_face_no == numbers::invalid_unsigned_int)
      {
        internal_fe_face_values.reinit(cell, face_no);
        fe_face_values = &internal_fe_face_values;
      }
    else
      {
        internal_fe_subface_values.reinit(cell, face_no, sub_face_no);
        fe_face_values = &internal_fe_subface_values;
      }
    if (sub_face_no_neighbor == numbers::invalid_unsigned_int)
      {
        internal_fe_face_values_neighbor.reinit(cell_neighbor,
                                                face_no_neighbor);
        fe_face_values_neighbor = &internal_fe_face_values_neighbor;
      }
    else
      {
        internal_fe_subface_values_neighbor.reinit(cell_neighbor,
                                                   face_no_neighbor,
                                                   sub_face_no_neighbor);
        fe_face_values_neighbor = &internal_fe_subface_values_neighbor;
      }

    std::vector<types::global_dof_index> v(n_dofs_fe);
    cell->get_dof_indices(v);
    std::vector<types::global_dof_index> v2(n_dofs_fe_neighbor);
    cell_neighbor->get_dof_indices(v2);

    interface_dof_indices.clear();
    interface_dof_indices.insert(interface_dof_indices.end(),
                                 v.begin(),
                                 v.end());
    interface_dof_indices.insert(interface_dof_indices.end(),
                                 v2.begin(),
                                 v2.end());
  }

  /**
   * Re-initialize this object to be used on a interface given by a single face
   * of a cell. This is useful to use FEInterfaceValues to work on boundaries of
   * the domain.
   *
   * As a consequence, queries like jump() will assume a value of zero for the
   * values on the "other" side.
   */
  template <class Iterator>
  void
  reinit(const Iterator &cell, const unsigned int &face_no)
  {
    internal_fe_face_values.reinit(cell, face_no);
    fe_face_values          = &internal_fe_face_values;
    fe_face_values_neighbor = nullptr;

    interface_dof_indices.resize(n_dofs_fe);
    cell->get_dof_indices(interface_dof_indices);
  }

  /**
   * Return the vector of JxW values for each quadrature point.
   */
  const std::vector<double> &
  get_JxW_values() const
  {
    return internal_fe_face_values.get_JxW_values();
  }

  /**
   * Return the normal vector of the interface in each quadrature point.
   */
  const std::vector<Tensor<1, spacedim>> &
  get_normal_vectors() const
  {
    return internal_fe_face_values.get_normal_vectors();
  }

  /**
   * Return a reference to the quadrature object in use.
   */
  const Quadrature<dim - 1> &
  get_quadrature() const
  {
    return internal_fe_face_values.get_quadrature();
  }

  /**
   * Return a reference to the quadrature object in use.
   */
  const std::vector<Point<spacedim>> &
  get_quadrature_points() const
  {
    return internal_fe_face_values.get_quadrature_points();
  }

  /**
   * Return the update flags set.
   */
  const UpdateFlags
  get_update_flags() const
  {
    return internal_fe_face_values.get_update_flags();
  }


  /**
   * Return the number of DoFs on the current interface. This is the sum of the
   * number of DoFs on the two adjacent cells.
   */
  unsigned
  n_interface_dofs() const
  {
    return n_dofs_fe + (at_boundary() ? 0 : n_dofs_fe_neighbor);
  }

  /**
   * Return if the current interface, which is set by reinit(), is an internal
   * face with two adjacent cells or a boundary face.
   */
  bool
  at_boundary() const
  {
    return fe_face_values_neighbor == nullptr;
  }

  /**
   * Return the set of joint DoF indices. This includes indices from both cells.
   */
  std::vector<types::global_dof_index>
  get_interface_dof_indices() const
  {
    Assert(interface_dof_indices.size() == n_interface_dofs(),
           ExcInternalError());
    return interface_dof_indices;
  }

  /**
   * Translate a local interface_dof_index to a local dof_index of the
   * underlying FEFaceValues.
   */
  unsigned int
  interface_dof_index_to_fe_dof_index(unsigned int interface_dof_index) const
  {
    AssertIndexRange(interface_dof_index, n_interface_dofs());
    return (interface_dof_index >= n_dofs_fe) ?
             (interface_dof_index - n_dofs_fe) :
             interface_dof_index;
  }

  /**
   * For a given interface dof index return whether it belongs to the cell (0)
   * or to the neighbor (1).
   */
  unsigned int
  interface_dof_index_to_fe_index(const unsigned int interface_dof_index) const
  {
    AssertIndexRange(interface_dof_index, n_interface_dofs());
    return (interface_dof_index < n_dofs_fe) ? 0 : 1;
  }

  /**
   * Convert an FiniteElement index @p fe_index (0=cell, 1=neighbor) and dof
   * index @p face_dof_index into an interface DoF index.
   *
   * The inverse of this operation is done with
   * interface_dof_index_to_fe_index() and interface_dof_idx_to_fe_dof_idx().
   */
  unsigned int
  interface_dof_index(const unsigned int fe_index,
                      const unsigned int face_dof_index) const
  {
    Assert(fe_index <= 1, ExcMessage("fe_index should be 0 or 1"));
    if (at_boundary())
      Assert(fe_index == 0,
             ExcMessage("A boundary facet only has FE index 0."));

    if (fe_index == 0)
      {
        Assert(face_dof_index < n_dofs_fe, ExcMessage("invalid face_dof_idx"));
        return face_dof_index;
      }
    else if (fe_index == 1)
      {
        Assert(face_dof_index < n_dofs_fe_neighbor,
               ExcMessage("invalid face_dof_idx"));
        return face_dof_index + n_dofs_fe;
      }
    // return something invalid:
    return -1;
  }

  /**
   * Return the FEFaceValue object of the cell or the neighboring cell
   *
   * If @p cell_or_neighbor is 0, return the cell's, if it is 1, return the
   * neighboring cell's FEFaceValue object.
   *
   * @note The argument @p cell_or_neighbor is returned by
   * interface_dof_index_to_fe_index().
   */
  const FEFaceValuesBase<dim, spacedim> &
  get_fe_values(const unsigned int cell_or_neighbor) const
  {
    AssertIndexRange(cell_or_neighbor, 2);
    return (cell_or_neighbor == 0) ? get_fe_values() : get_fe_values_neighbor();
  }

  /**
   * Return the FEFaceValue object of the current cell
   */
  const FEFaceValuesBase<dim, spacedim> &
  get_fe_values() const
  {
    return *fe_face_values;
  }

  /**
   * Return the FEFaceValue object of the current neighboring cell
   */
  const FEFaceValuesBase<dim, spacedim> &
  get_fe_values_neighbor() const
  {
    Assert(!at_boundary(), ExcMessage("Not possible for boundary Facet."));
    return *fe_face_values_neighbor;
  }

  /**
   * Return the normal in a given quadrature point.
   */
  Tensor<1, spacedim>
  normal(const unsigned int q_point_index) const
  {
    return fe_face_values->normal_vector(q_point_index);
  }

  /**
   * Return a view of a given scalar FE given by the extractor @p scalar.
   *
   * The concept of views is explained in the documentation of
   * the namespace FEValuesViews and in particular in the
   * @ref vector_valued
   * module.
   */
  const FEInterfaceViews::Scalar<dim, spacedim> &
  operator[](const FEValuesExtractors::Scalar &scalar) const
  {
    Assert(scalar.component < fe_face_values->get_fe().n_components(),
           ExcMessage("Invalid FEValuesExtractors::Scalar!"));

    return cached_views_scalar[scalar.component];
  }

  /**
   * Create a view of the current FEValues object that represents a set of
   * <code>dim</code> scalar components (i.e. a vector) of the vector-valued
   * finite element. The concept of views is explained in the documentation of
   * the namespace FEValuesViews and in particular in the
   * @ref vector_valued
   * module.
   */
  const FEInterfaceViews::Vector<dim, spacedim> &
  operator[](const FEValuesExtractors::Vector &vector) const
  {
    Assert(vector.first_vector_component + spacedim <=
             fe_face_values->get_fe().n_components(),
           ExcMessage("Invalid FEValuesExtractors::Vector!"));

    return cached_views_vector[vector.first_vector_component];
  }

  /**
   * Return the current (if @p current_cell is true) or neighboring (if @p
   * current_cell is false) value on the interface for the shape function @p
   * interface_dof_index in the quadrature point @p q_point of component @p
   * component.
   *
   * @note This function is typically used to pick the upstream or downstream
   * value based on a direction. This can be achieved by using
   * <code>(direction * normal)>0</code> as the first argument of this
   * function.
   */
  double
  choose(const bool         current_cell,
         const unsigned int interface_dof_index,
         const unsigned int q_point,
         const unsigned int component = 0) const
  {
    const unsigned int shape_fct =
      interface_dof_index_to_fe_dof_index(interface_dof_index);
    const unsigned int fe_idx =
      interface_dof_index_to_fe_index(interface_dof_index);

    if (current_cell && fe_idx == 0)
      return get_fe_values().shape_value_component(shape_fct,
                                                   q_point,
                                                   component);
    if (!current_cell && fe_idx == 1)
      return get_fe_values_neighbor().shape_value_component(shape_fct,
                                                            q_point,
                                                            component);

    return 0.0;
  }

  /**
   * Return the jump $[u]=u_1 - u_2$ on the interface for the shape function
   * @p interface_dof_index in the quadrature point @p q_point of component @p
   * component.
   */
  double
  jump(const unsigned int interface_dof_index,
       const unsigned int q_point,
       const unsigned int component = 0) const
  {
    const unsigned int shape_fct =
      interface_dof_index_to_fe_dof_index(interface_dof_index);
    const unsigned int fe_idx =
      interface_dof_index_to_fe_index(interface_dof_index);

    if (fe_idx == 0)
      return get_fe_values().shape_value_component(shape_fct,
                                                   q_point,
                                                   component);
    else
      {
        if (at_boundary())
          return 0.0;
        else
          return -get_fe_values_neighbor().shape_value_component(shape_fct,
                                                                 q_point,
                                                                 component);
      }
  }

private:
  /**
   * Update the internal data structures for the view objects.
   */
  void
  update_view_cache()
  {
    auto &fe = internal_fe_face_values.get_fe();

    cached_views_scalar.clear();
    cached_views_scalar.reserve(fe.n_components());
    for (unsigned int c = 0; c < fe.n_components(); ++c)
      cached_views_scalar.emplace_back(*this, c);

    const unsigned int n_vec_views =
      (fe.n_components() < spacedim) ? 0 : (fe.n_components() - spacedim);
    cached_views_vector.clear();

    cached_views_vector.reserve(n_vec_views);
    for (unsigned int c = 0; c < n_vec_views; ++c)
      cached_views_vector.emplace_back(*this, c);
  }

  /**
   * The list of DoF indices for the current interface, filled in reinit().
   */
  std::vector<types::global_dof_index> interface_dof_indices;

  /**
   * The number of DoFs belonging to the current cell
   */
  unsigned int n_dofs_fe;

  /**
   * The number of DoFs belonging to the neighboring cell
   */
  unsigned int n_dofs_fe_neighbor;

  /**
   * Cache for the scalar views, filled in update_view_cache().
   */
  std::vector<FEInterfaceViews::Scalar<dim, spacedim>> cached_views_scalar;

  /**
   * Cache for the vector views, filled in update_view_cache().
   */
  std::vector<FEInterfaceViews::Vector<dim, spacedim>> cached_views_vector;

  /**
   * Pointer to internal_fe_face_values or internal_fe_subface_values,
   * respectively.
   */
  FEFaceValuesBase<dim> *fe_face_values = nullptr;
  /**
   * Pointer to internal_fe_face_values_neighbor,
   * internal_fe_subface_values_neighbor, or nullptr, respectively.
   */
  FEFaceValuesBase<dim> *fe_face_values_neighbor = nullptr;

  /**
   * The FEFaceValues object for the current cell.
   */
  FEFaceValues<dim> internal_fe_face_values;

  /**
   * The FEFaceValues object for the current cell if the cell is refined.
   */
  FESubfaceValues<dim> internal_fe_subface_values;

  /**
   * The FEFaceValues object for the neighboring cell.
   */
  FEFaceValues<dim> internal_fe_face_values_neighbor;

  /**
   * The FEFaceValues object for the neighboring cell if the cell is refined.
   */
  FESubfaceValues<dim> internal_fe_subface_values_neighbor;
};



// example scratch/copydata
template <int dim>
struct ScratchData
{
  ScratchData(const Mapping<dim> &      mapping,
              const FiniteElement<dim> &fe,
              const unsigned int        quadrature_degree,
              const UpdateFlags         update_flags = update_values |
                                               update_gradients |
                                               update_quadrature_points |
                                               update_JxW_values,
              const UpdateFlags interface_update_flags =
                update_values | update_gradients | update_quadrature_points |
                update_JxW_values | update_normal_vectors)
    : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags)
    , fe_interface_values(mapping,
                          fe,
                          QGauss<dim - 1>(quadrature_degree),
                          interface_update_flags)
  {}


  ScratchData(const ScratchData<dim> &scratch_data)
    : fe_values(scratch_data.fe_values.get_mapping(),
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                scratch_data.fe_values.get_update_flags())
    , fe_interface_values(
        scratch_data.fe_values
          .get_mapping(), // TODO: implement for fe_interface_values
        scratch_data.fe_values.get_fe(),
        scratch_data.fe_interface_values.get_quadrature(),
        scratch_data.fe_interface_values.get_update_flags())
  {}

  FEValues<dim>          fe_values;
  FEInterfaceValues<dim> fe_interface_values;
};



struct CopyDataFace
{
  FullMatrix<double>                   cell_matrix;
  std::vector<types::global_dof_index> joint_dof_indices;
};

struct CopyData
{
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<CopyDataFace>            face_data;

  template <class Iterator>
  void
  reinit(const Iterator &cell, unsigned int dofs_per_cell)
  {
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);

    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
  // CopyDataFace &reinit_face();
};


template <class MatrixType, class VectorType>
inline void
copy(const CopyData &                 c,
     const AffineConstraints<double> &constraints,
     MatrixType &                     system_matrix,
     VectorType &                     system_rhs)
{
  constraints.distribute_local_to_global(
    c.cell_matrix, c.cell_rhs, c.local_dof_indices, system_matrix, system_rhs);
  for (auto &cdf : c.face_data)
    {
      // TODO: use constraints.distribute(...)
      const unsigned int dofs_per_cell = cdf.joint_dof_indices.size();
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int k = 0; k < dofs_per_cell; ++k)
          system_matrix.add(cdf.joint_dof_indices[i],
                            cdf.joint_dof_indices[k],
                            cdf.cell_matrix(i, k));
    }
}

DEAL_II_NAMESPACE_CLOSE

#endif
