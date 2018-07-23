#ifndef loop_stuff_h
#define loop_stuff_h

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
class FEFacetValues;

/**
 * Namespace for views you get from accessing FEFacetValues using an Extractor.
 */
namespace FEFacetViews
{
    /**
   * Base class for the views.
   */
  template <int dim, int spacedim = dim>
  class Base
  {
  public:
    typedef double value_type;

    Base(const FEFacetValues<dim, spacedim> &fe_facet)
      : fe_facet(&fe_facet)
    {}

  protected:
    const FEFacetValues<dim, spacedim> *fe_facet;
  };

  /**
   * The view of a scalar variable.
   */
  template <int dim, int spacedim = dim>
  class Scalar : public Base<dim, spacedim>
  {
  public:
    // needed?   Scalar ();

      /**
     * This is the type returned for values.
     */
    typedef double value_type;

      /**
     * This is the type returned for gradients, for example from gradient_avg
     */
    typedef dealii::Tensor<1, spacedim> gradient_type;

    /**
     * Constructor for an object that represents a single scalar component
     */
    Scalar(const FEFacetValues<dim, spacedim> &fefacet,
           const unsigned int                  component)
      : Base<dim, spacedim>(fefacet),
      component(component),
        extractor(component)
    {}

    /**
     * Return the jump $[u]=u_1 - u_2$ on the facet for the shape function @p idx
     * in the quadrature point @p q_point.
     */
    value_type
    jump(const unsigned int idx, const unsigned int q_point) const
    {
      const unsigned int shape_fct =
        this->fe_facet->facet_dof_idx_to_fe_dof_idx(idx);

      if (shape_fct == idx)
        return this->fe_facet->get_fe_values().shape_value_component(shape_fct,
                                                                     q_point,
                                                                     component);
      else
        {
          if (this->fe_facet->is_boundary_facet())
            return 0.0;
          else
            return -this->fe_facet->get_fe_values_neighbor()
                      .shape_value_component(shape_fct, q_point, component);
        }
    }

    /**
     * Return the average value $\{u\}=\frac{1}{2}(u_1 + u_2)$ on the facet for the shape
     * function @p idx in the quadrature point @p q_point.
     */
    value_type
    avg(const unsigned int idx, const unsigned int q_point) const
    {
      const unsigned int shape_fct = this->fe_facet->facet_dof_idx_to_fe_dof_idx(idx);
      const unsigned int fe_idx = this->fe_facet->facet_dof_idx_fe(idx);
      return 0.5 * this->fe_facet->get_fe_values(fe_idx).shape_value_component(
                       shape_fct, q_point, component);
    }

    /**
     * Return the average of the gradient $\{\nabla u\}$ on the facet for the shape
     * function @p idx in the quadrature point @p q_point.
     */
    gradient_type
    gradient_avg(const unsigned int idx, const unsigned int q_point) const
    {
        const unsigned int shape_fct =
          this->fe_facet->facet_dof_idx_to_fe_dof_idx(idx);
        const unsigned int fe_idx = this->fe_facet->facet_dof_idx_fe(idx);

        return 0.5 * this->fe_facet->get_fe_values(fe_idx)[extractor].gradient(shape_fct, q_point);
    }

    /**
     * Return the left or the right (if @p left is false) value on the facet for the shape
     * function @p idx in the quadrature point @p q_point.
     */
    value_type
    choose(const bool         left,
           const unsigned int idx,
           const unsigned int q_point) const
    {
      const unsigned int shape_fct =
        this->fe_facet->facet_dof_idx_to_fe_dof_idx(idx);
      const unsigned int fe_idx = this->fe_facet->facet_dof_idx_fe(idx);

      if (left && fe_idx==0)
        return this->fe_facet->get_fe_values().shape_value_component(shape_fct,
                                                                     q_point,
                                                                     component);
      if (!left && fe_idx==1)
        return this->fe_facet->get_fe_values_neighbor().shape_value_component(
          shape_fct, q_point, component);

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
    Vector(const FEFacetValues<dim, spacedim> &fefacet,
           const unsigned int                  first_vector_component)
      : Base<dim, spacedim>(fefacet),
      first_vector_component(first_vector_component),
      extractor(first_vector_component)
    {}

    value_type
    jump(const unsigned int idx, const unsigned int q_point) const
    {
        const unsigned int shape_fct =
          this->fe_facet->facet_dof_idx_to_fe_dof_idx(idx);
        if (shape_fct == idx)
            return this->fe_facet->get_fe_values()[extractor].value(shape_fct, q_point);
        else
          {
            if (this->fe_facet->is_boundary_facet())
              return value_type(); // return 0 tensor
            else
              return -this->fe_facet->get_fe_values_neighbor()[extractor].value(shape_fct, q_point);
          }
    }

    value_type
    avg(const unsigned int idx, const unsigned int q_point) const;

    gradient_type
    gradient_avg(const unsigned int idx, const unsigned int q_point) const;

    value_type
    gradient_dot_n_avg(const unsigned int idx, const unsigned int q_point) const
    {
        const unsigned int shape_fct =
          this->fe_facet->facet_dof_idx_to_fe_dof_idx(idx);
        const unsigned int fe_idx = this->fe_facet->facet_dof_idx_fe(idx);
        return 0.5 * this->fe_facet->get_fe_values(fe_idx)[extractor].gradient(shape_fct, q_point)
                    * this->fe_facet->get_fe_values().normal_vector(q_point);
    }

    value_type
    choose(const bool         left,
           const unsigned int idx,
           const unsigned int q_point) const;

  private:
    const unsigned int first_vector_component;
    FEValuesExtractors::Vector extractor;
  };

} // namespace FEFacetViews

template <int dim, int spacedim = dim>
class FEFacetValues
{
public:
  FEFacetValues(const FEFaceValues<dim> &   fe,
                const FESubfaceValues<dim> &fe_sub,
                const FEFaceValues<dim> &   fe_neighbor,
                const FESubfaceValues<dim> &fe_sub_neighbor)
    : internal_fe_face_values(fe),
     internal_fe_subface_values(fe_sub),
     internal_fe_face_values_neighbor(fe_neighbor),
     internal_fe_subface_values_neighbor(fe_sub_neighbor),
     fe_face_values(nullptr),
     fe_face_values_neighbor(nullptr)
  {
    update_view_cache();
    n_dofs_fe          = fe.get_fe().n_dofs_per_cell();
    n_dofs_fe_neighbor = fe_neighbor.get_fe().n_dofs_per_cell();
  }

  /**
   * Construct the FeFacetValues with a single FiniteElement (same on both
   * sides of the facet)
   */
  FEFacetValues(const Mapping<dim, spacedim> &      mapping,
                const FiniteElement<dim, spacedim> &fe,
                const Quadrature<dim - 1> &         quadrature,
                const UpdateFlags                   update_flags)
    :
    internal_fe_face_values(mapping, fe, quadrature, update_flags),
     internal_fe_subface_values(mapping, fe, quadrature, update_flags),
     internal_fe_face_values_neighbor(mapping, fe, quadrature, update_flags),
     internal_fe_subface_values_neighbor(mapping, fe, quadrature, update_flags),
     fe_face_values(nullptr),
     fe_face_values_neighbor(nullptr)
  {
    update_view_cache();
    n_dofs_fe          = fe.n_dofs_per_cell();
    n_dofs_fe_neighbor = fe.n_dofs_per_cell();
  }

  const std::vector<double> &
  get_JxW_values() const
  {
    return internal_fe_face_values.get_JxW_values();
  }

  const std::vector<Tensor<1, spacedim>> &
  get_normal_vectors() const
  {
    return internal_fe_face_values.get_normal_vectors();
  }

  const Quadrature<dim - 1> &
  get_quadrature() const
  {
    return internal_fe_face_values.get_quadrature();
  }

  const UpdateFlags
  get_update_flags() const
  {
    return internal_fe_face_values.get_update_flags();
  }


  /**
   * Re-initialize this object to be used on a new facet given by two faces of
   * two neighboring cells.
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

    facet_dof_indices.clear();
    facet_dof_indices.insert(facet_dof_indices.end(), v.begin(), v.end());
    facet_dof_indices.insert(facet_dof_indices.end(), v2.begin(), v2.end());
  }

  /**
   * Re-initialize this object to be used on a facet given by a single face of
   * a cell. This is useful to use FEFacetValues to work on boundaries of the
   * domain.
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

    facet_dof_indices.resize(n_dofs_fe);
    cell->get_dof_indices(facet_dof_indices);
  }

  /**
   * Return the number of DoFs on this facet.
   */
  unsigned
  n_facet_dofs() const
  {
    return n_dofs_fe + is_boundary_facet() ? n_dofs_fe_neighbor : 0;
  }

  /**
   *
   */
  bool
  is_boundary_facet() const
  {
    return fe_face_values_neighbor == nullptr;
  }

  /**
   * Return the set of joint DoF indices (includes indices from both cells)
   */
  std::vector<types::global_dof_index>
  get_facet_dof_indices() const
  {
      Assert(facet_dof_indices.size() == n_facet_dofs(), ExcInternalError());
    return facet_dof_indices;
  }

  unsigned int
  facet_dof_idx_to_fe_dof_idx(unsigned int facet_dof_idx) const
  {
    AssertIndexRange(facet_dof_idx, n_facet_dofs());
    return (facet_dof_idx >= n_dofs_fe) ? (facet_dof_idx - n_dofs_fe) :
                                          facet_dof_idx;
  }

  /**
   * For a given facet dof index return whether it belongs to the cell (0) or
   * the neighbor (1).
   */
  unsigned int
  facet_dof_idx_fe(const unsigned int facet_dof_idx) const
  {
    AssertIndexRange(facet_dof_idx, n_facet_dofs());
    return (facet_dof_idx < n_dofs_fe) ? 0 : 1;
  }

  /**
   * Return the FEFaceValue object of the cell or the neighboring cell
   */
  const FEFaceValuesBase<dim, spacedim> &
  get_fe_values(const unsigned int cell_or_neighbor) const
  {
    AssertIndexRange(cell_or_neighbor, is_boundary_facet() ? 1 : 2);
    return (cell_or_neighbor == 0) ? get_fe_values() : get_fe_values_neighbor();
  }

  const FEFaceValuesBase<dim, spacedim> &
  get_fe_values() const
  {
    return *fe_face_values;
  }

  const FEFaceValuesBase<dim, spacedim> &
  get_fe_values_neighbor() const
  {
    return *fe_face_values_neighbor;
  }

  Tensor<1, spacedim>
  normal(const unsigned int q_point_index) const
  {
    return fe_face_values->normal_vector(q_point_index);
  }

  const FEFacetViews::Scalar<dim, spacedim> &
  operator[](const FEValuesExtractors::Scalar &scalar) const
  {
    Assert(scalar.component < fe_face_values->get_fe().n_components(),
           ExcMessage("Invalid FEValuesExtractors::Scalar!"));

    return cached_views_scalar[scalar.component];
  }

  /**
   * Return *this[FEValuesExtractors::Scalar(0)] if this FE is scalar (has
   * only a single component).
   *
   * Fail otherwise.
   */
  const FEFacetViews::Scalar<dim, spacedim> &
  scalar() const
  {
    Assert(
      fe_face_values->get_fe().n_components() == 1,
      ExcMessage(
        "FEFacet::scalar() is only valid if the FE has exactly one component!"));
    return cached_views_scalar[0];
  }

  /**
   * Create a view of the current FEValues object that represents a set of
   * <code>dim</code> scalar components (i.e. a vector) of the vector-valued
   * finite element. The concept of views is explained in the documentation of
   * the namespace FEValuesViews and in particular in the
   * @ref vector_valued
   * module.
   */
  const FEFacetViews::Vector<dim, spacedim> &
  operator[](const FEValuesExtractors::Vector &vector) const
  {
    Assert(vector.first_vector_component + spacedim <=
             fe_face_values->get_fe().n_components(),
           ExcMessage("Invalid FEValuesExtractors::Vector!"));

    return cached_views_vector[vector.first_vector_component];
  }

private:

  /**
   * update the internal data structures for the view objects.
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

  std::vector<types::global_dof_index> facet_dof_indices;
  unsigned int                         n_dofs_fe;
  unsigned int                         n_dofs_fe_neighbor;

  std::vector<FEFacetViews::Scalar<dim, spacedim>> cached_views_scalar;
  std::vector<FEFacetViews::Vector<dim, spacedim>> cached_views_vector;
  FEFaceValuesBase<dim> *                          fe_face_values = nullptr;
  FEFaceValuesBase<dim> *fe_face_values_neighbor                  = nullptr;

  FEFaceValues<dim>    internal_fe_face_values;
  FESubfaceValues<dim> internal_fe_subface_values;
  FEFaceValues<dim>    internal_fe_face_values_neighbor;
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
              const UpdateFlags facet_update_flags = update_values |
                                                     update_gradients |
                                                     update_quadrature_points |
                                                     update_JxW_values |
                                                     update_normal_vectors)
    :
    fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags)
    , fe_facet_values(mapping,
                      fe,
                      QGauss<dim - 1>(quadrature_degree),
                      facet_update_flags)
  {}


  ScratchData(const ScratchData<dim> &scratch_data)
    : fe_values(scratch_data.fe_values.get_mapping(),
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                scratch_data.fe_values.get_update_flags())
    , fe_facet_values(scratch_data.fe_values
                        .get_mapping(), // TODO: implement for fe_facet_values
                      scratch_data.fe_values.get_fe(),
                      scratch_data.fe_facet_values.get_quadrature(),
                      scratch_data.fe_facet_values.get_update_flags())
  {}

  FEValues<dim>      fe_values;
  FEFacetValues<dim> fe_facet_values;
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
copy(const CopyData &        c,
     const ConstraintMatrix &constraints,
     MatrixType &            system_matrix,
     VectorType &            system_rhs)
{
  constraints.distribute_local_to_global(c.cell_matrix,
                                         c.cell_rhs,
                                         c.local_dof_indices,
                                         system_matrix,
                                         system_rhs);
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
