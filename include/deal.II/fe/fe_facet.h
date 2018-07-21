#ifndef loop_stuff_h
#define loop_stuff_h

#include <deal.II/meshworker/mesh_loop.h>

namespace {
    using namespace dealii;


    template <int dim, int spacedim>
    class FEFacetValues;

    namespace FEFacetViews
    {
        template <int dim, int spacedim=dim>
        class Base
        {
        public:
             typedef double value_type;

            Base(const FEFacetValues<dim,spacedim> &fe_facet)
                : fe_facet (&fe_facet)
            {}

        protected:
            const SmartPointer<const FEFacetValues<dim,spacedim> > fe_facet;
        };

        template <int dim, int spacedim=dim>
        class Scalar: public Base<dim,spacedim>
        {
        public:
         // needed?   Scalar ();

            typedef double          value_type;

            typedef dealii::Tensor<1,spacedim> gradient_type;


            /**
             * Constructor for an object that represents a single scalar component
             */
            Scalar (const FEFacetValues<dim,spacedim> &fefacet,
                    const unsigned int component)
                : Base<dim,spacedim>(fefacet),
                  component (component)
            {}

            value_type jump (const unsigned int idx, const unsigned int q_point) const
            {
                const unsigned int shape_fct = this->fe_facet->facet_to_fe_dof_idx(idx);
                if (shape_fct == idx)
                    return this->fe_facet->get_fe_values().shape_value_component(shape_fct, q_point, component);
                else
                    return -this->fe_facet->get_fe_values_neighbor().shape_value_component(shape_fct, q_point, component);
            }

            value_type avg (const unsigned int idx, const unsigned int q_point) const
            {
                const unsigned int shape_fct = this->fe_facet->facet_to_fe_dof_idx(idx);
                if (shape_fct == idx)
                    return 0.5*this->fe_facet->get_fe_values().shape_value_component(shape_fct, q_point, component);
                else
                    return 0.5*this->fe_facet->get_fe_values_neighbor().shape_value_component(shape_fct, q_point, component);
            }

            gradient_type gradient_avg (const unsigned int idx, const unsigned int q_point) const
            {
                Assert(false ,ExcNotImplemented());
            }

            value_type choose (const bool left, const unsigned int idx, const unsigned int q_point) const
            {
                const unsigned int shape_fct = this->fe_facet->facet_to_fe_dof_idx(idx);
                if (left && shape_fct == idx)
                    return this->fe_facet->get_fe_values().shape_value_component(shape_fct, q_point, component);
                if (!left && shape_fct != idx)
                    return this->fe_facet->get_fe_values_neighbor().shape_value_component(shape_fct, q_point, component);
                return 0.0;
            }

        private:
            const unsigned int component;
        };

        template <int dim, int spacedim=dim>
        class Vector: public Base<dim,spacedim>
        {
        public:
            typedef dealii::Tensor<1,spacedim>          value_type;

            typedef dealii::Tensor<2,spacedim>          gradient_type;

            /**
             * Constructor for an object that represents a vector component
             */
            Vector (const FEFacetValues<dim,spacedim> &fefacet,
                    const unsigned int first_vector_component)
                : Base<dim,spacedim>(fefacet),
                  first_vector_component (first_vector_component)
            {}

            value_type jump (const unsigned int idx, const unsigned int q_point) const;

            value_type avg (const unsigned int idx, const unsigned int q_point) const;

            gradient_type gradient_avg (const unsigned int idx, const unsigned int q_point) const;

            value_type choose (const bool left, const unsigned int idx, const unsigned int q_point) const;

        private:
            const unsigned int first_vector_component;
        };

    }

    template <int dim, int spacedim=dim>
    // TODO: if we derive from FEFacetViews Scalar we can make scalar assembly easier. This would mirror FeValues::value() etc.
    class FEFacetValues// : public FEFacetViews::Scalar<dim, spacedim>
    {
    public:
        FEFacetValues(const FEFaceValues<dim> &fe, const FESubfaceValues<dim> &fe_sub,
                      const FEFaceValues<dim> &fe_neighbor, const FESubfaceValues<dim> &fe_sub_neighbor)
            : //FEFacetViews::Scalar<dim, spacedim>(),

              internal_fe_face_values (fe),
              internal_fe_subface_values (fe_sub),
              internal_fe_face_values_neighbor (fe_neighbor),
              internal_fe_subface_values_neighbor (fe_sub_neighbor),
              fe_face_values (nullptr),
              fe_face_values_neighbor (nullptr)
        {
            update_view_cache();
            n_dofs_fe = fe.get_fe().n_dofs_per_cell();
            n_dofs_fe_neighbor = fe_neighbor.get_fe().n_dofs_per_cell();
        }

        /**
         * Construct the FeFacetValues with a single FiniteElement (same on both sides of the facet)
         */
        FEFacetValues(const Mapping<dim,spacedim>       &mapping,
                      const FiniteElement<dim,spacedim> &fe,
                      const Quadrature<dim-1>           &quadrature,
                      const UpdateFlags                  update_flags)
            : //FEFacetViews::Scalar<dim, spacedim>(),
              internal_fe_face_values (mapping, fe, quadrature, update_flags),
              internal_fe_subface_values (mapping, fe, quadrature, update_flags),
              internal_fe_face_values_neighbor (mapping, fe, quadrature, update_flags),
              internal_fe_subface_values_neighbor (mapping, fe, quadrature, update_flags),
              fe_face_values (nullptr),
              fe_face_values_neighbor (nullptr)
        {
            update_view_cache();
            n_dofs_fe = fe.get_fe().n_dofs_per_cell();
            n_dofs_fe_neighbor = fe.get_fe().n_dofs_per_cell();
        }

        /**
         * Re-initialize this object to be used on a new facet given by two faces of two neighboring cells.
         */
        template <class Iterator>
        void reinit(const Iterator &cell, const unsigned int &face_no, const unsigned int &sub_face_no,
                    const Iterator &cell_neighbor, const unsigned int &face_no_neighbor, const unsigned int &sub_face_no_neighbor)

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
                internal_fe_face_values_neighbor.reinit(cell_neighbor, face_no_neighbor);
                fe_face_values_neighbor = &internal_fe_face_values_neighbor;
              }
            else
              {
                internal_fe_subface_values_neighbor.reinit(cell_neighbor, face_no_neighbor, sub_face_no_neighbor);
                fe_face_values_neighbor = &internal_fe_subface_values_neighbor;
              }
        }

        /**
         * Re-initialize this object to be used on a facet given by a single face of a cell. This is useful to use FeFacet to work on boundaries of the domain.
         *
         * As a consequence, queries like jump() will assume a value of zero for the values on the "other" side.
         */
        template <class Iterator>
        void reinit(const Iterator &cell, const unsigned int &face_no)
        {
            internal_fe_face_values.reinit(cell, face_no);
            fe_face_values = &internal_fe_face_values;
            fe_face_values_neighbor = nullptr;
        }

        /**
         * Return the number of DoFs on this facet.
         */
        unsigned n_facet_dofs () const
        {
            return n_dofs_fe + n_dofs_fe_neighbor;
        }

        /**
         *
         */
        bool is_boundary_facet() const
        {
            return fe_face_values_neighbor==nullptr;
        }

        /**
         * Return the set of joint DoF indices (includes indices from both cells)
         */
        std::vector<types::global_dof_index> get_facet_dof_indices() const;

        unsigned int facet_dof_idx_to_fe_dof_idx(unsigned int facet_dof_idx) const
        {
            AssertIndexRange(facet_dof_idx, n_facet_dofs());
            return (facet_dof_idx >= n_dofs_fe) ?  (facet_dof_idx-n_dofs_fe) : facet_dof_idx;
        }

        /**
         * For a given facet dof index return whether it belongs to the cell (0) or the neighbor (1).
         */
        unsigned int facet_dof_idx_fe(const unsigned int facet_dof_idx) const
        {
            AssertIndexRange(facet_dof_idx, n_facet_dofs());
            return (facet_dof_idx<n_dofs_fe) ? 0 : 1;
        }

        /**
         * Return the FeFaceValue object of the cell or the neighboring cell
         */
        const FEFaceValuesBase<dim,spacedim> & get_fe_values(const unsigned int cell_or_neighbor) const
        {
            AssertIndexRange(cell_or_neighbor, is_boundary_facet() ? 1 : 2);
            return (cell_or_neighbor==0) ? get_fe_values() : get_fe_values_neighbor();
        }

        const FEFaceValuesBase<dim,spacedim> & get_fe_values() const
        {
            return *fe_face_values;
        }

        const FEFaceValuesBase<dim,spacedim> & get_fe_values_neighbor() const
        {
            return *fe_face_values_neighbor;
        }

        Tensor<1,spacedim> normal(const unsigned int q_point_index) const
        {
            return fe_face_values->normal_vector(q_point_index);
        }

        const FEFacetViews::Scalar<dim,spacedim> &
        operator[] (const FEValuesExtractors::Scalar &scalar) const
        {
            Assert(scalar.component < fe_face_values->get_fe()->n_components(),
                   ExcMessage("Invalid FEValuesExtractors::Scalar!"));

            return cached_views_scalar[scalar.component];
        }

        /**
         * Return *this[FEValuesExtractors::Scalar(0)] if this FE is scalar (has only a single component).
         *
         * Fail otherwise.
         */
        const FEFacetViews::Scalar<dim,spacedim> &
        scalar () const
        {
            Assert(fe_face_values->get_fe()->n_components()==1,
                   ExcMessage("FEFacet::scalar() is only valid if the FE has exactly one component!"));
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
        const FEFacetViews::Vector<dim,spacedim> &
        operator[] (const FEValuesExtractors::Vector &vector) const
        {
            Assert(vector.first_vector_component
                   + spacedim <= fe_face_values->get_fe()->n_components(),
                   ExcMessage("Invalid FEValuesExtractors::Vector!"));

            return cached_views_vector[vector.first_vector_component];
        }


//        /**
//         * Return the vector of normal vectors
//         */
//        const std::vector<Tensor<1,spacedim> > &get_normal_vectors () const
//        {
//            return fe_face_values->get_normal_vectors();
//        }

//        /**
//         * @brief JxW
//         * @param quadrature_point
//         * @return
//         */
//        double JxW (const unsigned int quadrature_point) const
//        {
//            return fe_face_values->JxW(quadrature_point);
//        }

//        /**
//         * Return a reference to the array holding the values returned by JxW().
//         */
//        const std::vector<double> &get_JxW_values () const
//        {
//            return fe_face_values->get_JxW_values();
//        }

//        const std::vector<Point<spacedim> > &
//        FEValuesBase<dim,spacedim>::get_quadrature_points () const
//        {
//            return fe_face_values->get_quadrature_points();
//        }

    private:
        void update_view_cache ()
        {
            auto & fe = internal_fe_face_values.get_fe();

            cached_views_scalar.clear();
            cached_views_scalar.reserve(fe.n_components());
            for (unsigned int c=0; c<fe.n_components(); ++c)
                cached_views_scalar.emplace_back(*this, c);

            const unsigned int n_vec_views = (fe.n_components() < spacedim)
                    ? 0 : (fe.n_components()-spacedim);
            cached_views_vector.clear();
            cached_views_vector.reserve(n_vec_views);
            for (unsigned int c=0; c<n_vec_views; ++c)
                cached_views_vector.emplace_back(*this, c);
        }

        unsigned int n_dofs_fe;
        unsigned int n_dofs_fe_neighbor;

        std::vector<FEFacetViews::Scalar<dim,spacedim> > cached_views_scalar;
        std::vector<FEFacetViews::Vector<dim,spacedim> > cached_views_vector;
        FEFaceValuesBase<dim> *fe_face_values = nullptr;
        FEFaceValuesBase<dim> *fe_face_values_neighbor = nullptr;

        FEFaceValues<dim> internal_fe_face_values;
        FESubfaceValues<dim> internal_fe_subface_values;
        FEFaceValues<dim> internal_fe_face_values_neighbor;
        FESubfaceValues<dim> internal_fe_subface_values_neighbor;
    };





// example scratch/copydata

template <int dim>
struct ScratchData
{
  ScratchData (const Mapping<dim> &mapping,
               const FiniteElement<dim> &fe,
               const unsigned int quadrature_degree,
               const UpdateFlags cell_flags =
      update_values | update_gradients | update_quadrature_points | update_JxW_values,
      const UpdateFlags face_flags =
      update_values | update_gradients | update_quadrature_points | update_JxW_values | update_normal_vectors)
    :
    fe_values (mapping,
               fe,
               QGauss<dim>(quadrature_degree),
               cell_flags),
    internal_fe_face_values (mapping,
                             fe,
                             QGauss<dim-1>(quadrature_degree),
                             face_flags),
    internal_fe_subface_values (mapping,
                                fe,
                                QGauss<dim-1>(quadrature_degree),
                                face_flags),
    internal_fe_face_values_neighbor (mapping,
                                      fe,
                                      QGauss<dim-1>(quadrature_degree),
                                      face_flags),
    internal_fe_subface_values_neighbor (mapping,
                                         fe,
                                         QGauss<dim-1>(quadrature_degree),
                                         face_flags)

  {}

  ScratchData (const ScratchData<dim> &scratch_data)
    :
    fe_values (scratch_data.fe_values.get_mapping(),
               scratch_data.fe_values.get_fe(),
               scratch_data.fe_values.get_quadrature(),
               scratch_data.fe_values.get_update_flags()),
    internal_fe_face_values (scratch_data.internal_fe_face_values.get_mapping(),
                             scratch_data.internal_fe_face_values.get_fe(),
                             scratch_data.internal_fe_face_values.get_quadrature(),
                             scratch_data.internal_fe_face_values.get_update_flags()),
    internal_fe_subface_values (scratch_data.internal_fe_face_values.get_mapping(),
                                scratch_data.internal_fe_face_values.get_fe(),
                                scratch_data.internal_fe_face_values.get_quadrature(),
                                scratch_data.internal_fe_face_values.get_update_flags()),
    internal_fe_face_values_neighbor (scratch_data.internal_fe_face_values.get_mapping(),
                                      scratch_data.internal_fe_face_values.get_fe(),
                                      scratch_data.internal_fe_face_values.get_quadrature(),
                                      scratch_data.internal_fe_face_values.get_update_flags()),
    internal_fe_subface_values_neighbor (scratch_data.internal_fe_face_values.get_mapping(),
                                         scratch_data.internal_fe_face_values.get_fe(),
                                         scratch_data.internal_fe_face_values.get_quadrature(),
                                         scratch_data.internal_fe_face_values.get_update_flags())
  {}

  template <class Iterator>
  void reinit(const Iterator &cell, const unsigned int &face_no, const unsigned int &sub_face_no,
              const Iterator &cell_neighbor, const unsigned int &face_no_neighbor, const unsigned int &sub_face_no_neighbor)
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
        internal_fe_face_values_neighbor.reinit(cell_neighbor, face_no_neighbor);
        fe_face_values_neighbor = &internal_fe_face_values_neighbor;
      }
    else
      {
        internal_fe_subface_values_neighbor.reinit(cell_neighbor, face_no_neighbor, sub_face_no_neighbor);
        fe_face_values_neighbor = &internal_fe_subface_values_neighbor;
      }
  }

  template <class Iterator>
  void reinit(const Iterator &cell)
  {
      fe_values.reinit (cell);
      fe_face_values = nullptr;
      fe_face_values_neighbor = nullptr;
  }

  double jump (const unsigned int idx, const unsigned int q_point) const
  {
      const FEFaceValuesBase<dim> &fe_v = *fe_face_values;
      const FEFaceValuesBase<dim> &fe_v_neighbor = *fe_face_values_neighbor;

      if (idx<fe_v.dofs_per_cell)
        return fe_v.shape_value(idx, q_point);
      else
        return -fe_v_neighbor.shape_value(idx-fe_v.dofs_per_cell, q_point);
    }

  Tensor<1,dim> jump (const unsigned int idx, const unsigned int q_point,
               const FEValuesExtractors::Vector &extractor) const
  {
      const FEFaceValuesBase<dim> &fe_v = *fe_face_values;
      const FEFaceValuesBase<dim> &fe_v_neighbor = *fe_face_values_neighbor;

      if (idx<fe_v.dofs_per_cell)
        return fe_v[extractor].value(idx, q_point);
      else
        return -fe_v_neighbor[extractor].value(idx-fe_v.dofs_per_cell, q_point);
    }


  double avg (const unsigned int idx, const unsigned int q_point) const
  {
      const FEFaceValuesBase<dim> &fe_v = *fe_face_values;
      const FEFaceValuesBase<dim> &fe_v_neighbor = *fe_face_values_neighbor;

      if (idx<fe_v.dofs_per_cell)
        return 0.5*fe_v.shape_value(idx, q_point);
      else
        return 0.5*fe_v_neighbor.shape_value(idx-fe_v.dofs_per_cell, q_point);
    }

  double avg (const unsigned int idx, const unsigned int q_point,
              const FEValuesExtractors::Scalar &extractor) const
  {
      const FEFaceValuesBase<dim> &fe_v = *fe_face_values;
      const FEFaceValuesBase<dim> &fe_v_neighbor = *fe_face_values_neighbor;

      if (idx<fe_v.dofs_per_cell)
        return 0.5*fe_v[extractor].value(idx, q_point);
      else
        return 0.5*fe_v_neighbor[extractor].value(idx-fe_v.dofs_per_cell, q_point);
    }

  Tensor<1,dim> gradient_avg (const unsigned int idx, const unsigned int q_point) const
  {
      const FEFaceValuesBase<dim> &fe_v = *fe_face_values;
      const FEFaceValuesBase<dim> &fe_v_neighbor = *fe_face_values_neighbor;

      if (idx<fe_v.dofs_per_cell)
        return 0.5*fe_v.shape_grad(idx, q_point);
      else
        return 0.5*fe_v_neighbor.shape_grad(idx-fe_v.dofs_per_cell, q_point);
    }

  Tensor<1,dim> gradient_n_avg (const unsigned int idx,
                                const unsigned int q_point,
                                const FEValuesExtractors::Vector &extractor
                                ) const
  {
      const FEFaceValuesBase<dim> &fe_v = *fe_face_values;
      const FEFaceValuesBase<dim> &fe_v_neighbor = *fe_face_values_neighbor;

      if (idx<fe_v.dofs_per_cell)
        return 0.5
                *fe_v[extractor].gradient(idx, q_point)
                *fe_v.normal_vector(q_point);
      else
        return 0.5
                *fe_v_neighbor[extractor].gradient(idx-fe_v.dofs_per_cell, q_point)
                *fe_v.normal_vector(q_point);
    }

  double choose (const bool left, const unsigned int idx, const unsigned int q_point) const
  {
      const FEFaceValuesBase<dim> &fe_v = *fe_face_values;
      const FEFaceValuesBase<dim> &fe_v_neighbor = *fe_face_values_neighbor;

      if (left && idx<fe_v.dofs_per_cell)
        return fe_v.shape_value(idx, q_point);
      if (!left && idx>=fe_v.dofs_per_cell)
        return fe_v_neighbor.shape_value(idx-fe_v.dofs_per_cell, q_point);
      return 0.0;
    };

  FEValues<dim>     fe_values;

  FEFaceValuesBase<dim> *fe_face_values;
  FEFaceValuesBase<dim> *fe_face_values_neighbor;

//  protected:
  FEFaceValues<dim> internal_fe_face_values;
  FESubfaceValues<dim> internal_fe_subface_values;
  FEFaceValues<dim> internal_fe_face_values_neighbor;
  FESubfaceValues<dim> internal_fe_subface_values_neighbor;

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
  std::vector<CopyDataFace> face_data;

};

template <class MatrixType, class VectorType>
inline void copy(const CopyData &c,
          const ConstraintMatrix &constraints,
          MatrixType & system_matrix,
          VectorType & system_rhs)
{
  constraints.distribute_local_to_global(c.cell_matrix,
                                         c.cell_rhs,
                                         c.local_dof_indices,
                                         system_matrix,
                                         system_rhs);
  for (auto &cdf : c.face_data)
    {
          // TODO: constraints?
      const unsigned int dofs_per_cell   = cdf.joint_dof_indices.size();
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          system_matrix.add(cdf.joint_dof_indices[i], cdf.joint_dof_indices[k],
                            cdf.cell_matrix(i,k));
    }
}



}
#endif
