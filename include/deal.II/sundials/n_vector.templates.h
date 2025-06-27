// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_sundials_n_vector_templates_h
#define dealii_sundials_n_vector_templates_h

#include <deal.II/base/config.h>

#include <deal.II/sundials/n_vector.h>
#include <deal.II/sundials/sundials_types.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/mpi.templates.h>

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/la_parallel_block_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/petsc_block_vector.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/vector_memory.h>

#  include <sundials/sundials_nvector.h>

#  include <limits>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  namespace internal
  {
    /**
     * An internal class that stores a pointer to a vector and manages the
     * memory if necessary. Objects of this class are used as the `content`
     * field in the SUNDIALS N_Vector module. In addition, this class has a
     * flag to store whether the stored vector should be treated as const. When
     * get() is called on a non-const object of this class the flag is checked
     * and an exception is thrown if the vector is actually const. Thus, we
     * retain a kind of "runtime const correctness" although the static const
     * correctness is lost because SUNDIALS N_Vector does not support constness.
     */
    template <typename VectorType>
    class NVectorContent
    {
    public:
      /**
       * Create a non-owning content with an existing @p vector.
       * @param vector The underlying vector to wrap in this object.
       */
      NVectorContent(VectorType *vector);

      /**
       * Create a non-owning content with an existing const @p vector. If this
       * constructor is used, access is only allowed via the get() const method.
       * @param vector The underlying vector to wrap in this object.
       */
      NVectorContent(const VectorType *vector);

      /**
       * Allocate a new (non-const) vector wrapped in a new content object. The
       * vector will be deallocated automatically when this object is destroyed.
       *
       * @note This constructor is intended for the N_VClone() call of SUNDIALS.
       */
      NVectorContent(const MPI_Comm comm);

      /**
       * Non-const access to the stored vector. Only allowed if a constructor
       * different than NVectorContent(const VectorType *vector) was used.
       * @return
       */
      VectorType *
      get();

      /**
       * Const access to the stored vector. Always allowed.
       */
      const VectorType *
      get() const;

      /**
       * Return a reference to a copy of the communicator the vector uses.
       * This function exists because the N_Vector
       * interface requires a function that returns a `void*` pointing
       * to the communicator object -- so somewhere, we need to have an
       * address to point to. The issue is that our vectors typically
       * return a *copy* of the communicator, rather than a reference to
       * the communicator they use, and so there is only a temporary
       * object and no address we can point to. To work around this
       * requirement, this class stores a copy of the communicator,
       * and this function here returns a reference to this copy.
       */
      const MPI_Comm &
      get_mpi_communicator() const;

    private:
      using PointerType =
        std::unique_ptr<VectorType, std::function<void(VectorType *)>>;

      /**
       * Vector memory which might be used to allocate storage if
       * NVectorContent() is called.
       */
      GrowingVectorMemory<VectorType> mem;

      /**
       * Actually stored vector content.
       */
      PointerType vector;

      /**
       * A copy of the communicator the vector uses, initialized in the
       * constructor of this class. We store this because the N_Vector
       * interface requires a function that returns a `void*` pointing
       * to the communicator object -- so somewhere, we need to have an
       * address to point to. The issue is that our vectors typically
       * return a *copy* of the communicator, rather than a reference to
       * the communicator they use, and so there is only a temporary
       * object and no address we can point to.
       */
      MPI_Comm comm;

      /**
       * Flag storing whether the stored pointer is to be treated as const. If
       * the pointer passed in the constructor was indeed const, it is cast away
       * but this flag will be set to true. Access to the pointer must then
       * check that this flag is set correctly.
       */
      const bool is_const;
    };

    /**
     * Helper to create a vector with all operations and the given @p content.
     *
     * @param content The vector content to attach to the N_Vector.
     * @return A new N_Vector
     */
    template <typename VectorType>
    N_Vector
    create_nvector(NVectorContent<VectorType> *content
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                   ,
                   SUNContext nvector_context
#  endif
    );

    /**
     * Helper to create an empty vector with all operations but no content.
     * @return A new N_Vector
     */
    template <typename VectorType>
    N_Vector
    create_empty_nvector(
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
      SUNContext nvector_context
#  endif
    );


    /**
     * Variable template that is true for distributed deal.II vectors and false
     * otherwise.
     */
    template <typename VectorType>
    constexpr bool is_dealii_compatible_distributed_vector =
      std::is_same_v<
        VectorType,
        LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                           MemorySpace::Host>> ||
      std::is_same_v<VectorType,
                     LinearAlgebra::distributed::BlockVector<
                       typename VectorType::value_type>>;

    template <typename VectorType,
              std::enable_if_t<!IsBlockVector<VectorType>::value> * = nullptr>
    unsigned int
    n_blocks(const VectorType &)
    {
      return 1;
    }



    template <typename VectorType,
              std::enable_if_t<IsBlockVector<VectorType>::value> * = nullptr>
    unsigned int
    n_blocks(const VectorType &vector)
    {
      return vector.n_blocks();
    }



    template <typename VectorType,
              std::enable_if_t<!IsBlockVector<VectorType>::value> * = nullptr>
    VectorType &
    block(VectorType &vector, const unsigned int b)
    {
      AssertDimension(b, 0);
      (void)b;
      return vector;
    }



    template <typename VectorType,
              std::enable_if_t<!IsBlockVector<VectorType>::value> * = nullptr>
    const VectorType &
    block(const VectorType &vector, const unsigned int b)
    {
      AssertDimension(b, 0);
      (void)b;
      return vector;
    }



    template <typename VectorType,
              std::enable_if_t<IsBlockVector<VectorType>::value> * = nullptr>
    typename VectorType::BlockType &
    block(VectorType &vector, const unsigned int b)
    {
      return vector.block(b);
    }



    template <typename VectorType,
              std::enable_if_t<IsBlockVector<VectorType>::value> * = nullptr>
    const typename VectorType::BlockType &
    block(const VectorType &vector, const unsigned int b)
    {
      return vector.block(b);
    }


    /**
     * Collection of all operations specified by SUNDIALS N_Vector
     * documentation. These functions are attached to the generic N_Vector
     * structure.
     */
    namespace NVectorOperations
    {
      N_Vector_ID
      get_vector_id(N_Vector v);

      N_Vector
      clone_empty(N_Vector w);

      template <typename VectorType>
      N_Vector
      clone(N_Vector w);

      template <typename VectorType>
      void
      destroy(N_Vector v);

      template <typename VectorType>
      sunindextype
      get_global_length(N_Vector v);

      template <typename VectorType>
      void
      linear_sum(SUNDIALS::realtype a,
                 N_Vector           x,
                 SUNDIALS::realtype b,
                 N_Vector           y,
                 N_Vector           z);

      template <
        typename VectorType,
        std::enable_if_t<is_dealii_compatible_distributed_vector<VectorType>>
          * = nullptr>
      int
      linear_combination(int nv, realtype *c, N_Vector *x, N_Vector z);

      template <typename VectorType>
      SUNDIALS::realtype
      dot_product(N_Vector x, N_Vector y);

      template <
        typename VectorType,
        std::enable_if_t<is_dealii_compatible_distributed_vector<VectorType>>
          * = nullptr>
      int
      dot_product_multi(int nv, N_Vector x, N_Vector *Y, realtype *d);

      template <
        typename VectorType,
        std::enable_if_t<is_dealii_compatible_distributed_vector<VectorType>>
          * = nullptr>
      int
      dot_product_multi_local(int nv, N_Vector x, N_Vector *y, realtype *d);


      template <
        typename VectorType,
        std::enable_if_t<is_dealii_compatible_distributed_vector<VectorType>>
          * = nullptr>
      int
      dot_product_multi_all_reduce(int nv, N_Vector x, realtype *d);

      template <typename VectorType>
      SUNDIALS::realtype
      weighted_l2_norm(N_Vector x, N_Vector y);

      template <typename VectorType>
      SUNDIALS::realtype
      l1_norm(N_Vector x);

      template <typename VectorType>
      void
      elementwise_product(N_Vector x, N_Vector y, N_Vector z);

      template <typename VectorType,
                std::enable_if_t<!IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_div(N_Vector x, N_Vector y, N_Vector z);

      template <typename VectorType,
                std::enable_if_t<IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_div(N_Vector x, N_Vector y, N_Vector z);

      template <typename VectorType,
                std::enable_if_t<!IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_inv(N_Vector x, N_Vector z);

      template <typename VectorType,
                std::enable_if_t<IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_inv(N_Vector x, N_Vector z);

      template <typename VectorType,
                std::enable_if_t<!IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_abs(N_Vector x, N_Vector z);

      template <typename VectorType,
                std::enable_if_t<IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_abs(N_Vector x, N_Vector z);

      template <typename VectorType>
      SUNDIALS::realtype
      weighted_rms_norm(N_Vector x, N_Vector w);

      template <typename VectorType>
      SUNDIALS::realtype
      weighted_rms_norm_mask(N_Vector x, N_Vector w, N_Vector mask);

      template <typename VectorType>
      SUNDIALS::realtype
      max_norm(N_Vector x);

      template <typename VectorType,
                std::enable_if_t<is_serial_vector<VectorType>::value, int> = 0>
      SUNDIALS::realtype
      min_element(N_Vector x);

      template <typename VectorType,
                std::enable_if_t<!is_serial_vector<VectorType>::value &&
                                   !IsBlockVector<VectorType>::value,
                                 int> = 0>
      SUNDIALS::realtype
      min_element(N_Vector x);

      template <typename VectorType,
                std::enable_if_t<!is_serial_vector<VectorType>::value &&
                                   IsBlockVector<VectorType>::value,
                                 int> = 0>
      SUNDIALS::realtype
      min_element(N_Vector x);

      template <typename VectorType>
      void
      scale(SUNDIALS::realtype c, N_Vector x, N_Vector z);

      template <typename VectorType>
      void
      set_constant(SUNDIALS::realtype c, N_Vector v);

      template <typename VectorType>
      void
      add_constant(N_Vector x, SUNDIALS::realtype b, N_Vector z);

      template <typename VectorType>
      const MPI_Comm &
      get_mpi_communicator(N_Vector v);

#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 0, 0)
      /**
       * Sundials likes their own communicator type by value.
       */
      template <typename VectorType>
      inline SUNComm
      get_mpi_communicator_by_value(N_Vector v);
#  else
      /**
       * Sundials likes a void* but we want to use the above functions
       * internally with a safe type.
       */
      template <typename VectorType>
      inline void *
      get_mpi_communicator_as_void_ptr(N_Vector v);
#  endif
    } // namespace NVectorOperations
  }   // namespace internal
} // namespace SUNDIALS



// -------------- Implementation of the functions declared above ---------

namespace SUNDIALS
{
  namespace internal
  {
    template <typename VectorType,
              std::enable_if_t<is_serial_vector<VectorType>::value, int> = 0>
    MPI_Comm
    get_mpi_communicator_from_vector(const VectorType &)
    {
      return MPI_COMM_SELF;
    }



    template <typename VectorType,
              std::enable_if_t<!is_serial_vector<VectorType>::value &&
                                 !IsBlockVector<VectorType>::value,
                               int> = 0>
    MPI_Comm
    get_mpi_communicator_from_vector(const VectorType &v)
    {
#  ifndef DEAL_II_WITH_MPI
      (void)v;
      return MPI_COMM_SELF;
#  else
      return v.get_mpi_communicator();
#  endif
    }



    template <typename VectorType,
              std::enable_if_t<!is_serial_vector<VectorType>::value &&
                                 IsBlockVector<VectorType>::value,
                               int> = 0>
    MPI_Comm
    get_mpi_communicator_from_vector(const VectorType &v)
    {
#  ifndef DEAL_II_WITH_MPI
      (void)v;
      return MPI_COMM_SELF;
#  else
      Assert(v.n_blocks() > 0,
             ExcMessage("You cannot ask a block vector without blocks "
                        "for its MPI communicator."));
      return v.block(0).get_mpi_communicator();
#  endif
    }


    template <typename VectorType>
    NVectorContent<VectorType>::NVectorContent(const MPI_Comm comm)
      : vector(typename VectorMemory<VectorType>::Pointer(mem))
      , comm(comm)
      , is_const(false)
    {}


    template <typename VectorType>
    NVectorContent<VectorType>::NVectorContent(VectorType *vector)
      : vector(vector,
               [](VectorType *) { /* not owning memory -> don't free*/ })
      , comm(get_mpi_communicator_from_vector(*vector))
      , is_const(false)
    {}



    template <typename VectorType>
    NVectorContent<VectorType>::NVectorContent(const VectorType *vector)
      : vector(const_cast<VectorType *>(vector),
               [](VectorType *) { /* not owning memory -> don't free*/ })
      , comm(get_mpi_communicator_from_vector(*vector))
      , is_const(true)
    {}



    template <typename VectorType>
    VectorType *
    NVectorContent<VectorType>::get()
    {
      AssertThrow(
        !is_const,
        ExcMessage(
          "Tried to access a constant vector content as non-const."
          " This most likely happened because a vector that is passed to a"
          " NVectorView() call should not be const.\n"
          "Alternatively, if you tried to access the vector, use"
          " unwrap_nvector_const() for this vector."));
      return vector.get();
    }



    template <typename VectorType>
    const VectorType *
    NVectorContent<VectorType>::get() const
    {
      return vector.get();
    }



    template <typename VectorType>
    const MPI_Comm &
    NVectorContent<VectorType>::get_mpi_communicator() const
    {
      return comm;
    }


#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    template <typename VectorType>
    NVectorView<VectorType>
    make_nvector_view(VectorType &vector, SUNContext nvector_context)
    {
      return NVectorView<VectorType>(vector, nvector_context);
    }
#  else
    template <typename VectorType>
    NVectorView<VectorType>
    make_nvector_view(VectorType &vector)
    {
      return NVectorView<VectorType>(vector);
    }
#  endif



    template <typename VectorType>
    VectorType *
    unwrap_nvector(N_Vector v)
    {
      Assert(v != nullptr, ExcInternalError());
      Assert(v->content != nullptr, ExcInternalError());
      auto *pContent =
        reinterpret_cast<NVectorContent<VectorType> *>(v->content);
      return pContent->get();
    }



    template <typename VectorType>
    const VectorType *
    unwrap_nvector_const(N_Vector v)
    {
      Assert(v != nullptr, ExcInternalError());
      Assert(v->content != nullptr, ExcInternalError());
      const auto *pContent =
        reinterpret_cast<const NVectorContent<VectorType> *>(v->content);
      return pContent->get();
    }



    template <typename VectorType>
    NVectorView<VectorType>::NVectorView(VectorType &vector
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                         ,
                                         SUNContext nvector_context
#  endif
                                         )
      : vector_ptr(
          create_nvector(new NVectorContent<std::remove_const_t<VectorType>>(
                           &vector)
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                           ,
                         nvector_context
#  endif
                         ),
          [](N_Vector v) { N_VDestroy(v); })
    {}



    template <typename VectorType>
    NVectorView<VectorType>::operator N_Vector() const
    {
      Assert(vector_ptr != nullptr, ExcNotInitialized());
      return vector_ptr.get();
    }



    template <typename VectorType>
    N_Vector
    NVectorView<VectorType>::operator->() const
    {
      Assert(vector_ptr != nullptr, ExcNotInitialized());
      return vector_ptr.get();
    }



    template <typename VectorType>
    N_Vector
    create_nvector(NVectorContent<VectorType> *content
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                   ,
                   SUNContext nvector_context
#  endif
    )
    {
      // Create an N_Vector with operators attached and empty content
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
      N_Vector v = create_empty_nvector<VectorType>();
#  else
      N_Vector v = create_empty_nvector<VectorType>(nvector_context);
#  endif
      Assert(v != nullptr, ExcInternalError());

      v->content = content;
      Assert(v->content != nullptr, ExcInternalError());
      return (v);
    }



    namespace NVectorOperations
    {
      N_Vector_ID
      get_vector_id(N_Vector)
      {
        return SUNDIALS_NVEC_CUSTOM;
      }



      N_Vector
      clone_empty(N_Vector w)
      {
        Assert(w != nullptr, ExcInternalError());
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
        N_Vector v = N_VNewEmpty();
#  else
        N_Vector v = N_VNewEmpty(w->sunctx);
#  endif
        Assert(v != nullptr, ExcInternalError());

        int status = N_VCopyOps(w, v);
        Assert(status == 0, ExcInternalError());
        (void)status;

        return v;
      }



      template <typename VectorType>
      N_Vector
      clone(N_Vector w)
      {
        N_Vector v = clone_empty(w);

        auto *w_dealii = unwrap_nvector_const<VectorType>(w);

        // Create the vector; the corresponding delete is called in destroy()
        auto cloned = new NVectorContent<VectorType>(
          get_mpi_communicator_from_vector(*w_dealii));

        // Then also copy the structure and values:
        *cloned->get() = *w_dealii;

        // Finally set the cloned object in 'v':
        v->content = cloned;
        return v;
      }



      template <typename VectorType>
      void
      destroy(N_Vector v)
      {
        // support destroying a nullptr because SUNDIALS vectors also do it
        if (v == nullptr)
          return;

        if (v->content != nullptr)
          {
            auto *content =
              reinterpret_cast<NVectorContent<VectorType> *>(v->content);
            // the NVectorContent knows if it owns the memory and will free
            // correctly
            delete content;
            v->content = nullptr;
          }

        N_VFreeEmpty(v);
      }



      template <typename VectorType>
      const MPI_Comm &
      get_mpi_communicator(N_Vector v)
      {
        Assert(v != nullptr, ExcInternalError());
        Assert(v->content != nullptr, ExcInternalError());
        auto *pContent =
          reinterpret_cast<NVectorContent<VectorType> *>(v->content);
        return pContent->get_mpi_communicator();
      }



#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 0, 0)
      template <typename VectorType>
      SUNComm
      get_mpi_communicator_by_value(N_Vector v)
      {
#    ifndef DEAL_II_WITH_MPI
        (void)v;
        return SUN_COMM_NULL;
#    else
        if (is_serial_vector<VectorType>::value == false)
          // SUNDIALS asks for a `SUNComm` object, which is a MPI-aware typedef
          // for `MPI_Comm` or `int`. To be clear, we adopt the SUNDIALS
          // interface here.
          //
          // Further, we need to cast away const here, as SUNDIALS demands the
          // communicator by value.
          return const_cast<SUNComm &>(get_mpi_communicator<VectorType>(v));
        else
          return SUN_COMM_NULL;
#    endif
      }
#  else
      template <typename VectorType>
      void *
      get_mpi_communicator_as_void_ptr(N_Vector v)
      {
#    ifndef DEAL_II_WITH_MPI
        (void)v;
        return nullptr;
#    else
        if (is_serial_vector<VectorType>::value == false)
          // We need to cast away const here, as SUNDIALS demands a pure
          // `void*`.
          return &(const_cast<MPI_Comm &>(get_mpi_communicator<VectorType>(v)));
        else
          return nullptr;
#    endif
      }
#  endif



      template <typename VectorType>
      sunindextype
      get_global_length(N_Vector v)
      {
        return unwrap_nvector_const<VectorType>(v)->size();
      }



      template <typename VectorType>
      void
      linear_sum(SUNDIALS::realtype a,
                 N_Vector           x,
                 SUNDIALS::realtype b,
                 N_Vector           y,
                 N_Vector           z)
      {
        auto *x_dealii = unwrap_nvector_const<VectorType>(x);
        auto *y_dealii = unwrap_nvector_const<VectorType>(y);
        auto *z_dealii = unwrap_nvector<VectorType>(z);

        if (z_dealii == x_dealii)
          z_dealii->sadd(a, b, *y_dealii);
        else if (z_dealii == y_dealii)
          z_dealii->sadd(b, a, *x_dealii);
        else
          {
            *z_dealii = 0;
            z_dealii->add(a, *x_dealii, b, *y_dealii);
          }
      }



      template <
        typename VectorType,
        std::enable_if_t<is_dealii_compatible_distributed_vector<VectorType>> *>
      int
      linear_combination(int nv, realtype *c, N_Vector *x, N_Vector z)
      {
        std::vector<const VectorType *> unwrapped_vectors(nv);
        for (int i = 0; i < nv; ++i)
          unwrapped_vectors[i] = unwrap_nvector_const<VectorType>(x[i]);
        auto *z_dealii = unwrap_nvector<VectorType>(z);

        // N.B. The first pointer may alias with z.
        for (unsigned int b = 0; b < n_blocks(*z_dealii); ++b)
          for (unsigned int j = 0; j < block(*z_dealii, b).locally_owned_size();
               ++j)
            {
              double temp = 0.;
              for (int i = 0; i < nv; ++i)
                temp += block(*unwrapped_vectors[i], b).local_element(j) * c[i];
              block(*z_dealii, b).local_element(j) = temp;
            }

        return 0;
      }



      template <typename VectorType>
      SUNDIALS::realtype
      dot_product(N_Vector x, N_Vector y)
      {
        return *unwrap_nvector_const<VectorType>(x) *
               *unwrap_nvector_const<VectorType>(y);
      }



      template <
        typename VectorType,
        std::enable_if_t<is_dealii_compatible_distributed_vector<VectorType>> *>
      int
      dot_product_multi(int nv, N_Vector x, N_Vector *y, realtype *d)
      {
        const int status = dot_product_multi_local<VectorType>(nv, x, y, d);
        if (status != 0)
          return status;

        return dot_product_multi_all_reduce<VectorType>(nv, x, d);
      }



      template <
        typename VectorType,
        std::enable_if_t<is_dealii_compatible_distributed_vector<VectorType>> *>
      int
      dot_product_multi_local(int nv, N_Vector x, N_Vector *y, realtype *d)
      {
        std::vector<const VectorType *> unwrapped_vectors(nv);
        for (int i = 0; i < nv; ++i)
          unwrapped_vectors[i] = unwrap_nvector_const<VectorType>(y[i]);
        const VectorType *x_dealii = unwrap_nvector_const<VectorType>(x);

        std::fill(d, d + nv, 0.);

        for (unsigned int b = 0; b < n_blocks(*x_dealii); ++b)
          for (unsigned int j = 0; j < block(*x_dealii, b).locally_owned_size();
               ++j)
            {
              for (int i = 0; i < nv; ++i)
                d[i] += block(*x_dealii, b).local_element(j) *
                        block(*unwrapped_vectors[i], b).local_element(j);
            }
        return 0;
      }



      template <
        typename VectorType,
        std::enable_if_t<is_dealii_compatible_distributed_vector<VectorType>> *>
      int
      dot_product_multi_all_reduce(int nv, N_Vector x, realtype *d)
      {
        ArrayView<realtype> products(d, nv);
        Utilities::MPI::sum(products,
                            get_mpi_communicator<VectorType>(x),
                            products);
        return 0;
      }



      template <typename VectorType>
      SUNDIALS::realtype
      weighted_l2_norm(N_Vector x, N_Vector w)
      {
        // TODO copy can be avoided by a custom kernel
        VectorType tmp      = *unwrap_nvector_const<VectorType>(x);
        auto      *w_dealii = unwrap_nvector_const<VectorType>(w);
        tmp.scale(*w_dealii);
        return tmp.l2_norm();
      }



      template <typename VectorType>
      SUNDIALS::realtype
      l1_norm(N_Vector x)
      {
        return unwrap_nvector_const<VectorType>(x)->l1_norm();
      }



      template <typename VectorType>
      void
      set_constant(SUNDIALS::realtype c, N_Vector v)
      {
        auto *v_dealii = unwrap_nvector<VectorType>(v);
        *v_dealii      = c;
      }



      template <typename VectorType>
      void
      add_constant(N_Vector x, SUNDIALS::realtype b, N_Vector z)
      {
        auto *x_dealii = unwrap_nvector_const<VectorType>(x);
        auto *z_dealii = unwrap_nvector<VectorType>(z);

        if (x_dealii == z_dealii)
          z_dealii->add(b);
        else
          {
            *z_dealii = *x_dealii;
            z_dealii->add(b);
          }
      }



      template <typename VectorType>
      void
      elementwise_product(N_Vector x, N_Vector y, N_Vector z)
      {
        auto *x_dealii = unwrap_nvector_const<VectorType>(x);
        auto *y_dealii = unwrap_nvector_const<VectorType>(y);
        auto *z_dealii = unwrap_nvector<VectorType>(z);
        if (z_dealii == x_dealii)
          z_dealii->scale(*y_dealii);
        else if (z_dealii == y_dealii)
          z_dealii->scale(*x_dealii);
        else
          {
            *z_dealii = *y_dealii;
            z_dealii->scale(*x_dealii);
          }
      }



      template <typename VectorType>
      SUNDIALS::realtype
      weighted_rms_norm(N_Vector x, N_Vector w)
      {
        // TODO copy can be avoided by a custom kernel
        VectorType tmp      = *unwrap_nvector_const<VectorType>(x);
        auto      *w_dealii = unwrap_nvector_const<VectorType>(w);
        const auto n        = tmp.size();
        tmp.scale(*w_dealii);
        return tmp.l2_norm() / std::sqrt(n);
      }



      template <typename VectorType>
      SUNDIALS::realtype
      weighted_rms_norm_mask(N_Vector x, N_Vector w, N_Vector mask)
      {
        // TODO copy can be avoided by a custom kernel
        VectorType tmp         = *unwrap_nvector_const<VectorType>(x);
        auto      *w_dealii    = unwrap_nvector_const<VectorType>(w);
        auto      *mask_dealii = unwrap_nvector_const<VectorType>(mask);
        const auto n           = tmp.size();
        tmp.scale(*w_dealii);
        tmp.scale(*mask_dealii);
        return tmp.l2_norm() / std::sqrt(n);
      }



      template <typename VectorType>
      SUNDIALS::realtype
      max_norm(N_Vector x)
      {
        return unwrap_nvector_const<VectorType>(x)->linfty_norm();
      }



      template <typename VectorType,
                std::enable_if_t<is_serial_vector<VectorType>::value, int>>
      SUNDIALS::realtype
      min_element(N_Vector x)
      {
        auto *vector = unwrap_nvector_const<VectorType>(x);
        return *std::min_element(vector->begin(), vector->end());
      }



      template <typename VectorType,
                std::enable_if_t<!is_serial_vector<VectorType>::value &&
                                   !IsBlockVector<VectorType>::value,
                                 int>>
      SUNDIALS::realtype
      min_element(N_Vector x)
      {
        auto *vector = unwrap_nvector_const<VectorType>(x);


        const auto indexed_less_than = [&](const IndexSet::size_type idxa,
                                           const IndexSet::size_type idxb) {
          return (*vector)[idxa] < (*vector)[idxb];
        };

        auto local_elements = vector->locally_owned_elements();

        const auto local_min = *std::min_element(local_elements.begin(),
                                                 local_elements.end(),
                                                 indexed_less_than);
        return Utilities::MPI::min((*vector)[local_min],
                                   get_mpi_communicator<VectorType>(x));
      }



      template <typename VectorType,
                std::enable_if_t<!is_serial_vector<VectorType>::value &&
                                   IsBlockVector<VectorType>::value,
                                 int>>
      SUNDIALS::realtype
      min_element(N_Vector x)
      {
        auto *vector = unwrap_nvector_const<VectorType>(x);

        // initialize local minimum to the largest possible value
        auto proc_local_min =
          std::numeric_limits<typename VectorType::value_type>::max();

        for (unsigned i = 0; i < vector->n_blocks(); ++i)
          {
            const auto indexed_less_than = [&](const IndexSet::size_type idxa,
                                               const IndexSet::size_type idxb) {
              return vector->block(i)[idxa] < vector->block(i)[idxb];
            };

            auto local_elements = vector->block(i).locally_owned_elements();

            const auto block_local_min_element =
              std::min_element(local_elements.begin(),
                               local_elements.end(),
                               indexed_less_than);

            // guard against empty blocks on this processor
            if (block_local_min_element != local_elements.end())
              proc_local_min =
                std::min(proc_local_min,
                         vector->block(i)[*block_local_min_element]);
          }

        return Utilities::MPI::min(proc_local_min,
                                   get_mpi_communicator<VectorType>(x));
      }



      template <typename VectorType>
      void
      scale(SUNDIALS::realtype c, N_Vector x, N_Vector z)
      {
        auto *x_dealii = unwrap_nvector_const<VectorType>(x);
        auto *z_dealii = unwrap_nvector<VectorType>(z);

        if (x_dealii == z_dealii)
          (*z_dealii) *= c;
        else
          z_dealii->sadd(0.0, c, *x_dealii);
      }



      template <typename VectorType,
                std::enable_if_t<!IsBlockVector<VectorType>::value, int>>
      void
      elementwise_div(N_Vector x, N_Vector y, N_Vector z)
      {
        auto *x_dealii = unwrap_nvector_const<VectorType>(x);
        auto *y_dealii = unwrap_nvector_const<VectorType>(y);
        auto *z_dealii = unwrap_nvector<VectorType>(z);

        AssertDimension(x_dealii->size(), z_dealii->size());
        AssertDimension(x_dealii->size(), y_dealii->size());

        auto x_ele = x_dealii->locally_owned_elements();
        for (const auto idx : x_ele)
          {
            (*z_dealii)[idx] = (*x_dealii)[idx] / (*y_dealii)[idx];
          }
        z_dealii->compress(VectorOperation::insert);
      }



      template <typename VectorType,
                std::enable_if_t<IsBlockVector<VectorType>::value, int>>
      void
      elementwise_div(N_Vector x, N_Vector y, N_Vector z)
      {
        auto *x_dealii = unwrap_nvector_const<VectorType>(x);
        auto *y_dealii = unwrap_nvector_const<VectorType>(y);
        auto *z_dealii = unwrap_nvector<VectorType>(z);

        AssertDimension(x_dealii->size(), z_dealii->size());
        AssertDimension(x_dealii->size(), y_dealii->size());
        AssertDimension(x_dealii->n_blocks(), z_dealii->n_blocks());
        AssertDimension(x_dealii->n_blocks(), y_dealii->n_blocks());

        for (unsigned i = 0; i < x_dealii->n_blocks(); ++i)
          {
            auto x_ele = x_dealii->block(i).locally_owned_elements();
            for (const auto idx : x_ele)
              {
                z_dealii->block(i)[idx] =
                  x_dealii->block(i)[idx] / y_dealii->block(i)[idx];
              }
          }
        z_dealii->compress(VectorOperation::insert);
      }



      template <typename VectorType,
                std::enable_if_t<!IsBlockVector<VectorType>::value, int>>
      void
      elementwise_inv(N_Vector x, N_Vector z)
      {
        auto *x_dealii = unwrap_nvector_const<VectorType>(x);
        auto *z_dealii = unwrap_nvector<VectorType>(z);

        AssertDimension(x_dealii->size(), z_dealii->size());

        auto x_ele = x_dealii->locally_owned_elements();
        for (const auto idx : x_ele)
          {
            (*z_dealii)[idx] = 1.0 / (*x_dealii)[idx];
          }
        z_dealii->compress(VectorOperation::insert);
      }



      template <typename VectorType,
                std::enable_if_t<IsBlockVector<VectorType>::value, int>>
      void
      elementwise_inv(N_Vector x, N_Vector z)
      {
        auto *x_dealii = unwrap_nvector_const<VectorType>(x);
        auto *z_dealii = unwrap_nvector<VectorType>(z);

        AssertDimension(x_dealii->size(), z_dealii->size());
        AssertDimension(x_dealii->n_blocks(), z_dealii->n_blocks());

        for (unsigned i = 0; i < x_dealii->n_blocks(); ++i)
          {
            auto x_ele = x_dealii->block(i).locally_owned_elements();
            for (const auto idx : x_ele)
              {
                z_dealii->block(i)[idx] = 1.0 / x_dealii->block(i)[idx];
              }
          }

        z_dealii->compress(VectorOperation::insert);
      }



      template <typename VectorType,
                std::enable_if_t<!IsBlockVector<VectorType>::value, int>>
      void
      elementwise_abs(N_Vector x, N_Vector z)
      {
        auto *x_dealii = unwrap_nvector_const<VectorType>(x);
        auto *z_dealii = unwrap_nvector<VectorType>(z);

        AssertDimension(x_dealii->size(), z_dealii->size());

        auto x_ele = x_dealii->locally_owned_elements();
        for (const auto idx : x_ele)
          {
            (*z_dealii)[idx] = std::fabs((*x_dealii)[idx]);
          }

        z_dealii->compress(VectorOperation::insert);
      }



      template <typename VectorType,
                std::enable_if_t<IsBlockVector<VectorType>::value, int>>
      void
      elementwise_abs(N_Vector x, N_Vector z)
      {
        auto *x_dealii = unwrap_nvector_const<VectorType>(x);
        auto *z_dealii = unwrap_nvector<VectorType>(z);

        AssertDimension(x_dealii->size(), z_dealii->size());
        AssertDimension(x_dealii->n_blocks(), z_dealii->n_blocks());

        for (unsigned i = 0; i < x_dealii->n_blocks(); ++i)
          {
            auto x_ele = x_dealii->block(i).locally_owned_elements();
            for (const auto idx : x_ele)
              {
                z_dealii->block(i)[idx] = std::fabs(x_dealii->block(i)[idx]);
              }
          }

        z_dealii->compress(VectorOperation::insert);
      }

    } // namespace NVectorOperations


    template <typename VectorType>
    N_Vector
    create_empty_nvector(
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
      SUNContext nvector_context
#  endif
    )
    {
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
      N_Vector v = N_VNewEmpty();
#  else
      N_Vector v = N_VNewEmpty(nvector_context);
#  endif
      Assert(v != nullptr, ExcInternalError());

      /* constructors, destructors, and utility operations */
      v->ops->nvgetvectorid = &NVectorOperations::get_vector_id;
      v->ops->nvclone       = &NVectorOperations::clone<VectorType>;
      v->ops->nvcloneempty  = &NVectorOperations::clone_empty;
      v->ops->nvdestroy     = &NVectorOperations::destroy<VectorType>;
      //  v->ops->nvspace           = undef;
#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 0, 0)
      v->ops->nvgetcommunicator =
        &NVectorOperations::get_mpi_communicator_by_value<VectorType>;
#  else
      v->ops->nvgetcommunicator =
        &NVectorOperations::get_mpi_communicator_as_void_ptr<VectorType>;
#  endif
      v->ops->nvgetlength = &NVectorOperations::get_global_length<VectorType>;

      /* standard vector operations */
      v->ops->nvlinearsum = &NVectorOperations::linear_sum<VectorType>;
      v->ops->nvconst     = &NVectorOperations::set_constant<VectorType>;
      v->ops->nvprod      = &NVectorOperations::elementwise_product<VectorType>;
      v->ops->nvdiv       = &NVectorOperations::elementwise_div<VectorType>;
      v->ops->nvscale     = &NVectorOperations::scale<VectorType>;
      v->ops->nvabs       = &NVectorOperations::elementwise_abs<VectorType>;
      v->ops->nvinv       = &NVectorOperations::elementwise_inv<VectorType>;
      v->ops->nvaddconst  = &NVectorOperations::add_constant<VectorType>;
      v->ops->nvdotprod   = &NVectorOperations::dot_product<VectorType>;
      v->ops->nvmaxnorm   = &NVectorOperations::max_norm<VectorType>;
      v->ops->nvwrmsnorm  = &NVectorOperations::weighted_rms_norm<VectorType>;
      v->ops->nvmin       = &NVectorOperations::min_element<VectorType>;
      v->ops->nvwl2norm   = &NVectorOperations::weighted_l2_norm<VectorType>;
      v->ops->nvl1norm    = &NVectorOperations::l1_norm<VectorType>;
      v->ops->nvwrmsnormmask =
        &NVectorOperations::weighted_rms_norm_mask<VectorType>;

      /* The following are declared as standard by SUNDIALS but were not
       * necessary so far.
       */
      v->ops->nvcompare     = nullptr;
      v->ops->nvinvtest     = nullptr;
      v->ops->nvconstrmask  = nullptr;
      v->ops->nvminquotient = nullptr;

      /* fused and vector array operations are disabled by default */

      /* OPTIONAL fused vector operations */

      // These operations are only implemented for deal.II vectors. SUNDIALS
      // will automatically fall back to a less optimized version implemented
      // through the standard functions above for other VectorTypes.
      if constexpr (is_dealii_compatible_distributed_vector<VectorType>)
        {
          v->ops->nvlinearcombination =
            &NVectorOperations::linear_combination<VectorType>;
          v->ops->nvdotprodmulti =
            &NVectorOperations::dot_product_multi<VectorType>;
        }
      else
        {
          v->ops->nvlinearcombination = nullptr;
          v->ops->nvdotprodmulti      = nullptr;
        }

      v->ops->nvscaleaddmulti = nullptr;

      /* OPTIONAL vector array operations */
      v->ops->nvlinearsumvectorarray         = nullptr;
      v->ops->nvscalevectorarray             = nullptr;
      v->ops->nvconstvectorarray             = nullptr;
      v->ops->nvwrmsnormvectorarray          = nullptr;
      v->ops->nvwrmsnormmaskvectorarray      = nullptr;
      v->ops->nvscaleaddmultivectorarray     = nullptr;
      v->ops->nvlinearcombinationvectorarray = nullptr;

      /* Local reduction kernels (no parallel communication) */
      v->ops->nvdotprodlocal     = nullptr;
      v->ops->nvmaxnormlocal     = nullptr;
      v->ops->nvminlocal         = nullptr;
      v->ops->nvl1normlocal      = nullptr;
      v->ops->nvinvtestlocal     = nullptr;
      v->ops->nvconstrmasklocal  = nullptr;
      v->ops->nvminquotientlocal = nullptr;
      v->ops->nvwsqrsumlocal     = nullptr;
      v->ops->nvwsqrsummasklocal = nullptr;

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
      /* Single buffer reduction operations */

      // These operations are only implemented for deal.II vectors. SUNDIALS
      // will automatically fall back to a less optimized version implemented
      // through the standard functions above for other VectorTypes.
      if constexpr (is_dealii_compatible_distributed_vector<VectorType>)
        {
          v->ops->nvdotprodmultilocal =
            &NVectorOperations::dot_product_multi_local<VectorType>;
          v->ops->nvdotprodmultiallreduce =
            &NVectorOperations::dot_product_multi_all_reduce<VectorType>;
        }
      else
        {
          v->ops->nvdotprodmultilocal     = nullptr;
          v->ops->nvdotprodmultiallreduce = nullptr;
        }
#  endif

      /* XBraid interface operations */
      v->ops->nvbufsize   = nullptr;
      v->ops->nvbufpack   = nullptr;
      v->ops->nvbufunpack = nullptr;

      /* Debugging functions (called when SUNDIALS_DEBUG_PRINTVEC is defined).
       */
      v->ops->nvprint     = nullptr;
      v->ops->nvprintfile = nullptr;

      return v;
    }

  } // namespace internal
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif
#endif
