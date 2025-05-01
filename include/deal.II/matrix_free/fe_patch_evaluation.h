#ifndef FE_PATCH_EVALUATION_H
#define FE_PATCH_EVALUATION_H

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/patch_storage.h>

#include <algorithm>
#include <vector>


DEAL_II_NAMESPACE_OPEN

using namespace dealii;

enum VectorizationType
{
  over_patches,
  within_patch
};

namespace internal
{
  // Taken from:
  // https://stackoverflow.com/questions/17923683/why-does-stdarray-not-have-an-constructor-that-takes-a-value-for-the-array-to
  template <typename T, std::size_t... Is>
  std::array<T, sizeof...(Is)>
  make_array(const T &value, std::index_sequence<Is...>)
  {
    return {{(static_cast<void>(Is), value)...}};
  }

  template <std::size_t N, typename T>
  std::array<T, N>
  make_array(const T &value)
  {
    return make_array(value, std::make_index_sequence<N>());
  }


  template <template <int, int, int, int, typename, typename> typename T,
            int dim,
            int degree,
            int n_components,
            int n_q_points,
            typename NumberType,
            typename VectorizedNumberType>
  constexpr int
  extract_degree(const T<dim,
                         degree,
                         n_components,
                         n_q_points,
                         NumberType,
                         VectorizedNumberType> *)
  {
    return degree;
  }


} // namespace internal

template <typename FEEval,
          typename Distributor,
          VectorizationType vectorizaton = within_patch>
class FEPatchEvaluation
{
public:
  using FEEvaluationType            = FEEval;
  using VectorizedNumber            = typename FEEvaluationType::NumberType;
  using Number                      = typename FEEvaluationType::number_type;
  const static unsigned int dim     = FEEvaluationType::dimension;
  const static unsigned int n_lanes = FEEvaluationType::n_lanes;

  const static unsigned int n_cells = 1 << dim;

  static constexpr int fe_degree =
    internal::extract_degree(static_cast<FEEvaluationType *>(nullptr));
  const static unsigned int n_dofs_per_cell = std::pow(fe_degree + 1, dim);

  const static unsigned int n_evaluators = 1;

  using PatchStorageType = PatchStorage<MatrixFree<dim, Number>>;
  using CellDofsViewRaw  = StridedArrayView<Number, n_lanes>;

  FEPatchEvaluation(const PatchStorageType &storage, const FEEval &fe_eval);

  void
  reinit(const std::size_t patch_index);

  constexpr unsigned int
  n_patch_dofs() const;

  const unsigned int &
  get_current_patch_index() const;

  template <typename NumberType>
  void
  distribute_patch_to_local(const ArrayView<const NumberType> &patch_vector,
                            const bool                         copy_duplicates);


  template <typename NumberType>
  void
  gather_local_to_patch(const ArrayView<NumberType> &patch_vector,
                        const bool                   sum_overlapping) const;

  template <typename VECTOR>
  void
  read_dof_values(const VECTOR &src);

  template <typename VECTOR>
  void
  distribute_local_to_global(VECTOR &dst) const;

  auto &
  get_dof_value(const unsigned int &cell, const unsigned int &dof);
  const auto &
  get_dof_value(const unsigned int &cell, const unsigned int &dof) const;


  std::array<FEEvaluationType, n_evaluators> fe_evaluations;

private:
  std::array<CellDofsViewRaw, n_cells> cell_dofs_view_raw;

  const PatchStorageType &storage;
  unsigned int            current_patch_index;

  Distributor distributor;
};



namespace internal
{
  // Helper using index_sequence to initialize std::array elements in-place
  template <typename FEEval,
            typename Number,
            unsigned int NLanes,
            std::size_t... Is>
  std::array<StridedArrayView<Number, NLanes>, sizeof...(Is)>
  create_cell_dofs_view_array_impl(
    FEEval      &fe_eval,
    unsigned int NDofsPerCell, // Pass non-const reference
    std::index_sequence<Is...>)
  {
    // Assumes StridedArrayView constructor takes (pointer, count)
    // Uses logic inspired by the commented-out code in reinit.
    return {
      {StridedArrayView<Number, NLanes>(&(fe_eval.begin_dof_values()[0][Is]),
                                        NDofsPerCell)...}};
  }

  // Helper function to create the array of StridedArrayViews
  template <typename FEEval,
            unsigned int NCells,
            typename Number,
            unsigned int NLanes>
  std::array<StridedArrayView<Number, NLanes>, NCells>
  create_cell_dofs_view_array(
    FEEval      &fe_eval,
    unsigned int NDofsPerCell) // Pass non-const reference
  {
    return create_cell_dofs_view_array_impl<FEEval, Number, NLanes>(
      fe_eval, NDofsPerCell, std::make_index_sequence<NCells>());
  }

} // namespace internal


template <typename FEEval, typename Distributor, VectorizationType vectorizaton>
FEPatchEvaluation<FEEval, Distributor, vectorizaton>::FEPatchEvaluation(
  const PatchStorageType &storage,
  const FEEval           &fe_eval)
  : fe_evaluations(internal::make_array<n_evaluators>(fe_eval))
  , cell_dofs_view_raw(
      internal::
        create_cell_dofs_view_array<FEEvaluationType, n_cells, Number, n_lanes>(
          fe_evaluations[0],
          n_dofs_per_cell))
  , storage(storage)
  , current_patch_index(numbers::invalid_unsigned_int)
{
  static_assert(n_evaluators == 1, "Only one evaluator is supported for now");
  Assert(storage.get_matrix_free().get() == &fe_eval.get_matrix_free(),
         ExcMessage("MatrixFree objects do not match!"));
}


template <typename FEEval, typename Distributor, VectorizationType vectorizaton>
void
FEPatchEvaluation<FEEval, Distributor, vectorizaton>::reinit(
  const std::size_t patch_index)
{
  static_assert(n_evaluators == 1, "Only one evaluator is supported for now");
  fe_evaluations[0].reinit(storage.get_regular_patch(patch_index).get_cells());
  current_patch_index = patch_index;

  for (unsigned int i = 0; i < n_cells; ++i)
    {
      // Assert that the view points to the correct memory location by
      // checking the address of the first element accessible through the
      // view's operator[].
      Assert(
        &cell_dofs_view_raw[i][0] ==
          &(fe_evaluations[0].begin_dof_values()[0][i]),
        ExcInternalError(
          "cell_dofs_view_raw does not point to the expected memory location after reinit."));
    }
  // The lexicographic view uses the updated cell_dofs_view_raw, so it should
  // reflect the changes.
}


template <typename FEEval, typename Distributor, VectorizationType vectorizaton>
inline constexpr unsigned int
FEPatchEvaluation<FEEval, Distributor, vectorizaton>::n_patch_dofs() const
{
  return distributor.n_patch_dofs();
}

template <typename FEEval, typename Distributor, VectorizationType vectorizaton>
const unsigned int &
FEPatchEvaluation<FEEval, Distributor, vectorizaton>::get_current_patch_index()
  const
{
  return current_patch_index;
}


template <typename FEEval, typename Distributor, VectorizationType vectorizaton>
auto &
FEPatchEvaluation<FEEval, Distributor, vectorizaton>::get_dof_value(
  const unsigned int &cell,
  const unsigned int &dof)
{
  static_assert(n_evaluators == 1, "Only one evaluator is supported for now");
  static_assert(vectorizaton == within_patch,
                "Only one evaluator is supported for now");

  return fe_evaluations[0].begin_dof_values()[dof][cell];
}


template <typename FEEval, typename Distributor, VectorizationType vectorizaton>
const auto &
FEPatchEvaluation<FEEval, Distributor, vectorizaton>::get_dof_value(
  const unsigned int &cell,
  const unsigned int &dof) const
{
  static_assert(n_evaluators == 1, "Only one evaluator is supported for now");
  static_assert(vectorizaton == within_patch,
                "Only one evaluator is supported for now");

  return fe_evaluations[0].begin_dof_values()[dof][cell];
}

template <typename FEEval, typename Distributor, VectorizationType vectorizaton>
template <typename NumberType>
inline void
FEPatchEvaluation<FEEval, Distributor, vectorizaton>::distribute_patch_to_local(
  const ArrayView<const NumberType> &patch_vector,
  const bool                         copy_duplicates)
{
  AssertDimension(patch_vector.size(), distributor.n_patch_dofs());

  distributor.distribute_patch_to_local(cell_dofs_view_raw,
                                        patch_vector,
                                        copy_duplicates);
}

template <typename FEEval, typename Distributor, VectorizationType vectorizaton>
template <typename NumberType>
inline void
FEPatchEvaluation<FEEval, Distributor, vectorizaton>::gather_local_to_patch(
  const ArrayView<NumberType> &patch_vector,
  const bool                   sum_overlapping) const
{
  AssertDimension(patch_vector.size(), distributor.n_patch_dofs());

  distributor.gather_local_to_patch(cell_dofs_view_raw,
                                    patch_vector,
                                    sum_overlapping);
}

template <typename FEEval, typename Distributor, VectorizationType vectorizaton>
template <typename VECTOR>
void
FEPatchEvaluation<FEEval, Distributor, vectorizaton>::read_dof_values(
  const VECTOR &src)
{
  Assert(current_patch_index != numbers::invalid_unsigned_int,
         ExcNotInitialized());
  for (auto &fe_eval : fe_evaluations)
    fe_eval.read_dof_values(src);
}

template <typename FEEval, typename Distributor, VectorizationType vectorizaton>
template <typename VECTOR>
void
FEPatchEvaluation<FEEval, Distributor, vectorizaton>::
  distribute_local_to_global(VECTOR &dst) const
{
  Assert(current_patch_index != numbers::invalid_unsigned_int,
         ExcNotInitialized());
  for (auto &fe_eval : fe_evaluations)
    fe_eval.distribute_local_to_global(dst);
}

DEAL_II_NAMESPACE_CLOSE

#endif // FE_PATCH_EVALUATION_H