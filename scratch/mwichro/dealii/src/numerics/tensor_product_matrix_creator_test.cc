#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/numerics/tensor_product_matrix_creator.h>

#include <numeric>


DEAL_II_NAMESPACE_OPEN


namespace TensorProductMatrixCreator
{


  FullMatrix<double>
  create_1d_cell_mass_matrix(const FiniteElement<1>     &fe,
                             const double               &h,
                             const std::pair<bool, bool> include_endpoints,
                             std::vector<unsigned int>   numbering)
  {
    if (dynamic_cast<const FE_DGQ<1> *>(&fe) == nullptr &&
        numbering.size() == 0)
      {
        Assert(
          include_endpoints.first == true && include_endpoints.second == true,
          ExcMessage(
            "You tried to generate a 1D mass matrix with excluding boundary "
            "dofs for a non-DGQ element without providing a numbering."));
      }


    if (numbering.size() == 0)
      {
        numbering.resize(fe.dofs_per_cell);
        std::iota(numbering.begin(), numbering.end(), 0);
      }
    const unsigned int degree          = fe.degree;
    const unsigned int n_dofs_per_cell = fe.dofs_per_cell;
    const double      &JxW             = h;
    QGauss<1>          quadrature(degree + 1);

    FullMatrix<double> cell_mass_matrix(n_dofs_per_cell, n_dofs_per_cell);
    cell_mass_matrix = 0;

    unsigned int start_dof = include_endpoints.first ? 0 : 1;
    unsigned int end_dof =
      include_endpoints.second ? n_dofs_per_cell : n_dofs_per_cell - 1;
    const unsigned int shift = include_endpoints.first ? 0 : 1;

    for (unsigned int i = start_dof; i < end_dof; ++i)
      for (unsigned int j = start_dof; j < end_dof; ++j)
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          cell_mass_matrix(i - shift, j - shift) +=
            (fe.shape_value(numbering[i], quadrature.point(q)) *
             fe.shape_value(numbering[j], quadrature.point(q))) *
            JxW * quadrature.weight(q);

    return cell_mass_matrix;
  }


} // namespace TensorProductMatrixCreator


DEAL_II_NAMESPACE_CLOSE