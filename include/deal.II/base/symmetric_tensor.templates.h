// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#ifndef dealii_symmetric_tensor_templates_h
#define dealii_symmetric_tensor_templates_h


#include <deal.II/base/symmetric_tensor.h>

#include <array>

DEAL_II_NAMESPACE_OPEN



template <typename Number>
std::array<Number,1>
eigenvalues (const SymmetricTensor<2,1,Number> &T)
{
  return { {T[0][0]} };
}



template <typename Number>
std::array<Number,2>
eigenvalues (const SymmetricTensor<2,2,Number> &T)
{
  const Number upp_tri_sq = T[0][1]*T[0][1];
  if (upp_tri_sq == Number(0.0))
    {
      // The tensor is diagonal
      std::array<Number,2> eig_vals =
      {
        {T[0][0], T[1][1]}
      };

      // Sort from largest to smallest.
      std::sort(eig_vals.begin(), eig_vals.end(), std::greater<Number>());
      return eig_vals;
    }
  else
    {
      const Number tr_T = trace(T);
      const Number det_T = determinant(T);
      const Number descrim = tr_T*tr_T - 4.0*det_T;
      Assert(descrim > Number(0.0), ExcMessage("The roots of the characteristic polynomial are complex valued."));
      const Number sqrt_desc = std::sqrt(descrim);

      const std::array<Number,2> eig_vals =
      {
        {
          static_cast<Number>(0.5*(tr_T + sqrt_desc)),
          static_cast<Number>(0.5*(tr_T - sqrt_desc))
        }
      };
      Assert(eig_vals[0] >= eig_vals[1], ExcMessage("The eigenvalue ordering is incorrect."));
      return eig_vals;
    }
}



template <typename Number>
std::array<Number,3>
eigenvalues (const SymmetricTensor<2,3,Number> &T)
{
  const Number upp_tri_sq = T[0][1]*T[0][1] + T[0][2]*T[0][2] + T[1][2]*T[1][2];
  if (upp_tri_sq == Number(0.0))
    {
      // The tensor is diagonal
      std::array<Number,3> eig_vals
      = { {T[0][0], T[1][1], T[2][2]} };

      // Sort from largest to smallest.
      std::sort(eig_vals.begin(), eig_vals.end(), std::greater<Number>());
      return eig_vals;
    }
  else
    {
      // Perform an affine change to T, and solve a different
      // characteristic equation that has a trigonometric solution.
      // Decompose T = p*B + q*I , and set q = tr(T)/3
      // and p = (tr((T - q.I)^{2})/6)^{1/2} . Then solve the equation
      // 0 = det(\lambda*I - B) = \lambda^{3} - 3*\lambda - det(B)
      // which has the solution
      // \lambda = 2*cos(1/3 * acos(det(B)/2) +2/3*pi*k ) ; k = 0,1,2
      // when substituting  \lambda = 2.cos(theta) and using trig identities.
      const Number tr_T = trace(T);
      const Number q = tr_T/3.0;
      const Number tmp1 = (  T[0][0] - q)*(T[0][0] - q)
                          + (T[1][1] - q)*(T[1][1] - q)
                          + (T[2][2] - q)*(T[2][2] - q)
                          + 2.0 * upp_tri_sq;
      const Number p = std::sqrt(tmp1/6.0);
      const SymmetricTensor<2,3,Number> B = Number(1.0/p)*(T - q*unit_symmetric_tensor<3,Number>());
      const Number tmp_2 = determinant(B)/2.0;

      // The value of tmp_2 should be within [-1,1], however
      // floating point errors might place it slightly outside
      // this range. It is therefore necessary to correct for it
      const Number phi =
        (tmp_2 <= -1.0 ? Number(numbers::PI/3.0) :
         (tmp_2 >= 1.0 ? Number(0.0) :
          std::acos(tmp_2)/3.0));

      // Due to the trigonometric solution, the computed eigenvalues
      // should be predictably in the order eig1 >= eig2 >= eig3...
      std::array<Number,3> eig_vals
      = { {
          static_cast<Number>(q + 2.0*p*std::cos(phi)),
          static_cast<Number>(0.0),
          static_cast<Number>(q + 2.0*p*std::cos(phi + (2.0/3.0*numbers::PI)))
        }
      };
      // Use the identity tr(T) = eig1 + eig2 + eig3
      eig_vals[1] = tr_T - eig_vals[0] - eig_vals[2];

      // ... however, when equal roots exist then floating point
      // errors may make this no longer be the case.
      // Sort from largest to smallest.
      std::sort(eig_vals.begin(), eig_vals.end(), std::greater<Number>());

      return eig_vals;
    }
}



namespace internal
{
  namespace SymmetricTensor
  {
    template <int dim, typename Number>
    void
    tridiagonalize (const dealii::SymmetricTensor<2,dim,Number> &A,
                    dealii::Tensor<2,dim,Number>                &Q,
                    std::array<Number,dim>                      &d,
                    std::array<Number,dim-1>                    &e)
    {
      // Create some intermediate storage
      Number h,g,omega_inv,K,f;

      // Initialize the transformation matrix as the
      // identity tensor
      Q = dealii::unit_symmetric_tensor<dim,Number>();

      // Make the first row and column to be of the
      // desired form
      h = 0.0;
      for (int i=1; i < dim; i++)
        h += A[0][i]*A[0][i];

      g = 0.0;
      if (A[0][1] > 0.0)
        g = -std::sqrt(h);
      else
        g = std::sqrt(h);
      e[0] = g;

      std::array<Number,dim> u;
      for (int i=1; i < dim; i++)
        {
          u[i] = A[0][i];
          if (i == 1)
            u[i] -= g;
        }

      std::array<Number,dim> q;
      const Number omega = h - g * A[0][1];
      if (omega > 0.0)
        {
          omega_inv = 1.0 / omega;
          K = 0.0;
          for (int i=1; i < dim; i++)
            {
              f = 0.0;
              for (int j=1; j < dim; j++)
                f += A[i][j] * u[j];
              q[i] = omega_inv * f;
              K   += u[i] * f;
            }
          K *= 0.5*omega_inv*omega_inv;

          for (int i=1; i < dim; i++)
            q[i] = q[i] - K * u[i];

          d[0] = A[0][0];
          for (int i=1; i < dim; i++)
            d[i] = A[i][i] - 2.0*q[i]*u[i];

          // Store inverse Householder transformation
          // in Q
          for (int j=1; j < dim; j++)
            {
              f = omega_inv * u[j];
              for (int i=1; i < dim; i++)
                Q[i][j] = Q[i][j] - f*u[i];
            }

          // For dim = 3: Calculate updated A[1][2] and
          // store it in e[1]
          for (int i=1; i < dim-1; i++)
            e[i] = A[i][i+1] - q[i]*u[i+1] - u[i]*q[i+1];
        }
      else
        {
          for (int i=0; i < dim; i++)
            d[i] = A[i][i];

          // For dim = 3:
          for (int i=1; i < dim-1; i++)
            e[i] = A[i][i+1];
        }
    }



    template <int dim, typename Number>
    std::array<std::pair<Number, Tensor<1,dim,Number> >,dim>
    ql_implicit_shifts (const dealii::SymmetricTensor<2,dim,Number> &A)
    {
      static_assert(numbers::NumberTraits<Number>::is_complex == false,
                    "This implementation of the QL implicit shift algorithm does "
                    "not support complex numbers");

      // Transform A to real tridiagonal form by the Householder method:
      // The orthogonal matrix effecting the transformation
      // this will ultimately store the eigenvectors
      dealii::Tensor<2,dim,Number> Q;
      // The diagonal elements of the tridiagonal matrix;
      // this will ultimately store the eigenvalues
      std::array<Number,dim>   w;
      // The off-diagonal elements of the tridiagonal
      std::array<Number,dim-1> ee;
      tridiagonalize<dim,Number>(A, Q, w, ee);

      // Number of iterations
      const unsigned int max_n_it = 30;

      // Transfer the off-diagonal entries to an auxiliary array
      // The third element is used only as temporary workspace
      std::array<Number,dim> e;
      for (unsigned int i=0; i<dim-1; ++i)
        e[i] = ee[i];

      // Create some intermediate storage
      Number g, r, p, f, b, s, c, t;

      // Loop over all off-diagonal elements
      for (int l=0; l < dim-1; l++)
        {
          for (unsigned int it=0; it <= max_n_it; ++it)
            {
              // Check for convergence and exit iteration loop
              // if the off-diagonal element e[l] is zero
              int m = l;
              for (; m <= dim-2; m++)
                {
                  g = std::abs(w[m]) + std::abs(w[m+1]);
                  if (std::abs(e[m]) + g == g)
                    break;
                }
              if (m == l)
                break;

              // Throw if no convergence is achieved within a
              // stipulated number of iterations
              if (it == max_n_it)
                {
                  AssertThrow(false, ExcMessage("No convergence in iterative QL eigenvector algorithm."))
                  return std::array<std::pair<Number, Tensor<1,dim,Number> >,dim> ();
                }

              // Calculate the shift..
              g = (w[l+1] - w[l]) / (e[l] + e[l]);
              r = std::sqrt(g*g + 1.0);
              // .. and then compute g = d_m - k_s for the
              // plane rotation (Press2007a eq 11.4.22)
              if (g > 0.0)
                g = w[m] - w[l] + e[l]/(g + r);
              else
                g = w[m] - w[l] + e[l]/(g - r);

              // Perform plane rotation, as is done in the
              // standard QL algorithm, followed by Givens
              // rotations to recover the tridiagonal form
              s = c = 1.0;
              p = 0.0;
              for (int i=m-1; i >= l; i--)
                {
                  f = s * e[i];
                  b = c * e[i];

                  // Branch to recover from underflow
                  if (std::abs(f) > std::abs(g))
                    {
                      c      = g / f;
                      r      = std::sqrt(c*c + 1.0);
                      e[i+1] = f * r;
                      c     *= (s = 1.0/r);
                    }
                  else
                    {
                      s      = f / g;
                      r      = std::sqrt(s*s + 1.0);
                      e[i+1] = g * r;
                      s     *= (c = 1.0/r);
                    }

                  g = w[i+1] - p;
                  r = (w[i] - g)*s + 2.0*c*b;
                  p = s * r;
                  w[i+1] = g + p;
                  g = c*r - b;

                  // Form the eigenvectors
                  for (int k=0; k < dim; k++)
                    {
                      t = Q[k][i+1];
                      Q[k][i+1] = s*Q[k][i] + c*t;
                      Q[k][i]   = c*Q[k][i] - s*t;
                    }
                }
              w[l] -= p;
              e[l]  = g;
              e[m]  = 0.0;
            }
        }

      // Structure the data to be outputted
      std::array<std::pair<Number, Tensor<1,dim,Number> >,dim> eig_vals_vecs;
      for (unsigned int e=0; e<dim; ++e)
        {
          eig_vals_vecs[e].first = w[e];

          // The column "e" of Q contains the non-normalized
          // eigenvector associated with the eigenvalue "e"
          for (unsigned int a=0; a<dim; ++a)
            {
              eig_vals_vecs[e].second[a] = Q[a][e];
            }

          // Normalize
          Assert(eig_vals_vecs[e].second.norm() != 0.0, ExcDivideByZero());
          eig_vals_vecs[e].second /= eig_vals_vecs[e].second.norm();
        }
      return eig_vals_vecs;
    }



    template <int dim, typename Number>
    std::array<std::pair<Number, Tensor<1,dim,Number> >,dim>
    jacobi (dealii::SymmetricTensor<2,dim,Number> A)
    {
      static_assert(numbers::NumberTraits<Number>::is_complex == false,
                    "This implementation of the Jacobi algorithm does "
                    "not support complex numbers");

      // Sums of diagonal resp. off-diagonal elements
      Number sd, so;
      // sin(phi), cos(phi), tan(phi) and temporary storage
      Number s, c, t;
      // More temporary storage
      Number g, h, z, theta;
      // Threshold value
      Number thresh;

      // Initialize the transformation matrix as the
      // identity tensor
      dealii::Tensor<2,dim,Number> Q (dealii::unit_symmetric_tensor<dim,Number>());

      // The diagonal elements of the tridiagonal matrix;
      // this will ultimately store the eigenvalues
      std::array<Number,dim> w;
      for (int i=0; i < dim; i++)
        w[i] = A[i][i];

      // Calculate (tr(A))^{2}
      sd = trace(A);
      sd *= sd;

      // Number of iterations
      const unsigned int max_n_it = 50;
      for (unsigned int it=0; it <= max_n_it; it++)
        {
          // Test for convergence
          so = 0.0;
          for (int p=0; p < dim; p++)
            for (int q=p+1; q < dim; q++)
              so += std::abs(A[p][q]);
          if (so == 0.0)
            break;

          // Throw if no convergence is achieved within a
          // stipulated number of iterations
          if (it == max_n_it)
            {
              AssertThrow(false, ExcMessage("No convergence in iterative Jacobi eigenvector algorithm."))
              return std::array<std::pair<Number, Tensor<1,dim,Number> >,dim> ();
            }

          // Compute threshold value which dictates whether or
          // not a Jacobi rotation is performed
          const unsigned int n_it_skip = 4;
          if (it < n_it_skip)
            thresh = 0.2 * so / (dim*dim);
          else
            thresh = 0.0;

          // Perform sweep
          for (int p=0; p < dim; p++)
            for (int q=p+1; q < dim; q++)
              {
                g = 100.0 * std::abs(A[p][q]);

                // After a given number of iterations the
                // rotation is skipped if the off-diagonal
                // element is small
                if (it > n_it_skip  &&
                    std::abs(w[p]) + g == std::abs(w[p])  &&
                    std::abs(w[q]) + g == std::abs(w[q]))
                  {
                    A[p][q] = 0.0;
                  }
                else if (std::abs(A[p][q]) > thresh)
                  {
                    // Calculate Jacobi transformation
                    h = w[q] - w[p];

                    // Compute surrogate for angle theta resulting from
                    // angle transformation and subsequent smallest solution
                    // of quadratic equation
                    if (std::abs(h) + g == std::abs(h))
                      {
                        // Prevent overflow for large theta^2. This computation
                        // is the algebraic equivalent of t = 1/(2*theta).
                        t = A[p][q] / h;
                      }
                    else
                      {
                        theta = 0.5 * h / A[p][q];
                        if (theta < 0.0)
                          t = -1.0 / (std::sqrt(1.0 + theta*theta) - theta);
                        else
                          t = 1.0 / (std::sqrt(1.0 + theta*theta) + theta);
                      }

                    // Compute trigonometric functions for rotation
                    // in such a way as to prevent overflow for
                    // large theta.
                    c = 1.0/std::sqrt(1.0 + t*t);
                    s = t * c;
                    z = t * A[p][q];

                    // Apply Jacobi transformation...
                    A[p][q] = 0.0;
                    w[p] -= z;
                    w[q] += z;
                    // ... by executing the various rotations in sequence
                    for (int r=0; r < p; r++)
                      {
                        t = A[r][p];
                        A[r][p] = c*t - s*A[r][q];
                        A[r][q] = s*t + c*A[r][q];
                      }
                    for (int r=p+1; r < q; r++)
                      {
                        t = A[p][r];
                        A[p][r] = c*t - s*A[r][q];
                        A[r][q] = s*t + c*A[r][q];
                      }
                    for (int r=q+1; r < dim; r++)
                      {
                        t = A[p][r];
                        A[p][r] = c*t - s*A[q][r];
                        A[q][r] = s*t + c*A[q][r];
                      }

                    // Update the eigenvectors
                    for (int r=0; r < dim; r++)
                      {
                        t = Q[r][p];
                        Q[r][p] = c*t - s*Q[r][q];
                        Q[r][q] = s*t + c*Q[r][q];
                      }
                  }
              }
        }

      // Structure the data to be outputted
      std::array<std::pair<Number, Tensor<1,dim,Number> >,dim> eig_vals_vecs;
      for (unsigned int e=0; e<dim; ++e)
        {
          eig_vals_vecs[e].first = w[e];

          // The column "e" of Q contains the non-normalized
          // eigenvector associated with the eigenvalue "e"
          for (unsigned int a=0; a<dim; ++a)
            {
              eig_vals_vecs[e].second[a] = Q[a][e];
            }

          // Normalize
          Assert(eig_vals_vecs[e].second.norm() != 0.0, ExcDivideByZero());
          eig_vals_vecs[e].second /= eig_vals_vecs[e].second.norm();
        }
      return eig_vals_vecs;
    }



    template <typename Number>
    std::array<std::pair<Number, Tensor<1,2,Number> >,2>
    hybrid (const dealii::SymmetricTensor<2,2,Number> &A)
    {
      static_assert(numbers::NumberTraits<Number>::is_complex == false,
                    "This implementation of the 2d Hybrid algorithm does "
                    "not support complex numbers");

      const unsigned int dim = 2;

      // Calculate eigenvalues
      const std::array<Number,dim> w = eigenvalues(A);

      std::array<std::pair<Number, Tensor<1,dim,Number> >,dim> eig_vals_vecs;

      Number t, u;          // Intermediate storage
      t = std::abs(w[0]);
      for (unsigned int i=1; i<dim; ++i)
        {
          u = std::abs(w[i]);
          if (u > t)
            t = u;
        }

      if (t < 1.0)
        u = t;
      else
        u = t*t;

      // Estimated maximum roundoff error
      const Number error = 256.0 * std::numeric_limits<double>::epsilon() * u*u;

      // Store eigenvalues
      eig_vals_vecs[0].first = w[0];
      eig_vals_vecs[1].first = w[1];

      // Compute eigenvectors
      // http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/
      // https://math.stackexchange.com/a/1548616
      if (A[1][0] != 0.0)
        {
          // First eigenvector
          eig_vals_vecs[0].second[0] = w[0] - A[1][1];
          eig_vals_vecs[0].second[1] = A[1][0];

          // Second eigenvector
          eig_vals_vecs[1].second[0] = w[1] - A[1][1];
          eig_vals_vecs[1].second[1] = A[1][0];
        }
      else
        {
          // First eigenvector
          eig_vals_vecs[0].second[0] = w[0];
          eig_vals_vecs[0].second[1] = 0.0;

          // Second eigenvector
          eig_vals_vecs[1].second[0] = 0.0;
          eig_vals_vecs[1].second[1] = w[1];
        }
      // Normalize
      eig_vals_vecs[0].second /= eig_vals_vecs[0].second.norm();
      eig_vals_vecs[1].second /= eig_vals_vecs[1].second.norm();

      // If vectors are nearly linearly dependent, or if there might have
      // been large cancelations in the calculation of A[i][i] - w[0], fall
      // back to QL algorithm
      if (eig_vals_vecs[0].second * eig_vals_vecs[1].second > error)
        {
          return ql_implicit_shifts(A);
        }

      return eig_vals_vecs;
    }



    template <typename Number>
    std::array<std::pair<Number, Tensor<1,3,Number> >,3>
    hybrid (const dealii::SymmetricTensor<2,3,Number> &A)
    {
      static_assert(numbers::NumberTraits<Number>::is_complex == false,
                    "This implementation of the 3d Hybrid algorithm does "
                    "not support complex numbers");

      const unsigned int dim = 3;
      Number norm;          // Squared norm or inverse norm of current eigenvector
      Number t, u;          // Intermediate storage

      // Calculate eigenvalues
      const std::array<Number,dim> w = eigenvalues(A);

      t = std::abs(w[0]);
      for (unsigned int i=1; i<dim; ++i)
        {
          u = std::abs(w[i]);
          if (u > t)
            t = u;
        }

      if (t < 1.0)
        u = t;
      else
        u = t*t;

      // Estimated maximum roundoff error
      const Number error = 256.0 * std::numeric_limits<double>::epsilon() * u*u;

      // Initialize the transformation matrix as the
      // identity tensor
      dealii::Tensor<2,dim,Number> Q;
      Q[0][1] = A[0][1]*A[1][2] - A[0][2]*A[1][1];
      Q[1][1] = A[0][2]*A[0][1] - A[1][2]*A[0][0];
      Q[2][1] = A[0][1]*A[0][1];

      // Calculate first eigenvector by the formula
      //   v[0] = (A - w[0]).e1 x (A - w[0]).e2
      Q[0][0] = Q[0][1] + A[0][2]*w[0];
      Q[1][0] = Q[1][1] + A[1][2]*w[0];
      Q[2][0] = (A[0][0] - w[0]) * (A[1][1] - w[0]) - Q[2][1];
      norm    = Q[0][0]*Q[0][0] + Q[1][0]*Q[1][0] + Q[2][0]*Q[2][0];

      // If vectors are nearly linearly dependent, or if there might have
      // been large cancellations in the calculation of A[i][i] - w[0], fall
      // back to QL algorithm
      // Note that this simultaneously ensures that multiple eigenvalues do
      // not cause problems: If w[0] = w[1], then A - w[0] * I has rank 1,
      // i.e. all columns of A - w[0] * I are linearly dependent.
      if (norm <= error)
        {
          return ql_implicit_shifts(A);
        }
      else                      // This is the standard branch
        {
          norm = std::sqrt(1.0 / norm);
          for (unsigned j=0; j < dim; j++)
            Q[j][0] = Q[j][0] * norm;
        }

      // Calculate second eigenvector by the formula
      //   v[1] = (A - w[1]).e1 x (A - w[1]).e2
      Q[0][1]  = Q[0][1] + A[0][2]*w[1];
      Q[1][1]  = Q[1][1] + A[1][2]*w[1];
      Q[2][1]  = (A[0][0] - w[1]) * (A[1][1] - w[1]) - Q[2][1];
      norm     = Q[0][1]*Q[0][1] + Q[1][1]*Q[1][1] + Q[2][1]*Q[2][1];
      if (norm <= error)
        {
          return ql_implicit_shifts(A);
        }
      else
        {
          norm = std::sqrt(1.0 / norm);
          for (unsigned int j=0; j < dim; j++)
            Q[j][1] = Q[j][1] * norm;
        }

      // Calculate third eigenvector according to
      //   v[2] = v[0] x v[1]
      Q[0][2] = Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1];
      Q[1][2] = Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1];
      Q[2][2] = Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1];

      // Structure the data to be outputted
      std::array<std::pair<Number, Tensor<1,dim,Number> >,dim> eig_vals_vecs;
      for (unsigned int e=0; e<dim; ++e)
        {
          eig_vals_vecs[e].first = w[e];

          // The column "e" of Q contains the non-normalized
          // eigenvector associated with the eigenvalue "e"
          for (unsigned int a=0; a<dim; ++a)
            {
              eig_vals_vecs[e].second[a] = Q[a][e];
            }

          // Normalize
          Assert(eig_vals_vecs[e].second.norm() != 0.0, ExcDivideByZero());
          eig_vals_vecs[e].second /= eig_vals_vecs[e].second.norm();
        }
      return eig_vals_vecs;
    }

  } // namespace SymmetricTensor
} // namespace internal



DEAL_II_NAMESPACE_CLOSE

#endif
