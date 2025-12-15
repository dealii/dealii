// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_fast_transcendental_h
#define dealii_fast_transcendental_h

#include <deal.II/base/config.h>

#include <deal.II/base/numbers.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/vectorization.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <type_traits>

#ifdef _MSC_VER
#  include <intrin.h>
#elif defined(__ARM_NEON)
#  include <arm_neon.h>
#elif defined(__x86_64__)
#  include <x86intrin.h>
#endif

DEAL_II_NAMESPACE_OPEN

namespace fast_transcendental
{
  namespace internal
  {
    /*
     * These arrays contain the polynomial coefficients used to approximate
     * the correction function K(y_f) that appears in the vectorized exponential
     * approximation. Each specialization provides the coefficients for a
     * specific polynomial degree.
     *
     * The coefficients have been determined via a least-squares fit using 1e8
     * uniformly spaced sample points over the interval [0, 1].
     *
     * The coefficients are ordered from the lowest to the highest polynomial
     * term, i.e., exp_pol_coeff[i] corresponds to the coefficient of
     * (y_f)^i.
     *
     * The primary template triggers a static assertion to ensure that an
     * unsupported polynomial degree results in a clear, compile-time error.
     */
    template <int degree, typename Number>
    constexpr std::array<Number, degree + 1> exp_pol_coeff =
      []() -> std::array<Number, degree + 1> {
      static_assert(degree > 2 && degree < 12,
                    "Unsupported polynomial degree for exp_pol_coeff.");
      return std::array<Number, degree + 1>{};
    }();

    template <typename Number>
    constexpr std::array<Number, 4> exp_pol_coeff<3, Number> = {
      {1.8809206555322363e-04,
       3.0316116404944521e-01,
       -2.2412622972531135e-01,
       -7.9019886953785187e-02}};

    template <typename Number>
    constexpr std::array<Number, 5> exp_pol_coeff<4, Number> = {
      {-7.2868375805796964e-06,
       3.0706874220004554e-01,
       -2.4171033150558635e-01,
       -5.1666839694436362e-02,
       -1.3676523629674565e-02}};

    template <typename Number>
    constexpr std::array<Number, 6> exp_pol_coeff<5, Number> = {
      {2.3053723876408911e-07,
       3.0684322094757277e-01,
       -2.4013168272249241e-01,
       -5.5876569798473615e-02,
       -8.9405772567122607e-03,
       -1.8943785491846432e-03}};

    template <typename Number>
    constexpr std::array<Number, 7> exp_pol_coeff<6, Number> = {
      {-6.1646222105268035e-09,
       3.0685316242622318e-01,
       -2.4023109751048086e-01,
       -5.5478910644035311e-02,
       -9.6861881733345725e-03,
       -1.2382409419015978e-03,
       -2.1871253576104385e-04}};

    template <typename Number>
    constexpr std::array<Number, 8> exp_pol_coeff<7, Number> = {
      {1.4276280848886066e-10,
       3.0685280921264796e-01,
       -2.4022632912712596e-01,
       -5.5505401662914296e-02,
       -9.6133378710895074e-03,
       -1.3431453773612783e-03,
       -1.4294822119902167e-04,
       -2.1646947017882258e-05}};

    template <typename Number>
    constexpr std::array<Number, 9> exp_pol_coeff<8, Number> = {
      {-2.9146986289710858e-12,
       3.0685281970141953e-01,
       -2.4022651268062989e-01,
       -5.5504055603868258e-02,
       -9.6183855925465120e-03,
       -1.3326461166925901e-03,
       -1.5519735866800178e-04,
       -1.4147475092332136e-05,
       -1.8748679813608566e-06}};

    template <typename Number>
    constexpr std::array<Number, 10> exp_pol_coeff<9, Number> = {
      {5.4337475594757437e-14,
       3.0685281943421089e-01,
       -2.4022650680203564e-01,
       -5.5504110470759059e-02,
       -9.6181181164443925e-03,
       -1.3333950497836748e-03,
       -1.5394913684739714e-04,
       -1.5370223000068715e-05,
       -1.2252831550509094e-06,
       -1.4435218373241001e-07}};

    template <typename Number>
    constexpr std::array<Number, 11> exp_pol_coeff<10, Number> = {
      {7.6072080297756739e-12,
       3.0685281862960340e-01,
       -2.4022648563631835e-01,
       -5.5504349821691862e-02,
       -9.6166785909533817e-03,
       -1.3384970969245304e-03,
       -1.4276337493295606e-04,
       -3.0711077326624553e-05,
       1.1584454591775664e-05,
       -6.0984990187293866e-06,
       1.1810109241254786e-06}};

    template <typename Number>
    constexpr std::array<Number, 12> exp_pol_coeff<11, Number> = {
      {3.6083394569535254e-13,
       3.0685281939395403e-01,
       -2.4022650549710606e-01,
       -5.5504128716495800e-02,
       -9.6179813320517119e-03,
       -1.3340079120703319e-03,
       -1.5221190054222076e-04,
       -1.8563222423500779e-05,
       2.5694826992914108e-06,
       -2.9577981133179920e-06,
       1.1827521620524028e-06,
       -2.1525063398151332e-07}};

    /**
     * This is 2^52 for double and 2^23 for float.
     */
    template <typename Number>
    constexpr Number exp_coeff_a = numbers::signaling_nan<Number>();
    template <>
    inline constexpr double exp_coeff_a<double> = 4503599627370496;
    template <>
    inline constexpr float exp_coeff_a<float> = 8388608;

    /**
     * This is 2^52 * 1023 for double and 2^23 * 127 for float.
     */
    template <typename Number>
    inline constexpr Number exp_coeff_b = numbers::signaling_nan<Number>();
    template <>
    inline constexpr double exp_coeff_b<double> = 4.60718241880001741e+18;
    template <>
    inline constexpr float exp_coeff_b<float> = 1065353216;

    /**
     * This is the maximum exponent for a floating point representation shifted
     * to the exponent bits. In numbers this is 2^52 * 2047 for double and 2^23
     * * 255 for float.
     */
    template <typename Number>
    inline constexpr Number max_exponent = numbers::signaling_nan<Number>();
    template <>
    inline constexpr double max_exponent<double> = 9.218868437227405312e+18;
    template <>
    inline constexpr float max_exponent<float> = 2139095040;

    /**
     * The largest |x| for which exp(x) remains a finite, representable value of
     * type Number. For double this is approximately 308 / log10(e), and for
     * float approximately 38 / log10(e).
     */
    template <typename Number>
    inline constexpr Number exp_max_abs_x = numbers::signaling_nan<Number>();
    template <>
    inline constexpr double exp_max_abs_x<double> = 709.196208642;
    template <>
    inline constexpr float exp_max_abs_x<float> = 87.49823353f;

    /**
     * The functions below provide architecture-specific specializations used
     * as a helper by the `fast_exp()` implementation. Each specialization
     * performs a two-step operation:
     *
     * <ol>
     *
     * <li> Static cast of each floating-point value to its corresponding integer
     * type (`float` to `int32_t`, `double` to int64_t`).
     *
     * <li> Reinterpretation of the resulting integer vector as a floating-point
     * vector via bitwise casting.
     *
     * </ol>
     */
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512 && defined(__AVX512F__)
    DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 8>
    type_cast(const VectorizedArray<double, 8> &in)
    {
      VectorizedArray<double, 8> out{};
      const __m512i              integer = _mm512_cvtpd_epi64(in.data);
      out.data                           = _mm512_castsi512_pd(integer);
      return out;
    }

    DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 16>
    type_cast(const VectorizedArray<float, 16> &in)
    {
      VectorizedArray<float, 16> out{};
      const __m512i              integer = _mm512_cvtps_epi32(in.data);
      out.data                           = _mm512_castsi512_ps(integer);
      return out;
    }
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && defined(__AVX__)
    DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 4>
    type_cast(const VectorizedArray<double, 4> &in)
    {
      VectorizedArray<double, 4> out;
      int64_t                    int_values[4];
      for (int i = 0; i < 4; i++)
        {
          int_values[i] = static_cast<int64_t>(in[i]);
        }
      out.data = _mm256_castsi256_pd(
        _mm256_loadu_si256(reinterpret_cast<__m256i *>(int_values)));
      return out;
    }

    DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 8>
    type_cast(const VectorizedArray<float, 8> &in)
    {
      VectorizedArray<float, 8> out;
      int32_t                   int_values[8];
      for (int i = 0; i < 8; i++)
        {
          int_values[i] = static_cast<int32_t>(in[i]);
        }
      out.data = _mm256_castsi256_ps(
        _mm256_loadu_si256(reinterpret_cast<__m256i *>(int_values)));
      return out;
    }
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__SSE2__)
    DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 2>
    type_cast(const VectorizedArray<double, 2> &in)
    {
      VectorizedArray<double, 2> out;
      int64_t                    int_values[2];
      for (int i = 0; i < 2; i++)
        {
          int_values[i] = static_cast<int64_t>(in[i]);
        }
      out.data = _mm_castsi128_pd(
        _mm_loadu_si128(reinterpret_cast<__m128i *>(int_values)));
      return out;
    }

    DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 4>
    type_cast(const VectorizedArray<float, 4> &in)
    {
      VectorizedArray<float, 4> out;
      int32_t                   int_values[4];
      for (int i = 0; i < 4; i++)
        {
          int_values[i] = static_cast<int32_t>(in[i]);
        }
      out.data = _mm_castsi128_ps(
        _mm_loadu_si128(reinterpret_cast<__m128i *>(int_values)));
      return out;
    }
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__ARM_NEON)
    DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 2>
    type_cast(const VectorizedArray<double, 2> &in)
    {
      VectorizedArray<double, 2> out;
      const int64x2_t            int_values = vcvtq_s64_f64(in.data);
      out.data                              = vreinterpretq_f64_s64(int_values);
      return out;
    }

    DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 4>
    type_cast(const VectorizedArray<float, 4> &in)
    {
      VectorizedArray<float, 4> out;
      const int32x4_t           int_values = vcvtq_s32_f32(in.data);
      out.data                             = vreinterpretq_f32_s32(int_values);
      return out;
    }
#endif

    template <typename number>
    DEAL_II_ALWAYS_INLINE inline VectorizedArray<number, 1>
    type_cast(const VectorizedArray<number, 1> &in)
    {
      VectorizedArray<number, 1> out;
      const auto                 result = static_cast<
        std::conditional_t<std::is_same_v<number, double>, int64_t, int32_t>>(
        in.data);
      std::memcpy(&out.data, &result, sizeof(out.data));
      return out;
    }

    /**
     * This function computes the value of a polynomial
     * @f[
     *   P(x) = c_0 + c_1 x + c_2 x^2 + \dots + c_{\text{degree}}
     * x^{\text{degree}}
     * @f]
     * where the coefficients $c_i$ are provided in ascending order using
     * Horner's scheme.
     *
     * @tparam degree      Degree of the polynomial.
     * @tparam CoeffsType  Type of the stored coefficients.
     * @tparam ValueType   Type of the evaluation point and result.
     *
     * @param c  Array of polynomial coefficients in ascending order.
     * @param x  The point at which the polynomial should be evaluated.
     *
     * @return The value $P(x)$.
     */
    template <int degree, typename CoeffsType, typename ValueType>
    DEAL_II_ALWAYS_INLINE inline ValueType
    horner_scheme(const std::array<CoeffsType, degree + 1> &c,
                  const ValueType                          &x)
    {
      ValueType result = c.back();
      for (int i = c.size() - 2; i >= 0; --i)
        result = c[i] + x * result;

      return result;
    }

    /**
     * This function computes the value of a polynomial
     * @f[
     *   P(x) = c_0 + c_1 x + c_2 x^2 + \dots + c_{\text{degree}}
     * x^{\text{degree}}
     * @f]
     * where the coefficients $c_i$ are provided in ascending order using
     * Estrin's scheme.
     *
     * @tparam degree      Degree of the polynomial.
     * @tparam CoeffsType  Type of the stored coefficients.
     * @tparam ValueType   Type of the evaluation point and result.
     *
     * @param c  Array of polynomial coefficients in ascending order.
     * @param x  The point at which the polynomial should be evaluated.
     *
     * @return The value $P(x)$.
     */
    template <int degree, typename CoeffsType, typename ValueType>
    DEAL_II_ALWAYS_INLINE inline ValueType
    estrin_scheme(const std::array<CoeffsType, degree + 1> &c,
                  const ValueType                          &x)
    {
      constexpr int n_precomputed_x_powers = degree / 2;
      std::array<ValueType, n_precomputed_x_powers> x_powers;
      x_powers[0] = x * x;
      for (unsigned i = 1; i < n_precomputed_x_powers; ++i)
        x_powers[i] = x_powers[i - 1] * x_powers[0];

      ValueType t1, t2;
      t1 = (c[0] + c[1] * x);

      for (unsigned i = 2; i < degree; i += 2)
        {
          t2 = (c[i] + c[i + 1] * x);
          t1 += t2 * x_powers[(i - 2) / 2];
        }

      if constexpr (degree % 2 == 0)
        t1 += c.back() * x_powers.back();
      return t1;
    }
  } // namespace internal

  /**
   * Enumeration of available polynomial evaluation schemes for the correction
   * polynomial in fast_exp_approximation().
   */
  enum class PolynomialEvalScheme : int
  {
    estrin,
    horner
  };


  /**
   * This function provides a fully vectorizable algorithm for approximating
   * the exponential function. In contrast to the standard `std::exp()`
   * implementation, this algorithm is designed for SIMD execution and
   * therefore offers significant performance advantages in vectorized
   * computations while maintaining high accuracy. The main algorithm is
   * briefly summarized here. For more details see @cite Proell2024.
   *
   * The algorithm is based on the IEEE-754 floating-point representation
   * @f[
   *   (-1)^s 2^{p-b} (1+m),
   * @f]
   * where $s$ is the sign bit, $p$ is the exponent and $m$ is the
   * mantissa with $0 \le m < 1$. $b$ is the bias, a constant
   * depending only on the specific type.
   *
   * For the exponential function one can write
   * @f[
   *   \exp(x)
   *     = 2^{x \log_2(e)}
   *     = 2^y
   *     = 2^{y_i}2^{y_f},
   * @f]
   * where $y = x \log_2(e)$, $y_i = \lfloor y \rfloor$ is the integer
   * part, and $y_f = y - y_i$ is the fractional part. From the IEEE-754
   * expression one can identify
   * @f[
   *   s = 0,\qquad 2^{y_i} = 2^{p-b},\qquad 2^{y_f} = 1 + m.
   * @f]
   * The integer exponent $y_i$ can be computed directly. However,
   * although
   * $y_f$ is known, reconstructing the mantissa $m = 2^{y_f} - 1$
   * still requires evaluating $2^{y_f}$, which is transcendental and
   * therefore expensive.
   *
   * To avoid transcendental operations, a correction function
   * @f[
   *   K(y_f) = 1 + y_f - 2^{y_f},
   * @f]
   * which is approximated using a polynomial of configurable degree, is
   * introduced. With this approximation, the mantissa can be written as
   * @f[
   *   m = y_f - K(y_f).
   * @f]
   * The accuracy of the exponential approximation depends primarily on the
   * polynomial used for $K(y_f)$. The degree of the polynomial can be
   * chosen using the template parameter `degree`. Supported degrees are in
   * the range [3, 11]. Higher degrees typically improve accuracy at the
   * expense of additional compute cost. In many applications, a degree of 5
   * already yields an excellent balance between precision and performance,
   * but users may tune it depending on their specific requirements.
   *
   * Once the mantissa and exponent have been computed, the final
   * floating-point value is assembled through internal bit manipulation that
   * constructs the correct IEEE-754 representation. For further algorithmic
   * details, see the referenced literature.
   *
   * This function optionally performs lower- and upper-bound checks on the
   * input values. If enabled:
   * - Underflow check: Values too small to represent are set to `0`.
   * - Overflow check: Values too large to represent are set to `+∞`.
   * These checks ensure numerical safety but incur additional runtime
   * overhead. They may be disabled when the user can guarantee that inputs
   * lie within the representable range.
   *
   * Two schemes are available for evaluating the correction polynomial:
   * (Horner's method)[https://en.wikipedia.org/wiki/Horner%27s_method], which
   * minimizes the operation count, and (Estrin's
   * method)[https://en.wikipedia.org/wiki/Estrin%27s_scheme], which performs
   * slightly more operations but makes better use of modern hardware
   * parallelism. The preferable scheme can vary depending on the application
   * and system architecture and should be assessed separately for each use
   * case.
   *
   * @tparam degree  Degree of the polynomial used to approximate the correction function.
   * @tparam Number  The floating-point data type. Only `double` and `float` are supported.
   * @tparam width   The SIMD vector width.
   * @tparam ensure_correct_treatment_of_zero  Whether to perform an underflow check.
   * @tparam ensure_correct_treatment_of_infinity  Whether to perform an overflow check.
   * @tparam polynomial_evaluation_scheme The polynomial evaluation scheme to be used for
   * the evaluation of the correction function.
   *
   * @param x The vector of values at which the exponential function is approximated.
   *
   * @return Approximation of $\exp(x)$.
   */
  template <int degree,
            typename Number,
            bool                 ensure_correct_treatment_of_zero     = true,
            bool                 ensure_correct_treatment_of_infinity = true,
            PolynomialEvalScheme polynomial_evaluation_scheme =
              PolynomialEvalScheme::horner>
  DEAL_II_ALWAYS_INLINE inline Number
  exp(Number x)
  {
    if constexpr (std::is_floating_point_v<Number>)
      {
        auto r = exp<degree,
                     dealii::VectorizedArray<Number, 1>,
                     ensure_correct_treatment_of_zero,
                     ensure_correct_treatment_of_infinity,
                     polynomial_evaluation_scheme>(
          dealii::VectorizedArray<Number, 1>(x));
        return r.data;
      }
    else
      {
        using floating_type = typename Number::value_type;

        static_assert(
          std::is_same_v<floating_type, double> ||
            std::is_same_v<floating_type, float>,
          "The fast transcendental approximation of the exp() function only supports the fundamental types float and double.");

        Number r = x * static_cast<floating_type>(numbers::LOG2E);

        Number fractional_exponent = r - r.get_floor();

        if constexpr (polynomial_evaluation_scheme ==
                      PolynomialEvalScheme::horner)
          r -= internal::horner_scheme<degree>(
            internal::exp_pol_coeff<degree, floating_type>,
            fractional_exponent);
        else if constexpr (polynomial_evaluation_scheme ==
                           PolynomialEvalScheme::estrin)
          r -= internal::estrin_scheme<degree>(
            internal::exp_pol_coeff<degree, floating_type>,
            fractional_exponent);
        else
          static_assert(
            polynomial_evaluation_scheme == PolynomialEvalScheme::estrin ||
              polynomial_evaluation_scheme == PolynomialEvalScheme::horner,
            "The provided scheme for the evaluation of the correction function is not supported!");

        r = internal::exp_coeff_a<floating_type> * r +
            internal::exp_coeff_b<floating_type>;

        // To ensure correct handling of zero and infinity, the input to the
        // static_cast operation is adjusted so that values falling outside the
        // valid range produce exactly zero or infinity.
        if constexpr (ensure_correct_treatment_of_zero)
          {
            // Denormal values are excluded since the approximation method does
            // not work for them.
            r = compare_and_apply_mask<SIMDComparison::less_than>(
              x,
              Number(-internal::exp_max_abs_x<floating_type>),
              Number(0.),
              r);
          }
        if constexpr (ensure_correct_treatment_of_infinity)
          {
            r = compare_and_apply_mask<SIMDComparison::greater_than>(
              x,
              Number(internal::exp_max_abs_x<floating_type>),
              Number(internal::max_exponent<floating_type>),
              r);
          }

        return internal::type_cast(r);
      }
  }
} // namespace fast_transcendental

DEAL_II_NAMESPACE_CLOSE

#endif
