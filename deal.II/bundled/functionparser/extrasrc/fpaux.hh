/***************************************************************************\
|* Function Parser for C++ v4.5.1                                          *|
|*-------------------------------------------------------------------------*|
|* Copyright: Juha Nieminen, Joel Yliluoma                                 *|
|*                                                                         *|
|* This library is distributed under the terms of the                      *|
|* GNU Lesser General Public License version 3.                            *|
|* (See lgpl.txt and gpl.txt for the license text.)                        *|
\***************************************************************************/

// NOTE:
// This file contains only internal types for the function parser library.
// You don't need to include this file in your code. Include "fparser.hh"
// only.

#ifndef ONCE_FPARSER_AUX_H_
#define ONCE_FPARSER_AUX_H_

#include "fptypes.hh"

#include <cmath>

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
#include "mpfr/MpfrFloat.hh"
#endif

#ifdef FP_SUPPORT_GMP_INT_TYPE
#include "mpfr/GmpInt.hh"
#endif

#ifdef FP_SUPPORT_COMPLEX_NUMBERS
#include <complex>
#endif

#ifdef ONCE_FPARSER_H_
namespace FUNCTIONPARSERTYPES
{
    template<typename>
    struct IsIntType
    {
        enum { result = false };
    };
    template<>
    struct IsIntType<long>
    {
        enum { result = true };
    };
#ifdef FP_SUPPORT_GMP_INT_TYPE
    template<>
    struct IsIntType<GmpInt>
    {
        enum { result = true };
    };
#endif

    template<typename>
    struct IsComplexType
    {
        enum { result = false };
    };
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
    template<typename T>
    struct IsComplexType<std::complex<T> >
    {
        enum { result = true };
    };
#endif


//==========================================================================
// Constants
//==========================================================================
    template<typename Value_t>
    inline Value_t fp_const_pi() // CONSTANT_PI
    {
        return Value_t(3.1415926535897932384626433832795028841971693993751L);
    }

    template<typename Value_t>
    inline Value_t fp_const_e() // CONSTANT_E
    {
        return Value_t(2.7182818284590452353602874713526624977572L);
    }
    template<typename Value_t>
    inline Value_t fp_const_einv() // CONSTANT_EI
    {
        return Value_t(0.367879441171442321595523770161460867445811131L);
    }
    template<typename Value_t>
    inline Value_t fp_const_log2() // CONSTANT_L2, CONSTANT_L2EI
    {
        return Value_t(0.69314718055994530941723212145817656807550013436025525412L);
    }
    template<typename Value_t>
    inline Value_t fp_const_log10() // CONSTANT_L10, CONSTANT_L10EI
    {
        return Value_t(2.302585092994045684017991454684364207601101488628772976L);
    }
    template<typename Value_t>
    inline Value_t fp_const_log2inv() // CONSTANT_L2I, CONSTANT_L2E
    {
        return Value_t(1.442695040888963407359924681001892137426645954L);
    }
    template<typename Value_t>
    inline Value_t fp_const_log10inv() // CONSTANT_L10I, CONSTANT_L10E
    {
        return Value_t(0.434294481903251827651128918916605082294397L);
    }

    template<typename Value_t>
    inline const Value_t& fp_const_deg_to_rad() // CONSTANT_DR
    {
        static const Value_t factor = fp_const_pi<Value_t>() / Value_t(180); // to rad from deg
        return factor;
    }

    template<typename Value_t>
    inline const Value_t& fp_const_rad_to_deg() // CONSTANT_RD
    {
        static const Value_t factor = Value_t(180) / fp_const_pi<Value_t>(); // to deg from rad
        return factor;
    }

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
    template<>
    inline MpfrFloat fp_const_pi<MpfrFloat>() { return MpfrFloat::const_pi(); }

    template<>
    inline MpfrFloat fp_const_e<MpfrFloat>() { return MpfrFloat::const_e(); }

    template<>
    inline MpfrFloat fp_const_einv<MpfrFloat>() { return MpfrFloat(1) / MpfrFloat::const_e(); }

    template<>
    inline MpfrFloat fp_const_log2<MpfrFloat>() { return MpfrFloat::const_log2(); }

    /*
    template<>
    inline MpfrFloat fp_const_log10<MpfrFloat>() { return fp_log(MpfrFloat(10)); }

    template<>
    inline MpfrFloat fp_const_log2inv<MpfrFloat>() { return MpfrFloat(1) / MpfrFloat::const_log2(); }

    template<>
    inline MpfrFloat fp_const_log10inv<MpfrFloat>() { return fp_log10(MpfrFloat::const_e()); }
    */
#endif


//==========================================================================
// Generic math functions
//==========================================================================
    template<typename Value_t>
    inline Value_t fp_abs(const Value_t& x) { return std::fabs(x); }

    template<typename Value_t>
    inline Value_t fp_acos(const Value_t& x) { return std::acos(x); }

    template<typename Value_t>
    inline Value_t fp_asin(const Value_t& x) { return std::asin(x); }

    template<typename Value_t>
    inline Value_t fp_atan(const Value_t& x) { return std::atan(x); }

    template<typename Value_t>
    inline Value_t fp_atan2(const Value_t& x, const Value_t& y)
    { return std::atan2(x, y); }

    template<typename Value_t>
    inline Value_t fp_ceil(const Value_t& x) { return std::ceil(x); }

    template<typename Value_t>
    inline Value_t fp_cos(const Value_t& x) { return std::cos(x); }

    template<typename Value_t>
    inline Value_t fp_cosh(const Value_t& x) { return std::cosh(x); }

    template<typename Value_t>
    inline Value_t fp_exp(const Value_t& x) { return std::exp(x); }

    template<typename Value_t>
    inline Value_t fp_floor(const Value_t& x) { return std::floor(x); }

    template<typename Value_t>
    inline Value_t fp_log(const Value_t& x) { return std::log(x); }

    template<typename Value_t>
    inline Value_t fp_mod(const Value_t& x, const Value_t& y)
    { return std::fmod(x, y); }

    template<typename Value_t>
    inline Value_t fp_sin(const Value_t& x) { return std::sin(x); }

    template<typename Value_t>
    inline Value_t fp_sinh(const Value_t& x) { return std::sinh(x); }

    template<typename Value_t>
    inline Value_t fp_sqrt(const Value_t& x) { return std::sqrt(x); }

    template<typename Value_t>
    inline Value_t fp_tan(const Value_t& x) { return std::tan(x); }

    template<typename Value_t>
    inline Value_t fp_tanh(const Value_t& x) { return std::tanh(x); }

#ifdef FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS
    template<typename Value_t>
    inline Value_t fp_asinh(const Value_t& x) { return std::asinh(x); }

    template<typename Value_t>
    inline Value_t fp_acosh(const Value_t& x) { return std::acosh(x); }

    template<typename Value_t>
    inline Value_t fp_atanh(const Value_t& x) { return std::atanh(x); }
#else
    template<typename Value_t>
    inline Value_t fp_asinh(const Value_t& x)
    { return fp_log(x + fp_sqrt(x*x + Value_t(1))); }

    template<typename Value_t>
    inline Value_t fp_acosh(const Value_t& x)
    { return fp_log(x + fp_sqrt(x*x - Value_t(1))); }

    template<typename Value_t>
    inline Value_t fp_atanh(const Value_t& x)
    {
        return fp_log( (Value_t(1)+x) / (Value_t(1)-x)) * Value_t(0.5);
        // Note: x = +1 causes division by zero
        //       x = -1 causes log(0)
        // Thus, x must not be +-1
    }
#endif // FP_SUPPORT_ASINH

#ifdef FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS
    template<typename Value_t>
    inline Value_t fp_hypot(const Value_t& x, const Value_t& y)
    { return std::hypot(x,y); }
#else
    template<typename Value_t>
    inline Value_t fp_hypot(const Value_t& x, const Value_t& y)
    { return fp_sqrt(x*x + y*y); }
#endif

    template<typename Value_t>
    inline Value_t fp_pow_base(const Value_t& x, const Value_t& y)
    { return std::pow(x, y); }

#ifdef FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS
    template<typename Value_t>
    inline Value_t fp_log2(const Value_t& x) { return std::log2(x); }
#else
    template<typename Value_t>
    inline Value_t fp_log2(const Value_t& x)
    {
        return fp_log(x) * fp_const_log2inv<Value_t>();
    }
#endif // FP_SUPPORT_LOG2

    template<typename Value_t>
    inline Value_t fp_log10(const Value_t& x)
    {
        return fp_log(x) * fp_const_log10inv<Value_t>();
    }

    template<typename Value_t>
    inline Value_t fp_trunc(const Value_t& x)
    {
        return x < Value_t() ? fp_ceil(x) : fp_floor(x);
    }

    template<typename Value_t>
    inline Value_t fp_int(const Value_t& x)
    {
        return x < Value_t() ?
            fp_ceil(x - Value_t(0.5)) : fp_floor(x + Value_t(0.5));
    }

    template<typename Value_t>
    inline void fp_sinCos(Value_t& sinvalue, Value_t& cosvalue,
                          const Value_t& param)
    {
        // Assuming that "cosvalue" and "param" do not
        // overlap, but "sinvalue" and "param" may.
        cosvalue = fp_cos(param);
        sinvalue = fp_sin(param);
    }

    template<typename Value_t>
    inline void fp_sinhCosh(Value_t& sinhvalue, Value_t& coshvalue,
                            const Value_t& param)
    {
        const Value_t ex(fp_exp(param)), emx(fp_exp(-param));
        sinhvalue = Value_t(0.5)*(ex-emx);
        coshvalue = Value_t(0.5)*(ex+emx);
    }

    template<typename Value_t>
    struct Epsilon
    {
        static Value_t value;
        static Value_t defaultValue() { return 0; }
    };

    template<> inline double Epsilon<double>::defaultValue() { return 1E-12; }
    template<> inline float Epsilon<float>::defaultValue() { return 1E-5F; }
    template<> inline long double Epsilon<long double>::defaultValue() { return 1E-14L; }

    template<> inline std::complex<double>
    Epsilon<std::complex<double> >::defaultValue() { return 1E-12; }

    template<> inline std::complex<float>
    Epsilon<std::complex<float> >::defaultValue() { return 1E-5F; }

    template<> inline std::complex<long double>
    Epsilon<std::complex<long double> >::defaultValue() { return 1E-14L; }

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
    template<> inline MpfrFloat
    Epsilon<MpfrFloat>::defaultValue() { return MpfrFloat::someEpsilon(); }
#endif

    template<typename Value_t> Value_t Epsilon<Value_t>::value =
        Epsilon<Value_t>::defaultValue();


#ifdef _GNU_SOURCE
    inline void fp_sinCos(double& sin, double& cos, const double& a)
    {
        sincos(a, &sin, &cos);
    }
    inline void fp_sinCos(float& sin, float& cos, const float& a)
    {
        sincosf(a, &sin, &cos);
    }
    inline void fp_sinCos(long double& sin, long double& cos,
                          const long double& a)
    {
        sincosl(a, &sin, &cos);
    }
#endif


// -------------------------------------------------------------------------
// Long int
// -------------------------------------------------------------------------
    inline long fp_abs(const long& x) { return x < 0 ? -x : x; }
    inline long fp_acos(const long&) { return 0; }
    inline long fp_asin(const long&) { return 0; }
    inline long fp_atan(const long&) { return 0; }
    inline long fp_atan2(const long&, const long&) { return 0; }
    inline long fp_cbrt(const long&) { return 0; }
    inline long fp_ceil(const long& x) { return x; }
    inline long fp_cos(const long&) { return 0; }
    inline long fp_cosh(const long&) { return 0; }
    inline long fp_exp(const long&) { return 0; }
    inline long fp_exp2(const long&) { return 0; }
    inline long fp_floor(const long& x) { return x; }
    inline long fp_log(const long&) { return 0; }
    inline long fp_log2(const long&) { return 0; }
    inline long fp_log10(const long&) { return 0; }
    inline long fp_mod(const long& x, const long& y) { return x % y; }
    inline long fp_pow(const long&, const long&) { return 0; }
    inline long fp_sin(const long&) { return 0; }
    inline long fp_sinh(const long&) { return 0; }
    inline long fp_sqrt(const long&) { return 1; }
    inline long fp_tan(const long&) { return 0; }
    inline long fp_tanh(const long&) { return 0; }
    inline long fp_asinh(const long&) { return 0; }
    inline long fp_acosh(const long&) { return 0; }
    inline long fp_atanh(const long&) { return 0; }
    inline long fp_pow_base(const long&, const long&) { return 0; }
    inline void fp_sinCos(long&, long&, const long&) {}
    inline void fp_sinhCosh(long&, long&, const long&) {}

    //template<> inline long fp_epsilon<long>() { return 0; }


// -------------------------------------------------------------------------
// MpfrFloat
// -------------------------------------------------------------------------
#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
    inline MpfrFloat fp_abs(const MpfrFloat& x) { return MpfrFloat::abs(x); }
    inline MpfrFloat fp_acos(const MpfrFloat& x) { return MpfrFloat::acos(x); }
    inline MpfrFloat fp_acosh(const MpfrFloat& x) { return MpfrFloat::acosh(x); }
    inline MpfrFloat fp_asin(const MpfrFloat& x) { return MpfrFloat::asin(x); }
    inline MpfrFloat fp_asinh(const MpfrFloat& x) { return MpfrFloat::asinh(x); }
    inline MpfrFloat fp_atan(const MpfrFloat& x) { return MpfrFloat::atan(x); }
    inline MpfrFloat fp_atan2(const MpfrFloat& x, const MpfrFloat& y)
    { return MpfrFloat::atan2(x, y); }
    inline MpfrFloat fp_atanh(const MpfrFloat& x) { return MpfrFloat::atanh(x); }
    inline MpfrFloat fp_cbrt(const MpfrFloat& x) { return MpfrFloat::cbrt(x); }
    inline MpfrFloat fp_ceil(const MpfrFloat& x) { return MpfrFloat::ceil(x); }
    inline MpfrFloat fp_cos(const MpfrFloat& x) { return MpfrFloat::cos(x); }
    inline MpfrFloat fp_cosh(const MpfrFloat& x) { return MpfrFloat::cosh(x); }
    inline MpfrFloat fp_exp(const MpfrFloat& x) { return MpfrFloat::exp(x); }
    inline MpfrFloat fp_exp2(const MpfrFloat& x) { return MpfrFloat::exp2(x); }
    inline MpfrFloat fp_floor(const MpfrFloat& x) { return MpfrFloat::floor(x); }
    inline MpfrFloat fp_hypot(const MpfrFloat& x, const MpfrFloat& y)
    { return MpfrFloat::hypot(x, y); }
    inline MpfrFloat fp_int(const MpfrFloat& x) { return MpfrFloat::round(x); }
    inline MpfrFloat fp_log(const MpfrFloat& x) { return MpfrFloat::log(x); }
    inline MpfrFloat fp_log2(const MpfrFloat& x) { return MpfrFloat::log2(x); }
    inline MpfrFloat fp_log10(const MpfrFloat& x) { return MpfrFloat::log10(x); }
    inline MpfrFloat fp_mod(const MpfrFloat& x, const MpfrFloat& y) { return x % y; }
    inline MpfrFloat fp_sin(const MpfrFloat& x) { return MpfrFloat::sin(x); }
    inline MpfrFloat fp_sinh(const MpfrFloat& x) { return MpfrFloat::sinh(x); }
    inline MpfrFloat fp_sqrt(const MpfrFloat& x) { return MpfrFloat::sqrt(x); }
    inline MpfrFloat fp_tan(const MpfrFloat& x) { return MpfrFloat::tan(x); }
    inline MpfrFloat fp_tanh(const MpfrFloat& x) { return MpfrFloat::tanh(x); }
    inline MpfrFloat fp_trunc(const MpfrFloat& x) { return MpfrFloat::trunc(x); }

    inline MpfrFloat fp_pow(const MpfrFloat& x, const MpfrFloat& y) { return MpfrFloat::pow(x, y); }
    inline MpfrFloat fp_pow_base(const MpfrFloat& x, const MpfrFloat& y) { return MpfrFloat::pow(x, y); }


    inline void fp_sinCos(MpfrFloat& sin, MpfrFloat& cos, const MpfrFloat& a)
    {
        MpfrFloat::sincos(a, sin, cos);
    }

    inline void fp_sinhCosh(MpfrFloat& sinhvalue, MpfrFloat& coshvalue,
                            const MpfrFloat& param)
    {
        const MpfrFloat paramCopy = param;
        sinhvalue = fp_sinh(paramCopy);
        coshvalue = fp_cosh(paramCopy);
    }
#endif // FP_SUPPORT_MPFR_FLOAT_TYPE


// -------------------------------------------------------------------------
// GMP int
// -------------------------------------------------------------------------
#ifdef FP_SUPPORT_GMP_INT_TYPE
    inline GmpInt fp_abs(const GmpInt& x) { return GmpInt::abs(x); }
    inline GmpInt fp_acos(const GmpInt&) { return 0; }
    inline GmpInt fp_acosh(const GmpInt&) { return 0; }
    inline GmpInt fp_asin(const GmpInt&) { return 0; }
    inline GmpInt fp_asinh(const GmpInt&) { return 0; }
    inline GmpInt fp_atan(const GmpInt&) { return 0; }
    inline GmpInt fp_atan2(const GmpInt&, const GmpInt&) { return 0; }
    inline GmpInt fp_atanh(const GmpInt&) { return 0; }
    inline GmpInt fp_cbrt(const GmpInt&) { return 0; }
    inline GmpInt fp_ceil(const GmpInt& x) { return x; }
    inline GmpInt fp_cos(const GmpInt&) { return 0; }
    inline GmpInt fp_cosh(const GmpInt&) { return 0; }
    inline GmpInt fp_exp(const GmpInt&) { return 0; }
    inline GmpInt fp_exp2(const GmpInt&) { return 0; }
    inline GmpInt fp_floor(const GmpInt& x) { return x; }
    inline GmpInt fp_hypot(const GmpInt&, const GmpInt&) { return 0; }
    inline GmpInt fp_int(const GmpInt& x) { return x; }
    inline GmpInt fp_log(const GmpInt&) { return 0; }
    inline GmpInt fp_log2(const GmpInt&) { return 0; }
    inline GmpInt fp_log10(const GmpInt&) { return 0; }
    inline GmpInt fp_mod(const GmpInt& x, const GmpInt& y) { return x % y; }
    inline GmpInt fp_pow(const GmpInt&, const GmpInt&) { return 0; }
    inline GmpInt fp_sin(const GmpInt&) { return 0; }
    inline GmpInt fp_sinh(const GmpInt&) { return 0; }
    inline GmpInt fp_sqrt(const GmpInt&) { return 0; }
    inline GmpInt fp_tan(const GmpInt&) { return 0; }
    inline GmpInt fp_tanh(const GmpInt&) { return 0; }
    inline GmpInt fp_trunc(const GmpInt& x) { return x; }
    inline GmpInt fp_pow_base(const GmpInt&, const GmpInt&) { return 0; }
    inline void fp_sinCos(GmpInt&, GmpInt&, const GmpInt&) {}
    inline void fp_sinhCosh(GmpInt&, GmpInt&, const GmpInt&) {}
#endif // FP_SUPPORT_GMP_INT_TYPE


#ifdef FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS
    template<typename Value_t>
    inline Value_t fp_cbrt(const Value_t& x) { return std::cbrt(x); }
#else
    template<typename Value_t>
    inline Value_t fp_cbrt(const Value_t& x)
    {
        return (x > Value_t() ?  fp_exp(fp_log( x) / Value_t(3)) :
                x < Value_t() ? -fp_exp(fp_log(-x) / Value_t(3)) :
                Value_t());
    }
#endif

// -------------------------------------------------------------------------
// Synthetic functions and fallbacks for when an optimized
// implementation or a library function is not available
// -------------------------------------------------------------------------
    template<typename Value_t> inline Value_t fp_arg(const Value_t& x);
    template<typename Value_t> inline Value_t fp_exp2(const Value_t& x);
    template<typename Value_t> inline Value_t fp_int(const Value_t& x);
    template<typename Value_t> inline Value_t fp_trunc(const Value_t& x);
    template<typename Value_t>
    inline void fp_sinCos(Value_t& , Value_t& , const Value_t& );
    template<typename Value_t>
    inline void fp_sinhCosh(Value_t& , Value_t& , const Value_t& );

#ifdef FP_SUPPORT_COMPLEX_NUMBERS
    /* NOTE: Complex multiplication of a and b can be done with:
        tmp = b.real * (a.real + a.imag)
        result.real = tmp - a.imag * (b.real + b.imag)
        result.imag = tmp + a.real * (b.imag - b.real)
        This has fewer multiplications than the standard
        algorithm. Take note, if you support mpfr complex one day.
    */

    template<typename T>
    struct FP_ProbablyHasFastLibcComplex
    { enum { result = false }; };
    /* The generic sqrt() etc. implementations in libstdc++
     * are very plain and non-optimized; however, it contains
     * callbacks to libc complex math functions where possible,
     * and I suspect that those may actually be well optimized.
     * So we use std:: functions when we suspect they may be fast,
     * and otherwise we use our own optimized implementations.
     */
#ifdef __GNUC__
    template<> struct FP_ProbablyHasFastLibcComplex<float>
    { enum { result = true }; };
    template<> struct FP_ProbablyHasFastLibcComplex<double>
    { enum { result = true }; };
    template<> struct FP_ProbablyHasFastLibcComplex<long double>
    { enum { result = true }; };
#endif

    template<typename T>
    inline const std::complex<T> fp_make_imag(const std::complex<T>& v)
    {
        return std::complex<T> ( T(), v.real() );
    }

    template<typename T>
    inline std::complex<T> fp_real(const std::complex<T>& x)
    {
        return x.real();
    }
    template<typename T>
    inline std::complex<T> fp_imag(const std::complex<T>& x)
    {
        return x.imag();
    }
    template<typename T>
    inline std::complex<T> fp_arg(const std::complex<T>& x)
    {
        return std::arg(x);
    }
    template<typename T>
    inline std::complex<T> fp_conj(const std::complex<T>& x)
    {
        return std::conj(x);
    }
    template<typename T, bool>
    inline std::complex<T> fp_polar(const T& x, const T& y)
    {
        T si, co; fp_sinCos(si, co, y);
        return std::complex<T> (x*co, x*si);
    }
    template<typename T>
    inline std::complex<T> fp_polar(const std::complex<T>& x, const std::complex<T>& y)
    {
        // x * cos(y) + i * x * sin(y) -- arguments are supposed to be REAL numbers
        return fp_polar<T,true> (x.real(), y.real());
        //return std::polar(x.real(), y.real());
        //return x * (fp_cos(y) + (std::complex<T>(0,1) * fp_sin(y));
    }

    // These provide fallbacks in case there's no library function
    template<typename T>
    inline std::complex<T> fp_floor(const std::complex<T>& x)
    {
        return std::complex<T> (fp_floor(x.real()), fp_floor(x.imag()));
    }
    template<typename T>
    inline std::complex<T> fp_trunc(const std::complex<T>& x)
    {
        return std::complex<T> (fp_trunc(x.real()), fp_trunc(x.imag()));
    }
    template<typename T>
    inline std::complex<T> fp_int(const std::complex<T>& x)
    {
        return std::complex<T> (fp_int(x.real()), fp_int(x.imag()));
    }
    template<typename T>
    inline std::complex<T> fp_ceil(const std::complex<T>& x)
    {
        return std::complex<T> (fp_ceil(x.real()), fp_ceil(x.imag()));
    }
    template<typename T>
    inline std::complex<T> fp_abs(const std::complex<T>& x)
    {
        return std::abs(x);
        //T extent = fp_max(fp_abs(x.real()), fp_abs(x.imag()));
        //if(extent == T()) return x;
        //return extent * fp_hypot(x.real() / extent, x.imag() / extent);
    }
    template<typename T>
    inline std::complex<T> fp_exp(const std::complex<T>& x)
    {
        if(FP_ProbablyHasFastLibcComplex<T>::result)
            return std::exp(x);
        return fp_polar<T,true>(fp_exp(x.real()), x.imag());
    }
    template<typename T>
    inline std::complex<T> fp_log(const std::complex<T>& x)
    {
        if(FP_ProbablyHasFastLibcComplex<T>::result)
            return std::log(x);
        // log(abs(x))        + i*arg(x)
        // log(Xr^2+Xi^2)*0.5 + i*arg(x)
        if(x.imag()==T())
            return std::complex<T>( fp_log(fp_abs(x.real())),
                                    fp_arg(x.real()) ); // Note: Uses real-value fp_arg() here!
        return std::complex<T>(
            fp_log(std::norm(x)) * T(0.5),
            fp_arg(x).real() );
    }
    template<typename T>
    inline std::complex<T> fp_sqrt(const std::complex<T>& x)
    {
        if(FP_ProbablyHasFastLibcComplex<T>::result)
            return std::sqrt(x);
        return fp_polar<T,true> (fp_sqrt(fp_abs(x).real()),
                                 T(0.5)*fp_arg(x).real());
    }
    template<typename T>
    inline std::complex<T> fp_acos(const std::complex<T>& x)
    {
        // -i * log(x + i * sqrt(1 - x^2))
        const std::complex<T> i (T(), T(1));
        return -i *  fp_log(x + i * fp_sqrt(T(1) - x*x));
        // Note: Real version of acos() cannot handle |x| > 1,
        //       because it would cause sqrt(negative value).
    }
    template<typename T>
    inline std::complex<T> fp_asin(const std::complex<T>& x)
    {
        // -i * log(i*x + sqrt(1 - x^2))
        const std::complex<T> i (T(), T(1));
        return -i * fp_log(i*x + fp_sqrt(T(1) - x*x));
        // Note: Real version of asin() cannot handle |x| > 1,
        //       because it would cause sqrt(negative value).
    }
    template<typename T>
    inline std::complex<T> fp_atan(const std::complex<T>& x)
    {
        // 0.5i * (log(1-i*x) - log(1+i*x))
        // -0.5i * log( (1+i*x) / (1-i*x) )
        const std::complex<T> i (T(), T(1));
        return (T(-0.5)*i) * fp_log( (T(1)+i*x) / (T(1)-i*x) );
        // Note: x = -1i causes division by zero
        //       x = +1i causes log(0)
        // Thus, x must not be +-1i
    }
    template<typename T>
    inline std::complex<T> fp_cos(const std::complex<T>& x)
    {
        return std::cos(x);
        // // (exp(i*x) + exp(-i*x)) / (2)
        // //const std::complex<T> i (T(), T(1));
        // //return (fp_exp(i*x) + fp_exp(-i*x)) * T(0.5);
        // // Also: cos(Xr)*cosh(Xi) - i*sin(Xr)*sinh(Xi)
        // return std::complex<T> (
        //     fp_cos(x.real())*fp_cosh(x.imag()),
        //     -fp_sin(x.real())*fp_sinh(x.imag()));
    }
    template<typename T>
    inline std::complex<T> fp_sin(const std::complex<T>& x)
    {
        return std::sin(x);
        // // (exp(i*x) - exp(-i*x)) / (2i)
        // //const std::complex<T> i (T(), T(1));
        // //return (fp_exp(i*x) - fp_exp(-i*x)) * (T(-0.5)*i);
        // // Also: sin(Xr)*cosh(Xi) + cos(Xr)*sinh(Xi)
        // return std::complex<T> (
        //     fp_sin(x.real())*fp_cosh(x.imag()),
        //     fp_cos(x.real())*fp_sinh(x.imag()));
    }
    template<typename T>
    inline void fp_sinCos(
        std::complex<T>& sinvalue,
        std::complex<T>& cosvalue,
        const std::complex<T>& x)
    {
        //const std::complex<T> i (T(), T(1)), expix(fp_exp(i*x)), expmix(fp_exp((-i)*x));
        //cosvalue = (expix + expmix) * T(0.5);
        //sinvalue = (expix - expmix) * (i*T(-0.5));
        // The above expands to the following:
        T srx, crx; fp_sinCos(srx, crx, x.real());
        T six, cix; fp_sinhCosh(six, cix, x.imag());
        sinvalue = std::complex<T>(srx*cix,  crx*six);
        cosvalue = std::complex<T>(crx*cix, -srx*six);
    }
    template<typename T>
    inline void fp_sinhCosh(
        std::complex<T>& sinhvalue,
        std::complex<T>& coshvalue,
        const std::complex<T>& x)
    {
        T srx, crx; fp_sinhCosh(srx, crx, x.real());
        T six, cix; fp_sinCos(six, cix, x.imag());
        sinhvalue = std::complex<T>(srx*cix, crx*six);
        coshvalue = std::complex<T>(crx*cix, srx*six);
    }
    template<typename T>
    inline std::complex<T> fp_tan(const std::complex<T>& x)
    {
        return std::tan(x);
        //std::complex<T> si, co;
        //fp_sinCos(si, co, x);
        //return si/co;
        // // (i-i*exp(2i*x)) / (exp(2i*x)+1)
        // const std::complex<T> i (T(), T(1)), exp2ix=fp_exp((2*i)*x);
        // return (i-i*exp2ix) / (exp2ix+T(1));
        // // Also: sin(x)/cos(y)
        // // return fp_sin(x)/fp_cos(x);
    }
    template<typename T>
    inline std::complex<T> fp_cosh(const std::complex<T>& x)
    {
        return std::cosh(x);
        // // (exp(x) + exp(-x)) * 0.5
        // // Also: cosh(Xr)*cos(Xi) + i*sinh(Xr)*sin(Xi)
        // return std::complex<T> (
        //     fp_cosh(x.real())*fp_cos(x.imag()),
        //     fp_sinh(x.real())*fp_sin(x.imag()));
    }
    template<typename T>
    inline std::complex<T> fp_sinh(const std::complex<T>& x)
    {
        return std::sinh(x);
        // // (exp(x) - exp(-x)) * 0.5
        // // Also: sinh(Xr)*cos(Xi) + i*cosh(Xr)*sin(Xi)
        // return std::complex<T> (
        //     fp_sinh(x.real())*fp_cos(x.imag()),
        //     fp_cosh(x.real())*fp_sin(x.imag()));
    }
    template<typename T>
    inline std::complex<T> fp_tanh(const std::complex<T>& x)
    {
        return std::tanh(x);
        //std::complex<T> si, co;
        //fp_sinhCosh(si, co, x);
        //return si/co;
        // // (exp(2*x)-1) / (exp(2*x)+1)
        // // Also: sinh(x)/tanh(x)
        // const std::complex<T> exp2x=fp_exp(x+x);
        // return (exp2x-T(1)) / (exp2x+T(1));
    }

#ifdef FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS
    template<typename T>
    inline std::complex<T> fp_acosh(const std::complex<T>& x)
    { return fp_log(x + fp_sqrt(x*x - std::complex<T>(1))); }
    template<typename T>
    inline std::complex<T> fp_asinh(const std::complex<T>& x)
    { return fp_log(x + fp_sqrt(x*x + std::complex<T>(1))); }
    template<typename T>
    inline std::complex<T> fp_atanh(const std::complex<T>& x)
    { return fp_log( (std::complex<T>(1)+x) / (std::complex<T>(1)-x))
           * std::complex<T>(0.5); }
#endif
    template<typename T>
    inline std::complex<T> fp_pow(const std::complex<T>& x, const std::complex<T>& y)
    {
        // return std::pow(x,y);

        // With complex numbers, pow(x,y) can be solved with
        // the general formula: exp(y*log(x)). It handles
        // all special cases gracefully.
        // It expands to the following:
        // A)
        //     t1 = log(x)
        //     t2 = y * t1
        //     res = exp(t2)
        // B)
        //     t1.r = log(x.r * x.r + x.i * x.i) * 0.5  \ fp_log()
        //     t1.i = atan2(x.i, x.r)                   /
        //     t2.r = y.r*t1.r - y.i*t1.i          \ multiplication
        //     t2.i = y.r*t1.i + y.i*t1.r          /
        //     rho   = exp(t2.r)        \ fp_exp()
        //     theta = t2.i             /
        //     res.r = rho * cos(theta)   \ fp_polar(), called from
        //     res.i = rho * sin(theta)   / fp_exp(). Uses sincos().
        // Aside from the common "norm" calculation in atan2()
        // and in the log parameter, both of which are part of fp_log(),
        // there does not seem to be any incentive to break this
        // function down further; it would not help optimizing it.
        // However, we do handle the following special cases:
        //
        // When x is real (positive or negative):
        //     t1.r = log(abs(x.r))
        //     t1.i = x.r<0 ? -pi : 0
        // When y is real:
        //     t2.r = y.r * t1.r
        //     t2.i = y.r * t1.i
        const std::complex<T> t =
            (x.imag() != T())
            ? fp_log(x)
            : std::complex<T> (fp_log(fp_abs(x.real())),
                               fp_arg(x.real())); // Note: Uses real-value fp_arg() here!
        return y.imag() != T()
            ? fp_exp(y * t)
            : fp_polar<T,true> (fp_exp(y.real()*t.real()), y.real()*t.imag());
    }
    template<typename T>
    inline std::complex<T> fp_cbrt(const std::complex<T>& x)
    {
        // For real numbers, prefer giving a real solution
        // rather than a complex solution.
        // For example, cbrt(-3) has the following three solutions:
        //  A) 0.7211247966535 + 1.2490247864016i
        //  B) 0.7211247966535 - 1.2490247864016i
        //  C) -1.442249593307
        // exp(log(x)/3) gives A, but we prefer to give C.
        if(x.imag() == T()) return fp_cbrt(x.real());
        const std::complex<T> t(fp_log(x));
        return fp_polar<T,true> (fp_exp(t.real() / T(3)), t.imag() / T(3));
    }

    template<typename T>
    inline std::complex<T> fp_exp2(const std::complex<T>& x)
    {
        // pow(2, x)
        // polar(2^Xr, Xi*log(2))
        return fp_polar<T,true> (fp_exp2(x.real()), x.imag()*fp_const_log2<T>());
    }
    template<typename T>
    inline std::complex<T> fp_mod(const std::complex<T>& x, const std::complex<T>& y)
    {
        // Modulo function is probably not defined for complex numbers.
        // But we do our best to calculate it the same way as it is done
        // with real numbers, so that at least it is in some way "consistent".
        if(y.imag() == 0) return fp_mod(x.real(), y.real()); // optimization
        std::complex<T> n = fp_trunc(x / y);
        return x - n * y;
    }

    /* libstdc++ already defines a streaming operator for complex values,
     * but we redefine our own that it is compatible with the input
     * accepted by fparser. I.e. instead of (5,3) we now get (5+3i),
     * and instead of (-4,0) we now get -4.
     */
    template<typename T>
    inline std::ostream& operator<<(std::ostream& os, const std::complex<T>& value)
    {
        if(value.imag() == T()) return os << value.real();
        if(value.real() == T()) return os << value.imag() << 'i';
        if(value.imag() < T())
            return os << '(' << value.real() << "-" << -value.imag() << "i)";
        else
            return os << '(' << value.real() << "+" << value.imag() << "i)";
    }

    /* Less-than or greater-than operators are not technically defined
     * for Complex types. However, in fparser and its tool set, these
     * operators are widely required to be present.
     * Our implementation here is based on converting the complex number
     * into a scalar and the doing a scalar comparison on the value.
     * The means by which the number is changed into a scalar is based
     * on the following principles:
     * - Does not introduce unjustified amounts of extra inaccuracy
     * - Is transparent to purely real values
     *     (this disqualifies something like x.real() + x.imag())
     * - Does not ignore the imaginary value completely
     *     (this may be relevant especially in testbed)
     * - Is not so complicated that it would slow down a great deal
     *
     * Basically our formula here is the same as std::abs(),
     * except that it keeps the sign of the original real number,
     * and it does not do a sqrt() calculation that is not really
     * needed because we are only interested in the relative magnitudes.
     *
     * Equality and nonequality operators must not need to be overloaded.
     * They are already implemented in standard, and we must
     * not introduce flawed equality assumptions.
     */
    template<typename T>
    inline T fp_complexScalarize(const std::complex<T>& x)
    {
        T res(std::norm(x));
        if(x.real() < T()) res = -res;
        return res;
    }
    template<typename T>
    inline T fp_realComplexScalarize(const T& x)
    {
        T res(x*x);
        if(x < T()) res = -res;
        return res;
    }
    //    { return x.real() * (T(1.0) + fp_abs(x.imag())); }
    #define d(op) \
    template<typename T> \
    inline bool operator op (const std::complex<T>& x, T y) \
        { return fp_complexScalarize(x) op fp_realComplexScalarize(y); } \
    template<typename T> \
    inline bool operator op (const std::complex<T>& x, const std::complex<T>& y) \
        { return fp_complexScalarize(x) op \
                 fp_complexScalarize(y); } \
    template<typename T> \
    inline bool operator op (T x, const std::complex<T>& y) \
        { return fp_realComplexScalarize(x) op fp_complexScalarize(y); }
    d( < ) d( <= ) d( > ) d( >= )
    #undef d
#endif

    template<typename Value_t>
    inline Value_t fp_real(const Value_t& x) { return x; }
    template<typename Value_t>
    inline Value_t fp_imag(const Value_t& ) { return Value_t(); }
    template<typename Value_t>
    inline Value_t fp_arg(const Value_t& x)
        { return x < Value_t() ? -fp_const_pi<Value_t>() : Value_t(); }
    template<typename Value_t>
    inline Value_t fp_conj(const Value_t& x) { return x; }
    template<typename Value_t>
    inline Value_t fp_polar(const Value_t& x, const Value_t& y)
        { return x * fp_cos(y); } // This is of course a meaningless function.

    template<typename Value_t>
    inline std::complex<Value_t> fp_atan2(const std::complex<Value_t>& y,
                                          const std::complex<Value_t>& x)
    {
        if(y == Value_t()) return fp_arg(x);
        if(x == Value_t()) return fp_const_pi<Value_t>() * Value_t(-0.5);
        // 2*atan(y / (sqrt(x^2+y^2) + x)    )
        // 2*atan(    (sqrt(x^2+y^2) - x) / y)
        std::complex<Value_t> res( fp_atan(y / (fp_hypot(x,y) + x)) );
        return res+res;
    }

// -------------------------------------------------------------------------
// Comparison
// -------------------------------------------------------------------------
    template<typename Value_t>
    inline bool fp_equal(const Value_t& x, const Value_t& y)
    { return IsIntType<Value_t>::result
            ? (x == y)
            : (fp_abs(x - y) <= Epsilon<Value_t>::value); }

    template<typename Value_t>
    inline bool fp_nequal(const Value_t& x, const Value_t& y)
    { return IsIntType<Value_t>::result
            ? (x != y)
            : (fp_abs(x - y) > Epsilon<Value_t>::value); }

    template<typename Value_t>
    inline bool fp_less(const Value_t& x, const Value_t& y)
    { return IsIntType<Value_t>::result
            ? (x < y)
            : (x < y - Epsilon<Value_t>::value); }

    template<typename Value_t>
    inline bool fp_lessOrEq(const Value_t& x, const Value_t& y)
    { return IsIntType<Value_t>::result
            ? (x <= y)
            : (x <= y + Epsilon<Value_t>::value); }


    template<typename Value_t>
    inline bool fp_greater(const Value_t& x, const Value_t& y)
    { return fp_less(y, x); }

    template<typename Value_t>
    inline bool fp_greaterOrEq(const Value_t& x, const Value_t& y)
    { return fp_lessOrEq(y, x); }

    template<typename Value_t>
    inline bool fp_truth(const Value_t& d)
    {
        return IsIntType<Value_t>::result
                ? d != Value_t()
                : fp_abs(d) >= Value_t(0.5);
    }

    template<typename Value_t>
    inline bool fp_absTruth(const Value_t& abs_d)
    {
        return IsIntType<Value_t>::result
                ? abs_d > Value_t()
                : abs_d >= Value_t(0.5);
    }

    template<typename Value_t>
    inline const Value_t& fp_min(const Value_t& d1, const Value_t& d2)
        { return d1<d2 ? d1 : d2; }

    template<typename Value_t>
    inline const Value_t& fp_max(const Value_t& d1, const Value_t& d2)
        { return d1>d2 ? d1 : d2; }

    template<typename Value_t>
    inline const Value_t fp_not(const Value_t& b)
        { return Value_t(!fp_truth(b)); }

    template<typename Value_t>
    inline const Value_t fp_notNot(const Value_t& b)
        { return Value_t(fp_truth(b)); }

    template<typename Value_t>
    inline const Value_t fp_absNot(const Value_t& b)
        { return Value_t(!fp_absTruth(b)); }

    template<typename Value_t>
    inline const Value_t fp_absNotNot(const Value_t& b)
        { return Value_t(fp_absTruth(b)); }

    template<typename Value_t>
    inline const Value_t fp_and(const Value_t& a, const Value_t& b)
        { return Value_t(fp_truth(a) && fp_truth(b)); }

    template<typename Value_t>
    inline const Value_t fp_or(const Value_t& a, const Value_t& b)
        { return Value_t(fp_truth(a) || fp_truth(b)); }

    template<typename Value_t>
    inline const Value_t fp_absAnd(const Value_t& a, const Value_t& b)
        { return Value_t(fp_absTruth(a) && fp_absTruth(b)); }

    template<typename Value_t>
    inline const Value_t fp_absOr(const Value_t& a, const Value_t& b)
        { return Value_t(fp_absTruth(a) || fp_absTruth(b)); }

    template<typename Value_t>
    inline const Value_t fp_make_imag(const Value_t& ) // Imaginary 1. In real mode, always zero.
    {
        return Value_t();
    }

    /////////////
    /* Opcode analysis functions are used by fp_opcode_add.inc */
    /* Moved here from fparser.cc because fp_opcode_add.inc
     * is also now included by fpoptimizer.cc
     */
    bool IsLogicalOpcode(unsigned op);
    bool IsComparisonOpcode(unsigned op);
    unsigned OppositeComparisonOpcode(unsigned op);
    bool IsNeverNegativeValueOpcode(unsigned op);
    bool IsAlwaysIntegerOpcode(unsigned op);
    bool IsUnaryOpcode(unsigned op);
    bool IsBinaryOpcode(unsigned op);
    bool IsVarOpcode(unsigned op);
    bool IsCommutativeOrParamSwappableBinaryOpcode(unsigned op);
    unsigned GetParamSwappedBinaryOpcode(unsigned op);

    template<bool ComplexType>
    bool HasInvalidRangesOpcode(unsigned op);

    template<typename Value_t>
    inline Value_t DegreesToRadians(const Value_t& degrees)
    {
        return degrees * fp_const_deg_to_rad<Value_t>();
    }

    template<typename Value_t>
    inline Value_t RadiansToDegrees(const Value_t& radians)
    {
        return radians * fp_const_rad_to_deg<Value_t>();
    }

    template<typename Value_t>
    inline long makeLongInteger(const Value_t& value)
    {
        return (long) fp_int(value);
    }

#ifdef FP_SUPPORT_COMPLEX_NUMBERS
    template<typename T>
    inline long makeLongInteger(const std::complex<T>& value)
    {
        return (long) fp_int( std::abs(value) );
    }
#endif

    // Is value an integer that fits in "long" datatype?
    template<typename Value_t>
    inline bool isLongInteger(const Value_t& value)
    {
        return value == Value_t( makeLongInteger(value) );
    }

    template<typename Value_t>
    inline bool isOddInteger(const Value_t& value)
    {
        const Value_t halfValue = (value + Value_t(1)) * Value_t(0.5);
        return fp_equal(halfValue, fp_floor(halfValue));
    }

    template<typename Value_t>
    inline bool isEvenInteger(const Value_t& value)
    {
        const Value_t halfValue = value * Value_t(0.5);
        return fp_equal(halfValue, fp_floor(halfValue));
    }

    template<typename Value_t>
    inline bool isInteger(const Value_t& value)
    {
        return fp_equal(value, fp_floor(value));
    }

#ifdef FP_SUPPORT_LONG_INT_TYPE
    template<>
    inline bool isEvenInteger(const long& value)
    {
        return value%2 == 0;
    }

    template<>
    inline bool isInteger(const long&) { return true; }

    template<>
    inline bool isLongInteger(const long&) { return true; }

    template<>
    inline long makeLongInteger(const long& value)
    {
        return value;
    }
#endif

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
    template<>
    inline bool isInteger(const MpfrFloat& value) { return value.isInteger(); }

    template<>
    inline bool isEvenInteger(const MpfrFloat& value)
    {
        return isInteger(value) && value%2 == 0;
    }

    template<>
    inline long makeLongInteger(const MpfrFloat& value)
    {
        return (long) value.toInt();
    }
#endif

#ifdef FP_SUPPORT_GMP_INT_TYPE
    template<>
    inline bool isEvenInteger(const GmpInt& value)
    {
        return value%2 == 0;
    }

    template<>
    inline bool isInteger(const GmpInt&) { return true; }

    template<>
    inline long makeLongInteger(const GmpInt& value)
    {
        return (long) value.toInt();
    }
#endif

#ifdef FP_SUPPORT_LONG_INT_TYPE
    template<>
    inline bool isOddInteger(const long& value)
    {
        return value%2 != 0;
    }
#endif

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
    template<>
    inline bool isOddInteger(const MpfrFloat& value)
    {
        return value.isInteger() && value%2 != 0;
    }
#endif

#ifdef FP_SUPPORT_GMP_INT_TYPE
    template<>
    inline bool isOddInteger(const GmpInt& value)
    {
        return value%2 != 0;
    }
#endif


// -------------------------------------------------------------------------
// fp_pow
// -------------------------------------------------------------------------
    // Commented versions in fparser.cc
    template<typename Value_t>
    inline Value_t fp_pow_with_exp_log(const Value_t& x, const Value_t& y)
    {
        return fp_exp(fp_log(x) * y);
    }

    template<typename Value_t>
    inline Value_t fp_powi(Value_t x, unsigned long y)
    {
        Value_t result(1);
        while(y != 0)
        {
            if(y & 1) { result *= x; y -= 1; }
            else      { x *= x;      y /= 2; }
        }
        return result;
    }

    template<typename Value_t>
    Value_t fp_pow(const Value_t& x, const Value_t& y)
    {
        if(x == Value_t(1)) return Value_t(1);
        if(isLongInteger(y))
        {
            if(y >= Value_t(0))
                return fp_powi(x, makeLongInteger(y));
            else
                return Value_t(1) / fp_powi(x, -makeLongInteger(y));
        }
        if(y >= Value_t(0))
        {
            if(x > Value_t(0)) return fp_pow_with_exp_log(x, y);
            if(x == Value_t(0)) return Value_t(0);
            if(!isInteger(y*Value_t(16)))
                return -fp_pow_with_exp_log(-x, y);
        }
        else
        {
            if(x > Value_t(0)) return fp_pow_with_exp_log(Value_t(1) / x, -y);
            if(x < Value_t(0))
            {
                if(!isInteger(y*Value_t(-16)))
                    return -fp_pow_with_exp_log(Value_t(-1) / x, -y);
            }
        }
        return fp_pow_base(x, y);
    }

    template<typename Value_t>
    inline Value_t fp_exp2(const Value_t& x)
    {
        return fp_pow(Value_t(2), x);
    }
} // namespace FUNCTIONPARSERTYPES

#endif // ONCE_FPARSER_H_
#endif // ONCE_FPARSER_AUX_H_
