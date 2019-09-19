#ifndef _SCTL_VEC_WRAPPER_HPP_
#define _SCTL_VEC_WRAPPER_HPP_

#include SCTL_INCLUDE(math_utils.hpp)
#include SCTL_INCLUDE(common.hpp)
#include <cstdint>
#include <ostream>

#ifdef __SSE__
#include <xmmintrin.h>
#endif
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __SSE3__
#include <pmmintrin.h>
#endif
#ifdef __SSE4_2__
#include <smmintrin.h>
#endif
#ifdef __AVX__
#include <immintrin.h>
#endif
#if defined(__MIC__)
#include <immintrin.h>
#endif

// TODO: Implement AVX versions of floats, int32_t, int64_t

// TODO: Check alignment when SCTL_MEMDEBUG is defined
// TODO: Replace pointers with iterators

namespace SCTL_NAMESPACE {

  enum class DataType {
    Integer,
    Real,
    Bool
  };

  template <class ValueType> class TypeTraits {
    public:
      static constexpr DataType Type = DataType::Bool;
      static constexpr Integer Size = sizeof(ValueType);
      static constexpr Integer SigBits = 1;
  };
  template <> class TypeTraits<int32_t> {
    public:
      static constexpr DataType Type = DataType::Integer;
      static constexpr Integer Size = sizeof(int32_t);
      static constexpr Integer SigBits = Size * 8;
  };
  template <> class TypeTraits<int64_t> {
    public:
      static constexpr DataType Type = DataType::Integer;
      static constexpr Integer Size = sizeof(int64_t);
      static constexpr Integer SigBits = Size * 8;
  };
  template <> class TypeTraits<float> {
    public:
      static constexpr DataType Type = DataType::Real;
      static constexpr Integer Size = sizeof(float);
      static constexpr Integer SigBits = 23;
  };
  template <> class TypeTraits<double> {
    public:
      static constexpr DataType Type = DataType::Real;
      static constexpr Integer Size = sizeof(double);
      static constexpr Integer SigBits = 52;
  };

  template <DataType type, Integer size> class GetType {
    public:
      typedef bool ValueType;
  };
  template <> class GetType<DataType::Integer,4> {
    public:
      typedef int32_t ValueType;
  };
  template <> class GetType<DataType::Integer,8> {
    public:
      typedef int64_t ValueType;
  };
  template <> class GetType<DataType::Real,4> {
    public:
      typedef float ValueType;
  };
  template <> class GetType<DataType::Real,8> {
    public:
      typedef double ValueType;
  };

  template <class ValueType, Integer N> class alignas(sizeof(ValueType) * N) Vec {

    public:

      typedef typename GetType<DataType::Integer,TypeTraits<ValueType>::Size>::ValueType IntegerType;
      typedef typename GetType<DataType::Real,TypeTraits<ValueType>::Size>::ValueType RealType;
      typedef Vec<IntegerType,N> IntegerVec;
      typedef Vec<RealType,N> RealVec;
      typedef ValueType ScalarType;

      static constexpr Integer Size() {
        return N;
      }

      static Vec Zero() {
        Vec r;
        for (Integer i = 0; i < N; i++) r.v[i] = 0;
        return r;
      }

      static Vec Load1(ValueType const* p) {
        Vec r;
        for (Integer i = 0; i < N; i++) r.v[i] = p[0];
        return r;
      }
      static Vec Load(ValueType const* p) {
        Vec r;
        for (Integer i = 0; i < N; i++) r.v[i] = p[i];
        return r;
      }
      static Vec LoadAligned(ValueType const* p) {
        Vec r;
        for (Integer i = 0; i < N; i++) r.v[i] = p[i];
        return r;
      }

      Vec() = default;

      Vec(const ValueType& a) {
        for (Integer i = 0; i < N; i++) v[i] = a;
      }

      void Store(ValueType* p) const {
        for (Integer i = 0; i < N; i++) p[i] = v[i];
      }
      void StoreAligned(ValueType* p) const {
        for (Integer i = 0; i < N; i++) p[i] = v[i];
      }

      // Bitwise NOT
      Vec operator~() const {
        Vec r;
        char* vo = (char*)r.v;
        const char* vi = (const char*)this->v;
        for (Integer i = 0; i < (Integer)(N*sizeof(ValueType)); i++) vo[i] = ~vi[i];
        return r;
      }

      // Unary plus and minus
      Vec operator+() const {
        return *this;
      }
      Vec operator-() const {
        Vec r;
        for (Integer i = 0; i < N; i++) r.v[i] = -v[i];
        return r;
      }

      // C-style cast
      //template <class RetValueType> explicit operator Vec<RetValueType,N>() const {
      //  Vec<RetValueType,N> r;
      //  for (Integer i = 0; i < N; i++) r.v[i] = (RetValueType)v[i];
      //  return r;
      //}

      // Arithmetic operators
      friend Vec operator*(Vec lhs, const Vec& rhs) {
        for (Integer i = 0; i < N; i++) lhs.v[i] *= rhs.v[i];
        return lhs;
      }
      friend Vec operator+(Vec lhs, const Vec& rhs) {
        for (Integer i = 0; i < N; i++) lhs.v[i] += rhs.v[i];
        return lhs;
      }
      friend Vec operator-(Vec lhs, const Vec& rhs) {
        for (Integer i = 0; i < N; i++) lhs.v[i] -= rhs.v[i];
        return lhs;
      }
      friend Vec FMA(Vec a, const Vec& b, const Vec& c) {
        for (Integer i = 0; i < N; i++) a.v[i] = a.v[i] * b.v[i] + c.v[i];
        return a;
      }

      // Comparison operators
      friend Vec operator< (Vec lhs, const Vec& rhs) {
        static const ValueType value_zero = const_zero();
        static const ValueType value_one = const_one();
        for (Integer i = 0; i < N; i++) lhs.v[i] = (lhs.v[i] < rhs.v[i] ? value_one : value_zero);
        return lhs;
      }
      friend Vec operator<=(Vec lhs, const Vec& rhs) {
        static const ValueType value_zero = const_zero();
        static const ValueType value_one = const_one();
        for (Integer i = 0; i < N; i++) lhs.v[i] = (lhs.v[i] <= rhs.v[i] ? value_one : value_zero);
        return lhs;
      }
      friend Vec operator>=(Vec lhs, const Vec& rhs) {
        static const ValueType value_zero = const_zero();
        static const ValueType value_one = const_one();
        for (Integer i = 0; i < N; i++) lhs.v[i] = (lhs.v[i] >= rhs.v[i] ? value_one : value_zero);
        return lhs;
      }
      friend Vec operator> (Vec lhs, const Vec& rhs) {
        static const ValueType value_zero = const_zero();
        static const ValueType value_one = const_one();
        for (Integer i = 0; i < N; i++) lhs.v[i] = (lhs.v[i] > rhs.v[i] ? value_one : value_zero);
        return lhs;
      }
      friend Vec operator==(Vec lhs, const Vec& rhs) {
        static const ValueType value_zero = const_zero();
        static const ValueType value_one = const_one();
        for (Integer i = 0; i < N; i++) lhs.v[i] = (lhs.v[i] == rhs.v[i] ? value_one : value_zero);
        return lhs;
      }
      friend Vec operator!=(Vec lhs, const Vec& rhs) {
        static const ValueType value_zero = const_zero();
        static const ValueType value_one = const_one();
        for (Integer i = 0; i < N; i++) lhs.v[i] = (lhs.v[i] != rhs.v[i] ? value_one : value_zero);
        return lhs;
      }

      // Bitwise operators
      friend Vec operator&(Vec lhs, const Vec& rhs) {
        char* vo = (char*)lhs.v;
        const char* vi = (const char*)rhs.v;
        for (Integer i = 0; i < (Integer)sizeof(ValueType)*N; i++) vo[i] &= vi[i];
        return lhs;
      }
      friend Vec operator^(Vec lhs, const Vec& rhs) {
        char* vo = (char*)lhs.v;
        const char* vi = (const char*)rhs.v;
        for (Integer i = 0; i < (Integer)sizeof(ValueType)*N; i++) vo[i] ^= vi[i];
        return lhs;
      }
      friend Vec operator|(Vec lhs, const Vec& rhs) {
        char* vo = (char*)lhs.v;
        const char* vi = (const char*)rhs.v;
        for (Integer i = 0; i < (Integer)sizeof(ValueType)*N; i++) vo[i] |= vi[i];
        return lhs;
      }
      friend Vec AndNot(Vec lhs, const Vec& rhs) {
        return lhs & (~rhs);
      }

      // Bitshift
      friend IntegerVec operator<<(const Vec& lhs, const Integer& rhs) {
        IntegerVec r = IntegerVec::LoadAligned(&lhs.v[0]);
        for (Integer i = 0; i < N; i++) r.v[i] = r.v[i] << rhs;
        return r;
      }

      // Assignment operators
      Vec& operator+=(const Vec& rhs) {
        for (Integer i = 0; i < N; i++) v[i] += rhs.v[i];
        return *this;
      }
      Vec& operator-=(const Vec& rhs) {
        for (Integer i = 0; i < N; i++) v[i] -= rhs.v[i];
        return *this;
      }
      Vec& operator*=(const Vec& rhs) {
        for (Integer i = 0; i < N; i++) v[i] *= rhs.v[i];
        return *this;
      }
      Vec& operator&=(const Vec& rhs) {
        char* vo = (char*)this->v;
        const char* vi = (const char*)rhs.v;
        for (Integer i = 0; i < (Integer)sizeof(ValueType)*N; i++) vo[i] &= vi[i];
        return *this;
      }
      Vec& operator^=(const Vec& rhs) {
        char* vo = (char*)this->v;
        const char* vi = (const char*)rhs.v;
        for (Integer i = 0; i < (Integer)sizeof(ValueType)*N; i++) vo[i] ^= vi[i];
        return *this;
      }
      Vec& operator|=(const Vec& rhs) {
        char* vo = (char*)this->v;
        const char* vi = (const char*)rhs.v;
        for (Integer i = 0; i < (Integer)sizeof(ValueType)*N; i++) vo[i] |= vi[i];
        return *this;
      }

      // Conversion operators
      // /

      // Other operators
      friend Vec max(Vec lhs, const Vec& rhs) {
        for (Integer i = 0; i < N; i++) {
          if (lhs.v[i] < rhs.v[i]) lhs.v[i] = rhs.v[i];
        }
        return lhs;
      }
      friend Vec min(Vec lhs, const Vec& rhs) {
        for (Integer i = 0; i < N; i++) {
          if (lhs.v[i] > rhs.v[i]) lhs.v[i] = rhs.v[i];
        }
        return lhs;
      }

      friend std::ostream& operator<<(std::ostream& os, const Vec& in) {
        //for (Integer i = 0; i < (Integer)sizeof(ValueType)*8; i++) os << ((*(uint64_t*)in.v) & (1UL << i) ? '1' : '0');
        //os << '\n';
        for (Integer i = 0; i < N; i++) os << in.v[i] << ' ';
        return os;
      }
      friend Vec approx_rsqrt(const Vec& x) {
        Vec r;
        for (int i = 0; i < N; i++) r.v[i] = 1 / sqrt<ValueType>(x.v[i]);
        return r;
      }

      template <class Vec1, class Vec2> friend Vec1 reinterpret(const Vec2& x);

    private:

      static const ValueType const_zero() {
        union {
          ValueType value;
          unsigned char cvalue[sizeof(ValueType)];
        };
        for (Integer i = 0; i < (Integer)sizeof(ValueType); i++) cvalue[i] = 0;
        return value;
      }
      static const ValueType const_one() {
        union {
          ValueType value;
          unsigned char cvalue[sizeof(ValueType)];
        };
        for (Integer i = 0; i < (Integer)sizeof(ValueType); i++) cvalue[i] = ~(unsigned char)0;
        return value;
      }

      ValueType v[N];
  };

  // Other operators
  template <class RetVec, class Vec> RetVec reinterpret(const Vec& v){
    static_assert(sizeof(RetVec) == sizeof(Vec));
    RetVec& r = *(RetVec*)&v;
    return r;
  }
  template <class RealVec, class IntVec> RealVec ConvertInt2Real(const IntVec& x) {
    typedef typename RealVec::ScalarType Real;
    typedef typename IntVec::ScalarType Int;
    assert(sizeof(RealVec) == sizeof(IntVec));
    assert(sizeof(Real) == sizeof(Int));
    static constexpr Integer SigBits = TypeTraits<Real>::SigBits;
    union {
      Int Cint = (1UL << (SigBits - 1)) + ((SigBits + ((1UL<<(sizeof(Real)*8 - SigBits - 2))-1)) << SigBits);
      Real Creal;
    };
    IntVec l(x + IntVec(Cint));
    return *(RealVec*)&l - RealVec(Creal);
  }
  template <class Vec> typename Vec::IntegerVec RoundReal2Int(const Vec& x) {
    using IntegerType = typename Vec::IntegerType;
    using RealType = typename Vec::RealType;
    using IntegerVec = typename Vec::IntegerVec;
    using RealVec = typename Vec::RealVec;
    static_assert(std::is_same<RealVec,Vec>::value, "RoundReal2Int: expected real input argument!");

    static constexpr Integer SigBits = TypeTraits<RealType>::SigBits;
    union {
      IntegerType Cint = (1UL << (SigBits - 1)) + ((SigBits + ((1UL<<(sizeof(RealType)*8 - SigBits - 2))-1)) << SigBits);
      RealType Creal;
    };
    RealVec d = x + RealVec(Creal);
    return reinterpret<IntegerVec>(d) - IntegerVec(Cint);
  }
  template <class Vec> Vec RoundReal2Real(const Vec& x) {
    typedef typename Vec::ScalarType Real;
    static constexpr Integer SigBits = TypeTraits<Real>::SigBits;
    union {
      int64_t Cint = (1UL << (SigBits - 1)) + ((SigBits + ((1UL<<(sizeof(Real)*8 - SigBits - 2))-1)) << SigBits);
      Real Creal;
    };
    Vec Vreal(Creal);
    return (x + Vreal) - Vreal;
  }
  template <class Vec> void sincos_intrin(Vec& sinx, Vec& cosx, const Vec& x) {
    constexpr Integer ORDER = 13;
    // ORDER    ERROR
    //     1 8.81e-02
    //     3 2.45e-03
    //     5 3.63e-05
    //     7 3.11e-07
    //     9 1.75e-09
    //    11 6.93e-12
    //    13 2.09e-14
    //    15 6.66e-16
    //    17 6.66e-16

    using Real = typename Vec::ScalarType;
    static constexpr Integer SigBits = TypeTraits<Real>::SigBits;
    static constexpr Real coeff3  = -1/(((Real)2)*3);
    static constexpr Real coeff5  =  1/(((Real)2)*3*4*5);
    static constexpr Real coeff7  = -1/(((Real)2)*3*4*5*6*7);
    static constexpr Real coeff9  =  1/(((Real)2)*3*4*5*6*7*8*9);
    static constexpr Real coeff11 = -1/(((Real)2)*3*4*5*6*7*8*9*10*11);
    static constexpr Real coeff13 =  1/(((Real)2)*3*4*5*6*7*8*9*10*11*12*13);
    static constexpr Real coeff15 = -1/(((Real)2)*3*4*5*6*7*8*9*10*11*12*13*14*15);
    static constexpr Real coeff17 =  1/(((Real)2)*3*4*5*6*7*8*9*10*11*12*13*14*15*16*17);
    static constexpr Real coeff19 = -1/(((Real)2)*3*4*5*6*7*8*9*10*11*12*13*14*15*16*17*18*19);
    static constexpr Real x0 = (Real)1.570796326794896619231321691639l;
    static constexpr Real invx0 = 1 / x0;

    Vec x_ = RoundReal2Real(x * invx0); // 4.5 - cycles
    Vec x1 = x - x_ * x0; // 2 - cycles
    Vec x2, x3, x5, x7, x9, x11, x13, x15, x17, x19;

    Vec s1 = x1;
    if (ORDER >= 3) { // 5 - cycles
      x2 = x1 * x1;
      x3 = x1 * x2;
      s1 += x3 * coeff3;
    }
    if (ORDER >= 5) { // 3 - cycles
      x5 = x3 * x2;
      s1 += x5 * coeff5;
    }
    if (ORDER >= 7) {
      x7 = x5 * x2;
      s1 += x7 * coeff7;
    }
    if (ORDER >= 9) {
      x9 = x7 * x2;
      s1 += x9 * coeff9;
    }
    if (ORDER >= 11) {
      x11 = x9 * x2;
      s1 += x11 * coeff11;
    }
    if (ORDER >= 13) {
      x13 = x11 * x2;
      s1 += x13 * coeff13;
    }
    if (ORDER >= 15) {
      x15 = x13 * x2;
      s1 += x15 * coeff15;
    }
    if (ORDER >= 17) {
      x17 = x15 * x2;
      s1 += x17 * coeff17;
    }
    if (ORDER >= 19) {
      x19 = x17 * x2;
      s1 += x19 * coeff19;
    }

    Vec cos_squared = (Real)1.0 - s1 * s1;
    Vec inv_cos = approx_rsqrt(cos_squared); // 1.5 - cycles
    if (ORDER < 5) {
    } else if (ORDER < 9) {
      inv_cos *= ((3.0) - cos_squared * inv_cos * inv_cos) * 0.5; // 7 - cycles
    } else if (ORDER < 15) {
      inv_cos *= ((3.0) - cos_squared * inv_cos * inv_cos); // 7 - cycles
      inv_cos *= ((3.0 * pow<pow<0>(3)*3-1>(2.0)) - cos_squared * inv_cos * inv_cos) * (pow<(pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
    } else {
      inv_cos *= ((3.0) - cos_squared * inv_cos * inv_cos); // 7 - cycles
      inv_cos *= ((3.0 * pow<pow<0>(3)*3-1>(2.0)) - cos_squared * inv_cos * inv_cos); // 7 - cycles
      inv_cos *= ((3.0 * pow<pow<1>(3)*3-1>(2.0)) - cos_squared * inv_cos * inv_cos) * (pow<(pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
    }
    Vec c1 = cos_squared * inv_cos; // 1 - cycle

    union {
      int64_t int_zero = 0 + (1UL << (SigBits - 1)) + ((SigBits + ((1UL<<(sizeof(Real)*8 - SigBits - 2))-1)) << SigBits);
      Real real_zero;
    };
    union {
      int64_t int_one = 1 + (1UL << (SigBits - 1)) + ((SigBits + ((1UL<<(sizeof(Real)*8 - SigBits - 2))-1)) << SigBits);
      Real real_one;
    };
    union {
      int64_t int_two = 2 + (1UL << (SigBits - 1)) + ((SigBits + ((1UL<<(sizeof(Real)*8 - SigBits - 2))-1)) << SigBits);
      Real real_two;
    };
    Vec x_offset(real_zero);
    auto xAnd1 = (((x_+x_offset) & Vec(real_one)) == x_offset);
    auto xAnd2 = (((x_+x_offset) & Vec(real_two)) == x_offset);

    Vec s2 = AndNot( c1,xAnd1) | (s1 & xAnd1);
    Vec c2 = AndNot(-s1,xAnd1) | (c1 & xAnd1);
    Vec s3 = AndNot(-s2,xAnd2) | (s2 & xAnd2);
    Vec c3 = AndNot(-c2,xAnd2) | (c2 & xAnd2);

    sinx = s3;
    cosx = c3;
  }
  template <class Vec> void exp_intrin(Vec& expx, const Vec& x) {
    constexpr Integer ORDER = 10;
    using IntegerType = typename Vec::IntegerType;
    using RealType = typename Vec::RealType;
    using IntegerVec = typename Vec::IntegerVec;
    using RealVec = typename Vec::RealVec;
    static_assert(std::is_same<Vec,RealVec>::value, "exp_intrin: expected a real argument");

    using Real = typename RealVec::ScalarType;
    static constexpr Integer SigBits = TypeTraits<Real>::SigBits;
    static constexpr Real coeff2  = 1/(((Real)2));
    static constexpr Real coeff3  = 1/(((Real)2)*3);
    static constexpr Real coeff4  = 1/(((Real)2)*3*4);
    static constexpr Real coeff5  = 1/(((Real)2)*3*4*5);
    static constexpr Real coeff6  = 1/(((Real)2)*3*4*5*6);
    static constexpr Real coeff7  = 1/(((Real)2)*3*4*5*6*7);
    static constexpr Real coeff8  = 1/(((Real)2)*3*4*5*6*7*8);
    static constexpr Real coeff9  = 1/(((Real)2)*3*4*5*6*7*8*9);
    static constexpr Real coeff10 = 1/(((Real)2)*3*4*5*6*7*8*9*10);
    static constexpr Real x0 = (Real)0.693147180559945309417232121458l; // ln(2)
    static constexpr Real invx0 = 1 / x0;

    RealVec x_ = RoundReal2Real(x * invx0); // 4.5 - cycles
    IntegerVec int_x_ = RoundReal2Int<RealVec>(x_);
    RealVec x1 = x - x_ * x0; // 2 - cycles
    RealVec x2, x3, x4, x5, x6, x7, x8, x9, x10;

    RealVec e1 = 1.0 + x1;
    if (ORDER >= 2) {
      x2 = x1 * x1;
      e1 += x2 * coeff2;
    }
    if (ORDER >= 3) {
      x3 = x2 * x1;
      e1 += x3 * coeff3;
    }
    if (ORDER >= 4) {
      x4 = x2 * x2;
      e1 += x4 * coeff4;
    }
    if (ORDER >= 5) {
      x5 = x3 * x2;
      e1 += x5 * coeff5;
    }
    if (ORDER >= 6) {
      x6 = x3 * x3;
      e1 += x6 * coeff6;
    }
    if (ORDER >= 7) {
      x7 = x4 * x3;
      e1 += x7 * coeff7;
    }
    if (ORDER >= 8) {
      x8 = x4 * x4;
      e1 += x8 * coeff8;
    }
    if (ORDER >= 9) {
      x9 = x5 * x4;
      e1 += x9 * coeff9;
    }
    if (ORDER >= 10) {
      x10 = x5 * x5;
      e1 += x10 * coeff10;
    }

    RealVec e2;
    { // set e2 = 2 ^ x_
      union {
        RealType real_one = 1.0;
        IntegerType int_one;
      };
      //__m256i int_e2 = _mm256_add_epi64(
      //                                  _mm256_set1_epi64x(int_one),
      //                                  _mm256_slli_epi64(
      //                                                      _mm256_load_si256((__m256i const*)&int_x_),
      //                                                      SigBits
      //                                                    )
      //                                  ); // int_e2 = int_one + (int_x_ << SigBits);
      IntegerVec int_e2 = IntegerVec(int_one) + (int_x_ << SigBits);

      // Handle underflow
      static constexpr IntegerType max_exp = -(IntegerType)(1UL<<((sizeof(Real)*8-SigBits-2)));
      int_e2 &= (int_x_ > IntegerVec(max_exp));

      e2 = RealVec::LoadAligned((RealType*)&int_e2);
    }

    expx = e1 * e2;
  }

#ifdef __AVX__
  template <> class alignas(sizeof(double)*4) Vec<double,4> {
    typedef __m256d VecType;
    typedef double ValueType;
    static constexpr Integer N = 4;
    public:

      typedef typename GetType<DataType::Integer,TypeTraits<ValueType>::Size>::ValueType IntegerType;
      typedef typename GetType<DataType::Real,TypeTraits<ValueType>::Size>::ValueType RealType;
      typedef Vec<IntegerType,N> IntegerVec;
      typedef Vec<RealType,N> RealVec;
      typedef ValueType ScalarType;

      static constexpr Integer Size() {
        return N;
      }

      static Vec Zero() {
        Vec r;
        r.v = _mm256_setzero_pd();
        return r;
      }

      static Vec Load1(ValueType const* p) {
        Vec r;
        r.v = _mm256_broadcast_sd(p);
        return r;
      }
      static Vec Load(ValueType const* p) {
        Vec r;
        r.v = _mm256_loadu_pd(p);
        return r;
      }
      static Vec LoadAligned(ValueType const* p) {
        Vec r;
        r.v = _mm256_load_pd(p);
        return r;
      }

      Vec() = default;

      Vec(const ValueType& a) {
        v = _mm256_set1_pd(a);
      }

      void Store(ValueType* p) const {
        _mm256_storeu_pd(p, v);
      }
      void StoreAligned(ValueType* p) const {
        _mm256_store_pd(p, v);
      }

      // Bitwise NOT
      Vec operator~() const {
        Vec r;
        static constexpr ScalarType Creal = -1.0;
        r.v = _mm256_xor_pd(v, _mm256_set1_pd(Creal));
        return r;
      }

      // Unary plus and minus
      Vec operator+() const {
        return *this;
      }
      Vec operator-() const {
        return Zero() - (*this);
      }

      // C-style cast
      //template <class RetValueType> explicit operator Vec<RetValueType,N>() const {
      //}

      // Arithmetic operators
      friend Vec operator*(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_mul_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec operator+(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_add_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec operator-(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_sub_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec FMA(Vec a, const Vec& b, const Vec& c) {
        #ifdef __FMA__
        a.v = _mm256_fmadd_pd(a.v, b.v, c.v);
        #else
        a.v = _mm256_add_pd(_mm256_mul_pd(a.v, b.v), c.v);
        #endif
        return a;
      }

      // Comparison operators
      friend Vec operator< (Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_cmp_pd(lhs.v, rhs.v, _CMP_LT_OS);
        return lhs;
      }
      friend Vec operator<=(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_cmp_pd(lhs.v, rhs.v, _CMP_LE_OS);
        return lhs;
      }
      friend Vec operator>=(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_cmp_pd(lhs.v, rhs.v, _CMP_GE_OS);
        return lhs;
      }
      friend Vec operator> (Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_cmp_pd(lhs.v, rhs.v, _CMP_GT_OS);
        return lhs;
      }
      friend Vec operator==(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_cmp_pd(lhs.v, rhs.v, _CMP_EQ_OS);
        return lhs;
      }
      friend Vec operator!=(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_cmp_pd(lhs.v, rhs.v, _CMP_NEQ_OS);
        return lhs;
      }

      // Bitwise operators
      friend Vec operator&(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_and_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec operator^(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_xor_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec operator|(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_or_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec AndNot(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_andnot_pd(rhs.v, lhs.v);
        return lhs;
      }

      // Assignment operators
      Vec& operator*=(const Vec& rhs) {
        v = _mm256_mul_pd(v, rhs.v);
        return *this;
      }
      Vec& operator+=(const Vec& rhs) {
        v = _mm256_add_pd(v, rhs.v);
        return *this;
      }
      Vec& operator-=(const Vec& rhs) {
        v = _mm256_sub_pd(v, rhs.v);
        return *this;
      }
      Vec& operator&=(const Vec& rhs) {
        v = _mm256_and_pd(v, rhs.v);
        return *this;
      }
      Vec& operator^=(const Vec& rhs) {
        v = _mm256_xor_pd(v, rhs.v);
        return *this;
      }
      Vec& operator|=(const Vec& rhs) {
        v = _mm256_or_pd(v, rhs.v);
        return *this;
      }

      // Other operators
      friend Vec max(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_max_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec min(Vec lhs, const Vec& rhs) {
        lhs.v = _mm256_min_pd(lhs.v, rhs.v);
        return lhs;
      }

      friend std::ostream& operator<<(std::ostream& os, const Vec& in) {
        union {
          VecType vec;
          ValueType val[N];
        };
        vec = in.v;
        for (Integer i = 0; i < N; i++) os << val[i] << ' ';
        return os;
      }
      friend Vec approx_rsqrt(const Vec& x) {
        Vec r;
        r.v = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x.v)));
        return r;
      }

      template <class Vec1, class Vec2> friend Vec1 reinterpret(const Vec2& x);
      template <class Vec> friend Vec RoundReal2Real(const Vec& x);
      template <class Vec> friend void sincos_intrin(Vec& sinx, Vec& cosx, const Vec& x);
      template <class Vec> friend void exp_intrin(Vec& expx, const Vec& x);

    private:

      VecType v;
  };

  template <> inline Vec<int64_t,4> reinterpret<Vec<int64_t,4>,Vec<double,4>>(const Vec<double,4>& x){
    union {
      Vec<int64_t,4> r;
      __m256i y;
    };
    y = _mm256_castpd_si256(x.v);
    return r;
  }

  template <> inline Vec<double,4> RoundReal2Real(const Vec<double,4>& x) {
    Vec<double,4> r;
    r.v = _mm256_round_pd(x.v,_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);
    return r;
  }

  #ifdef SCTL_HAVE_SVML
  template <> inline void sincos_intrin(Vec<double,4>& sinx, Vec<double,4>& cosx, const Vec<double,4>& x) {
    sinx.v = _mm256_sin_pd(x.v);
    cosx.v = _mm256_cos_pd(x.v);
  }

  template <> inline void exp_intrin(Vec<double,4>& expx, const Vec<double,4>& x) {
    expx.v = _mm256_exp_pd(x.v);
  }
  #endif

#endif

#ifdef __AVX512F__
  template <> class alignas(sizeof(double)*8) Vec<double,8> {
    typedef __m512d VecType;
    typedef double ValueType;
    static constexpr Integer N = 8;
    public:

      typedef typename GetType<DataType::Integer,TypeTraits<ValueType>::Size>::ValueType IntegerType;
      typedef typename GetType<DataType::Real,TypeTraits<ValueType>::Size>::ValueType RealType;
      typedef Vec<IntegerType,N> IntegerVec;
      typedef Vec<RealType,N> RealVec;
      typedef ValueType ScalarType;

      static constexpr Integer Size() {
        return N;
      }

      static Vec Zero() {
        Vec r;
        r.v = _mm512_setzero_pd();
        return r;
      }

      static Vec Load1(ValueType const* p) {
        Vec r;
        // TODO: different from _m256d, could make it faster?
        // r.v = _mm512_broadcast_f64x4(_mm256_broadcast_sd(p));
        r.v = _mm512_set1_pd(*p);
        return r;
      }
      static Vec Load(ValueType const* p) {
        Vec r;
        r.v = _mm512_loadu_pd(p);
        return r;
      }
      static Vec LoadAligned(ValueType const* p) {
        Vec r;
        r.v = _mm512_load_pd(p);
        return r;
      }

      Vec() = default;

      Vec(const ValueType& a) {
        v = _mm512_set1_pd(a);
      }

      Vec(const __mmask8& a) = delete; // disallow implicit conversions

      void Store(ValueType* p) const {
        _mm512_storeu_pd(p, v);
      }
      void StoreAligned(ValueType* p) const {
        _mm512_store_pd(p, v);
      }

      // Bitwise NOT
      Vec operator~() const {
        Vec r;
        static constexpr ScalarType Creal = -1.0;
        r.v = _mm512_xor_pd(v, _mm512_set1_pd(Creal));
        return r;
      }

      // Unary plus and minus
      Vec operator+() const {
        return *this;
      }
      Vec operator-() const {
        return Zero() - (*this);
      }

      // C-style cast
      //template <class RetValueType> explicit operator Vec<RetValueType,N>() const {
      //}

      // Arithmetic operators
      friend Vec operator*(Vec lhs, const Vec& rhs) {
        lhs.v = _mm512_mul_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec operator+(Vec lhs, const Vec& rhs) {
        lhs.v = _mm512_add_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec operator-(Vec lhs, const Vec& rhs) {
        lhs.v = _mm512_sub_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec FMA(Vec a, const Vec& b, const Vec& c) {
        a.v = _mm512_fmadd_pd(a.v, b.v, c.v);
        //a.v = _mm512_add_pd(_mm512_mul_pd(a.v, b.v), c.v);
        return a;
      }

      // Comparison operators
      //friend Vec operator< (Vec lhs, const Vec& rhs) {
      //  lhs.v = _mm512_castsi512_pd(_mm512_movm_epi64(_mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_LT_OS)));
      //  return lhs;
      //}
      //friend Vec operator<=(Vec lhs, const Vec& rhs) {
      //  lhs.v = _mm512_castsi512_pd(_mm512_movm_epi64(_mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_LE_OS)));
      //  return lhs;
      //}
      //friend Vec operator>=(Vec lhs, const Vec& rhs) {
      //  lhs.v = _mm512_castsi512_pd(_mm512_movm_epi64(_mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_GE_OS)));
      //  return lhs;
      //}
      //friend Vec operator> (Vec lhs, const Vec& rhs) {
      //  lhs.v = _mm512_castsi512_pd(_mm512_movm_epi64(_mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_GT_OS)));
      //  return lhs;
      //}
      //friend Vec operator==(Vec lhs, const Vec& rhs) {
      //  lhs.v = _mm512_castsi512_pd(_mm512_movm_epi64(_mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_EQ_OS)));
      //  return lhs;
      //}
      //friend Vec operator!=(Vec lhs, const Vec& rhs) {
      //  lhs.v = _mm512_castsi512_pd(_mm512_movm_epi64(_mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_NEQ_OS)));
      //  return lhs;
      //}

      friend __mmask8 operator< (Vec lhs, const Vec& rhs) {
        return _mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_LT_OS);
      }
      friend __mmask8 operator<=(Vec lhs, const Vec& rhs) {
        return _mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_LE_OS);
      }
      friend __mmask8 operator>=(Vec lhs, const Vec& rhs) {
        return _mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_GE_OS);
      }
      friend __mmask8 operator> (Vec lhs, const Vec& rhs) {
        return _mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_GT_OS);
      }
      friend __mmask8 operator==(Vec lhs, const Vec& rhs) {
        return _mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_EQ_OS);
      }
      friend __mmask8 operator!=(Vec lhs, const Vec& rhs) {
        return _mm512_cmp_pd_mask(lhs.v, rhs.v, _CMP_NEQ_OS);
      }

      // Bitwise operators
      friend Vec operator&(Vec lhs, const Vec& rhs) {
        lhs.v = _mm512_and_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec operator^(Vec lhs, const Vec& rhs) {
        lhs.v = _mm512_xor_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec operator|(Vec lhs, const Vec& rhs) {
        lhs.v = _mm512_or_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec AndNot(Vec lhs, const Vec& rhs) {
        lhs.v = _mm512_andnot_pd(rhs.v, lhs.v);
        return lhs;
      }
      friend Vec operator&(Vec lhs, const __mmask8& rhs) {
        lhs.v = _mm512_maskz_mov_pd(rhs, lhs.v);
        return lhs;
      }
      friend Vec AndNot(Vec lhs, const __mmask8& rhs) {
        lhs.v = _mm512_mask_mov_pd(lhs.v, rhs, _mm512_setzero_pd());
        return lhs;
      }

      // Assignment operators
      Vec& operator*=(const Vec& rhs) {
        v = _mm512_mul_pd(v, rhs.v);
        return *this;
      }
      Vec& operator+=(const Vec& rhs) {
        v = _mm512_add_pd(v, rhs.v);
        return *this;
      }
      Vec& operator-=(const Vec& rhs) {
        v = _mm512_sub_pd(v, rhs.v);
        return *this;
      }
      Vec& operator&=(const Vec& rhs) {
        v = _mm512_and_pd(v, rhs.v);
        return *this;
      }
      Vec& operator^=(const Vec& rhs) {
        v = _mm512_xor_pd(v, rhs.v);
        return *this;
      }
      Vec& operator|=(const Vec& rhs) {
        v = _mm512_or_pd(v, rhs.v);
        return *this;
      }
      Vec& operator&=(const __mmask8& rhs) {
        v = _mm512_maskz_mov_pd(rhs, v);
        return *this;
      }

      // Other operators
      friend Vec max(Vec lhs, const Vec& rhs) {
        lhs.v = _mm512_max_pd(lhs.v, rhs.v);
        return lhs;
      }
      friend Vec min(Vec lhs, const Vec& rhs) {
        lhs.v = _mm512_min_pd(lhs.v, rhs.v);
        return lhs;
      }

      friend std::ostream& operator<<(std::ostream& os, const Vec& in) {
        union {
          VecType vec;
          ValueType val[N];
        };
        vec = in.v;
        for (Integer i = 0; i < N; i++) os << val[i] << ' ';
        return os;
      }
      friend Vec approx_rsqrt(const Vec& x) {
        Vec r;
        r.v = _mm512_cvtps_pd(_mm256_rsqrt_ps(_mm512_cvtpd_ps(x.v)));
        return r;
      }

      template <class Vec1, class Vec2> friend Vec1 reinterpret(const Vec2& x);
      template <class Vec> friend Vec RoundReal2Real(const Vec& x);
      template <class Vec> friend void sincos_intrin(Vec& sinx, Vec& cosx, const Vec& x);
      template <class Vec> friend void exp_intrin(Vec& expx, const Vec& x);

    private:

      VecType v;
  };

  template <> inline Vec<int64_t,8> reinterpret<Vec<int64_t,8>,Vec<double,8>>(const Vec<double,8>& x){
    union {
      Vec<int64_t,8> r;
      __m512i y;
    };
    y = _mm512_castpd_si512(x.v);
    return r;
  }

  template <> inline Vec<double,8> RoundReal2Real(const Vec<double,8>& x) {
    Vec<double,8> r;
    // TODO: need double check
    r.v = _mm512_roundscale_pd(x.v,_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);
    return r;
  }

  #ifdef SCTL_HAVE_SVML
  template <> inline void sincos_intrin(Vec<double,8>& sinx, Vec<double,8>& cosx, const Vec<double,8>& x) {
    sinx.v = _mm512_sin_pd(x.v);
    cosx.v = _mm512_cos_pd(x.v);
  }

  template <> inline void exp_intrin(Vec<double,8>& expx, const Vec<double,8>& x) {
    expx.v = _mm512_exp_pd(x.v);
  }
  #endif

#endif

}

#endif  //_SCTL_VEC_WRAPPER_HPP_
