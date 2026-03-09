#ifndef _SCTL_COMMON_HPP_
#define _SCTL_COMMON_HPP_

#ifndef SCTL_NAMESPACE
#define SCTL_NAMESPACE sctl
#endif
#define SCTL_QUOTEME(x) SCTL_QUOTEME_1(x)
#define SCTL_QUOTEME_1(x) #x
#define SCTL_INCLUDE(x) SCTL_QUOTEME(SCTL_NAMESPACE/x)

// Profiling parameters
#ifndef SCTL_PROFILE
#define SCTL_PROFILE -1 // Granularity level
#endif

#if defined(__AVX512__) || defined(__AVX512F__)
  #define SCTL_ALIGN_BYTES 64
#elif defined(__AVX__)
  #define SCTL_ALIGN_BYTES 32
#elif defined(__SSE__)
  #define SCTL_ALIGN_BYTES 16
#else
  #define SCTL_ALIGN_BYTES 8
#endif

// Parameters for memory manager
#ifndef SCTL_MEM_ALIGN
#define SCTL_MEM_ALIGN (64 > SCTL_ALIGN_BYTES ? 64 : SCTL_ALIGN_BYTES)
#endif
#ifndef SCTL_GLOBAL_MEM_BUFF
#define SCTL_GLOBAL_MEM_BUFF 1024LL * 0LL  // in MB
#endif

// Define NULL
#ifndef NULL
#define NULL 0
#endif

#include <cstddef>
#include <cstdint>
#if defined(_WIN32)
// MinGW lacks the POSIX drand48/srand48 API used by SCTL.
inline uint64_t& sctl_drand48_state() {
  static uint64_t state = 0x1234abcd330eULL;
  return state;
}

inline void srand48(long seedval) {
  sctl_drand48_state() = ((static_cast<uint64_t>(seedval) << 16) + 0x330eULL) & ((1ULL << 48) - 1);
}

inline double drand48() {
  uint64_t& state = sctl_drand48_state();
  state = (0x5deece66dULL * state + 0xbULL) & ((1ULL << 48) - 1);
  return static_cast<double>(state) / static_cast<double>(1ULL << 48);
}
#endif

namespace SCTL_NAMESPACE {
typedef long Integer;  // bounded numbers < 32k
typedef int64_t Long;  // problem size
}

#include <iostream>

#define SCTL_WARN(msg)                                         \
  do {                                                          \
    std::cerr << "\n\033[1;31mWarning:\033[0m " << msg << '\n'; \
  } while (0)

#define SCTL_ERROR(msg)                                      \
  do {                                                        \
    std::cerr << "\n\033[1;31mError:\033[0m " << msg << '\n'; \
    abort();                                                  \
  } while (0)

#define SCTL_ASSERT_MSG(cond, msg) \
  do {                              \
    if (!(cond)) SCTL_ERROR(msg);  \
  } while (0)

#define SCTL_ASSERT(cond)                                                                                      \
  do {                                                                                                          \
    if (!(cond)) {                                                                                              \
      fprintf(stderr, "\n%s:%d: %s: Assertion `%s' failed.\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond); \
      abort();                                                                                                  \
    }                                                                                                           \
  } while (0)

#define SCTL_UNUSED(x) (void)(x)  // to ignore unused variable warning.

namespace SCTL_NAMESPACE {
#ifdef SCTL_MEMDEBUG
template <class ValueType> class ConstIterator;
template <class ValueType> class Iterator;
template <class ValueType, Long DIM> class StaticArray;
#else
template <typename ValueType> using Iterator = ValueType*;
template <typename ValueType> using ConstIterator = const ValueType*;
template <typename ValueType, Long DIM> using StaticArray = ValueType[DIM];
#endif
}

#endif  //_SCTL_COMMON_HPP_
