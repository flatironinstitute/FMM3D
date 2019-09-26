#ifndef _SCTL_MAT_UTILS_
#define _SCTL_MAT_UTILS_

#include SCTL_INCLUDE(common.hpp)

namespace SCTL_NAMESPACE {
namespace mat {

template <class T> void gemm(char TransA, char TransB, int M, int N, int K, T alpha, T *A, int lda, T *B, int ldb, T beta, T *C, int ldc);

template <class T> void cublasgemm(char TransA, char TransB, int M, int N, int K, T alpha, T *A, int lda, T *B, int ldb, T beta, T *C, int ldc);

template <class T> void svd(char *JOBU, char *JOBVT, int *M, int *N, Iterator<T> A, int *LDA, Iterator<T> S, Iterator<T> U, int *LDU, Iterator<T> VT, int *LDVT, Iterator<T> WORK, int *LWORK, int *INFO);

/**
 * \brief Computes the pseudo inverse of matrix M(n1xn2) (in row major form)
 * and returns the output M_(n2xn1).
 */
template <class T> void pinv(Iterator<T> M, int n1, int n2, T eps, Iterator<T> M_);

}  // end namespace mat
}  // end namespace SCTL_NAMESPACE

#include SCTL_INCLUDE(mat_utils.txx)

#endif  //_SCTL_MAT_UTILS_
