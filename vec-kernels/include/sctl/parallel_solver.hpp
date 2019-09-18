#ifndef _SCTL_PARALLEL_SOLVER_HPP_
#define _SCTL_PARALLEL_SOLVER_HPP_

#include SCTL_INCLUDE(vector.hpp)
#include SCTL_INCLUDE(comm.hpp)

#include <functional>

namespace SCTL_NAMESPACE {

template <class Real> class ParallelSolver {

 public:
  using ParallelOp = std::function<void(Vector<Real>*, const Vector<Real>&)>;

  ParallelSolver(const Comm& comm = Comm::Self(), bool verbose = true) : comm_(comm), verbose_(verbose) {}

  void operator()(Vector<Real>* x, const ParallelOp& A, const Vector<Real>& b, Real tol, Integer max_iter = -1, bool use_abs_tol = false);

 private:
  Comm comm_;
  bool verbose_;
};

}  // end namespace

#ifdef SCTL_HAVE_PETSC

#include <petscksp.h>

namespace SCTL_NAMESPACE {

template <class Real> int ParallelSolverMatVec(Mat M_, Vec x_, Vec Mx_) {
  PetscErrorCode ierr;

  PetscInt N, N_;
  VecGetLocalSize(x_, &N);
  VecGetLocalSize(Mx_, &N_);
  SCTL_ASSERT(N == N_);

  void* data = nullptr;
  MatShellGetContext(M_, &data);
  auto& M = dynamic_cast<const typename ParallelSolver<Real>::ParallelOp&>(*(typename ParallelSolver<Real>::ParallelOp*)data);

  const PetscScalar* x_ptr;
  ierr = VecGetArrayRead(x_, &x_ptr);
  CHKERRQ(ierr);

  Vector<Real> x(N);
  for (Long i = 0; i < N; i++) x[i] = x_ptr[i];
  Vector<Real> Mx(N);
  M(&Mx, x);

  PetscScalar* Mx_ptr;
  ierr = VecGetArray(Mx_, &Mx_ptr);
  CHKERRQ(ierr);

  for (long i = 0; i < N; i++) Mx_ptr[i] = Mx[i];
  ierr = VecRestoreArray(Mx_, &Mx_ptr);
  CHKERRQ(ierr);

  return 0;
}

template <class Real> inline void ParallelSolver<Real>::operator()(Vector<Real>* x, const ParallelOp& A, const Vector<Real>& b, Real tol, Integer max_iter, bool use_abs_tol) {
  PetscInt N = b.Dim();
  if (max_iter < 0) max_iter = N;
  MPI_Comm comm = comm_.GetMPI_Comm();
  PetscErrorCode ierr;

  Mat PetscA;
  {  // Create Matrix. PetscA
    MatCreateShell(comm, N, N, PETSC_DETERMINE, PETSC_DETERMINE, (void*)&A, &PetscA);
    MatShellSetOperation(PetscA, MATOP_MULT, (void (*)(void))ParallelSolverMatVec<Real>);
  }

  Vec Petsc_x, Petsc_b;
  {  // Create vectors
    VecCreateMPI(comm, N, PETSC_DETERMINE, &Petsc_b);
    VecCreateMPI(comm, N, PETSC_DETERMINE, &Petsc_x);

    PetscScalar* b_ptr;
    ierr = VecGetArray(Petsc_b, &b_ptr);
    CHKERRABORT(comm, ierr);
    for (long i = 0; i < N; i++) b_ptr[i] = b[i];
    ierr = VecRestoreArray(Petsc_b, &b_ptr);
    CHKERRABORT(comm, ierr);
  }

  // Create linear solver context
  KSP ksp;
  ierr = KSPCreate(comm, &ksp);
  CHKERRABORT(comm, ierr);

  // Set operators. Here the matrix that defines the linear system
  // also serves as the preconditioning matrix.
  ierr = KSPSetOperators(ksp, PetscA, PetscA);
  CHKERRABORT(comm, ierr);

  // Set runtime options
  KSPSetType(ksp, KSPGMRES);
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
  if (use_abs_tol) KSPSetTolerances(ksp, PETSC_DEFAULT, tol, PETSC_DEFAULT, max_iter);
  else KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
  KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
  //if (verbose_) KSPMonitorSet(ksp, KSPMonitorDefault, nullptr, nullptr); // Doesn't work for some versions of PETSc!! WTH!!
  KSPGMRESSetRestart(ksp, max_iter);
  ierr = KSPSetFromOptions(ksp);
  CHKERRABORT(comm, ierr);

  // -------------------------------------------------------------------
  // Solve the linear system: Ax=b
  // -------------------------------------------------------------------
  ierr = KSPSolve(ksp, Petsc_b, Petsc_x);
  CHKERRABORT(comm, ierr);

  // View info about the solver
  // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(comm, ierr);

  // Iterations
  // PetscInt its;
  // ierr = KSPGetIterationNumber(ksp,&its); CHKERRABORT(comm, ierr);
  // ierr = PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its); CHKERRABORT(comm, ierr);

  {  // Set x
    const PetscScalar* x_ptr;
    ierr = VecGetArrayRead(Petsc_x, &x_ptr);
    CHKERRABORT(comm, ierr);

    if (x->Dim() != N) x->ReInit(N);
    for (long i = 0; i < N; i++) (*x)[i] = x_ptr[i];
  }

  ierr = KSPDestroy(&ksp);
  CHKERRABORT(comm, ierr);
  ierr = MatDestroy(&PetscA);
  CHKERRABORT(comm, ierr);
  ierr = VecDestroy(&Petsc_x);
  CHKERRABORT(comm, ierr);
  ierr = VecDestroy(&Petsc_b);
  CHKERRABORT(comm, ierr);
}

}  // end namespace

#else

namespace SCTL_NAMESPACE {

template <class Real> static Real inner_prod(const Vector<Real>& x, const Vector<Real>& y, const Comm& comm) {
  Real x_dot_y = 0;
  Long N = x.Dim();
  SCTL_ASSERT(y.Dim() == N);
  for (Long i = 0; i < N; i++) x_dot_y += x[i] * y[i];

  Real x_dot_y_glb = 0;
  comm.Allreduce(Ptr2ConstItr<Real>(&x_dot_y, 1), Ptr2Itr<Real>(&x_dot_y_glb, 1), 1, Comm::CommOp::SUM);

  return x_dot_y_glb;
}

template <class Real> inline void ParallelSolver<Real>::operator()(Vector<Real>* x, const ParallelOp& A, const Vector<Real>& b, Real tol, Integer max_iter, bool use_abs_tol) {
  Long N = b.Dim();
  if (max_iter < 0) max_iter = N;

  { // Initialize x
    if (x->Dim() != N) x->ReInit(N);
    x->SetZero();
  }
  Real b_norm = sqrt(inner_prod(b, b, comm_));
  if (b_norm == 0) return;

  Vector<Real> q(N);
  Vector<Vector<Real>> Q;
  Vector<Vector<Real>> H;
  {  // Initialize q, Q
    q = b;
    Real one_over_q_norm = 1.0 / sqrt<Real>(inner_prod(q, q, comm_));
    for (Long j = 0; j < N; j++) q[j] *= one_over_q_norm;
    Q.PushBack(q);
  }

  Matrix<Real> H_;
  Vector<Real> Aq(N), y, h, r = b;
  Real abs_tol = tol * (use_abs_tol ? 1 : b_norm);
  while (1) {
    Real r_norm = sqrt(inner_prod(r, r, comm_));
    if (verbose_ && !comm_.Rank()) printf("%3lld KSP Residual norm %.12e\n", (long long)H.Dim(), r_norm);
    if (r_norm < abs_tol || H.Dim() == max_iter) break;

    A(&Aq, q);
    q = Aq;

    h.ReInit(Q.Dim() + 1);
    for (Integer i = 0; i < Q.Dim(); i++) {  // Orthogonalized q
      h[i] = inner_prod(q, Q[i], comm_);
      for (Long j = 0; j < N; j++) q[j] -= h[i] * Q[i][j];
    }
    {  // Normalize q
      h[Q.Dim()] = sqrt<Real>(inner_prod(q, q, comm_));
      if (h[Q.Dim()] == 0) break;
      Real one_over_q_norm = 1.0 / h[Q.Dim()];
      for (Long j = 0; j < N; j++) q[j] *= one_over_q_norm;
    }
    Q.PushBack(q);
    H.PushBack(h);

    {  // Set y
      H_.ReInit(H.Dim(), Q.Dim());
      H_.SetZero();
      for (Integer i = 0; i < H.Dim(); i++) {
        for (Integer j = 0; j < H[i].Dim(); j++) {
          H_[i][j] = H[i][j];
        }
      }
      H_ = H_.pinv();

      y.ReInit(H_.Dim(1));
      for (Integer i = 0; i < y.Dim(); i++) {
        y[i] = H_[0][i] * b_norm;
      }
    }

    {  // Compute residual
      Vector<Real> Hy(Q.Dim());
      Hy.SetZero();
      for (Integer i = 0; i < H.Dim(); i++) {
        for (Integer j = 0; j < H[i].Dim(); j++) {
          Hy[j] += H[i][j] * y[i];
        }
      }

      Vector<Real> QHy(N);
      QHy.SetZero();
      for (Integer i = 0; i < Q.Dim(); i++) {
        for (Long j = 0; j < N; j++) {
          QHy[j] += Q[i][j] * Hy[i];
        }
      }

      for (Integer j = 0; j < N; j++) {  // Set r
        r[j] = b[j] - QHy[j];
      }
    }
  }

  {  // Set x
    for (Integer i = 0; i < y.Dim(); i++) {
      for (Integer j = 0; j < N; j++) {
        (*x)[j] += y[i] * Q[i][j];
      }
    }
  }
}

}  // end namespace

#endif

#endif  //_SCTL_PARALLEL_SOLVER_HPP_
