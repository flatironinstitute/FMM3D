#ifndef _VEC_KERNELS_HPP_
#define _VEC_KERNELS_HPP_

#define NDEBUG
#include <sctl.hpp>


template <class Real> void h3ddirectcp_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
  static constexpr int32_t COORD_DIM = 3;
  static constexpr int32_t KDIM0 = 2;
  static constexpr int32_t KDIM1 = 2;
  long Ntrg = nt[0];
  long Nsrc = ns[0];
  long nd_ = nd[0];

  Real thresh2 = thresh[0] * thresh[0];
  //#pragma omp parallel for schedule(static)
  //for (long t = 0; t < Ntrg; t++) {
    //for (long s = 0; s < Nsrc; s++) {
  for (long s = 0; s < Nsrc; s++) {
    #pragma omp parallel for simd simdlen(32) schedule(static)
    for (long t = 0; t < Ntrg; t++) {
       Real dX[3];
       dX[0] = ztarg[t*COORD_DIM+0] - sources[s*COORD_DIM+0];
       dX[1] = ztarg[t*COORD_DIM+1] - sources[s*COORD_DIM+1];
       dX[2] = ztarg[t*COORD_DIM+2] - sources[s*COORD_DIM+2];
       Real R2 = dX[0]*dX[0] + dX[1]*dX[1] + dX[2]*dX[2];
       Real Rinv = (R2 > thresh2 ? 1/sqrt(R2) : 0);

       Real R = R2 * Rinv;
       Real G0 = cos(zk[0]*R) * exp(-zk[1]*R) * Rinv;
       Real G1 = sin(zk[0]*R) * exp(-zk[1]*R) * Rinv;

       for (long k = 0; k < nd_; k++) {
         pot[(t*nd_+k)*KDIM1+0] += G0 * charge[(s*nd_+k)*KDIM0+0] - G1 * charge[(s*nd_+k)*KDIM0+1];
         pot[(t*nd_+k)*KDIM1+1] += G1 * charge[(s*nd_+k)*KDIM0+0] + G0 * charge[(s*nd_+k)*KDIM0+1];
       }
    }
  }
}


// charge, potential 
template <class Real, sctl::Integer MaxVecLen=4> void h3ddirectcp_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;
  static constexpr sctl::Integer KDIM0 = 2;
  static constexpr sctl::Integer KDIM1 = 2;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_*KDIM0, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_*KDIM1, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_*KDIM0, sctl::Ptr2Itr<Real>((Real*)charge , nd_*KDIM0*Nsrc), false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    transpose(Xt, Xt_);
    Vt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec zk_[2];
  zk_[0].Load1(zk+0);
  zk_[1].Load1(zk+1);
  Vec thresh2 = thresh[0] * thresh[0];
  // load charge and source
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  sctl::Matrix<Real> Vs_(Nsrc, nd_*KDIM0, sctl::Ptr2Itr<Real>((Real*)charge , nd_*KDIM0*Nsrc), false);
  //sctl::Vector<Vec> Xsrc(Nsrc*COORD_DIM);
  sctl::Vector<Vec> Vsrc(Nsrc*nd_*KDIM0);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    //for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      //Xsrc[s*COORD_DIM+k] = Vec::Load1(&Xs_[s][k]);
    //}
    for (long i = 0; i < nd_; i++) {
      Vsrc[s*nd_*KDIM0+i*KDIM0+0] = Vec::Load1(&Vs_[s][i*KDIM0+0]);
      Vsrc[s*nd_*KDIM0+i*KDIM0+1] = Vec::Load1(&Vs_[s][i*KDIM0+1]);
    }
  }
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential
    Vec Vtrg[nd_][KDIM1];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0] = Vec::LoadAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1] = Vec::LoadAligned(&Vt[i*KDIM1+1][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec R = R2 * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::sincos_intrin<Vec>(sin_izkR, cos_izkR, izkR[1]);
      sctl::exp_intrin(exp_izkR, izkR[0]);
      Vec G0 = cos_izkR * exp_izkR * Rinv;
      Vec G1 = sin_izkR * exp_izkR * Rinv;

      for (long i = 0; i < nd_; i++) {
        Vtrg[i][0] += G0 * Vsrc[s*nd_*KDIM0+i*KDIM0+0] - G1 * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        Vtrg[i][1] += G1 * Vsrc[s*nd_*KDIM0+i*KDIM0+0] + G0 * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
      }
    }
    // store
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0].StoreAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1].StoreAligned(&Vt[i*KDIM1+1][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_*KDIM1; j++) {
      pot[i*nd_*KDIM1+j] += Vt[j][i];
    }
  }
}


// charge, potential and gradient
template <class Real, sctl::Integer MaxVecLen=4> void h3ddirectcg_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;
  static constexpr sctl::Integer KDIM0 = 2;
  static constexpr sctl::Integer KDIM1 = 2;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_*KDIM0, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_*KDIM1, Ntrg_);
  sctl::Matrix<Real> Gt(nd_*KDIM1*COORD_DIM, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_*KDIM0, sctl::Ptr2Itr<Real>((Real*)charge , nd_*KDIM0*Nsrc), false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    transpose(Xt, Xt_);
    Vt = 0;
    Gt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec zk_[2];
  zk_[0].Load1(zk+0);
  zk_[1].Load1(zk+1);
  Vec thresh2 = thresh[0] * thresh[0];
  // load charge and source
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  sctl::Matrix<Real> Vs_(Nsrc, nd_*KDIM0, sctl::Ptr2Itr<Real>((Real*)charge , nd_*KDIM0*Nsrc), false);
  //sctl::Vector<Vec> Xsrc(Nsrc*COORD_DIM);
  sctl::Vector<Vec> Vsrc(Nsrc*nd_*KDIM0);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    //for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      //Xsrc[s*COORD_DIM+k] = Vec::Load1(&Xs_[s][k]);
    //}
    for (long i = 0; i < nd_; i++) {
      Vsrc[s*nd_*KDIM0+i*KDIM0+0] = Vec::Load1(&Vs_[s][i*KDIM0+0]);
      Vsrc[s*nd_*KDIM0+i*KDIM0+1] = Vec::Load1(&Vs_[s][i*KDIM0+1]);
    }
  }
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential and gradient
    Vec Vtrg[nd_][KDIM1], Gtrg[nd_][COORD_DIM][KDIM1];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0] = Vec::LoadAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1] = Vec::LoadAligned(&Vt[i*KDIM1+1][t]);
      Gtrg[i][0][0] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+0][t]);
      Gtrg[i][0][1] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+1][t]);
      Gtrg[i][1][0] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+0][t]);
      Gtrg[i][1][1] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+1][t]);
      Gtrg[i][2][0] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+0][t]);
      Gtrg[i][2][1] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+1][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::sincos_intrin<Vec>(sin_izkR, cos_izkR, izkR[1]);
      sctl::exp_intrin(exp_izkR, izkR[0]);
      // exp(ikr)/r
      Vec G0 = cos_izkR * exp_izkR * Rinv;
      Vec G1 = sin_izkR * exp_izkR * Rinv;
      // (ikr-1)*exp(ikr)/r^3
      Vec H0 = (izkR[0] - (1.0))*G0 - izkR[1]*G1;
      Vec H1 = izkR[1]*G0 + (izkR[0] - (1.0))*G1;
      H0 = H0 * Rinv2;
      H1 = H1 * Rinv2;
      // (ikr-1)*exp(ikr)/r^3 * (xt - xs)
      Vec Ztmp[COORD_DIM][KDIM0];
      for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          Ztmp[dimi][0] = H0*dX[dimi];
          Ztmp[dimi][1] = H1*dX[dimi];
      }

      for (long i = 0; i < nd_; i++) {
        Vtrg[i][0] += G0 * Vsrc[s*nd_*KDIM0+i*KDIM0+0] - G1 * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        Vtrg[i][1] += G1 * Vsrc[s*nd_*KDIM0+i*KDIM0+0] + G0 * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          Gtrg[i][dimi][0] += Ztmp[dimi][0] * Vsrc[s*nd_*KDIM0+i*KDIM0+0] - Ztmp[dimi][1] * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
          Gtrg[i][dimi][1] += Ztmp[dimi][1] * Vsrc[s*nd_*KDIM0+i*KDIM0+0] + Ztmp[dimi][0] * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        }
      }
    }
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0].StoreAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1].StoreAligned(&Vt[i*KDIM1+1][t]);
      Gtrg[i][0][0].StoreAligned(&Gt[(0*nd_+i)*KDIM1+0][t]);
      Gtrg[i][0][1].StoreAligned(&Gt[(0*nd_+i)*KDIM1+1][t]);
      Gtrg[i][1][0].StoreAligned(&Gt[(1*nd_+i)*KDIM1+0][t]);
      Gtrg[i][1][1].StoreAligned(&Gt[(1*nd_+i)*KDIM1+1][t]);
      Gtrg[i][2][0].StoreAligned(&Gt[(2*nd_+i)*KDIM1+0][t]);
      Gtrg[i][2][1].StoreAligned(&Gt[(2*nd_+i)*KDIM1+1][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_*KDIM1; j++) {
      pot[i*nd_*KDIM1+j] += Vt[j][i];
    }
    for (long j = 0; j < nd_*KDIM1*COORD_DIM; j++) {
      grad[i*nd_*KDIM1*COORD_DIM+j] += Gt[j][i];
    }
  }
}


// dipole, potential 
template <class Real, sctl::Integer MaxVecLen=4> void h3ddirectdp_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;
  static constexpr sctl::Integer KDIM0 = 2;
  static constexpr sctl::Integer KDIM1 = 2;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Gs(nd_*KDIM0*COORD_DIM, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_*KDIM1, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Gs_(Nsrc, nd_*KDIM0*COORD_DIM, sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*KDIM0*COORD_DIM*Nsrc), false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Gs, Gs_);
    transpose(Xt, Xt_);
    Vt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec zk_[2];
  zk_[0].Load1(zk+0);
  zk_[1].Load1(zk+1);
  Vec thresh2 = thresh[0] * thresh[0];
  // load dipvec and source
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  sctl::Matrix<Real> Gs_(Nsrc, nd_*KDIM0*COORD_DIM, sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*KDIM0*COORD_DIM*Nsrc), false);
  //sctl::Vector<Vec> Xsrc(Nsrc*COORD_DIM);
  sctl::Vector<Vec> Gsrc(Nsrc*nd_*COORD_DIM*KDIM0);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    //for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      //Xsrc[s*COORD_DIM+k] = Vec::Load1(&Xs_[s][k]);
    //}
    long s_ind = s*nd_*COORD_DIM*KDIM0;
    for (long i = 0; i < nd_; i++) {
      long si_ind = s_ind+i*COORD_DIM*KDIM0;
      Gsrc[si_ind+0] = Vec::Load1(&Gs_[s][(0*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+1] = Vec::Load1(&Gs_[s][(0*nd_+i)*KDIM0+1]);
      Gsrc[si_ind+2] = Vec::Load1(&Gs_[s][(1*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+3] = Vec::Load1(&Gs_[s][(1*nd_+i)*KDIM0+1]);
      Gsrc[si_ind+4] = Vec::Load1(&Gs_[s][(2*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+5] = Vec::Load1(&Gs_[s][(2*nd_+i)*KDIM0+1]);
    }
  }
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential
    Vec Vtrg[nd_][KDIM1];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0] = Vec::LoadAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1] = Vec::LoadAligned(&Vt[i*KDIM1+1][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::sincos_intrin<Vec>(sin_izkR, cos_izkR, izkR[1]);
      sctl::exp_intrin(exp_izkR, izkR[0]);
      // exp(ikr)/r
      Vec G0 = cos_izkR * exp_izkR * Rinv;
      Vec G1 = sin_izkR * exp_izkR * Rinv;
      Vec tmp0 = (1.0)-izkR[0];
      Vec tmp1 = -izkR[1];
      // (1-ikr)*exp(ikr)/r^3
      Vec H0 = (tmp0*G0 - tmp1*G1) * Rinv2;
      Vec H1 = (tmp1*G0 + tmp0*G1) * Rinv2;

      long s_ind = s*nd_*COORD_DIM*KDIM0;
      for (long i = 0; i < nd_; i++) {
        long si_ind = s_ind + i*COORD_DIM*KDIM0;
        Vec D0 = dX[0]*Gsrc[si_ind+0] + dX[1]*Gsrc[si_ind+2] + dX[2]*Gsrc[si_ind+4];
        Vec D1 = dX[0]*Gsrc[si_ind+1] + dX[1]*Gsrc[si_ind+3] + dX[2]*Gsrc[si_ind+5];
        Vtrg[i][0] += H0 * D0 - H1 * D1;
        Vtrg[i][1] += H1 * D0 + H0 * D1;
      }
    }
    // store
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0].StoreAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1].StoreAligned(&Vt[i*KDIM1+1][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_*KDIM1; j++) {
      pot[i*nd_*KDIM1+j] += Vt[j][i];
    }
  }
}


// dipole, potential and gradient
template <class Real, sctl::Integer MaxVecLen=4> void h3ddirectdg_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;
  static constexpr sctl::Integer KDIM0 = 2;
  static constexpr sctl::Integer KDIM1 = 2;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_*KDIM0*COORD_DIM, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_*KDIM1, Ntrg_);
  sctl::Matrix<Real> Gt(nd_*KDIM1*COORD_DIM, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_*KDIM0*COORD_DIM, sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*KDIM0*COORD_DIM*Nsrc), false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    transpose(Xt, Xt_);
    Vt = 0;
    Gt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec zk_[2];
  zk_[0].Load1(zk+0);
  zk_[1].Load1(zk+1);
  Vec thresh2 = thresh[0] * thresh[0];
  // load dipvec and source
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  sctl::Matrix<Real> Gs_(Nsrc, nd_*KDIM0*COORD_DIM, sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*KDIM0*COORD_DIM*Nsrc), false);
  //sctl::Vector<Vec> Xsrc(Nsrc*COORD_DIM);
  sctl::Vector<Vec> Gsrc(Nsrc*nd_*COORD_DIM*KDIM0);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    //for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      //Xsrc[s*COORD_DIM+k] = Vec::Load1(&Xs_[s][k]);
    //}
    long s_ind = s*nd_*COORD_DIM*KDIM0;
    for (long i = 0; i < nd_; i++) {
      long si_ind = s_ind+i*COORD_DIM*KDIM0;
      Gsrc[si_ind+0] = Vec::Load1(&Gs_[s][(0*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+1] = Vec::Load1(&Gs_[s][(0*nd_+i)*KDIM0+1]);
      Gsrc[si_ind+2] = Vec::Load1(&Gs_[s][(1*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+3] = Vec::Load1(&Gs_[s][(1*nd_+i)*KDIM0+1]);
      Gsrc[si_ind+4] = Vec::Load1(&Gs_[s][(2*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+5] = Vec::Load1(&Gs_[s][(2*nd_+i)*KDIM0+1]);
    }
  }
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential and gradient
    Vec Vtrg[nd_][KDIM1], Gtrg[nd_][COORD_DIM][KDIM1];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0] = Vec::LoadAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1] = Vec::LoadAligned(&Vt[i*KDIM1+1][t]);
      Gtrg[i][0][0] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+0][t]);
      Gtrg[i][0][1] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+1][t]);
      Gtrg[i][1][0] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+0][t]);
      Gtrg[i][1][1] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+1][t]);
      Gtrg[i][2][0] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+0][t]);
      Gtrg[i][2][1] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+1][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::sincos_intrin<Vec>(sin_izkR, cos_izkR, izkR[1]);
      sctl::exp_intrin(exp_izkR, izkR[0]);
      // exp(ikr)/r
      Vec G0 = cos_izkR * exp_izkR * Rinv;
      Vec G1 = sin_izkR * exp_izkR * Rinv;
      Vec tmp0 = (1.0) - izkR[0];
      Vec tmp1 = -izkR[1];
      // (1-ikr)*exp(ikr)/r^3
      Vec H0 = (tmp0*G0 - tmp1*G1) * Rinv2;
      Vec H1 = (tmp1*G0 + tmp0*G1) * Rinv2;
      tmp0 = zk[1]*zk[1] - zk[0]*zk[0];
      tmp1 = zk[0]*zk[1]*(-2.0);
      tmp0 = (-3.0)*(Rinv*zk[1]+Rinv2) - tmp0;
      tmp1 = (3.0)*Rinv*zk[0] - tmp1;
      Vec J0 = (G0 * tmp0 - G1 * tmp1)*Rinv2;
      Vec J1 = (G1 * tmp0 + G0 * tmp1)*Rinv2;

      long s_ind = s*nd_*COORD_DIM*KDIM0;
      for (long i = 0; i < nd_; i++) {
        long si_ind = s_ind+i*COORD_DIM*KDIM0;
        Vec Ztmp[KDIM0];
        Vec D0 = dX[0]*Gsrc[si_ind+0] + dX[1]*Gsrc[si_ind+2] + dX[2]*Gsrc[si_ind+4];
        Vec D1 = dX[0]*Gsrc[si_ind+1] + dX[1]*Gsrc[si_ind+3] + dX[2]*Gsrc[si_ind+5];
        Ztmp[0] = J0*D0 - J1*D1;
        Ztmp[1] = J1*D0 + J0*D1;
        Vtrg[i][0] += H0 * D0 - H1 * D1;
        Vtrg[i][1] += H1 * D0 + H0 * D1;
        for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          Gtrg[i][dimi][0] += Ztmp[0] * dX[dimi];
          Gtrg[i][dimi][1] += Ztmp[1] * dX[dimi];
          Gtrg[i][dimi][0] += H0*Gsrc[si_ind+dimi*KDIM0+0] - H1*Gsrc[si_ind+dimi*KDIM0+1];
          Gtrg[i][dimi][1] += H1*Gsrc[si_ind+dimi*KDIM0+0] + H0*Gsrc[si_ind+dimi*KDIM0+1];
        }
      }
    }
    // store
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0].StoreAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1].StoreAligned(&Vt[i*KDIM1+1][t]);
      Gtrg[i][0][0].StoreAligned(&Gt[(0*nd_+i)*KDIM1+0][t]);
      Gtrg[i][0][1].StoreAligned(&Gt[(0*nd_+i)*KDIM1+1][t]);
      Gtrg[i][1][0].StoreAligned(&Gt[(1*nd_+i)*KDIM1+0][t]);
      Gtrg[i][1][1].StoreAligned(&Gt[(1*nd_+i)*KDIM1+1][t]);
      Gtrg[i][2][0].StoreAligned(&Gt[(2*nd_+i)*KDIM1+0][t]);
      Gtrg[i][2][1].StoreAligned(&Gt[(2*nd_+i)*KDIM1+1][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_*KDIM1; j++) {
      pot[i*nd_*KDIM1+j] += Vt[j][i];
    }
    for (long j=0; j < nd_*KDIM1*COORD_DIM; j++) {
      grad[i*nd_*KDIM1*COORD_DIM+j] += Gt[j][i];
    }
  }
}


// charge and dipole, potential
template <class Real, sctl::Integer MaxVecLen=4> void h3ddirectcdp_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;
  static constexpr sctl::Integer KDIM0 = 2;
  static constexpr sctl::Integer KDIM1 = 2;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_*KDIM0, Nsrc);
  //sctl::Matrix<Real> Gs(nd_*KDIM0*COORD_DIM, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_*KDIM1, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_*KDIM0, sctl::Ptr2Itr<Real>((Real*)charge , nd_*KDIM0*Nsrc), false);
    //sctl::Matrix<Real> Gs_(Nsrc, nd_*KDIM0*COORD_DIM, sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*KDIM0*COORD_DIM*Nsrc), false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    //transpose(Gs, Gs_);
    transpose(Xt, Xt_);
    Vt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec zk_[2];
  zk_[0].Load1(zk+0);
  zk_[1].Load1(zk+1);
  Vec thresh2 = thresh[0] * thresh[0];
  // load source, charge and dipvec
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  sctl::Matrix<Real> Vs_(Nsrc, nd_*KDIM0, sctl::Ptr2Itr<Real>((Real*)charge , nd_*KDIM0*Nsrc), false);
  sctl::Matrix<Real> Gs_(Nsrc, nd_*KDIM0*COORD_DIM, sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*KDIM0*COORD_DIM*Nsrc), false);
  //sctl::Vector<Vec> Xsrc(Nsrc*COORD_DIM);
  sctl::Vector<Vec> Vsrc(Nsrc*nd_*KDIM0);
  sctl::Vector<Vec> Gsrc(Nsrc*nd_*COORD_DIM*KDIM0);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    //for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      //Xsrc[s*COORD_DIM+k] = Vec::Load1(&Xs_[s][k]);
    //}
    long s_ind = s*nd_*COORD_DIM*KDIM0;
    for (long i = 0; i < nd_; i++) {
      long si_ind = s_ind+i*COORD_DIM*KDIM0;
      Vsrc[s*nd_*KDIM0+i*KDIM0+0] = Vec::Load1(&Vs_[s][i*KDIM0+0]);
      Vsrc[s*nd_*KDIM0+i*KDIM0+1] = Vec::Load1(&Vs_[s][i*KDIM0+1]);
      Gsrc[si_ind+0] = Vec::Load1(&Gs_[s][(0*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+1] = Vec::Load1(&Gs_[s][(0*nd_+i)*KDIM0+1]);
      Gsrc[si_ind+2] = Vec::Load1(&Gs_[s][(1*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+3] = Vec::Load1(&Gs_[s][(1*nd_+i)*KDIM0+1]);
      Gsrc[si_ind+4] = Vec::Load1(&Gs_[s][(2*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+5] = Vec::Load1(&Gs_[s][(2*nd_+i)*KDIM0+1]);
    }
  }
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential
    Vec Vtrg[nd_][KDIM1];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0] = Vec::LoadAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1] = Vec::LoadAligned(&Vt[i*KDIM1+1][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::sincos_intrin<Vec>(sin_izkR, cos_izkR, izkR[1]);
      sctl::exp_intrin(exp_izkR, izkR[0]);
      // exp(ikr)/r
      Vec G0 = cos_izkR * exp_izkR * Rinv;
      Vec G1 = sin_izkR * exp_izkR * Rinv;
      Vec tmp0 = (1.0)-izkR[0];
      Vec tmp1 = -izkR[1];
      // (1-ikr)*exp(ikr)/r^3
      Vec H0 = (tmp0*G0 - tmp1*G1) * Rinv2;
      Vec H1 = (tmp1*G0 + tmp0*G1) * Rinv2;

      long s_ind = s*nd_*COORD_DIM*KDIM0;
      for (long i = 0; i < nd_; i++) {
        long si_ind = s_ind + i*COORD_DIM*KDIM0;
        Vtrg[i][0] += G0 * Vsrc[s*nd_*KDIM0+i*KDIM0+0] - G1 * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        Vtrg[i][1] += G1 * Vsrc[s*nd_*KDIM0+i*KDIM0+0] + G0 * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        Vec D0 = dX[0]*Gsrc[si_ind+0] + dX[1]*Gsrc[si_ind+2] + dX[2]*Gsrc[si_ind+4];
        Vec D1 = dX[0]*Gsrc[si_ind+1] + dX[1]*Gsrc[si_ind+3] + dX[2]*Gsrc[si_ind+5];
        Vtrg[i][0] += H0 * D0 - H1 * D1;
        Vtrg[i][1] += H1 * D0 + H0 * D1;
      }
    }
    // store
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0].StoreAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1].StoreAligned(&Vt[i*KDIM1+1][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_*KDIM1; j++) {
      pot[i*nd_*KDIM1+j] += Vt[j][i];
    }
  }
}


// charge and dipole, potential and gradient
template <class Real, sctl::Integer MaxVecLen=4> void h3ddirectcdg_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;
  static constexpr sctl::Integer KDIM0 = 2;
  static constexpr sctl::Integer KDIM1 = 2;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_*KDIM0, Nsrc);
  //sctl::Matrix<Real> Gs(nd_*KDIM0*COORD_DIM, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_*KDIM1, Ntrg_);
  sctl::Matrix<Real> Gt(nd_*KDIM1*COORD_DIM, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_*KDIM0, sctl::Ptr2Itr<Real>((Real*)charge , nd_*KDIM0*Nsrc), false);
    //sctl::Matrix<Real> Gs_(Nsrc, nd_*KDIM0*COORD_DIM, sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*KDIM0*COORD_DIM*Nsrc), false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    //transpose(Gs, Gs_);
    transpose(Xt, Xt_);
    Vt = 0;
    Gt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec zk_[2];
  zk_[0].Load1(zk+0);
  zk_[1].Load1(zk+1);
  Vec thresh2 = thresh[0] * thresh[0];
  // load source, charge and dipvec
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  sctl::Matrix<Real> Vs_(Nsrc, nd_*KDIM0, sctl::Ptr2Itr<Real>((Real*)charge , nd_*KDIM0*Nsrc), false);
  sctl::Matrix<Real> Gs_(Nsrc, nd_*KDIM0*COORD_DIM, sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*KDIM0*COORD_DIM*Nsrc), false);
  //sctl::Vector<Vec> Xsrc(Nsrc*COORD_DIM);
  sctl::Vector<Vec> Vsrc(Nsrc*nd_*KDIM0);
  sctl::Vector<Vec> Gsrc(Nsrc*nd_*COORD_DIM*KDIM0);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    //for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      //Xsrc[s*COORD_DIM+k] = Vec::Load1(&Xs_[s][k]);
    //}
    long s_ind = s*nd_*COORD_DIM*KDIM0;
    for (long i = 0; i < nd_; i++) {
      long si_ind = s_ind+i*COORD_DIM*KDIM0;
      Vsrc[s*nd_*KDIM0+i*KDIM0+0] = Vec::Load1(&Vs_[s][i*KDIM0+0]);
      Vsrc[s*nd_*KDIM0+i*KDIM0+1] = Vec::Load1(&Vs_[s][i*KDIM0+1]);
      Gsrc[si_ind+0] = Vec::Load1(&Gs_[s][(0*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+1] = Vec::Load1(&Gs_[s][(0*nd_+i)*KDIM0+1]);
      Gsrc[si_ind+2] = Vec::Load1(&Gs_[s][(1*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+3] = Vec::Load1(&Gs_[s][(1*nd_+i)*KDIM0+1]);
      Gsrc[si_ind+4] = Vec::Load1(&Gs_[s][(2*nd_+i)*KDIM0+0]);
      Gsrc[si_ind+5] = Vec::Load1(&Gs_[s][(2*nd_+i)*KDIM0+1]);
    }
  }
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential and gradient
    Vec Vtrg[nd_][KDIM1], Gtrg[nd_][COORD_DIM][KDIM1];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0] = Vec::LoadAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1] = Vec::LoadAligned(&Vt[i*KDIM1+1][t]);
      Gtrg[i][0][0] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+0][t]);
      Gtrg[i][0][1] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+1][t]);
      Gtrg[i][1][0] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+0][t]);
      Gtrg[i][1][1] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+1][t]);
      Gtrg[i][2][0] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+0][t]);
      Gtrg[i][2][1] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+1][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::sincos_intrin<Vec>(sin_izkR, cos_izkR, izkR[1]);
      sctl::exp_intrin(exp_izkR, izkR[0]);
      // exp(ikr)/r
      Vec G0 = cos_izkR * exp_izkR * Rinv;
      Vec G1 = sin_izkR * exp_izkR * Rinv;
      Vec tmp0 = (1.0)-izkR[0];
      Vec tmp1 = -izkR[1];
      // (1-ikr)*exp(ikr)/r^3
      Vec H0 = (tmp0*G0 - tmp1*G1) * Rinv2;
      Vec H1 = (tmp1*G0 + tmp0*G1) * Rinv2;
      tmp0 = zk[1]*zk[1] - zk[0]*zk[0];
      tmp1 = zk[0]*zk[1]*(-2.0);
      tmp0 = (-3.0)*(Rinv*zk[1]+Rinv2) - tmp0;
      tmp1 = (3.0)*Rinv*zk[0] - tmp1;
      Vec J0 = (G0 * tmp0 - G1 * tmp1)*Rinv2;
      Vec J1 = (G1 * tmp0 + G0 * tmp1)*Rinv2;
      // (ikr-1)*exp(ikr)/r^3 * (xt - xs)
      Vec Zctmp[COORD_DIM][KDIM0];
      for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          Zctmp[dimi][0] = -H0*dX[dimi];
          Zctmp[dimi][1] = -H1*dX[dimi];
      }

      Vec D0, D1;
      long s_ind = s*nd_*COORD_DIM*KDIM0;
      for (long i = 0; i < nd_; i++) {
        long si_ind = s_ind + i*COORD_DIM*KDIM0;

        // temp values
        D0 = dX[0]*Gsrc[si_ind+0] + dX[1]*Gsrc[si_ind+2] + dX[2]*Gsrc[si_ind+4];
        D1 = dX[0]*Gsrc[si_ind+1] + dX[1]*Gsrc[si_ind+3] + dX[2]*Gsrc[si_ind+5];
        tmp0 = J0*D0 - J1*D1;
        tmp1 = J1*D0 + J0*D1;

        // potential
        Vtrg[i][0] += G0 * Vsrc[s*nd_*KDIM0+i*KDIM0+0] - G1 * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        Vtrg[i][1] += G1 * Vsrc[s*nd_*KDIM0+i*KDIM0+0] + G0 * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        Vtrg[i][0] += H0 * D0 - H1 * D1;
        Vtrg[i][1] += H1 * D0 + H0 * D1;
        
        // gradient
        for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          Gtrg[i][dimi][0] += Zctmp[dimi][0] * Vsrc[s*nd_*KDIM0+i*KDIM0+0] - Zctmp[dimi][1] * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
          Gtrg[i][dimi][1] += Zctmp[dimi][1] * Vsrc[s*nd_*KDIM0+i*KDIM0+0] + Zctmp[dimi][0] * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
          
          Gtrg[i][dimi][0] += tmp0 * dX[dimi];
          Gtrg[i][dimi][1] += tmp1 * dX[dimi];

          Gtrg[i][dimi][0] += H0*Gsrc[si_ind+dimi*KDIM0+0] - H1*Gsrc[si_ind+dimi*KDIM0+1];
          Gtrg[i][dimi][1] += H1*Gsrc[si_ind+dimi*KDIM0+0] + H0*Gsrc[si_ind+dimi*KDIM0+1];
        }
      }
    }
    // store
    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0].StoreAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1].StoreAligned(&Vt[i*KDIM1+1][t]);
      Gtrg[i][0][0].StoreAligned(&Gt[(0*nd_+i)*KDIM1+0][t]);
      Gtrg[i][0][1].StoreAligned(&Gt[(0*nd_+i)*KDIM1+1][t]);
      Gtrg[i][1][0].StoreAligned(&Gt[(1*nd_+i)*KDIM1+0][t]);
      Gtrg[i][1][1].StoreAligned(&Gt[(1*nd_+i)*KDIM1+1][t]);
      Gtrg[i][2][0].StoreAligned(&Gt[(2*nd_+i)*KDIM1+0][t]);
      Gtrg[i][2][1].StoreAligned(&Gt[(2*nd_+i)*KDIM1+1][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_*KDIM1; j++) {
      pot[i*nd_*KDIM1+j] += Vt[j][i];
    }
    for (long j=0; j < nd_*KDIM1*COORD_DIM; j++) {
      grad[i*nd_*KDIM1*COORD_DIM+j] += Gt[j][i];
    }
  }
}


// Laplace kernels
// charge, potential 
template <class Real, sctl::Integer MaxVecLen=4> void l3ddirectcp_vec_cpp(const int32_t* nd, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_,       sctl::Ptr2Itr<Real>((Real*)charge , nd_*Nsrc),       false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    transpose(Xt, Xt_);
    Vt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec thresh2 = thresh[0] * thresh[0];
  // load charge
  sctl::Matrix<Real> Vs_(Nsrc, nd_,       sctl::Ptr2Itr<Real>((Real*)charge , nd_*Nsrc),       false);
  //Vec Vsrc[Nsrc][nd_];
  sctl::Vector<Vec> Vsrc(Nsrc*nd_);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (long i = 0; i < nd_; i++) {
      Vsrc[s*nd_+i] = Vec::Load1(&Vs_[s][i]);
    }
  }
  // load source
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  /*
  Vec Xss[Nsrc][COORD_DIM];
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xss[s][k] = Vec::Load1(&Xs_[s][k]);
    }
  }
  */
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential
    Vec Vtrg[nd_];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i] = Vec::LoadAligned(&Vt[i][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      for (long i = 0; i < nd_; i++) {
        Vtrg[i] += Vsrc[s*nd_+i]*Rinv;
      }
    }
    for (long i = 0; i < nd_; i++) {
      Vtrg[i].StoreAligned(&Vt[i][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_; j++) {
      pot[i*nd_+j] += Vt[j][i];
    }
  }
}


// charge, potential and gradient
template <class Real, sctl::Integer MaxVecLen=4> void l3ddirectcg_vec_cpp(const int32_t* nd, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_, Ntrg_);
  sctl::Matrix<Real> Gt(nd_*COORD_DIM, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_,       sctl::Ptr2Itr<Real>((Real*)charge , nd_*Nsrc),       false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    transpose(Xt, Xt_);
    Vt = 0;
    Gt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec thresh2 = thresh[0] * thresh[0];
  // load charge
  sctl::Matrix<Real> Vs_(Nsrc, nd_,       sctl::Ptr2Itr<Real>((Real*)charge , nd_*Nsrc),       false);
  //Vec Vsrc[Nsrc][nd_];
  sctl::Vector<Vec> Vsrc(Nsrc*nd_);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (long i = 0; i < nd_; i++) {
      Vsrc[s*nd_+i] = Vec::Load1(&Vs_[s][i]);
    }
  }
  // load source
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  /*
  Vec Xss[Nsrc][COORD_DIM];
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xss[s][k] = Vec::Load1(&Xs_[s][k]);
    }
  }
  */
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential and gradient
    Vec Vtrg[nd_], Gtrg[nd_][COORD_DIM];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i] = Vec::LoadAligned(&Vt[i][t]);
      Gtrg[i][0] = Vec::LoadAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1] = Vec::LoadAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2] = Vec::LoadAligned(&Gt[2*nd_+i][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec nRinv3 = -Rinv*Rinv*Rinv;
      Vec ztmp[COORD_DIM] = {nRinv3*dX[0], nRinv3*dX[1], nRinv3*dX[2]};
      for (long i = 0; i < nd_; i++) {
        Vtrg[i] += Vsrc[s*nd_+i]*Rinv;
        Gtrg[i][0] += Vsrc[s*nd_+i]*ztmp[0];
        Gtrg[i][1] += Vsrc[s*nd_+i]*ztmp[1];
        Gtrg[i][2] += Vsrc[s*nd_+i]*ztmp[2];
      }
    }
    for (long i = 0; i < nd_; i++) {
      Vtrg[i].StoreAligned(&Vt[i][t]);
      Gtrg[i][0].StoreAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1].StoreAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2].StoreAligned(&Gt[2*nd_+i][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_; j++) {
      pot[i*nd_+j] += Vt[j][i];
    }
    for (long j=0; j < nd_*COORD_DIM; j++) {
      grad[i*nd_*COORD_DIM+j] += Gt[j][i];
    }
  }
}


// dipole, potential
template <class Real, sctl::Integer MaxVecLen=4> void l3ddirectdp_vec_cpp(const int32_t* nd, const Real* sources, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_,       sctl::Ptr2Itr<Real>((Real*)charge , nd_*Nsrc),       false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    transpose(Xt, Xt_);
    Vt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec thresh2 = thresh[0] * thresh[0];
  // load dipole
  sctl::Matrix<Real> Gs_(Nsrc, nd_*COORD_DIM,       sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*Nsrc),       false);
  //Vec Gsrc[Nsrc][nd_][COORD_DIM];
  sctl::Vector<Vec> Gsrc(Nsrc*nd_*COORD_DIM);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (long i = 0; i < nd_; i++) {
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] = Vec::Load1(&Gs_[s][0*nd_+i]);
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] = Vec::Load1(&Gs_[s][1*nd_+i]);
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2] = Vec::Load1(&Gs_[s][2*nd_+i]);
    }
  }
  // load source
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  /*
  Vec Xss[Nsrc][COORD_DIM];
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xss[s][k] = Vec::Load1(&Xs_[s][k]);
    }
  }
  */
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential
    Vec Vtrg[nd_];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i] = Vec::LoadAligned(&Vt[i][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec Rinv3 = Rinv * Rinv * Rinv;
      // TODO: test move Dprod out, faster?
      for (long i = 0; i < nd_; i++) {
        Vec Dprod = dX[0]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] + dX[1]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] + dX[2]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2];
        Vtrg[i] += Dprod*Rinv3;
      }
    }
    for (long i = 0; i < nd_; i++) {
      Vtrg[i].StoreAligned(&Vt[i][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_; j++) {
      pot[i*nd_+j] += Vt[j][i];
    }
  }
}


// dipole, potential and gradient
template <class Real, sctl::Integer MaxVecLen=4> void l3ddirectdg_vec_cpp(const int32_t* nd, const Real* sources, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_, Ntrg_);
  sctl::Matrix<Real> Gt(nd_*COORD_DIM, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_,       sctl::Ptr2Itr<Real>((Real*)charge , nd_*Nsrc),       false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    transpose(Xt, Xt_);
    Vt = 0;
    Gt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec thresh2 = thresh[0] * thresh[0];
  // load dipole
  sctl::Matrix<Real> Gs_(Nsrc, nd_*COORD_DIM,       sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*Nsrc),       false);
  //Vec Gsrc[Nsrc][nd_][COORD_DIM];
  sctl::Vector<Vec> Gsrc(Nsrc*nd_*COORD_DIM);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (long i = 0; i < nd_; i++) {
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] = Vec::Load1(&Gs_[s][0*nd_+i]);
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] = Vec::Load1(&Gs_[s][1*nd_+i]);
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2] = Vec::Load1(&Gs_[s][2*nd_+i]);
    }
  }
  // load source
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  /*
  Vec Xss[Nsrc][COORD_DIM];
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xss[s][k] = Vec::Load1(&Xs_[s][k]);
    }
  }
  */
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential and gradient
    Vec Vtrg[nd_], Gtrg[nd_][COORD_DIM];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i] = Vec::LoadAligned(&Vt[i][t]);
      Gtrg[i][0] = Vec::LoadAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1] = Vec::LoadAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2] = Vec::LoadAligned(&Gt[2*nd_+i][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec Rinv2 = Rinv * Rinv;
      Vec Rinv3 = Rinv * Rinv2;
      Vec Rinv5 =  -3.0*Rinv2*Rinv3;
      for (long i = 0; i < nd_; i++) {
        Vec Dprod = dX[0]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] + dX[1]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] + dX[2]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2];
        Vec RinvDprod = Rinv5*Dprod;
        Vtrg[i] += Dprod*Rinv3;
        Gtrg[i][0] += RinvDprod*dX[0] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0];
        Gtrg[i][1] += RinvDprod*dX[1] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1];
        Gtrg[i][2] += RinvDprod*dX[2] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2];
      }
    }
    for (long i = 0; i < nd_; i++) {
      Vtrg[i].StoreAligned(&Vt[i][t]);
      Gtrg[i][0].StoreAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1].StoreAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2].StoreAligned(&Gt[2*nd_+i][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_; j++) {
      pot[i*nd_+j] += Vt[j][i];
    }
    for (long j=0; j < nd_*COORD_DIM; j++) {
      grad[i*nd_*COORD_DIM+j] += Gt[j][i];
    }
  }
}


// charge and dipole, potential
template <class Real, sctl::Integer MaxVecLen=4> void l3ddirectcdp_vec_cpp(const int32_t* nd, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_,       sctl::Ptr2Itr<Real>((Real*)charge , nd_*Nsrc),       false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    transpose(Xt, Xt_);
    Vt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec thresh2 = thresh[0] * thresh[0];
  // load charge and dipole
  sctl::Matrix<Real> Vs_(Nsrc, nd_,           sctl::Ptr2Itr<Real>((Real*)charge , nd_*Nsrc),       false);
  sctl::Matrix<Real> Gs_(Nsrc, nd_*COORD_DIM, sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*Nsrc),       false);
  //Vec Vsrc[Nsrc][nd_];
  //Vec Gsrc[Nsrc][nd_][COORD_DIM];
  sctl::Vector<Vec> Vsrc(Nsrc*nd_);
  sctl::Vector<Vec> Gsrc(Nsrc*nd_*COORD_DIM);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (long i = 0; i < nd_; i++) {
      Vsrc[s*nd_+i] = Vec::Load1(&Vs_[s][i]);
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] = Vec::Load1(&Gs_[s][0*nd_+i]);
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] = Vec::Load1(&Gs_[s][1*nd_+i]);
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2] = Vec::Load1(&Gs_[s][2*nd_+i]);
    }
  }
  // load source
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  /*
  Vec Xss[Nsrc][COORD_DIM];
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xss[s][k] = Vec::Load1(&Xs_[s][k]);
    }
  }
  */
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential
    Vec Vtrg[nd_];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i] = Vec::LoadAligned(&Vt[i][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec Rinv3 = Rinv * Rinv * Rinv;
      for (long i = 0; i < nd_; i++) {
        Vec Dprod = dX[0]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] + dX[1]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] + dX[2]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2];
        Vtrg[i] += Dprod*Rinv3;
        Vtrg[i] += Vsrc[s*nd_+i]*Rinv;
      }
    }
    for (long i = 0; i < nd_; i++) {
      Vtrg[i].StoreAligned(&Vt[i][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_; j++) {
      pot[i*nd_+j] += Vt[j][i];
    }
  }
}


// charge and dipole, potential and gradient
template <class Real, sctl::Integer MaxVecLen=4> void l3ddirectcdg_vec_cpp(const int32_t* nd, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
  static constexpr sctl::Integer COORD_DIM = 3;

  sctl::Long nd_ = nd[0];
  sctl::Long Nsrc = ns[0];
  sctl::Long Ntrg = nt[0];
  sctl::Long Ntrg_ = ((Ntrg+MaxVecLen-1)/MaxVecLen)*MaxVecLen;

  //sctl::Matrix<Real> Xs(COORD_DIM, Nsrc);
  //sctl::Matrix<Real> Vs(nd_, Nsrc);
  sctl::Matrix<Real> Xt(COORD_DIM, Ntrg_);
  sctl::Matrix<Real> Vt(nd_, Ntrg_);
  sctl::Matrix<Real> Gt(nd_*COORD_DIM, Ntrg_);
  { // Set Xs, Vs, Xt, Vt
    auto transpose = [](sctl::Matrix<Real>& A, const sctl::Matrix<Real>& B) {
      sctl::Long d0 = std::min(A.Dim(0), B.Dim(1));
      sctl::Long d1 = std::min(A.Dim(1), B.Dim(0));
      for (long i = 0; i < d0; i++) {
        for (long j = 0; j < d1; j++) {
          A[i][j] = B[j][i];
        }
      }
    };
    //sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
    //sctl::Matrix<Real> Vs_(Nsrc, nd_,       sctl::Ptr2Itr<Real>((Real*)charge , nd_*Nsrc),       false);
    sctl::Matrix<Real> Xt_(Ntrg, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)ztarg  , COORD_DIM*Ntrg), false);
    //transpose(Xs, Xs_);
    //transpose(Vs, Vs_);
    transpose(Xt, Xt_);
    Vt = 0;
    Gt = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec thresh2 = thresh[0] * thresh[0];
  // load charge and dipole
  sctl::Matrix<Real> Vs_(Nsrc, nd_,           sctl::Ptr2Itr<Real>((Real*)charge , nd_*Nsrc),       false);
  sctl::Matrix<Real> Gs_(Nsrc, nd_*COORD_DIM, sctl::Ptr2Itr<Real>((Real*)dipvec , nd_*Nsrc),       false);
  //Vec Vsrc[Nsrc][nd_];
  //Vec Gsrc[Nsrc][nd_][COORD_DIM];
  sctl::Vector<Vec> Vsrc(Nsrc*nd_);
  sctl::Vector<Vec> Gsrc(Nsrc*nd_*COORD_DIM);
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (long i = 0; i < nd_; i++) {
      Vsrc[s*nd_+i] = Vec::Load1(&Vs_[s][i]);
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] = Vec::Load1(&Gs_[s][0*nd_+i]);
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] = Vec::Load1(&Gs_[s][1*nd_+i]);
      Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2] = Vec::Load1(&Gs_[s][2*nd_+i]);
    }
  }
  // load source
  sctl::Matrix<Real> Xs_(Nsrc, COORD_DIM, sctl::Ptr2Itr<Real>((Real*)sources, COORD_DIM*Nsrc), false);
  /*
  Vec Xss[Nsrc][COORD_DIM];
  for (sctl::Long s = 0; s < Nsrc; s++) {
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xss[s][k] = Vec::Load1(&Xs_[s][k]);
    }
  }
  */
  #pragma omp parallel for schedule(static)
  for (sctl::Long t = 0; t < Ntrg_; t += VecLen) {
    Vec Xtrg[COORD_DIM];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      Xtrg[k] = Vec::LoadAligned(&Xt[k][t]);
    }
    // load potential and gradient
    Vec Vtrg[nd_], Gtrg[nd_][COORD_DIM];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i] = Vec::LoadAligned(&Vt[i][t]);
      Gtrg[i][0] = Vec::LoadAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1] = Vec::LoadAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2] = Vec::LoadAligned(&Gt[2*nd_+i][t]);
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = approx_rsqrt(R2);
      if (sizeof(Real) <= 2) {
      } else if (sizeof(Real) <= 4) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv) * 0.5; // 7 - cycles
      } else if (sizeof(Real) <= 8) {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<0>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      } else {
        Rinv *= ((3.0) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<0>(3)*3-1>(2.0)) - R2 * Rinv * Rinv); // 7 - cycles
        Rinv *= ((3.0 * sctl::pow<sctl::pow<1>(3)*3-1>(2.0)) - R2 * Rinv * Rinv) * (sctl::pow<(sctl::pow<1>(3)*3-1)*3/2+1>(0.5)); // 8 - cycles
      }
      Rinv &= (R2 > thresh2);

      Vec Rinv2 = Rinv * Rinv;
      Vec Rinv3 = Rinv * Rinv2;
      Vec nRinv3 = -Rinv*Rinv*Rinv;
      Vec Rinv5 =  -3.0*Rinv2*Rinv3;
      Vec ztmp[COORD_DIM] = {nRinv3*dX[0], nRinv3*dX[1], nRinv3*dX[2]};
      for (long i = 0; i < nd_; i++) {
        Vec Dprod = dX[0]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] + dX[1]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] + dX[2]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2];
        Vtrg[i] += Dprod*Rinv3;
        Vtrg[i] += Vsrc[s*nd_+i]*Rinv;
        Vec RinvDprod = Rinv5*Dprod;
        Gtrg[i][0] += RinvDprod*dX[0] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] + Vsrc[s*nd_+i]*ztmp[0];
        Gtrg[i][1] += RinvDprod*dX[1] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] + Vsrc[s*nd_+i]*ztmp[1];
        Gtrg[i][2] += RinvDprod*dX[2] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2] + Vsrc[s*nd_+i]*ztmp[2];
      }
    }
    for (long i = 0; i < nd_; i++) {
      Vtrg[i].StoreAligned(&Vt[i][t]);
      Gtrg[i][0].StoreAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1].StoreAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2].StoreAligned(&Gt[2*nd_+i][t]);
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_; j++) {
      pot[i*nd_+j] += Vt[j][i];
    }
    for (long j=0; j < nd_*COORD_DIM; j++) {
      grad[i*nd_*COORD_DIM+j] += Gt[j][i];
    }
  }
}


#endif //_VEC_KERNELS_HPP_
