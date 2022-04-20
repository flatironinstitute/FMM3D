#ifndef _VEC_KERNELS_HPP_
#define _VEC_KERNELS_HPP_

#define NDEBUG
#include <sctl.hpp>

template <class Real, class VecType, sctl::Integer DIM, sctl::Integer KDIM0, sctl::Integer KDIM1, sctl::Integer SCDIM, class uKernel, sctl::Integer digits> struct uKerHelper {
  template <class CtxType> static inline void Eval(VecType* vt, const VecType (&dX)[DIM], const Real* vs, const sctl::Integer nd, const CtxType& ctx) {
    VecType M[KDIM0][KDIM1][SCDIM];
    uKernel::template uKerMatrix<digits>(M, dX, ctx);
    for (sctl::Integer i = 0; i < nd; i++) {
      const Real* vs_ = vs+i*SCDIM;
      for (sctl::Integer k1 = 0; k1 < KDIM1; k1++) {
        VecType* vt_ = vt+(k1*nd+i)*SCDIM;
        for (sctl::Integer k0 = 0; k0 < KDIM0; k0++) {
          const VecType vs0(vs_[(k0*nd)*SCDIM+0]);
          vt_[0] = FMA(M[k0][k1][0], vs0, vt_[0]);
          if (SCDIM == 2) {
            const VecType vs1(vs_[(k0*nd)*SCDIM+1]);
            vt_[0] = FMA(M[k0][k1][1],-vs1, vt_[0]);
            vt_[1] = FMA(M[k0][k1][1], vs0, vt_[1]);
            vt_[1] = FMA(M[k0][k1][0], vs1, vt_[1]);
          }
        }
      }
    }
  }
  template <sctl::Integer nd, class CtxType> static inline void EvalND(VecType* vt, const VecType (&dX)[DIM], const Real* vs, const CtxType& ctx) {
    VecType M[KDIM0][KDIM1][SCDIM];
    uKernel::template uKerMatrix<digits>(M, dX, ctx);
    for (sctl::Integer i = 0; i < nd; i++) {
      const Real* vs_ = vs+i*SCDIM;
      for (sctl::Integer k1 = 0; k1 < KDIM1; k1++) {
        VecType* vt_ = vt+(k1*nd+i)*SCDIM;
        for (sctl::Integer k0 = 0; k0 < KDIM0; k0++) {
          const VecType vs0(vs_[(k0*nd)*SCDIM+0]);
          vt_[0] = FMA(M[k0][k1][0], vs0, vt_[0]);
          if (SCDIM == 2) {
            const VecType vs1(vs_[(k0*nd)*SCDIM+1]);
            vt_[0] = FMA(M[k0][k1][1],-vs1, vt_[0]);
            vt_[1] = FMA(M[k0][k1][1], vs0, vt_[1]);
            vt_[1] = FMA(M[k0][k1][0], vs1, vt_[1]);
          }
        }
      }
    }
  }
};

template <class uKernel> class GenericKernel : public uKernel {
    static constexpr sctl::Integer VecLen = uKernel::VecLen;
    using VecType = typename uKernel::VecType;
    using Real = typename uKernel::RealType;

    template <sctl::Integer K0, sctl::Integer K1, sctl::Integer Q, sctl::Integer D, class ...T> static constexpr sctl::Integer get_DIM  (void (*uKer)(VecType (&M)[K0][K1][Q], const VecType (&r)[D], T... args)) { return D;  }
    template <sctl::Integer K0, sctl::Integer K1, sctl::Integer Q, sctl::Integer D, class ...T> static constexpr sctl::Integer get_SCDIM(void (*uKer)(VecType (&M)[K0][K1][Q], const VecType (&r)[D], T... args)) { return Q;  }
    template <sctl::Integer K0, sctl::Integer K1, sctl::Integer Q, sctl::Integer D, class ...T> static constexpr sctl::Integer get_KDIM0(void (*uKer)(VecType (&M)[K0][K1][Q], const VecType (&r)[D], T... args)) { return K0; }
    template <sctl::Integer K0, sctl::Integer K1, sctl::Integer Q, sctl::Integer D, class ...T> static constexpr sctl::Integer get_KDIM1(void (*uKer)(VecType (&M)[K0][K1][Q], const VecType (&r)[D], T... args)) { return K1; }

    static constexpr sctl::Integer DIM   = get_DIM  (uKernel::template uKerMatrix<0,GenericKernel>);
    static constexpr sctl::Integer SCDIM = get_SCDIM(uKernel::template uKerMatrix<0,GenericKernel>);
    static constexpr sctl::Integer KDIM0 = get_KDIM0(uKernel::template uKerMatrix<0,GenericKernel>);
    static constexpr sctl::Integer KDIM1 = get_KDIM1(uKernel::template uKerMatrix<0,GenericKernel>);

  public:

    GenericKernel() : ctx_ptr(this) {}

    static constexpr sctl::Integer CoordDim() {
      return DIM;
    }
    static constexpr sctl::Integer SrcDim() {
      return KDIM0*SCDIM;
    }
    static constexpr sctl::Integer TrgDim() {
      return KDIM1*SCDIM;
    }

    template <bool enable_openmp=false, sctl::Integer digits=-1> void Eval(sctl::Vector<sctl::Vector<Real>>& v_trg_, const sctl::Vector<Real>& r_trg, const sctl::Vector<Real>& r_src, const sctl::Vector<sctl::Vector<Real>>& v_src_, const sctl::Integer nd) const {
      if      (nd == 1) EvalHelper<enable_openmp, digits, 1>(v_trg_, r_trg, r_src, v_src_, nd);
      else if (nd == 2) EvalHelper<enable_openmp, digits, 2>(v_trg_, r_trg, r_src, v_src_, nd);
      else if (nd == 3) EvalHelper<enable_openmp, digits, 3>(v_trg_, r_trg, r_src, v_src_, nd);
      else if (nd == 4) EvalHelper<enable_openmp, digits, 4>(v_trg_, r_trg, r_src, v_src_, nd);
      else if (nd == 5) EvalHelper<enable_openmp, digits, 5>(v_trg_, r_trg, r_src, v_src_, nd);
      else if (nd == 6) EvalHelper<enable_openmp, digits, 6>(v_trg_, r_trg, r_src, v_src_, nd);
      else if (nd == 7) EvalHelper<enable_openmp, digits, 7>(v_trg_, r_trg, r_src, v_src_, nd);
      else if (nd == 8) EvalHelper<enable_openmp, digits, 8>(v_trg_, r_trg, r_src, v_src_, nd);
      else  EvalHelper<enable_openmp, digits, 0>(v_trg_, r_trg, r_src, v_src_, nd);
    }

  private:

    template <bool enable_openmp=false, sctl::Integer digits=-1, sctl::Integer ND=0> void EvalHelper(sctl::Vector<sctl::Vector<Real>>& v_trg_, const sctl::Vector<Real>& r_trg, const sctl::Vector<Real>& r_src, const sctl::Vector<sctl::Vector<Real>>& v_src_, const sctl::Integer nd) const {
      static constexpr sctl::Integer digits_ = (digits==-1 ? (sctl::Integer)(sctl::TypeTraits<Real>::SigBits*0.3010299957) : digits);
      auto uKerEval = [this](VecType* vt, const VecType (&dX)[DIM], const Real* vs, const sctl::Integer nd) {
        if (ND > 0) uKerHelper<Real,VecType,DIM,KDIM0,KDIM1,SCDIM,uKernel,digits_>::template EvalND<ND>(vt, dX, vs, *this);
        else uKerHelper<Real,VecType,DIM,KDIM0,KDIM1,SCDIM,uKernel,digits_>::Eval(vt, dX, vs, nd, *this);
      };

      const sctl::Long Ns = r_src.Dim() / DIM;
      const sctl::Long Nt = r_trg.Dim() / DIM;
      SCTL_ASSERT(r_trg.Dim() == Nt*DIM);
      SCTL_ASSERT(r_src.Dim() == Ns*DIM);

      sctl::Vector<sctl::Long> src_cnt(v_src_.Dim()), src_dsp(v_src_.Dim()); src_dsp = 0;
      sctl::Vector<sctl::Long> trg_cnt(v_trg_.Dim()), trg_dsp(v_trg_.Dim()); trg_dsp = 0;
      for (sctl::Integer i = 0; i < trg_cnt.Dim(); i++) {
        trg_cnt[i] = v_trg_[i].Dim()/Nt;
        trg_dsp[i] = (i ? trg_dsp[i-1]+trg_cnt[i-1] : 0);
      }
      for (sctl::Integer i = 0; i < src_cnt.Dim(); i++) {
        src_cnt[i] = v_src_[i].Dim()/Ns;
        src_dsp[i] = (i ? src_dsp[i-1]+src_cnt[i-1] : 0);
      }
      SCTL_ASSERT(src_cnt[src_cnt.Dim()-1] + src_dsp[src_dsp.Dim()-1] == SrcDim()*nd);
      SCTL_ASSERT(trg_cnt[trg_cnt.Dim()-1] + trg_dsp[trg_dsp.Dim()-1] == TrgDim()*nd);

      sctl::Vector<Real> v_src(Ns*SrcDim()*nd);
      for (sctl::Integer j = 0; j < src_cnt.Dim(); j++) {
        const sctl::Integer src_cnt_ = src_cnt[j];
        const sctl::Integer src_dsp_ = src_dsp[j];
        for (sctl::Integer k = 0; k < src_cnt_; k++) {
          for (sctl::Long i = 0; i < Ns; i++) {
            v_src[i*SrcDim()*nd+src_dsp_+k] = v_src_[j][i*src_cnt_+k];
          }
        }
      }

      const sctl::Long NNt = ((Nt + VecLen - 1) / VecLen) * VecLen;
      //if (NNt == VecLen) {
      //  VecType xt[DIM], vt[KDIM1], xs[DIM];
      //  Real vs[KDIM0];
      //  for (sctl::Integer k = 0; k < KDIM1; k++) vt[k] = VecType::Zero();
      //  for (sctl::Integer k = 0; k < DIM; k++) {
      //    alignas(sizeof(VecType)) sctl::StaticArray<Real,VecLen> Xt;
      //    VecType::Zero().StoreAligned(&Xt[0]);
      //    for (sctl::Integer i = 0; i < Nt; i++) Xt[i] = r_trg[i*DIM+k];
      //    xt[k] = VecType::LoadAligned(&Xt[0]);
      //  }
      //  for (sctl::Long s = 0; s < Ns; s++) {
      //    for (sctl::Integer k = 0; k < DIM; k++) xs[k] = VecType::Load1(&r_src[s*DIM+k]);
      //    for (sctl::Integer k = 0; k < KDIM0; k++) vs[k] = v_src[s*KDIM0+k];
      //    uKerEval(vt, xt, xs, vs, nd);
      //  }
      //  for (sctl::Integer k = 0; k < KDIM1; k++) {
      //    alignas(sizeof(VecType)) sctl::StaticArray<Real,VecLen> out;
      //    vt[k].StoreAligned(&out[0]);
      //    for (sctl::Long t = 0; t < Nt; t++) {
      //      v_trg[t*KDIM1+k] += out[t] * uKernel::uKerScaleFactor();
      //    }
      //  }
      //} else
      {
        const sctl::Matrix<Real> Xs_(Ns, DIM, (sctl::Iterator<Real>)r_src.begin(), false);
        const sctl::Matrix<Real> Vs_(Ns, SrcDim()*nd, (sctl::Iterator<Real>)v_src.begin(), false);

        sctl::Matrix<Real> Xt_(DIM, NNt), Vt_(TrgDim()*nd, NNt);
        for (sctl::Long k = 0; k < DIM; k++) { // Set Xt_
          for (sctl::Long i = 0; i < Nt; i++) {
            Xt_[k][i] = r_trg[i*DIM+k];
          }
          for (sctl::Long i = Nt; i < NNt; i++) {
            Xt_[k][i] = 0;
          }
        }
        if (enable_openmp) { // Compute Vt_
          #pragma omp parallel for schedule(static)
          for (sctl::Long t = 0; t < NNt; t += VecLen) {
            VecType xt[DIM], vt[TrgDim()*nd];
            for (sctl::Integer k = 0; k < TrgDim()*nd; k++) vt[k] = VecType::Zero();
            for (sctl::Integer k = 0; k < DIM; k++) xt[k] = VecType::LoadAligned(&Xt_[k][t]);

            for (sctl::Long s = 0; s < Ns; s++) {
              VecType dX[DIM];
              for (sctl::Integer k = 0; k < DIM; k++) dX[k] = xt[k] - Xs_[s][k];
              uKerEval(vt, dX, &Vs_[s][0], nd);
            }
            for (sctl::Integer k = 0; k < TrgDim()*nd; k++) vt[k].StoreAligned(&Vt_[k][t]);
          }
        } else {
          for (sctl::Long t = 0; t < NNt; t += VecLen) {
            VecType xt[DIM], vt[TrgDim()*nd];
            for (sctl::Integer k = 0; k < TrgDim()*nd; k++) vt[k] = VecType::Zero();
            for (sctl::Integer k = 0; k < DIM; k++) xt[k] = VecType::LoadAligned(&Xt_[k][t]);

            for (sctl::Long s = 0; s < Ns; s++) {
              VecType dX[DIM];
              for (sctl::Integer k = 0; k < DIM; k++) dX[k] = xt[k] - Xs_[s][k];
              uKerEval(vt, dX, &Vs_[s][0], nd);
            }
            for (sctl::Integer k = 0; k < TrgDim()*nd; k++) vt[k].StoreAligned(&Vt_[k][t]);
          }
        }

        for (sctl::Integer j = 0; j < trg_cnt.Dim(); j++) {
          const sctl::Integer trg_cnt_ = trg_cnt[j];
          const sctl::Integer trg_dsp_ = trg_dsp[j];
          for (sctl::Long i = 0; i < Nt; i++) {
            for (sctl::Integer k = 0; k < trg_cnt_; k++) {
              v_trg_[j][i*trg_cnt_+k] += Vt_[trg_dsp_+k][i] * uKernel::uKerScaleFactor();
            }
          }
        }
      }

    }

    void* ctx_ptr;
};


template <class Real, sctl::Integer VecLen_, sctl::Integer chrg, sctl::Integer dipo, sctl::Integer poten, sctl::Integer grad> struct Helmholtz3D {
  static constexpr sctl::Integer VecLen = VecLen_;
  using VecType = sctl::Vec<Real, VecLen>;
  using RealType = Real;

  VecType thresh2;
  VecType zk[2];

  static constexpr Real uKerScaleFactor() {
    return 1;
  }
  template <sctl::Integer digits, class CtxType> static inline void uKerMatrix(VecType (&M)[chrg+dipo][poten+grad][2], const VecType (&dX)[3], const CtxType& ctx) {
    using RealType = typename VecType::ScalarType;
    static constexpr sctl::Integer COORD_DIM = 3;

    const VecType& thresh2 = ctx.thresh2;
    const VecType (&zk)[2] = ctx.zk;

    const VecType R2 = dX[0]*dX[0]+dX[1]*dX[1]+dX[2]*dX[2];
    const VecType Rinv = sctl::approx_rsqrt<digits>(R2, (R2 > thresh2));
    const VecType Rinv2 = Rinv * Rinv;

    const VecType R = R2 * Rinv;
    const VecType izkR[2] = {-zk[1]*R, zk[0]*R};

    VecType sin_izkR, cos_izkR;
    sctl::approx_sincos<digits>(sin_izkR, cos_izkR, izkR[1]);
    const VecType exp_izkR = sctl::approx_exp<digits>(izkR[0]);

    // exp(ikr)/r
    const VecType G0 = cos_izkR * exp_izkR * Rinv;
    const VecType G1 = sin_izkR * exp_izkR * Rinv;

    // (1-ikr)*exp(ikr)/r^3
    const VecType H0 = ((izkR[0]-(RealType)1)*G0 - izkR[1]*G1) * Rinv2;
    const VecType H1 = ((izkR[0]-(RealType)1)*G1 + izkR[1]*G0) * Rinv2;

    const VecType tmp0 = (-3.0)*(Rinv*zk[1]+Rinv2) - zk[1]*zk[1] + zk[0]*zk[0];
    const VecType tmp1 = (3.0)*Rinv*zk[0] - zk[0]*zk[1]*(-2.0);
    const VecType J0 = (G0 * tmp0 - G1 * tmp1)*Rinv2;
    const VecType J1 = (G1 * tmp0 + G0 * tmp1)*Rinv2;

    if (chrg && poten) { // charge potential
      M[0][0][0] =  G0;
      M[0][0][1] =  G1;
    }

    if (chrg && grad) { // charge gradient
      for (sctl::Integer i = 0; i < COORD_DIM; i++){
        M[0][poten+i][0] = H0*dX[i];
        M[0][poten+i][1] = H1*dX[i];
      }
    }

    if (dipo && poten) { // dipole potential
      for (sctl::Integer i = 0; i < COORD_DIM; i++){
        M[chrg+i][0][0] = -H0*dX[i];
        M[chrg+i][0][1] = -H1*dX[i];
      }
    }

    if (dipo && grad) { // dipole gradient
      for (sctl::Integer i = 0; i < COORD_DIM; i++){
        const VecType J0_dXi = J0*dX[i];
        const VecType J1_dXi = J1*dX[i];
        for (sctl::Integer j = 0; j < COORD_DIM; j++){
          M[chrg+i][poten+j][0] = (i==j ? J0_dXi*dX[j]-H0 : J0_dXi*dX[j]);
          M[chrg+i][poten+j][1] = (i==j ? J1_dXi*dX[j]-H1 : J1_dXi*dX[j]);
        }
      }
    }
  }
};

template <class Real, sctl::Integer VecLen, sctl::Integer chrg, sctl::Integer dipo, sctl::Integer poten, sctl::Integer grad> static void EvalHelmholtz(sctl::Vector<sctl::Vector<Real>>& v_trg, const sctl::Vector<Real>& r_trg, const sctl::Vector<Real>& r_src, const sctl::Vector<sctl::Vector<Real>>& v_src, const sctl::Integer nd, const Real* zk, const Real thresh, const sctl::Integer digits) {
  GenericKernel<Helmholtz3D<Real,VecLen, chrg,dipo,poten,grad>> ker;
  ker.thresh2 = thresh*thresh;
  ker.zk[0] = zk[0];
  ker.zk[1] = zk[1];

  if (digits < 0) ker.template Eval<true, -1>(v_trg, r_trg, r_src, v_src, nd);
  else if (digits <= 3) ker.template Eval<true, 3>(v_trg, r_trg, r_src, v_src, nd);
  else if (digits <= 6) ker.template Eval<true, 6>(v_trg, r_trg, r_src, v_src, nd);
  else if (digits <= 9) ker.template Eval<true, 9>(v_trg, r_trg, r_src, v_src, nd);
  else if (digits <=12) ker.template Eval<true,12>(v_trg, r_trg, r_src, v_src, nd);
  else if (digits <=15) ker.template Eval<true,15>(v_trg, r_trg, r_src, v_src, nd);
  else ker.template Eval<true,-1>(v_trg, r_trg, r_src, v_src, nd);
}

template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void h3ddirectcdg_new_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh, const sctl::Integer digits) {
  const sctl::Vector<Real> r_trg(nt[0]*3, (sctl::Iterator<Real>)sctl::Ptr2ConstItr<Real>(ztarg,nt[0]*3), false);
  const sctl::Vector<Real> r_src(ns[0]*3, (sctl::Iterator<Real>)sctl::Ptr2ConstItr<Real>(sources,ns[0]*3), false);

  sctl::Integer offset = 0;
  sctl::Vector<sctl::Vector<Real>> v_src((charge?1:0)+(dipvec?1:0));
  if (charge) {
    v_src[offset].ReInit(ns[0]*nd[0]*2, (sctl::Iterator<Real>)sctl::Ptr2ConstItr<Real>(charge,ns[0]*nd[0]*2), false);
    offset++;
  }
  if (dipvec) {
    v_src[offset].ReInit(ns[0]*3*nd[0]*2, (sctl::Iterator<Real>)sctl::Ptr2ConstItr<Real>(dipvec,ns[0]*3*nd[0]*2), false);
  }

  offset = 0;
  sctl::Vector<sctl::Vector<Real>> v_trg((pot?1:0)+(grad?1:0));
  if (pot) {
    v_trg[offset].ReInit(nt[0]*nd[0]*2, sctl::Ptr2Itr<Real>(pot,nt[0]*nd[0]*2), false);
    offset++;
  }
  if (grad) {
    v_trg[offset].ReInit(nt[0]*3*nd[0]*2, sctl::Ptr2Itr<Real>(grad,nt[0]*3*nd[0]*2), false);
  }

  if ((charge) && (!dipvec) && (pot) && (!grad)) EvalHelmholtz<Real, MaxVecLen, 1,0,1,0>(v_trg, r_trg, r_src, v_src, nd[0], zk, thresh[0], digits);
  else if ((charge) && (!dipvec) && (pot) && (grad)) EvalHelmholtz<Real, MaxVecLen, 1,0,1,3>(v_trg, r_trg, r_src, v_src, nd[0], zk, thresh[0], digits);
  else if ((!charge) && (dipvec) && (pot) && (!grad)) EvalHelmholtz<Real, MaxVecLen, 0,3,1,0>(v_trg, r_trg, r_src, v_src, nd[0], zk, thresh[0], digits);
  else if ((!charge) && (dipvec) && (pot) && (grad)) EvalHelmholtz<Real, MaxVecLen, 0,3,1,3>(v_trg, r_trg, r_src, v_src, nd[0], zk, thresh[0], digits);
  else if ((charge) && (dipvec) && (pot) && (!grad)) EvalHelmholtz<Real, MaxVecLen, 1,3,1,0>(v_trg, r_trg, r_src, v_src, nd[0], zk, thresh[0], digits);
  else if ((charge) && (dipvec) && (pot) && (grad)) EvalHelmholtz<Real, MaxVecLen, 1,3,1,3>(v_trg, r_trg, r_src, v_src, nd[0], zk, thresh[0], digits);
}




template <class Real, sctl::Integer VecLen_, sctl::Integer chrg, sctl::Integer dipo, sctl::Integer poten, sctl::Integer grad> struct Laplace3D {
  static constexpr sctl::Integer VecLen = VecLen_;
  using VecType = sctl::Vec<Real, VecLen>;
  using RealType = Real;

  VecType thresh2;

  static constexpr Real uKerScaleFactor() {
    return 1;
  }
  template <sctl::Integer digits, class CtxType> static inline void uKerMatrix(VecType (&M)[chrg+dipo][poten+grad][1], const VecType (&dX)[3], const CtxType& ctx) {
    using RealType = typename VecType::ScalarType;
    static constexpr sctl::Integer COORD_DIM = 3;

    const VecType& thresh2 = ctx.thresh2;

    const VecType R2 = dX[0]*dX[0]+dX[1]*dX[1]+dX[2]*dX[2];
    const VecType Rinv = sctl::approx_rsqrt<digits>(R2, (R2 > thresh2));
    const VecType Rinv2 = Rinv * Rinv;
    const VecType Rinv3 = Rinv * Rinv2;

    if (chrg && poten) { // charge potential
      M[0][0][0] = Rinv;
    }

    if (chrg && grad) { // charge gradient
      for (sctl::Integer i = 0; i < COORD_DIM; i++){
        M[0][poten+i][0] = -Rinv3*dX[i];
      }
    }

    if (dipo && poten) { // dipole potential
      for (sctl::Integer i = 0; i < COORD_DIM; i++){
        M[chrg+i][0][0] = Rinv3*dX[i];
      }
    }

    if (dipo && grad) { // dipole gradient
      const VecType J0 = Rinv3 * Rinv2 * (RealType)(-3);
      for (sctl::Integer i = 0; i < COORD_DIM; i++){
        const VecType J0_dXi = J0*dX[i];
        for (sctl::Integer j = 0; j < COORD_DIM; j++){
          M[chrg+i][poten+j][0] = (i==j ? J0_dXi*dX[j]+Rinv3 : J0_dXi*dX[j]);
        }
      }
    }
  }
};

template <class Real, sctl::Integer VecLen, sctl::Integer chrg, sctl::Integer dipo, sctl::Integer poten, sctl::Integer grad> static void EvalLaplace(sctl::Vector<sctl::Vector<Real>>& v_trg, const sctl::Vector<Real>& r_trg, const sctl::Vector<Real>& r_src, const sctl::Vector<sctl::Vector<Real>>& v_src, const sctl::Integer nd, const Real thresh, const sctl::Integer digits) {
  GenericKernel<Laplace3D<Real,VecLen, chrg,dipo,poten,grad>> ker;
  ker.thresh2 = thresh*thresh;

  if (digits < 0) ker.template Eval<true, -1>(v_trg, r_trg, r_src, v_src, nd);
  else if (digits <= 3) ker.template Eval<true, 3>(v_trg, r_trg, r_src, v_src, nd);
  else if (digits <= 6) ker.template Eval<true, 6>(v_trg, r_trg, r_src, v_src, nd);
  else if (digits <= 9) ker.template Eval<true, 9>(v_trg, r_trg, r_src, v_src, nd);
  else if (digits <=12) ker.template Eval<true,12>(v_trg, r_trg, r_src, v_src, nd);
  else if (digits <=15) ker.template Eval<true,15>(v_trg, r_trg, r_src, v_src, nd);
  else ker.template Eval<true,-1>(v_trg, r_trg, r_src, v_src, nd);
}

template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void l3ddirectcdg_new_cpp(const int32_t* nd, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh, const sctl::Integer digits) {
  const sctl::Vector<Real> r_trg(nt[0]*3, (sctl::Iterator<Real>)sctl::Ptr2ConstItr<Real>(ztarg,nt[0]*3), false);
  const sctl::Vector<Real> r_src(ns[0]*3, (sctl::Iterator<Real>)sctl::Ptr2ConstItr<Real>(sources,ns[0]*3), false);

  sctl::Integer offset = 0;
  sctl::Vector<sctl::Vector<Real>> v_src((charge?1:0)+(dipvec?1:0));
  if (charge) {
    v_src[offset].ReInit(ns[0]*nd[0], (sctl::Iterator<Real>)sctl::Ptr2ConstItr<Real>(charge,ns[0]*nd[0]), false);
    offset++;
  }
  if (dipvec) {
    v_src[offset].ReInit(ns[0]*3*nd[0], (sctl::Iterator<Real>)sctl::Ptr2ConstItr<Real>(dipvec,ns[0]*3*nd[0]), false);
  }

  offset = 0;
  sctl::Vector<sctl::Vector<Real>> v_trg((pot?1:0)+(grad?1:0));
  if (pot) {
    v_trg[offset].ReInit(nt[0]*nd[0], sctl::Ptr2Itr<Real>(pot,nt[0]*nd[0]), false);
    offset++;
  }
  if (grad) {
    v_trg[offset].ReInit(nt[0]*3*nd[0], sctl::Ptr2Itr<Real>(grad,nt[0]*3*nd[0]), false);
  }

  if ((charge) && (!dipvec) && (pot) && (!grad)) EvalLaplace<Real, MaxVecLen, 1,0,1,0>(v_trg, r_trg, r_src, v_src, nd[0], thresh[0], digits);
  else if ((charge) && (!dipvec) && (pot) && (grad)) EvalLaplace<Real, MaxVecLen, 1,0,1,3>(v_trg, r_trg, r_src, v_src, nd[0], thresh[0], digits);
  else if ((!charge) && (dipvec) && (pot) && (!grad)) EvalLaplace<Real, MaxVecLen, 0,3,1,0>(v_trg, r_trg, r_src, v_src, nd[0], thresh[0], digits);
  else if ((!charge) && (dipvec) && (pot) && (grad)) EvalLaplace<Real, MaxVecLen, 0,3,1,3>(v_trg, r_trg, r_src, v_src, nd[0], thresh[0], digits);
  else if ((charge) && (dipvec) && (pot) && (!grad)) EvalLaplace<Real, MaxVecLen, 1,3,1,0>(v_trg, r_trg, r_src, v_src, nd[0], thresh[0], digits);
  else if ((charge) && (dipvec) && (pot) && (grad)) EvalLaplace<Real, MaxVecLen, 1,3,1,3>(v_trg, r_trg, r_src, v_src, nd[0], thresh[0], digits);
}


















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
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void h3ddirectcp_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
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
  zk_[0] = Vec::Load1(zk+0);
  zk_[1] = Vec::Load1(zk+1);
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec R = R2 * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::approx_sincos<-1>(sin_izkR, cos_izkR, izkR[1]);
      exp_izkR = sctl::approx_exp<-1>(izkR[0]);
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
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void h3ddirectcg_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
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
  zk_[0] = Vec::Load1(zk+0);
  zk_[1] = Vec::Load1(zk+1);
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::approx_sincos<-1>(sin_izkR, cos_izkR, izkR[1]);
      exp_izkR = sctl::approx_exp<-1>(izkR[0]);
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


// charge, potential, gradient and hessian
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void h3ddirectch_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, Real* hess, const Real* thresh) {
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
  sctl::Matrix<Real> Ht(nd_*KDIM1*COORD_DIM*2, Ntrg_);
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
    Ht = 0;
  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec zk_[2];
  zk_[0] = Vec::Load1(zk+0);
  zk_[1] = Vec::Load1(zk+1);

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
    // load potential and gradient and hessian
    Vec Vtrg[nd_][KDIM1], Gtrg[nd_][COORD_DIM][KDIM1], Htrg[nd_][COORD_DIM*2][KDIM1];

    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0] = Vec::LoadAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1] = Vec::LoadAligned(&Vt[i*KDIM1+1][t]);
      Gtrg[i][0][0] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+0][t]);
      Gtrg[i][0][1] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+1][t]);
      Gtrg[i][1][0] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+0][t]);
      Gtrg[i][1][1] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+1][t]);
      Gtrg[i][2][0] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+0][t]);
      Gtrg[i][2][1] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+1][t]);
      for (long j = 0; j < 6; j++) {
        Htrg[i][j][0] = Vec::LoadAligned(&Ht[(j*nd_+i)*KDIM1+0][t]);
        Htrg[i][j][1] = Vec::LoadAligned(&Ht[(j*nd_+i)*KDIM1+1][t]);
      }

    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::approx_sincos<-1>(sin_izkR, cos_izkR, izkR[1]);
      exp_izkR = sctl::approx_exp<-1>(izkR[0]);
      // exp(ikr)/r
      Vec G0 = cos_izkR * exp_izkR * Rinv;
      Vec G1 = sin_izkR * exp_izkR * Rinv;
      // (ikr-1)*exp(ikr)/r^3
      Vec H0 = (izkR[0] - (1.0))*G0 - izkR[1]*G1;
      Vec H1 = izkR[1]*G0 + (izkR[0] - (1.0))*G1;
      H0 = H0 * Rinv2;
      H1 = H1 * Rinv2;

      Vec tmp0 = zk_[1]*zk_[1] - zk_[0]*zk_[0];
      Vec tmp1 = zk_[0]*zk_[1]*(-2.0);
      tmp0 = (-3.0)*(Rinv*zk_[1]+Rinv2) - tmp0;
      tmp1 = (3.0)*Rinv*zk_[0] - tmp1;
      // exp(ikr)*(-(ik)^2-3/r^2+3ik/r)/r^3
      Vec J0 = (G0 * tmp0 - G1 * tmp1)*Rinv2;
      Vec J1 = (G1 * tmp0 + G0 * tmp1)*Rinv2;

      // (ikr-1)*exp(ikr)/r^3 * (xt - xs)
      Vec Ztmp[COORD_DIM][KDIM0];
      for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          Ztmp[dimi][0] = H0*dX[dimi];
          Ztmp[dimi][1] = H1*dX[dimi];
      }
      Vec Hctmp[COORD_DIM*2][KDIM0];
      Hctmp[0][0] = -J0*dX[0]*dX[0]+H0;
      Hctmp[0][1] = -J1*dX[0]*dX[0]+H1;
      Hctmp[1][0] = -J0*dX[1]*dX[1]+H0;
      Hctmp[1][1] = -J1*dX[1]*dX[1]+H1;
      Hctmp[2][0] = -J0*dX[2]*dX[2]+H0;
      Hctmp[2][1] = -J1*dX[2]*dX[2]+H1;
      Hctmp[3][0] = -J0*dX[0]*dX[1];
      Hctmp[3][1] = -J1*dX[0]*dX[1];
      Hctmp[4][0] = -J0*dX[0]*dX[2];
      Hctmp[4][1] = -J1*dX[0]*dX[2];
      Hctmp[5][0] = -J0*dX[1]*dX[2];
      Hctmp[5][1] = -J1*dX[1]*dX[2];


      for (long i = 0; i < nd_; i++) {
        Vtrg[i][0] += G0 * Vsrc[s*nd_*KDIM0+i*KDIM0+0] - G1 * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        Vtrg[i][1] += G1 * Vsrc[s*nd_*KDIM0+i*KDIM0+0] + G0 * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          Gtrg[i][dimi][0] += Ztmp[dimi][0] * Vsrc[s*nd_*KDIM0+i*KDIM0+0] - Ztmp[dimi][1] * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
          Gtrg[i][dimi][1] += Ztmp[dimi][1] * Vsrc[s*nd_*KDIM0+i*KDIM0+0] + Ztmp[dimi][0] * Vsrc[s*nd_*KDIM0+i*KDIM0+1];
        }
        // hessian
        for (sctl::Integer dimi = 0; dimi < COORD_DIM*2; dimi++){
          Htrg[i][dimi][0] += Hctmp[dimi][0]*Vsrc[s*nd_*KDIM0+i*KDIM0+0] - Hctmp[dimi][1]*Vsrc[s*nd_*KDIM0+i*KDIM0+1];
          Htrg[i][dimi][1] += Hctmp[dimi][0]*Vsrc[s*nd_*KDIM0+i*KDIM0+1] + Hctmp[dimi][1]*Vsrc[s*nd_*KDIM0+i*KDIM0+0];
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
      for (long j = 0; j < 6; j++) {
        Htrg[i][j][0].StoreAligned(&Ht[(j*nd_+i)*KDIM1+0][t]);
        Htrg[i][j][1].StoreAligned(&Ht[(j*nd_+i)*KDIM1+1][t]);
      }

    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_*KDIM1; j++) {
      pot[i*nd_*KDIM1+j] += Vt[j][i];
    }
    for (long j = 0; j < nd_*KDIM1*COORD_DIM; j++) {
      grad[i*nd_*KDIM1*COORD_DIM+j] += Gt[j][i];
    }
    for (long j=0; j < nd_*KDIM1*COORD_DIM*2; j++) {
      hess[i*nd_*KDIM1*COORD_DIM*2+j] += Ht[j][i];
    }
  }
}


// dipole, potential
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void h3ddirectdp_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
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
  zk_[0] = Vec::Load1(zk+0);
  zk_[1] = Vec::Load1(zk+1);
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::approx_sincos<-1>(sin_izkR, cos_izkR, izkR[1]);
      exp_izkR = sctl::approx_exp<-1>(izkR[0]);
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
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void h3ddirectdg_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
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
  zk_[0] = Vec::Load1(zk+0);
  zk_[1] = Vec::Load1(zk+1);
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::approx_sincos<-1>(sin_izkR, cos_izkR, izkR[1]);
      exp_izkR = sctl::approx_exp<-1>(izkR[0]);
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


// dipole, potential, gradient and hessian
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void h3ddirectdh_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, Real* hess, const Real* thresh) {
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
  sctl::Matrix<Real> Ht(nd_*KDIM1*COORD_DIM*2, Ntrg_);

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
    Ht = 0;

  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec zk_[2];
  zk_[0] = Vec::Load1(zk+0);
  zk_[1] = Vec::Load1(zk+1);
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
    // load potential and gradient and hessian
    Vec Vtrg[nd_][KDIM1], Gtrg[nd_][COORD_DIM][KDIM1], Htrg[nd_][COORD_DIM*2][KDIM1];

    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0] = Vec::LoadAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1] = Vec::LoadAligned(&Vt[i*KDIM1+1][t]);
      Gtrg[i][0][0] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+0][t]);
      Gtrg[i][0][1] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+1][t]);
      Gtrg[i][1][0] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+0][t]);
      Gtrg[i][1][1] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+1][t]);
      Gtrg[i][2][0] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+0][t]);
      Gtrg[i][2][1] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+1][t]);
      for (long j = 0; j < 6; j++) {
        Htrg[i][j][0] = Vec::LoadAligned(&Ht[(j*nd_+i)*KDIM1+0][t]);
        Htrg[i][j][1] = Vec::LoadAligned(&Ht[(j*nd_+i)*KDIM1+1][t]);
      }

    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec Rinv4 = Rinv2 * Rinv2;
      Vec Rinv5 = Rinv4 * Rinv;
      Vec Rinv6 = Rinv4 * Rinv2;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::approx_sincos<-1>(sin_izkR, cos_izkR, izkR[1]);
      exp_izkR = sctl::approx_exp<-1>(izkR[0]);
      // exp(ikr)/r
      Vec G0 = cos_izkR * exp_izkR * Rinv;
      Vec G1 = sin_izkR * exp_izkR * Rinv;
      Vec zf1[2];
      zf1[0] = izkR[0]-(1.0);
      zf1[1] = izkR[1];
      // (1-ikr)*exp(ikr)/r^3
      Vec H0 = ( zf1[1]*G1 - zf1[0]*G0) * Rinv2;
      Vec H1 = (-zf1[1]*G0 - zf1[0]*G1) * Rinv2;
      Vec tmp0 = zk[1]*zk[1] - zk[0]*zk[0];
      Vec tmp1 = zk[0]*zk[1]*(-2.0);
      tmp0 = (-3.0)*(Rinv*zk[1]+Rinv2) - tmp0;
      tmp1 = (3.0)*Rinv*zk[0] - tmp1;
      // exp(ikr)*(-(ik)^2-3/r^2+3ik/r)/r^3

      Vec J0 = (G0 * tmp0 - G1 * tmp1)*Rinv2;
      Vec J1 = (G1 * tmp0 + G0 * tmp1)*Rinv2;

      Vec cdc[2], cdc2[2], cdc3[2], cdc4[2];
      Vec izkR2[2];
      izkR2[0] = izkR[0]*izkR[0] - izkR[1]*izkR[1];
      izkR2[1] = (2.0)*izkR[0]*izkR[1];
      Vec zk2[2];
      zk2[0] = zk_[0]*zk_[0] - zk_[1]*zk_[1];
      zk2[1] = (2.0)*zk_[0]*zk_[1];
      tmp0 = (3.0)*zf1[0] - izkR2[0];
      tmp1 = (3.0)*zf1[1] - izkR2[1];
      cdc[0] = (G0*tmp0 - G1*tmp1)*Rinv4;
      cdc[1] = (G1*tmp0 + G0*tmp1)*Rinv4;
      Vec tmp2,tmp3,tmp4,tmp5;
      tmp2 = zf1[0]+(2.0);
      tmp3 = zf1[1];
      tmp0 = (zk2[0]*tmp2-zk2[1]*tmp3)*Rinv4;
      tmp1 = (zk2[0]*tmp3+zk2[1]*tmp2)*Rinv4;
      tmp2 = (-zk_[1]*izkR[0] - zk_[0]*izkR[1])*(7.0)*Rinv5;
      tmp3 = (-zk_[1]*izkR[1] + zk_[0]*izkR[0])*(7.0)*Rinv5;
      tmp0 = tmp0 - (15.0)*Rinv6*zf1[0] + tmp2;
      tmp1 = tmp1 - (15.0)*Rinv6*zf1[1] + tmp3;
      cdc2[0] = G0*tmp0 - G1*tmp1;
      cdc2[1] = G1*tmp0 + G0*tmp1;
      tmp0 = (6.0)*zf1[0] - (2.0)*izkR2[0];
      tmp1 = (6.0)*zf1[1] - (2.0)*izkR2[1];
      cdc3[0] = (G0*tmp0 - G1*tmp1)*Rinv4;
      cdc3[1] = (G0*tmp1 + G1*tmp0)*Rinv4;
      tmp2 = zf1[0]+(2.0);
      tmp3 = zf1[1];
      tmp0 = zk2[0]*tmp2 - zk2[1]*tmp3;
      tmp1 = zk2[0]*tmp3 + zk2[1]*tmp2;
      tmp0 += - (15.0)*Rinv2*zf1[0] + (7.0)*Rinv2*izkR2[0];
      tmp1 += - (15.0)*Rinv2*zf1[1] + (7.0)*Rinv2*izkR2[1];
      cdc4[0] = (G0*tmp0-G1*tmp1)*Rinv4;
      cdc4[1] = (G0*tmp1+G1*tmp0)*Rinv4;

      Vec D0, D1;
      Vec crossval[3][2];
      long s_ind = s*nd_*COORD_DIM*KDIM0;
      for (long i = 0; i < nd_; i++) {
        long si_ind = s_ind+i*COORD_DIM*KDIM0;
        D0 = dX[0]*Gsrc[si_ind+0] + dX[1]*Gsrc[si_ind+2] + dX[2]*Gsrc[si_ind+4];
        D1 = dX[0]*Gsrc[si_ind+1] + dX[1]*Gsrc[si_ind+3] + dX[2]*Gsrc[si_ind+5];
        tmp0 = J0*D0 - J1*D1;
        tmp1 = J1*D0 + J0*D1;
        crossval[0][0] = Gsrc[si_ind+0]*dX[1] + Gsrc[si_ind+2]*dX[0];
        crossval[0][1] = Gsrc[si_ind+1]*dX[1] + Gsrc[si_ind+3]*dX[0];
        crossval[1][0] = Gsrc[si_ind+0]*dX[2] + Gsrc[si_ind+4]*dX[0];
        crossval[1][1] = Gsrc[si_ind+1]*dX[2] + Gsrc[si_ind+5]*dX[0];
        crossval[2][0] = Gsrc[si_ind+4]*dX[1] + Gsrc[si_ind+2]*dX[2];
        crossval[2][1] = Gsrc[si_ind+5]*dX[1] + Gsrc[si_ind+3]*dX[2];

        Vtrg[i][0] += H0 * D0 - H1 * D1;
        Vtrg[i][1] += H1 * D0 + H0 * D1;
        for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          Gtrg[i][dimi][0] += tmp0 * dX[dimi];
          Gtrg[i][dimi][1] += tmp1 * dX[dimi];
          Gtrg[i][dimi][0] += H0*Gsrc[si_ind+dimi*KDIM0+0] - H1*Gsrc[si_ind+dimi*KDIM0+1];
          Gtrg[i][dimi][1] += H1*Gsrc[si_ind+dimi*KDIM0+0] + H0*Gsrc[si_ind+dimi*KDIM0+1];
        }

        // hessian
        // dipole part
        tmp0 = cdc[0]*D0 - cdc[1]*D1;
        tmp1 = cdc[0]*D1 + cdc[1]*D0;
        tmp2 = cdc2[0]*D0 - cdc2[1]*D1;
        tmp3 = cdc2[0]*D1 + cdc2[1]*D0;
        // hess component 1-3
        for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          tmp4 = cdc3[0]*Gsrc[si_ind+dimi*KDIM0+0] - cdc3[1]*Gsrc[si_ind+dimi*KDIM0+1];
          tmp5 = cdc3[0]*Gsrc[si_ind+dimi*KDIM0+1] + cdc3[1]*Gsrc[si_ind+dimi*KDIM0+0];

          Htrg[i][dimi][0] += tmp0 + tmp2*dX[dimi]*dX[dimi] + tmp4*dX[dimi];
          Htrg[i][dimi][1] += tmp1 + tmp3*dX[dimi]*dX[dimi] + tmp5*dX[dimi];
        }
        tmp0 = cdc4[0]*D0 - cdc4[1]*D1;
        tmp1 = cdc4[0]*D1 + cdc4[1]*D0;
        // hess component 4-5
        for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          if(dimi==0)
            tmp5 = dX[0]*dX[1];
          else if(dimi==1)
            tmp5 = dX[0]*dX[2];
          else
            tmp5 = dX[2]*dX[1];
          tmp3 = cdc[0]*crossval[dimi][0] - cdc[1]*crossval[dimi][1];
          tmp4 = cdc[0]*crossval[dimi][1] + cdc[1]*crossval[dimi][0];
          Htrg[i][3+dimi][0] += tmp3 + tmp0*tmp5;
          Htrg[i][3+dimi][1] += tmp4 + tmp1*tmp5;
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
      for (long j = 0; j < 6; j++) {
        Htrg[i][j][0].StoreAligned(&Ht[(j*nd_+i)*KDIM1+0][t]);
        Htrg[i][j][1].StoreAligned(&Ht[(j*nd_+i)*KDIM1+1][t]);
      }

    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_*KDIM1; j++) {
      pot[i*nd_*KDIM1+j] += Vt[j][i];
    }
    for (long j=0; j < nd_*KDIM1*COORD_DIM; j++) {
      grad[i*nd_*KDIM1*COORD_DIM+j] += Gt[j][i];
    }
    for (long j=0; j < nd_*KDIM1*COORD_DIM*2; j++) {
      hess[i*nd_*KDIM1*COORD_DIM*2+j] += Ht[j][i];
    }

  }
}


// charge and dipole, potential
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void h3ddirectcdp_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
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
  zk_[0] = Vec::Load1(zk+0);
  zk_[1] = Vec::Load1(zk+1);
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::approx_sincos<-1>(sin_izkR, cos_izkR, izkR[1]);
      exp_izkR = sctl::approx_exp<-1>(izkR[0]);
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
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void h3ddirectcdg_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
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
  zk_[0] = Vec::Load1(zk+0);
  zk_[1] = Vec::Load1(zk+1);

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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec izkR[2] = {-zk[1]*R, zk[0]*R};
      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::approx_sincos<-1>(sin_izkR, cos_izkR, izkR[1]);
      exp_izkR = sctl::approx_exp<-1>(izkR[0]);
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


// charge and dipole, potential, gradient and hessian
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void h3ddirectcdh_vec_cpp(const int32_t* nd, const Real* zk, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, Real* hess, const Real* thresh) {
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
  sctl::Matrix<Real> Ht(nd_*KDIM1*COORD_DIM*2, Ntrg_);

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
    Ht = 0;

  }

  static constexpr sctl::Integer VecLen = MaxVecLen;
  using Vec = sctl::Vec<Real,VecLen>;
  Vec zk_[2];
  zk_[0] = Vec::Load1(zk+0);
  zk_[1] = Vec::Load1(zk+1);
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
    // load potential and gradient and hessian
    Vec Vtrg[nd_][KDIM1], Gtrg[nd_][COORD_DIM][KDIM1], Htrg[nd_][COORD_DIM*2][KDIM1];

    for (long i = 0; i < nd_; i++) {
      Vtrg[i][0] = Vec::LoadAligned(&Vt[i*KDIM1+0][t]);
      Vtrg[i][1] = Vec::LoadAligned(&Vt[i*KDIM1+1][t]);
      Gtrg[i][0][0] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+0][t]);
      Gtrg[i][0][1] = Vec::LoadAligned(&Gt[(0*nd_+i)*KDIM1+1][t]);
      Gtrg[i][1][0] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+0][t]);
      Gtrg[i][1][1] = Vec::LoadAligned(&Gt[(1*nd_+i)*KDIM1+1][t]);
      Gtrg[i][2][0] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+0][t]);
      Gtrg[i][2][1] = Vec::LoadAligned(&Gt[(2*nd_+i)*KDIM1+1][t]);
      for (long j = 0; j < 6; j++) {
        Htrg[i][j][0] = Vec::LoadAligned(&Ht[(j*nd_+i)*KDIM1+0][t]);
        Htrg[i][j][1] = Vec::LoadAligned(&Ht[(j*nd_+i)*KDIM1+1][t]);
      }

    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec R = R2 * Rinv;
      Vec Rinv2 = Rinv * Rinv;
      Vec Rinv4 = Rinv2 * Rinv2;
      Vec Rinv5 = Rinv4 * Rinv;
      Vec Rinv6 = Rinv4 * Rinv2;
      Vec izkR[2] = {-zk_[1]*R, zk_[0]*R};

      Vec sin_izkR, cos_izkR, exp_izkR;
      sctl::approx_sincos<-1>(sin_izkR, cos_izkR, izkR[1]);
      exp_izkR = sctl::approx_exp<-1>(izkR[0]);
      // exp(ikr)/r
      Vec G0 = cos_izkR * exp_izkR * Rinv;
      Vec G1 = sin_izkR * exp_izkR * Rinv;
      Vec zf1[2];
      zf1[0] = izkR[0]-(1.0);
      zf1[1] = izkR[1];
      // (1-ikr)*exp(ikr)/r^3
      Vec H0 = ( zf1[1]*G1 - zf1[0]*G0) * Rinv2;
      Vec H1 = (-zf1[1]*G0 - zf1[0]*G1) * Rinv2;
      Vec tmp0 = zk_[1]*zk_[1] - zk_[0]*zk_[0];
      Vec tmp1 = zk_[0]*zk_[1]*(-2.0);
      tmp0 = (-3.0)*(Rinv*zk_[1]+Rinv2) - tmp0;
      tmp1 = (3.0)*Rinv*zk_[0] - tmp1;
      // exp(ikr)*(-(ik)^2-3/r^2+3ik/r)/r^3

      Vec J0 = (G0 * tmp0 - G1 * tmp1)*Rinv2;
      Vec J1 = (G1 * tmp0 + G0 * tmp1)*Rinv2;
      // (ikr-1)*exp(ikr)/r^3 * (xt - xs)
      Vec Zctmp[COORD_DIM][KDIM0];
      for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          Zctmp[dimi][0] = -H0*dX[dimi];
          Zctmp[dimi][1] = -H1*dX[dimi];
      }
      Vec Hctmp[COORD_DIM*2][KDIM0];
      Hctmp[0][0] = -J0*dX[0]*dX[0]-H0;
      Hctmp[0][1] = -J1*dX[0]*dX[0]-H1;
      Hctmp[1][0] = -J0*dX[1]*dX[1]-H0;
      Hctmp[1][1] = -J1*dX[1]*dX[1]-H1;
      Hctmp[2][0] = -J0*dX[2]*dX[2]-H0;
      Hctmp[2][1] = -J1*dX[2]*dX[2]-H1;
      Hctmp[3][0] = -J0*dX[0]*dX[1];
      Hctmp[3][1] = -J1*dX[0]*dX[1];
      Hctmp[4][0] = -J0*dX[0]*dX[2];
      Hctmp[4][1] = -J1*dX[0]*dX[2];
      Hctmp[5][0] = -J0*dX[1]*dX[2];
      Hctmp[5][1] = -J1*dX[1]*dX[2];

      Vec cdc[2], cdc2[2], cdc3[2], cdc4[2];
      Vec izkR2[2];
      izkR2[0] = izkR[0]*izkR[0] - izkR[1]*izkR[1];
      izkR2[1] = (2.0)*izkR[0]*izkR[1];
      Vec zk2[2];
      zk2[0] = zk_[0]*zk_[0] - zk_[1]*zk_[1];
      zk2[1] = (2.0)*zk_[0]*zk_[1];
      tmp0 = (3.0)*zf1[0] - izkR2[0];
      tmp1 = (3.0)*zf1[1] - izkR2[1];
      cdc[0] = (G0*tmp0 - G1*tmp1)*Rinv4;
      cdc[1] = (G1*tmp0 + G0*tmp1)*Rinv4;
      Vec tmp2,tmp3,tmp4,tmp5;
      tmp2 = zf1[0]+(2.0);
      tmp3 = zf1[1];
      tmp0 = (zk2[0]*tmp2-zk2[1]*tmp3)*Rinv4;
      tmp1 = (zk2[0]*tmp3+zk2[1]*tmp2)*Rinv4;
      tmp2 = (-zk_[1]*izkR[0] - zk_[0]*izkR[1])*(7.0)*Rinv5;
      tmp3 = (-zk_[1]*izkR[1] + zk_[0]*izkR[0])*(7.0)*Rinv5;
      tmp0 = tmp0 - (15.0)*Rinv6*zf1[0] + tmp2;
      tmp1 = tmp1 - (15.0)*Rinv6*zf1[1] + tmp3;
      cdc2[0] = G0*tmp0 - G1*tmp1;
      cdc2[1] = G1*tmp0 + G0*tmp1;
      tmp0 = (6.0)*zf1[0] - (2.0)*izkR2[0];
      tmp1 = (6.0)*zf1[1] - (2.0)*izkR2[1];
      cdc3[0] = (G0*tmp0 - G1*tmp1)*Rinv4;
      cdc3[1] = (G0*tmp1 + G1*tmp0)*Rinv4;
      tmp2 = zf1[0]+(2.0);
      tmp3 = zf1[1];
      tmp0 = zk2[0]*tmp2 - zk2[1]*tmp3;
      tmp1 = zk2[0]*tmp3 + zk2[1]*tmp2;
      tmp0 += - (15.0)*Rinv2*zf1[0] + (7.0)*Rinv2*izkR2[0];
      tmp1 += - (15.0)*Rinv2*zf1[1] + (7.0)*Rinv2*izkR2[1];
      cdc4[0] = (G0*tmp0-G1*tmp1)*Rinv4;
      cdc4[1] = (G0*tmp1+G1*tmp0)*Rinv4;

      Vec D0, D1;
      Vec crossval[3][2];

      long s_ind = s*nd_*COORD_DIM*KDIM0;
      for (long i = 0; i < nd_; i++) {
        long si_ind = s_ind + i*COORD_DIM*KDIM0;

        // temp values
        D0 = dX[0]*Gsrc[si_ind+0] + dX[1]*Gsrc[si_ind+2] + dX[2]*Gsrc[si_ind+4];
        D1 = dX[0]*Gsrc[si_ind+1] + dX[1]*Gsrc[si_ind+3] + dX[2]*Gsrc[si_ind+5];
        tmp0 = J0*D0 - J1*D1;
        tmp1 = J1*D0 + J0*D1;
        crossval[0][0] = Gsrc[si_ind+0]*dX[1] + Gsrc[si_ind+2]*dX[0];
        crossval[0][1] = Gsrc[si_ind+1]*dX[1] + Gsrc[si_ind+3]*dX[0];
        crossval[1][0] = Gsrc[si_ind+0]*dX[2] + Gsrc[si_ind+4]*dX[0];
        crossval[1][1] = Gsrc[si_ind+1]*dX[2] + Gsrc[si_ind+5]*dX[0];
        crossval[2][0] = Gsrc[si_ind+4]*dX[1] + Gsrc[si_ind+2]*dX[2];
        crossval[2][1] = Gsrc[si_ind+5]*dX[1] + Gsrc[si_ind+3]*dX[2];


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
        // hessian
        // charge part
        for (sctl::Integer dimi = 0; dimi < COORD_DIM*2; dimi++){
          Htrg[i][dimi][0] += Hctmp[dimi][0]*Vsrc[s*nd_*KDIM0+i*KDIM0+0] - Hctmp[dimi][1]*Vsrc[s*nd_*KDIM0+i*KDIM0+1];
          Htrg[i][dimi][1] += Hctmp[dimi][0]*Vsrc[s*nd_*KDIM0+i*KDIM0+1] + Hctmp[dimi][1]*Vsrc[s*nd_*KDIM0+i*KDIM0+0];
        }

        // dipole part
        tmp0 = cdc[0]*D0 - cdc[1]*D1;
        tmp1 = cdc[0]*D1 + cdc[1]*D0;
        tmp2 = cdc2[0]*D0 - cdc2[1]*D1;
        tmp3 = cdc2[0]*D1 + cdc2[1]*D0;
        // hess component 1-3
        for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          tmp4 = cdc3[0]*Gsrc[si_ind+dimi*KDIM0+0] - cdc3[1]*Gsrc[si_ind+dimi*KDIM0+1];
          tmp5 = cdc3[0]*Gsrc[si_ind+dimi*KDIM0+1] + cdc3[1]*Gsrc[si_ind+dimi*KDIM0+0];

          Htrg[i][dimi][0] += tmp0 + tmp2*dX[dimi]*dX[dimi] + tmp4*dX[dimi];
          Htrg[i][dimi][1] += tmp1 + tmp3*dX[dimi]*dX[dimi] + tmp5*dX[dimi];
        }
        tmp0 = cdc4[0]*D0 - cdc4[1]*D1;
        tmp1 = cdc4[0]*D1 + cdc4[1]*D0;
        // hess component 4-5
        for (sctl::Integer dimi = 0; dimi < COORD_DIM; dimi++){
          if(dimi==0)
            tmp5 = dX[0]*dX[1];
          else if(dimi==1)
            tmp5 = dX[0]*dX[2];
          else
            tmp5 = dX[2]*dX[1];
          tmp3 = cdc[0]*crossval[dimi][0] - cdc[1]*crossval[dimi][1];
          tmp4 = cdc[0]*crossval[dimi][1] + cdc[1]*crossval[dimi][0];
          Htrg[i][3+dimi][0] += tmp3 + tmp0*tmp5;
          Htrg[i][3+dimi][1] += tmp4 + tmp1*tmp5;
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

      for (long j = 0; j < 6; j++) {
        Htrg[i][j][0].StoreAligned(&Ht[(j*nd_+i)*KDIM1+0][t]);
        Htrg[i][j][1].StoreAligned(&Ht[(j*nd_+i)*KDIM1+1][t]);
      }


    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_*KDIM1; j++) {
      pot[i*nd_*KDIM1+j] += Vt[j][i];
    }
    for (long j=0; j < nd_*KDIM1*COORD_DIM; j++) {
      grad[i*nd_*KDIM1*COORD_DIM+j] += Gt[j][i];
    }
    for (long j=0; j < nd_*KDIM1*COORD_DIM*2; j++) {
      hess[i*nd_*KDIM1*COORD_DIM*2+j] += Ht[j][i];
    }

  }
}


// Laplace kernels
// charge, potential
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void l3ddirectcp_vec_cpp(const int32_t* nd, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

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
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void l3ddirectcg_vec_cpp(const int32_t* nd, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

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


// charge, potential, gradient and hessian
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void l3ddirectch_vec_cpp(const int32_t* nd, const Real* sources, const Real* charge, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, Real* hess, const Real* thresh) {
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
  sctl::Matrix<Real> Ht(nd_*COORD_DIM*2, Ntrg_);
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
    Ht = 0;
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
    // load potential and gradient and hessian
    Vec Vtrg[nd_], Gtrg[nd_][COORD_DIM], Htrg[nd_][COORD_DIM*2];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i] = Vec::LoadAligned(&Vt[i][t]);
      Gtrg[i][0] = Vec::LoadAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1] = Vec::LoadAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2] = Vec::LoadAligned(&Gt[2*nd_+i][t]);
      for (long j = 0; j < 6; j++) {
        Htrg[i][j] = Vec::LoadAligned(&Ht[j*nd_+i][t]);
      }
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec nRinv3 = -Rinv*Rinv*Rinv;
      Vec Rinv5 = -nRinv3*Rinv*Rinv;
      Vec htmp1 = Rinv5*(3.0*dX[0]*dX[0]-R2);
      Vec htmp2 = Rinv5*(3.0*dX[1]*dX[1]-R2);
      Vec htmp3 = Rinv5*(3.0*dX[2]*dX[2]-R2);
      Vec htmp4 = 3.0*Rinv5*dX[0]*dX[1];
      Vec htmp5 = 3.0*Rinv5*dX[0]*dX[2];
      Vec htmp6 = 3.0*Rinv5*dX[1]*dX[2];
      Vec ztmp[COORD_DIM] = {nRinv3*dX[0], nRinv3*dX[1], nRinv3*dX[2]};
      for (long i = 0; i < nd_; i++) {
        Vtrg[i] += Vsrc[s*nd_+i]*Rinv;
        Gtrg[i][0] += Vsrc[s*nd_+i]*ztmp[0];
        Gtrg[i][1] += Vsrc[s*nd_+i]*ztmp[1];
        Gtrg[i][2] += Vsrc[s*nd_+i]*ztmp[2];
        Htrg[i][0] += htmp1*Vsrc[s*nd_+i];
        Htrg[i][1] += htmp2*Vsrc[s*nd_+i];
        Htrg[i][2] += htmp3*Vsrc[s*nd_+i];
        Htrg[i][3] += htmp4*Vsrc[s*nd_+i];
        Htrg[i][4] += htmp5*Vsrc[s*nd_+i];
        Htrg[i][5] += htmp6*Vsrc[s*nd_+i];
      }
    }
    for (long i = 0; i < nd_; i++) {
      Vtrg[i].StoreAligned(&Vt[i][t]);
      Gtrg[i][0].StoreAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1].StoreAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2].StoreAligned(&Gt[2*nd_+i][t]);
      for (long j = 0; j < 6; j++) {
        Htrg[i][j].StoreAligned(&Ht[j*nd_+i][t]);
      }
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_; j++) {
      pot[i*nd_+j] += Vt[j][i];
    }
    for (long j=0; j < nd_*COORD_DIM; j++) {
      grad[i*nd_*COORD_DIM+j] += Gt[j][i];
    }
    for (long j = 0; j < nd_*COORD_DIM*2; j++) {
      hess[i*nd_*COORD_DIM*2+j] += Ht[j][i];
    }
  }
}


// dipole, potential
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void l3ddirectdp_vec_cpp(const int32_t* nd, const Real* sources, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

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
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void l3ddirectdg_vec_cpp(const int32_t* nd, const Real* sources, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

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


// dipole, potential, gradient and hessian
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void l3ddirectdh_vec_cpp(const int32_t* nd, const Real* sources, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, Real* hess, const Real* thresh) {
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
  sctl::Matrix<Real> Ht(nd_*COORD_DIM*2, Ntrg_);
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
    Ht = 0;
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
    // load potential and gradient and hessian
    Vec Vtrg[nd_], Gtrg[nd_][COORD_DIM], Htrg[nd_][COORD_DIM*2];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i] = Vec::LoadAligned(&Vt[i][t]);
      Gtrg[i][0] = Vec::LoadAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1] = Vec::LoadAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2] = Vec::LoadAligned(&Gt[2*nd_+i][t]);
      for (long j = 0; j < 6; j++) {
        Htrg[i][j] = Vec::LoadAligned(&Ht[j*nd_+i][t]);
      }
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec Rinv2 = Rinv * Rinv;
      Vec Rinv3 = Rinv * Rinv2;
      Vec Rinv5 = Rinv2 * Rinv3;
      Vec nRinv5 =  -3.0*Rinv5;
      Vec dx = dX[0]*Rinv;
      Vec dy = dX[1]*Rinv;
      Vec dz = dX[2]*Rinv;
      for (long i = 0; i < nd_; i++) {
        Vec Dprod = dX[0]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] + dX[1]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] + dX[2]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2];
        Vec RinvDprod = nRinv5*Dprod;
        Vtrg[i] += Dprod*Rinv3;
        Gtrg[i][0] += RinvDprod*dX[0] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0];
        Gtrg[i][1] += RinvDprod*dX[1] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1];
        Gtrg[i][2] += RinvDprod*dX[2] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2];
        Htrg[i][0] -= (Dprod*(5.0*dx*dx-1.0) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0]*dX[0]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0]*dX[0]))*nRinv5;
        Htrg[i][1] -= (Dprod*(5.0*dy*dy-1.0) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1]*dX[1]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1]*dX[1]))*nRinv5;
        Htrg[i][2] -= (Dprod*(5.0*dz*dz-1.0) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2]*dX[2]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2]*dX[2]))*nRinv5;
        Htrg[i][3] -= (Dprod*(5.0*dx*dy) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1]*dX[0]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0]*dX[1]))*nRinv5;
        Htrg[i][4] -= (Dprod*(5.0*dx*dz) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2]*dX[0]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0]*dX[2]))*nRinv5;
        Htrg[i][5] -= (Dprod*(5.0*dy*dz) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2]*dX[1]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1]*dX[2]))*nRinv5;
      }
    }
    for (long i = 0; i < nd_; i++) {
      Vtrg[i].StoreAligned(&Vt[i][t]);
      Gtrg[i][0].StoreAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1].StoreAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2].StoreAligned(&Gt[2*nd_+i][t]);
      for (long j = 0; j < 6; j++) {
        Htrg[i][j].StoreAligned(&Ht[j*nd_+i][t]);
      }
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_; j++) {
      pot[i*nd_+j] += Vt[j][i];
    }
    for (long j=0; j < nd_*COORD_DIM; j++) {
      grad[i*nd_*COORD_DIM+j] += Gt[j][i];
    }
    for (long j = 0; j < nd_*COORD_DIM*2; j++) {
      hess[i*nd_*COORD_DIM*2+j] += Ht[j][i];
    }
  }
}


// charge and dipole, potential
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void l3ddirectcdp_vec_cpp(const int32_t* nd, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, const Real* thresh) {
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

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
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void l3ddirectcdg_vec_cpp(const int32_t* nd, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, const Real* thresh) {
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

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

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


// charge and dipole, potential, gradient and hessian
template <class Real, sctl::Integer MaxVecLen=sctl::DefaultVecLen<Real>()> void l3ddirectcdh_vec_cpp(const int32_t* nd, const Real* sources, const Real* charge, const Real* dipvec, const int32_t* ns, const Real* ztarg, const int32_t* nt, Real* pot, Real* grad, Real* hess, const Real* thresh) {
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
  sctl::Matrix<Real> Ht(nd_*COORD_DIM*2, Ntrg_);
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
    Ht = 0;
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
    // load potential and gradient and hessian
    Vec Vtrg[nd_], Gtrg[nd_][COORD_DIM], Htrg[nd_][2*COORD_DIM];
    for (long i = 0; i < nd_; i++) {
      Vtrg[i] = Vec::LoadAligned(&Vt[i][t]);
      Gtrg[i][0] = Vec::LoadAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1] = Vec::LoadAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2] = Vec::LoadAligned(&Gt[2*nd_+i][t]);
      for (long j = 0; j < 6; j++) {
        Htrg[i][j] = Vec::LoadAligned(&Ht[j*nd_+i][t]);
      }
    }
    for (sctl::Long s = 0; s < Nsrc; s++) {
      Vec dX[COORD_DIM], R2 = Vec::Zero();
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        dX[k] = Xtrg[k] - Vec::Load1(&Xs_[s][k]);
        R2 += dX[k]*dX[k];
      }

      Vec Rinv = sctl::approx_rsqrt<-1>(R2, (R2 > thresh2));

      Vec Rinv2 = Rinv * Rinv;
      Vec Rinv3 = Rinv * Rinv2;
      Vec nRinv3 = -Rinv*Rinv*Rinv;
      Vec Rinv5 = Rinv2*Rinv3;
      Vec nRinv5 =  -3.0*Rinv5;
      Vec ztmp[COORD_DIM] = {nRinv3*dX[0], nRinv3*dX[1], nRinv3*dX[2]};
      Vec dx = dX[0]*Rinv;
      Vec dy = dX[1]*Rinv;
      Vec dz = dX[2]*Rinv;
      Vec htmp1 = Rinv5*(3.0*dX[0]*dX[0]-R2);
      Vec htmp2 = Rinv5*(3.0*dX[1]*dX[1]-R2);
      Vec htmp3 = Rinv5*(3.0*dX[2]*dX[2]-R2);
      Vec htmp4 = 3.0*Rinv5*dX[0]*dX[1];
      Vec htmp5 = 3.0*Rinv5*dX[0]*dX[2];
      Vec htmp6 = 3.0*Rinv5*dX[1]*dX[2];
      for (long i = 0; i < nd_; i++) {
        Vec Dprod = dX[0]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] + dX[1]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] + dX[2]*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2];
        Vtrg[i] += Dprod*Rinv3;
        Vtrg[i] += Vsrc[s*nd_+i]*Rinv;
        Vec RinvDprod = nRinv5*Dprod;
        Gtrg[i][0] += RinvDprod*dX[0] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0] + Vsrc[s*nd_+i]*ztmp[0];
        Gtrg[i][1] += RinvDprod*dX[1] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1] + Vsrc[s*nd_+i]*ztmp[1];
        Gtrg[i][2] += RinvDprod*dX[2] + Rinv3*Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2] + Vsrc[s*nd_+i]*ztmp[2];
        Htrg[i][0] += htmp1*Vsrc[s*nd_+i] - (Dprod*(5.0*dx*dx-1.0) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0]*dX[0]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0]*dX[0]))*nRinv5;
        Htrg[i][1] += htmp2*Vsrc[s*nd_+i] - (Dprod*(5.0*dy*dy-1.0) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1]*dX[1]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1]*dX[1]))*nRinv5;
        Htrg[i][2] += htmp3*Vsrc[s*nd_+i] - (Dprod*(5.0*dz*dz-1.0) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2]*dX[2]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2]*dX[2]))*nRinv5;
        Htrg[i][3] += htmp4*Vsrc[s*nd_+i] - (Dprod*(5.0*dx*dy) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1]*dX[0]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0]*dX[1]))*nRinv5;
        Htrg[i][4] += htmp5*Vsrc[s*nd_+i] - (Dprod*(5.0*dx*dz) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2]*dX[0]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+0]*dX[2]))*nRinv5;
        Htrg[i][5] += htmp6*Vsrc[s*nd_+i] - (Dprod*(5.0*dy*dz) -
                      (Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+2]*dX[1]+Gsrc[s*nd_*COORD_DIM+i*COORD_DIM+1]*dX[2]))*nRinv5;
      }
    }
    for (long i = 0; i < nd_; i++) {
      Vtrg[i].StoreAligned(&Vt[i][t]);
      Gtrg[i][0].StoreAligned(&Gt[0*nd_+i][t]);
      Gtrg[i][1].StoreAligned(&Gt[1*nd_+i][t]);
      Gtrg[i][2].StoreAligned(&Gt[2*nd_+i][t]);
      for (long j = 0; j < 6; j++) {
        Htrg[i][j].StoreAligned(&Ht[j*nd_+i][t]);
      }
    }
  }

  for (long i = 0; i < Ntrg; i++) {
    for (long j = 0; j < nd_; j++) {
      pot[i*nd_+j] += Vt[j][i];
    }
    for (long j = 0; j < nd_*COORD_DIM; j++) {
      grad[i*nd_*COORD_DIM+j] += Gt[j][i];
    }
    for (long j = 0; j < nd_*COORD_DIM*2; j++) {
      hess[i*nd_*COORD_DIM*2+j] += Ht[j][i];
    }
  }
}


#endif //_VEC_KERNELS_HPP_
