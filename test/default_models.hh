#ifndef Y3_CLUSTER_DEFAULT_MODELS_HH
#define Y3_CLUSTER_DEFAULT_MODELS_HH

#include <gamma_t.hh>

#include <test/mor_t.hh>
#include <test/lo_lc_t.hh>
#include <test/lc_lt_t.hh>
#include <test/zo_zt_t.hh>
#include <test/roffset_t.hh>
#include <test/t_cen_t.hh>
#include <test/t_mis_t.hh>
#include <test/a_cen_t.hh>
#include <test/a_mis_t.hh>
#include <test/hmf_t.hh>
#include <test/del_sig_t.hh>
#include <test/dv_do_dz_t.hh>
#include <test/omega_z_sdss.hh>

namespace y3_cluster
{
    template<typename MOR_,
             typename LO_LC_,
             typename LC_LT_,
             typename ZO_ZT_,
             typename ROFFSET_,
             typename T_CEN_,
             typename T_MIS_,
             typename A_CEN_,
             typename A_MIS_,
             typename HMF_,
             typename DEL_SIG_,
             typename DV_DO_DZ_,
             typename OMEGA_Z_>
    struct Models {
        using MOR = MOR_;
        using LO_LC = LO_LC_;
        using LC_LT = LC_LT_;
        using ZO_ZT = ZO_ZT_;
        using ROFFSET = ROFFSET_;
        using T_CEN = T_CEN_;
        using T_MIS = T_MIS_;
        using A_CEN = A_CEN_;
        using A_MIS = A_MIS_;
        using HMF = HMF_;
        using DEL_SIG = DEL_SIG_;
        using DV_DO_DZ = DV_DO_DZ_;
        using OMEGA_Z = OMEGA_Z_;
    };

    using DefaultModels = Models<MOR_t, LO_LC_t, LC_LT_t, ZO_ZT_t, ROFFSET_t, T_CEN_t, T_MIS_t, A_CEN_t, A_MIS_t, HMF_t, DEL_SIG_t, DV_DO_DZ_t, OMEGA_Z_SDSS>;
}

#endif
