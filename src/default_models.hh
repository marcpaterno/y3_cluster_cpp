#ifndef Y3_CLUSTER_DEFAULT_MODELS_HH
#define Y3_CLUSTER_DEFAULT_MODELS_HH

#include <gamma_t.hh>

#include <a_cen_t.hh>
#include <a_mis_t.hh>
#include <del_sig_y1.hh>
#include <dv_do_dz_t.hh>
#include <hmb_t.hh>
#include <hmf_t.hh>
#include <lc_lt_t.hh>
#include <lo_lc_t.hh>
#include <mor_t2.hh>
#include <roffset_t.hh>
#include <sample_variance.hh>
#include <t_cen_t.hh>
#include <t_mis_t.hh>
#include <omega_z_sdss.hh>
#include <zo_zt_t.hh>

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
             typename HMB_,
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
        // Halo Mass Bias and Halo Mass Function
        using HMB = HMB_;
        using HMF = HMF_;
        using DEL_SIG = DEL_SIG_;
        using DV_DO_DZ = DV_DO_DZ_;
        using OMEGA_Z = OMEGA_Z_;
    };

    using DefaultModels = Models<MOR_t2,
                                 LO_LC_t,
                                 LC_LT_t,
                                 ZO_ZT_t,
                                 ROFFSET_t,
                                 T_CEN_t,
                                 T_MIS_t,
                                 A_CEN_t,
                                 A_MIS_t,
                                 HMB_t,
                                 HMF_t,
                                 DEL_SIG_y1,
                                 DV_DO_DZ_t,
                                 OMEGA_Z_SDSS>;
}

#endif
