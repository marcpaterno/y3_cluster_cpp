#ifndef Y3_CLUSTER_DEFAULT_MODELS_HH
#define Y3_CLUSTER_DEFAULT_MODELS_HH

#include <mor_t.hh>
#include <lo_lc_t.hh>
#include <lc_lt_t.hh>
#include <zo_zt_t.hh>
#include <roffset_t.hh>
#include <t_cen_t.hh>
#include <t_mis_t.hh>
#include <a_cen_t.hh>
#include <a_mis_t.hh>
#include <hmf_t.hh>
#include <del_sig_cen_t.hh>
#include <del_sig_mis_t.hh>
#include <dv_do_dz_t.hh>
#include <omega_z_t.hh>

namespace y3_cluster
{
  struct DefaultModels
  {
    using MOR = MOR_t,
    using LO_LC = LO_LC_t,
    using LC_LT = LC_LT_t,
    using ZO_ZT = ZO_ZT_t,
    using ROFFSET = ROFFSET_t,
    using T_CEN = T_CEN_t,
    using T_MIS = T_MIS_t,
    using A_CEN = A_CEN_t,
    using A_MIS = A_MIS_t,
    using HMF = HMF_t,
    using DEL_SIG_CEN = DEL_SIG_CEN_t,
    using DEL_SIG_MIS = DEL_SIG_MIS_t,
    using DV_DO_DZ = DV_DO_DZ_t,
    using OMEGA_Z = OMEGA_Z_t,
  };
}

#endif
