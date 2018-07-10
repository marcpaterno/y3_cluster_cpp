#ifndef Y3_CLUSTER_DEFAULT_MODELS_HH
#define Y3_CLUSTER_DEFAULT_MODELS_HH

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
  struct DefaultModels
  {
    using MOR = MOR_t;
    using LO_LC = LO_LC_t;
    using LC_LT = LC_LT_t;
    using ZO_ZT = ZO_ZT_t;
    using ROFFSET = ROFFSET_t;
    using T_CEN = T_CEN_t;
    using T_MIS = T_MIS_t;
    using A_CEN = A_CEN_t;
    using A_MIS = A_MIS_t;
    using HMF = HMF_t;
    using DEL_SIG = DEL_SIG_t;
    using DV_DO_DZ = DV_DO_DZ_t;
    using OMEGA_Z = OMEGA_Z_SDSS;
  };
}

#endif
