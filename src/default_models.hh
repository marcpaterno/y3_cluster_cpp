#ifndef Y3_CLUSTER_DEFAULT_MODELS_HH
#define Y3_CLUSTER_DEFAULT_MODELS_HH


#include <a_cen_t.hh>
#include <a_mis_t.hh>
#include <del_sig_tom.hh>
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
#include <pzsource_t.hh>
#include <pzsource_gaussian_t.hh>
#include "models.hh"

namespace y3_cluster {
    using DefaultModels = Models<MOR_t2,
                                 LO_LC_t,
                                 LC_LT_t,
                                 ZO_ZT_t,
                                 PZSOURCE_GAUSSIAN_t,
                                 ROFFSET_t,
                                 T_CEN_t,
                                 T_MIS_t,
                                 A_CEN_t,
                                 A_MIS_t,
                                 HMB_t,
                                 HMF_t,
                                 DEL_SIG_TOM,
                                 DV_DO_DZ_t,
                                 OMEGA_Z_SDSS>;
}

#endif
