"""Load plob_ltr_params EMG coefficients into the datablock.

This module reads the pre-baked EMG spline coefficients from
data/prj_params/plob_ltr_params.npz and writes them to the
datablock section 'plob_ltr_params'.

The EMG kernel P(lambda_ob | lambda_tr, z) uses these coefficients
to compute (mu, sigma, tau, fprj) as functions of (lambda_tr, z).

Coefficients:
  - a_mu, b_mu: mu(ltr, z) = a_mu(z) + b_mu(z) * ltr
  - a_sig, b_sig: sigma(ltr, z) = b_sig(z) * ltr^a_sig(z)
  - a_tau, b_tau: tau(ltr, z) = b_tau(z) / ltr^a_tau(z)
  - a_fprj, b_fprj: fprj(ltr, z) = min(1, b_fprj(z) / (1+exp(-ltr))^a_fprj(z))
"""
import os
import numpy as np
from cosmosis.datablock import option_section


_PLOB_NPZ_DEFAULT = os.path.join(
    os.environ.get("Y3_CLUSTER_CPP_DIR",
                   "/pscratch/sd/j/jesteves/y3_cluster_cpp"),
    "data", "prj_params", "plob_ltr_params.npz",
)


def setup(options):
    cfg = {}
    try:
        cfg['npz_path'] = options.get_string(option_section, 'npz_path')
    except Exception:
        cfg['npz_path'] = os.environ.get("RICHNESS_SELECTION_PLOB_NPZ",
                                         _PLOB_NPZ_DEFAULT)
    return cfg


def execute(block, config):
    path = config['npz_path']
    print(f"[plob_ltr_loader] Loading EMG coefficients from {path}", flush=True)

    data = np.load(path)

    keys = ('z', 'a_tau', 'b_tau', 'a_mu', 'b_mu',
            'a_sig', 'b_sig', 'a_fprj', 'b_fprj')

    for k in keys:
        arr = np.asarray(data[k], dtype=np.float64)
        block['plob_ltr_params', k] = arr

    print(f"[plob_ltr_loader] Loaded {len(keys)} arrays, "
          f"z grid has {len(block['plob_ltr_params', 'z'])} points", flush=True)
    return 0


def cleanup(config):
    return 0
