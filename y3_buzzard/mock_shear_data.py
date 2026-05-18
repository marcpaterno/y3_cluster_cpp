# Mock Shear Data Module - creates fake shear data for testing
# this is not accurate, just a simple power-law with amplitude depending on lambda and z bins
# created by Arwa to do test with out running the mock pipeline b/c running that with GPU is giving errors for now
import numpy as np
from cosmosis.datablock import option_section

def setup(options):
    n_lambda_bins = options.get_int(option_section, "n_lambda_bins", default=4)
    n_z_bins = options.get_int(option_section, "n_z_bins", default=3)
    n_r_bins = options.get_int(option_section, "n_r_bins", default=10)
    r_min = options.get_double(option_section, "r_min", default=0.3)
    r_max = options.get_double(option_section, "r_max", default=5.0)
    print(f"MockShearData: {n_lambda_bins} lambda x {n_z_bins} z x {n_r_bins} r bins")
    return n_lambda_bins, n_z_bins, n_r_bins, r_min, r_max

def execute(block, config):
    n_lambda_bins, n_z_bins, n_r_bins, r_min, r_max = config
    r = np.logspace(np.log10(r_min), np.log10(r_max), n_r_bins)
    n_bins = n_lambda_bins * n_z_bins
    shear_cen = np.zeros((n_bins, n_r_bins))
    for ij in range(n_bins):
        l_bin = ij // n_z_bins
        z_bin = ij % n_z_bins
        amp = 50.0 * (1 + 0.5 * l_bin) * (1 - 0.1 * z_bin)
        shear_cen[ij] = amp * (r / 0.5)**(-0.8)
    block["shear", "r"] = r
    block["shear", "shear_cen"] = shear_cen
    print(f"MockShearData: wrote shear shape {shear_cen.shape}, r=[{r.min():.2f}, {r.max():.2f}]")
    return 0

def cleanup(config):
    pass
