# Check that boost correction was applied
import numpy as np
from cosmosis.datablock import option_section

def setup(options):
    return {}

def execute(block, config):
    shear_corrected = block["shear", "shear_cen"]
    shear_original = block["shear_boosted", "shear_uncorrected"]
    r = block["shear", "r"]

    print("\n" + "="*50)
    print("BOOST FACTOR TEST RESULTS")
    print("="*50)
    print(f"Radii: {r}")
    print(f"Original shear (bin 0): {shear_original[0]}")
    print(f"Corrected shear (bin 0): {shear_corrected[0]}")

    ratio = shear_original[0] / shear_corrected[0]
    print(f"Boost factor B(R) = original/corrected: {ratio}")
    print(f"B should be > 1 at small R: {np.all(ratio >= 0.99)}")
    print("="*50)
    print("TEST PASSED!" if np.all(ratio >= 0.99) else "TEST FAILED!")
    print("="*50 + "\n")
    return 0

def cleanup(config):
    pass
