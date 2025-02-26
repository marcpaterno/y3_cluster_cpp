import numpy as np
import os
from cosmosis.datablock import option_section, names

def setup(options):
    section = option_section
    fname_z = options.get_string(section, "fname_z")
    fname_lt = options.get_string(section, "fname_l")
    fname_bin = options.get_string(section, "fname_bin")
    return fname_z, fname_lt, fname_bin

def execute(block, config):
    fname_z, fname_lt, fname_bin = config

    basepath = "%s/data/" % os.environ["Y3_CLUSTER_CPP_DIR"]
    zt = np.array(np.genfromtxt(basepath+fname_z))
    lt = np.array(np.genfromtxt(basepath+fname_lt))
    bins = [np.array(np.genfromtxt(basepath+fname_bin.format(i))) for i in range(4)]

    check(zt, "zt")
    check(lt, "lt")
    check(bins[0], "bin1")

    block["grabLoLt", "zt"] = zt
    block["grabLoLt", "lt"] = lt
    block["grabLoLt", "bin1"] = bins[0]
    block["grabLoLt", "bin2"] = bins[1]
    block["grabLoLt", "bin3"] = bins[2]
    block["grabLoLt", "bin4"] = bins[3]

    return 0


def check(x,label):
    print('Check Variable: %s'%label,x)
    print('Shape of vector %s is: '%label, x.shape)

def cleanup(config):
    pass
