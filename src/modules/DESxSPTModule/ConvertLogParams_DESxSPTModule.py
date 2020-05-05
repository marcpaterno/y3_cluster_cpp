from cosmosis.datablock import names as section_names


def setup(options):
    """Nothing to read in or setup to perform"""
    return 0


def execute(block, config):
    logM1 = block["abundance_integral", "log10M1"]
    block.put_double("abundance_integral", "M1", 10.**logM1)
    logMmin = block["abundance_integral", "log10Mmin"]
    block.put_double("abundance_integral", "Mmin", 10.**logMmin)
    return 0


def cleanup(config):
    pass
