import time
from cosmosis.datablock import option_section, names

def setup(options):
    section = option_section
    delay = options[section, "delay"]
    return int(delay)

def execute(block, config):
    delay = config
    
    # Wait for one second
    time.sleep(delay)

    # The program will pause for one second before continuing
    print("%i second has passed."%delay)

    return 0

def cleanup(config):
    pass

