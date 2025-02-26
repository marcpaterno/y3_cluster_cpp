def find_cosmosis_sections(block, regex="numbercounts"):
    """
    Search for sections in the CosmoSIS datablock that start with regex, e.g. 'numberCounts'.

    Parameters:
    block (cosmosis.datablock.DataBlock): The CosmoSIS datablock.

    Returns:
    list: List of section names that start with regex, e.g. 'numberCounts'.
    """
    unique_sections = list(set(section for section, _ in block.keys()))  # Extract only section names
    return [section for section in unique_sections if section.startswith(regex.lower())]
