__author__ = 'Paulo A. D. Moraes'
__email__ = 'paulo.a.d.moraes@ufsc.br'
__version__ = '1.2'
__copyright__ = "MIT License, Copyright (c) 2022 Paulo A. D. Moraes"


from itertools import combinations_with_replacement as combine

# Single-letter abbreviation of the natural aminoacids
AMINONAMES = ['G', 'A', 'L', 'V', 'I', 'P', 'F', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'R', 'K', 'H', 'W', 'M']

# Monoisotopic masses of the aminoacids above-mentioned, in the same order
AMINOMASS = [57.02146, 71.03711, 113.08406, 99.06841, 113.08406, 97.05276, 147.06841, 87.03203, 101.04768, 103.00919,
             163.06333, 114.04293, 128.05858, 115.02694, 129.04259, 156.10111, 128.09496, 137.05891, 186.07931, 131.04049]

WATER_MONOISOTOPE = 18.0106
CARBON_MONOISOTOPE = 12.0000
HYDROGEN_MONOISOTOPE = 1.00782
OXYGEN_MONOISOTOPE = 15.99492


def lipophylic_chains_generator(chain: int) -> tuple:
    """Calculates monoisotopic masses of beta-hydroxy carboxylic acids from a given number of carbons.
    The output was restricted to saturated, mono- and di- unsaturated chains, avoiding a lengthy output but considering
    that some unsaturated fatty acids are found in vegetable oils.
    """

    saturated_mass = (CARBON_MONOISOTOPE * chain) + (HYDROGEN_MONOISOTOPE * chain * 2) + (OXYGEN_MONOISOTOPE * 3)
    monounsaturated_mass = saturated_mass - (2 * HYDROGEN_MONOISOTOPE)
    diunsaturated_mass = saturated_mass - (4 * HYDROGEN_MONOISOTOPE)

    return saturated_mass, monounsaturated_mass, diunsaturated_mass


def combine_aminoacids(chain_size: int, macrocycle_check: bool, is_surfactin: bool) -> (list, list):
    """Combines n aminoacids in tuples(single-letter/ masses) through the imported itertools function. Then, the values
    are concatenated/ added (removing the water condensation mass).

    NOTE: as expected from the global variables, this function results are exclusively related to natural amino acids.
    """
    if macrocycle_check:
        water_mass = (WATER_MONOISOTOPE * chain_size)
    if is_surfactin:
        water_mass = water_mass + WATER_MONOISOTOPE  # Due to surfactin lactone moiety
        chain_size -= 4
    else:
        water_mass = WATER_MONOISOTOPE * (chain_size - 1)  # Normal peptide condensation, w/ free terminal NH2/CO2H
    name_combinations = list(combine(AMINONAMES, chain_size))
    mass_combinations = list(combine(AMINOMASS, chain_size))
    name_strings = [''.join(name) for name in name_combinations]
    mass_values = [(sum(mass) - water_mass) for mass in mass_combinations]
    return name_strings, mass_values


def macromolecule_mass_filter(amino_list, mass_list, search_value, restriction: float) -> (list, list):
    """After inputting a given mass ±tolerance, creates a list of filtered results by list comprehension. This list is
    then used to index-search the respective single-letter amino-acid strings.
    """

    chain_mass = [mass for mass in mass_list if (((round(mass, 1)) >= float(search_value) - restriction) and
                                                 (round(mass, 1) <= float(search_value) + restriction))]

    mass_index = (mass_list.index(i) for i in chain_mass)
    chain_names = (amino_list[i] for i in mass_index)
    return chain_mass, chain_names


def main_loop():
    """"Disclaimer: while the search function itself is quite universal (searches for masses in range), presumptions  
    were made to define the searched mass:
        - Only one fatty acid moiety is present;
        - No other derivatization (e.g. esterification) occurs in the aminoacids.
    MINOR BUG: sometimes, there's a repeated amino acid string in the results, which cause might be in the itertool func.
    These dupes must be manually removed during analysis.
    """

    def default_input() -> (bool, bool, int, float):
        """Simplifies user input for surfactin search, with an arbitrary search tolerance of 0.5 Da"""

        defaulting = bool(int(input('Use default settings? (1/0):')))
        if defaulting:
            return True, True, 7, 0.5
        else:
            is_cyclic = bool(int(input("It's a macrocyle? (1/0): ")))
            is_surfactin = bool(int(input("It's a surfactin? (1/0): ")))
            total_aminoacids = int(input("Number of aminoacids residues: "))
            restrain = float(input("Tolerance range: "))

        return is_cyclic, is_surfactin, total_aminoacids, restrain

    # Determining characteristics of the polypeptide
    cycle_bool, surfactin_bool, aminoacid_count, restriction = default_input()
    searched_mass = float(input('Peak mass (number style: 0.0000): '))

    # Setting the fatty acid moiety
    hydrophobic_chain = int(input('Alkylic chain size: '))
    carbon_chain_list = lipophylic_chains_generator(hydrophobic_chain)

    # By definition, surfactin has 7 amino acids, which 4 are almost always present.
    # Thus, instead of performing calculations on 7 A.A., these 4 are subtracted in the next steps.
    if surfactin_bool:
        searched_mass = searched_mass - (129.04259 + (2 * 113.08406) + 115.02694)

    # Stores all possible values for n-sized amino acid strings
    total_names, total_mass = combine_aminoacids(aminoacid_count, cycle_bool, surfactin_bool)

    # Iters 03 times (previously calculated alkylic chains)
    for chain_mass, chain_type in zip(carbon_chain_list, ('Saturated', 'Monounsaturated', 'Diunsaturated')):
        neat_aminoacid = searched_mass - chain_mass

        mass_list, peptide_chain = macromolecule_mass_filter(total_names, total_mass, neat_aminoacid, restriction)
        print(f'Saturation: {chain_type}')
        for mass, chain in zip(mass_list, peptide_chain):
            print(f'{mass},{chain}')
        print(f'Total of {len(mass_list)} result(s)')
    return None


if __name__ == '__main__':
    looper = True
    while looper:
        main_loop()
        looper = bool(int(input('Continue? (1/0):')))
