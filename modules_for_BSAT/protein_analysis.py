from typing import List, Union

RESIDUES_NAMES = {'ALA': 'A',
                  'ARG': 'R',
                  'ASN': 'N',
                  'ASP': 'D',
                  'CYS': 'C',
                  'GLN': 'Q',
                  'GLU': 'E',
                  'GLY': 'G',
                  'HIS': 'H',
                  'ILE': 'I',
                  'LEU': 'L',
                  'LYS': 'K',
                  'MET': 'M',
                  'PHE': 'F',
                  'PRO': 'P',
                  'SER': 'S',
                  'THR': 'T',
                  'TRP': 'W',
                  'TYR': 'Y',
                  'VAL': 'V'
                  }

RESIDUES_NAMES_THREE = {'A': 'ALA',
                      'R': 'ARG',
                      'N': 'ASN',
                      'D': 'ASP',
                      'C': 'CYS',
                      'Q': 'GLN',
                      'E': 'GLU',
                      'G': 'GLY',
                      'H': 'HIS',
                      'I': 'ILE',
                      'L': 'LEU',
                      'K': 'LYS',
                      'M': 'MET',
                      'F': 'PHE',
                      'P': 'PRO',
                      'S': 'SER',
                      'T': 'THR',
                      'W': 'TRP',
                      'Y': 'TYR',
                      'V': 'VAL'
                      }

# first value is hydrophobicity index, second is pKa (pKa1, pKa2, pKa3 respectively), third is molecular mass in Da
RESIDUES_CHARACTERISTICS = {'A': [1.8, [2.34, 9.69, 0], 89],
                            'R': [-4.5, [2.17, 9.04, 12.48], 174],
                            'N': [-3.5, [2.02, 8.80, 0], 132],
                            'D': [-3.5, [1.88, 9.60, 3.65], 133],
                            'C': [2.5, [1.96, 10.28, 8.18], 121],
                            'Q': [-3.5, [2.17, 9.13, 0], 146],
                            'E': [-3.5, [2.19, 9.67, 4.25], 147],
                            'G': [-0.4, [2.34, 9.60, 0], 75],
                            'H': [-3.2, [1.82, 9.17, 6.00], 155],
                            'I': [4.5, [2.36, 9.60, 0], 131],
                            'L': [3.8, [2.36, 9.60, 0], 131],
                            'K': [-3.9, [2.18, 8.95, 10.53], 146],
                            'M': [1.9, [2.28, 9.21, 0], 149],
                            'F': [2.8, [1.83, 9.13, 0], 165],
                            'P': [-1.6, [1.99, 10.60, 0], 115],
                            'S': [-0.8, [2.21, 9.15, 0], 105],
                            'T': [-0.7, [2.09, 9.10, 0], 119],
                            'W': [-0.9, [2.83, 9.39, 0], 204],
                            'Y': [-1.3, [2.20, 9.11, 0], 181],
                            'V': [4.2, [2.32, 9.62, 0], 117]}

# amino acid with corresponding degenerate codon/codons
AMINO_ACID_TO_MRNA = {'A': 'GCN',
                      'R': '(CGN/AGR)',
                      'N': 'AAY',
                      'D': 'GAY',
                      'C': 'UGY',
                      'Q': 'CAR',
                      'E': 'GAR',
                      'G': 'GGN',
                      'H': 'CAY',
                      'I': 'AUH',
                      'L': '(CUN/UUR)',
                      'K': 'AAR',
                      'M': 'AUG',
                      'F': 'UUY',
                      'P': 'CCN',
                      'S': '(UCN/AGY)',
                      'T': 'ACN',
                      'W': 'UGG',
                      'Y': 'UAY',
                      'V': 'GUN'}

def save_register(seq: str) -> List[int]:
    """
    Additional function for correct change_residues_encoding return
    :param seq: protein seq (str)
    :return: List of type of register character in protein seq (List[num])    
    """
    register = []
    sep_seq = []
    if ' ' in seq:
        sep_seq = seq.split()
    else:
        sep_seq = list(seq)
    for residue in sep_seq:
        if residue.isupper():
            register.append(1)
        else:
            register.append(0)
    return register


def change_residues_encoding(seq: str, encoding = 'one') -> str:
    """
    Transfer amino acids from 3-letter to 1-letter code. By default, converts all seq into 1-letter
    :param seq: protein seq (str)
    :param encoding: specify target encoding (str)
    :return: same protein seq in another encoding (str)
    """
    encode_seq = []
    encode_seq_registered = []
    if ' ' in seq:
        sep_seq = seq.upper().split()
    else:
        sep_seq = list(seq)
    residue_code_length = len(sep_seq[0])
    for residue in sep_seq:
        if len(residue) != residue_code_length:
            raise ValueError (f'Wrong sequence format in {seq}!')
        if encoding == 'one':
            if len(sep_seq[0]) == 1:
                encode_seq.append(residue)
            if len(sep_seq[0]) == 3:
                encode_seq.append(RESIDUES_NAMES[residue])
        if encoding == 'three':
            if len(sep_seq[0]) == 3:
                encode_seq.append(residue)
            if len(sep_seq[0]) == 1:
                encode_seq.append(RESIDUES_NAMES_THREE[residue])
    for residue, reg in zip(encode_seq, save_register(seq)):
        if (reg == 1):
            encode_seq_registered += residue.upper()
        elif (reg == 0):
            encode_seq_registered += residue.lower()
    if encoding == 'one':
        return ''.join(encode_seq_registered)
    if encoding == 'three':
        fin_seq = []
        for i, residue in enumerate(encode_seq_registered):
            fin_seq.append(residue)
            if ((i+1) % 3 == 0):
                fin_seq.append(' ')
        return ' '.join(fin_seq)


def is_protein(seq: str) -> bool:
    """
    Check if sequence is protein or not by identify invalid seq elements, which are not presented in dicts above.
    :param seq: protein seq in 1-letter encoding (str)
    :return: if seq is correct protein seq or not (bool)
    """
    for residue in seq.upper():
        if residue not in RESIDUES_NAMES.values():
            return False
    return True


def get_seq_characteristic(seq: str) -> dict:
    """
    Count entry of each residue type in your seq. Get description of amino acid composition.
    
    :param seq: protein seq in 1-letter encoding (str)
    :return: each residue type in seq in 3-letter code and its amount in current seq (dict)
    """
    seq = seq.upper()
    residue_count = {}
    for residue in seq:
        residue_count[[tl_code for tl_code in RESIDUES_NAMES if RESIDUES_NAMES[tl_code] == residue][0]] = 0
    for residue in seq:
        residue_count[[tl_code for tl_code in RESIDUES_NAMES if RESIDUES_NAMES[tl_code] == residue][0]] += 1
    return residue_count


def find_residue(seq: str, residue_of_interest: str) -> str:
    """
    Find all positions of certain residue in your seq.
    
    :param seq: protein seq in 1-letter encoding (str)
    :param residue_of_interest: specify the residue of interest (str)
    :return: positions of specified residue in your seq (str)
    """
    residue_of_interest = residue_of_interest.upper()
    seq = seq.upper()
    if len(residue_of_interest) == 3:
        residue_of_interest = RESIDUES_NAMES[residue_of_interest]
    residue_of_interest_position = []
    for ind, residue in enumerate(seq, 1):
        if residue == residue_of_interest:
            residue_of_interest_position.append(ind)
    return f'{residue_of_interest} positions: {residue_of_interest_position}'


def find_site(seq: str, site: str) -> str:
    """
    Find if seq contains certain site and get positions of its site.
    
    :param seq: protein seq in 1-letter encoding (str)
    :param site: specify site of interest (str)
    :return: positions of residues for each certain site in seq (str)
    """
    site = change_residues_encoding(site).upper()
    seq = seq.upper()
    if not is_protein(site):
        return f'Site {site} is not a protein!'
    if site in seq:
        site_full_position = []
        site_count = seq.count(site)
        site_start_position = [(coordinate + 1) for coordinate in range(len(seq)) if seq.startswith(site, coordinate)]
        site_end_position = [(coordinate + len(site)) for coordinate in site_start_position]
        for counter in range(len(site_start_position)):
            site_full_position.append(f'{site_start_position[counter]}:{site_end_position[counter]}')
        return f'Site entry in sequence = {site_count}. Site residues can be found at positions: {site_full_position}'
    else:
        return f'{site} site is not in sequence!'


def calculate_protein_mass(seq: str) -> float:
    """
    Get mass of residues in your seq in Da.
    
    :param seq: protein seq in 1-letter encoding (str)
    :return: mass in Da (float)
    """
    total_mass = 0
    for res in seq.upper():
        total_mass += RESIDUES_CHARACTERISTICS[res][2]
    return total_mass


def calculate_average_hydrophobicity(seq: str) -> float:
    """
    Get hydrophobicity index for protein seq as sum of index for each residue 
    in your seq divided by its length.
    
    :param seq: protein seq in 1-letter encoding (str)
    :return: average hydrophobicity (float)
    """
    sum_hydrophobicity_ind = 0
    for res in seq.upper():
        sum_hydrophobicity_ind += RESIDUES_CHARACTERISTICS[res][0]
    average_hydrophobicity = round(sum_hydrophobicity_ind / len(seq),2)
    return average_hydrophobicity


def get_mrna(seq: str) -> List[str]:
    """
    Get encoding mRNA nucleotides for your seq.
    
    :param seq: protein seq in 1-letter encoding (str)
    :return: potential encoding mRNA sequence with multiple choice for some positions (str)
    """
    mrna_seq = []
    for residue in seq.upper():
        mrna_seq.append(AMINO_ACID_TO_MRNA[residue])
    return mrna_seq


def calculate_isoelectric_point(seq: str) -> float:
    """
    Find isoelectrinc point as sum of known pI for residues in your seq.
    
    :param seq: protein seq in 1-letter encoding (str)
    :return: isoelectric point (float)
    """
    sum_pka = 0
    pka_amount = 0
    for ind, res in enumerate(seq.upper(), 1):
        if ind == 1:
            sum_pka += RESIDUES_CHARACTERISTICS[res][1][1]
            pka_amount += 1
        elif RESIDUES_CHARACTERISTICS[res][1][2] != 0:
            sum_pka += RESIDUES_CHARACTERISTICS[res][1][2]
            pka_amount += 1
        elif ind == len(seq):
            sum_pka += RESIDUES_CHARACTERISTICS[res][1][0]
            pka_amount += 1
    pi = round(sum_pka / pka_amount, 2)
    return pi


def analyze_secondary_structure(seq: str) -> List[float]:
    """
    Calculate the percentage of amino acids found in the three main
    types of protein secondary structure: beta-turn, beta-sheet and alpha-helix.
    
    :param seq: protein seq in 1-letter encoding (str)
    :return: percentage of amino acids belonging to three types of secondary structure (list[str])    
    """
    b_turn_set = {'G', 'P', 'N', 'D'}
    b_sheet_set = {'F', 'Y', 'I', 'V', 'C', 'W'}
    alpha_helix_set = {'M', 'A', 'L', 'E', 'K'}
    protein_length = len(seq)
    result = []
    counter_b_turn = 0
    counter_b_sheet = 0
    counter_alpha_helix = 0
    for residue in seq.upper():
        if residue in b_turn_set:
            counter_b_turn += 1
        elif residue in b_sheet_set:
            counter_b_sheet += 1
        elif residue in alpha_helix_set:
            counter_alpha_helix += 1
    result.append(round(counter_b_turn / protein_length * 100,2))
    result.append(round(counter_b_sheet / protein_length * 100,2))
    result.append(round(counter_alpha_helix / protein_length * 100,2))
    return result

