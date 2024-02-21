import os
from typing import Union, Tuple, List, Dict
from Biopython 


def dna_rna_analysis(*args: str, operation: str) -> Union[List[float], List[str]]:
    """
    This function performs a number of operations on DNA or RNA.
    
    Operations supported by this functions:
    -transcribe - return transcribed sequence
    -reverse - return reverse sequence
    -complement - return complement sequence
    -reverse_complement - return reverse complement sequence
    -gc_calculate - return sequence GC-content in percent
    
    :param args: nucleic acid sequence
    :type param seqs: str
    
    :param operation: type of operation required
    :type param operation: str
    
    :return: Analysis of nucleic acid sequence
    :rtype: List[str]

    :raises ValueError: if sequence not RNA or DNA, also if the operation value out of OPERATION_DICT
    """
    OPERATION_DICT = {'transcribe': na.transcribe, 
                'reverse': na.reverse,
                'reverse_complement': na.reverse_complement,
                'complement': na.complement,
                'gc_calculate': na.gc_calculate}
    analysis = []
    for seq in args:
        na.is_na(seq)
        if operation in OPERATION_DICT.keys():
            analysis.append(OPERATION_DICT[operation](seq))
        else:
            raise ValueError(f'Wrong operation!')
    return analysis


def filter_fastq(input_path: str, 
                  gc_bounds: Union[int, float, Tuple [int], Tuple [float]] = (0, 100), 
                  length_bounds: Union[int, Tuple [int]] = (0, 2**32),
                  quality_threshold: float = 0.0, filtered_file_name: Union[None, str] = None) -> Dict[str,str]:
    """
    This function help analyze a set of reads obtained from next-generation sequencing. 
    
    The function allow to filter the desired reads according to three parameters:
    GC-content, length and reading quality.
    
    :param seqs: 
    Path to the file with FASTQ-sequences in the format. 
    :type seqs: str
    
    :param gc_bounds: 
    Boundary parameters for filtering sequences by GC-content. Save only reads with a GC-content between boundaries 
    or lower than one boundary. Lower boundary cannot be less than 0 and upper boundary cannot be greater than 100. 
    gc_bounds default value is (0,100).
    :type param gc_bounds: Union[int, float, Tuple [int], Tuple [float]]
    
    :param length_bounds: 
    Boundary parameters for filtering sequences by length. Works the same as gc_bounds. Lower boundary cannot be less 
    than 0 and upper boundary cannot be greater than 2^32. length_bounds default value is (0,2^32)
    :type param length_bounds: Union[int, Tuple [int]
    
    :param quality_threshold: 
    Threshold for quality of each nucleotide in read. Quality incodes by ASCII codes. The threshold cannot be more 
    than 40. quality_threshold default value is 0 
    :type param quality_threshold: float
    
    :return: 
    New dictionaries with fastq sequence.The first one consisting of filtered fastq sequences and the other one with 
    sequences that did not pass filters.
    :rtype: Dict[str]
    
    :raises ValueError: if sequence not RNA or DNA, also if the argument values are outside the allowed ones
    """
    seqs, path_to_file = fq.read_fastq(input_path)
    file_name = path_to_file.split("/")[-1]
    if type(gc_bounds) == float or type(gc_bounds) == int:
        gc_bounds = (0,gc_bounds)
        if gc_bounds[0] < 0 or gc_bounds[1] > 100:
            raise ValueError(f'Wrong boundaries!')
    if type(length_bounds) == int:
        length_bounds = (0,length_bounds)
        if length_bounds[0] < 0 or length_bounds[1] > 2**32:
            raise ValueError(f'Wrong boundaries!')        
    if quality_threshold > 40:
        raise ValueError(f'Wrong quality threshold!')
    analysed_seq = {}
    error_seq = {}
    for seq in seqs.items():
        if fq.analyse_gc(seq[1][0]) >= gc_bounds[0] and fq.analyse_gc(seq[1][0]) <= gc_bounds[1]:
            if fq.analyse_length(seq[1][0]) >= length_bounds[0] and fq.analyse_length(seq[1][0]) <= length_bounds[1]:
                if fq.analyse_quality(seq[1][1]) > quality_threshold:
                    analysed_seq[seq[0]] = (seq[1])
                else:
                    error_seq[seq[0]] = (seq[1])
            else:
                error_seq[seq[0]] = (seq[1])
        else:
            error_seq[seq[0]] = (seq[1])
    if filtered_file_name == None:
        new_file_name = file_name 
    else:
        new_file_name = filtered_file_name
    fq.write_fastq(analysed_seq, path_to_file, new_file_name)
    return error_seq


def run_protein_analysis(*args: str) -> Union[List[str], List[float]]:
    """
    Launch desired operation with proteins sequences. Pass comma-separated sequences,
    additional argument (if certain function requires it) and specify function name you want to apply to all
    sequences.
    Pass arguments strictly in this order, otherwise it won't be parsed.
    If the input is a sequence of amino acids written in a three-letter code, then the amino acids must be separated 
    by a space. If the input is a sequence of amino acids written in a single-letter code, then the amino acids may 
    not be separated by a space.

    :param args:
    - seq (str): amino acids sequences for analysis in 1-letter or 3-letter code (as many as you wish)
    - additional arg (str): necessary parameter for certain functions (for example, specify target protein site)
    - operation name (str): specify procedure you want to apply
    :type param args: str

    :return: the result of procedure in list or str format
    :rtype: Union[List[str], List[float]
    :raises ValueError: if sequence not protein or sequence in wrong format    
    """
    # first value is function name, second is real function, third is number of function arguments
    function_names = {'change_residues_encoding': [pa.change_residues_encoding, 2],
                      'is_protein': [pa.is_protein, 1],
                      'get_seq_characteristic': [pa.get_seq_characteristic, 1],
                      'find_residue': [pa.find_residue, 2],
                      'find_site': [pa.find_site, 2],
                      'calculate_protein_mass': [pa.calculate_protein_mass, 1],
                      'calculate_average_hydrophobicity': [pa.calculate_average_hydrophobicity, 1],
                      'get_mrna': [pa.get_mrna, 1],
                      'calculate_isoelectric_point': [pa.calculate_isoelectric_point, 1],
                      'analyze_secondary_structure': [pa.analyze_secondary_structure, 1]}

    procedure = args[-1]

    processed_result = []

    seqs = [pa.change_residues_encoding(seq) for seq in args[:-1 * (function_names[procedure][1])]]
    for idx, seq in enumerate(seqs):
        if not pa.is_protein(seq):
            raise ValueError(f'Sequence {seq} is not protein!')
            continue
        if function_names[procedure][1] == 1:
            processed_result.append(function_names[procedure][0](seq))
        elif function_names[procedure][1] == 2:
            add_arg = args[-2]
            processed_result.append(function_names[procedure][0](seq, add_arg))
    if len(processed_result) == 1:
        return processed_result[0]
    return processed_result


def read_fastq(path: str) -> Union[dict, str]:
    """
    The function takes as input a file with DNA sequences in FASTQ format and creates a dictionary. The key
    is the sequence name, the value is the nucleotide sequence and the quality of reading
    :param path: path to the file with sequences in FASTQ format
    :type path: str
    :rtype: Union [dict, Str]
    :return: Dictionary with FASTQ-sequences and path to the file with sequences in FASTQ format. The last one is
    necessary for the tool to work correctly.
    """
    seq_alphabet = {'A', 'T', 'G', 'C'}
    personal_path = path
    fastq_dict = {}
    name = ''
    seq_and_quality = []
    with open (path) as seq_fastq:
        for line in seq_fastq:
            if line.startswith('@SRX079804'):
                name = line.split(' ')[0]
            elif set(line.strip()) <= seq_alphabet:
                seq_and_quality += [line.strip()]
            elif line.startswith('+'):
                continue
            else:
                seq_and_quality += [line.strip()]
                fastq_dict[name] = seq_and_quality
                name = ''
                seq_and_quality = []
    return fastq_dict, personal_path


def analyse_gc(seq: str) -> float:
    """
    Return GC-content of DNA/RNA sequence.
    
    :param seq: DNA/RNA sequence
    :type seq: str
    :rtype: float
    :return: User sequence GC-content
    """
    length = len(seq)
    gc_content = 0.0
    seq_up = seq.upper()
    c = seq_up.count("C")
    g = seq_up.count("G")
    gc_content = round(((c+g)/length*100),2)
    return gc_content


def analyse_length(seq: str) -> int:
    """
    Return length of DNA/RNA sequence.
    
    :param seq: DNA/RNA sequence
    :type seq: str
    :rtype: int
    :return: User sequence length
    """
    length = len(seq)
    return length


def analyse_quality(quality: str) -> float:
    """
    Return quality score of read, that coding by ASCII code.
    
    :param seq: quality symbols for each nucleotide
    :type seq: str
    :rtype: float
    :return: User sequence quality
    """    
    q_score = 0
    for char in quality:
        q_score += ord(char) - 33
    q_score_sum = q_score/len(quality)
    return round(q_score_sum,2)


def write_fastq(fastq_dict, personal_path, new_file_name):
    """
    The function writes filtered FASTQ reads into new file and save it into folder "fastq_filtrator_resuls".
    :param seqs: dictionary with filtered sequences
    :type seqs: dict
    :param personal_path: path to the file with original sequences in FASTQ format
    :type personal_path: str
    :param new_file_name: name for file with filtered sequences
    :type new_file_name: str
    :rtype: None
    :return: None
    """
    os.mkdir("fastq_filtrator_resuls")
    with open (personal_path) as seq_fastq:
        res = []
        count = 0
        for line in seq_fastq:
            if line.split(' ')[0] in fastq_dict:
                count += 1
                res += [line]
                continue
            count += 1
            res += [line]
            if count == 4:
                count = 0
    with open (os.path.join('fastq_filtrator_resuls', new_file_name), mode='w') as new_seq_fastq:
        for el in res:
            new_seq_fastq.write(el)


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


def transcribe(seq: str) -> str:
    """
    Function return return transcribed sequence.
    
    :param seq: DNA sequence
    :type seq: str
    :rtype: str
    :return: transcribed sequence       
    """
    nucleotide_type = set(seq.upper())
    if 'U' in nucleotide_type:
        raise ValueError(f'Sequence {seq} is not DNA!')
    rna_seq = seq.replace('T', 'U').replace('t', 'u') 
    return rna_seq


def reverse(seq: str) -> str:
    """
    Function return return reversed sequence.
    
    :param seq: DNA or RNA sequence
    :type seq: str
    :rtype: str
    :return: reversed sequence    
    """
    reverse_seq = seq[::-1]
    return reverse_seq


def complement(seq: str) -> str:
    """
    Function return return complement sequence.
    
    :param seq: DNA or RNA sequence
    :type seq: str
    :rtype: str
    :return: complement sequence   
    """
    complement_dict = {'A': 'T', 'C': 'G', 
                   'G': 'C', 'T': 'A', 'U': 'A', 'a': 't',
                   'c': 'g', 'g': 'c', 't': 'a', 'u': 'a'}
    complement_seq = []
    length = len(seq)
    for i in range (length):
        if seq[i] in complement_dict:
            complement_seq.append(complement_dict[seq[i]])
    return ''.join(complement_seq)


def reverse_complement(seq: str) -> str:
    """
    Function return return reverse complement sequence.
    
    :param seq: DNA or RNA sequence
    :type seq: str
    :rtype: str
    :return: reverse complement sequence  
    """
    seq = complement(seq)
    reverse_complement_seq = reverse(seq)
    return reverse_complement_seq


def gc_calculate(seq: str) -> float:
    """
    Function return sequence GC-content in percent.
    
    :param seq: DNA or RNA sequence
    :type seq: str
    :rtype: float
    :return: GC-contentn percent 
    """
    length = len(seq)
    gc_content = 0.0
    seq_up = seq.upper()
    c = seq_up.count("C")
    g = seq_up.count("G")
    gc_content = round(((c+g)/length*100),2)
    return gc_content


def is_na(seq: str) -> None:
    """
    Function return None if sequence is DNA or RNA.
    
    :param seq: some sequence
    :type seq: str
    :rtype: None
    :return: None   
    """
    nucleotide_set = {'A', 'G', 'C', 'T', 'U'}
    seq_nucleotide_type = set(seq.upper())
    if seq_nucleotide_type.issubset(nucleotide_set) != True:
        raise ValueError(f'Sequence {seq} is not nucleic acid!')
    elif 'T' in seq_nucleotide_type and 'U' in seq_nucleotide_type:
        raise ValueError(f'Sequence {seq} is not nucleic acid!')
