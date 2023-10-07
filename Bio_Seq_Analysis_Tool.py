from typing import Union, Tuple, List, Set, Dict
from modules_for_BSAT import fastq_analysis as fq
from modules_for_BSAT  import dna_rna_analysis as na
from modules_for_BSAT import protein_analysis as pa


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


def analyse_fastq(seqs: dict, 
                  gc_bounds: Union[int, float, Tuple [int], Tuple [float]] = (0, 100), 
                  length_bounds: Union[int, Tuple [int]] = (0, 2**32),
                  quality_threshold: float = 0.0) -> Dict[str,str]:
    """
    This function help analyze a set of reads obtained from next-generation sequencing. 
    
    The function allow to filter the desired reads according to three parameters:
    GC-content, length and reading quality.
    
    :param seqs: 
    A dictionary consisting of fastq sequences. The structure is as follows: Key - string, sequence name. 
    The value is a tuple of two strings: sequence and quality. The sequence is RNA or DNA.  
    :type seqs: Dict[str]
    
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
        if fq.is_nucleotide(seq[1][0]) != True:
            raise TypeError(f'Wrong sequence format')        
        if fq.analyse_gc(seq[1][0]) > gc_bounds[0] and fq.analyse_gc(seq[1][0]) < gc_bounds[1]:
            if fq.analyse_length(seq[1][0]) > length_bounds[0] and fq.analyse_length(seq[1][0]) < length_bounds[1]:
                if fq.analyse_quality(seq[1][1]) > quality_threshold:
                    analysed_seq[seq[0]] = (seq[1])
                else:
                    error_seq[seq[0]] = (seq[1])
            else:
                error_seq[seq[0]] = (seq[1])
        else:
            error_seq[seq[0]] = (seq[1])
    return analysed_seq, error_seq


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

