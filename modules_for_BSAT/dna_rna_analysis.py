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
