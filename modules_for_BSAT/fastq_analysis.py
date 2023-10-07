def is_nucleotide(seq: str) -> bool:
    """
    Function return 'True' if sequence is DNA or RNA.
    
    :param seq: some sequence
    :type seq: str
    :rtype: bool
    :return: 'True' if sequence is DNA or RNA   
    """
    seq_alphabet = {'A', 'T', 'U', 'G', 'C'}
    if set(seq.upper()) <= seq_alphabet:
        return True


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
