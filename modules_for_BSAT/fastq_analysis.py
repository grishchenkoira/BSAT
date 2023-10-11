from typing import Union
import os

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
    personal_path = path
    fastq_dict = {}
    name = ''
    seq_and_quality = []
    with open (path) as seq_fastq:
        for line in seq_fastq:
            if line.startswith('@SRX079804'):
                name = line.split(' ')[0]
            elif set(line.strip()) <= SEQ_ALPHABET:
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


def write_fastq(seqs, personal_path, new_file_name):
    """
    The function writes filtered FASTQ reads into new file and save it into folder "fastq_filtrator_resuls".
    :param seqs: path to the file with sequences in FASTQ format
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

