import os
from typing import List

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> str:
    """
    Function conver multiline fasta file into fasta with name line marked by '>' and the other line with sequence.
    Results save into new folder "Converted_data". If the folder already exists, new data is written to it.
    :param input_fasta: Path to fasta-file with your seqs with. 
    It is necessary to indicate the name along with extensions (.fasta)
    :type input_fasta: str
    :param output_fasta: Name of new fasta-file with your seqs in one line, default value is None
    It is necessary to indicate the name along with extensions (.fasta)
    :type output_fasta: str
    :rtype: str
    :return: Script completion message 
    """
    if os.path.exists(os.path.join('.', 'Converted_data')) == False:
        os.mkdir(os.path.join('.', 'Converted_data'))
    if output_fasta == None:
        out_name = 'one_line_' + input_fasta
    else:
        out_path = os.path.abspath(output_fasta)
    counter = 0
    with open (input_fasta) as seq_fasta:
        for line in seq_fasta:
            if counter == 0:
                if line.startswith('>'):
                    with open (os.path.join('.', 'Converted_data', out_name), mode = 'w') as new_seq_fasta:
                        new_seq_fasta.write(line)
                    counter += 1
                    continue
                else:
                    raise ValueError('Wrong start of FASTA')
            if counter != 0:
                if line.startswith('>'):
                    with open (os.path.join('.', 'Converted_data', out_name), mode = 'a') as new_seq_fasta:
                        new_seq_fasta.write('\n' + line)
                else:
                    with open (os.path.join('.', 'Converted_data', out_name), mode = 'a') as new_seq_fasta:
                        new_seq_fasta.write(line.strip('\n'))
                    continue
        return 'All sequences processed!'