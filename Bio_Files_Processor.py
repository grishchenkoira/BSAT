import os
from typing import List

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> str:
    """
    Function conver multiline fasta file into fasta with name line marked by '>' and the other line with sequence.
    Results save into new folder "Converted_data". If the folder already exists, new data is written to it.
    
    :param input_fasta: Path to fasta-file with your seqs with. 
     !! It is necessary to indicate the name along with extensions (.fasta) !!
    :type input_fasta: str
    :param output_fasta: Name of new fasta-file with your seqs in one line, default value = None
    It is necessary to indicate the name along with extensions (.fasta)
    :type output_fasta: str
    :rtype: str
    :return: Script completion message 
    """
    if input_fasta.find('.fasta') == 0:
        raise ValueError(f'Wrong file format in input!')
    if os.path.exists(os.path.join('.', 'Converted_data')) == False:
        os.mkdir(os.path.join('.', 'Converted_data'))
    if output_fasta == None:
        out_name = 'one_line_' + input_fasta.strip("/")[-1]
    else:
        out_name = output_fasta
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


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: List[str], n_before: int = 1, n_after: int = 1, 
                                   output_fasta: str = None) -> str:
    '''
    Function help to search neighbours of GOI (gene of interest). Function writes neighbours of GOI in new FASTA-file 
    as: name of gene, protein sequence. Results save into new folder "Analyzed_data". If the folder already exists, 
    new data is written to it.
    
    :input_gbk: Path to gbk-file with your seqs with.
    !! It is necessary to indicate the name along with extensions (.gbk) !!
    :type input_gbk: str
    :param genes: Gene of interest names
    :type genes: List[str]
    :param n_before: number of genes before GOI (>0), default value = 1
    :type n_before: int
    :param n_after: number of genes after GOI (>0), default value = 1
    :type n_after: int
    :output_fasta: Name of FASTA-file with neighbours of GOI (names and seqs), default value = None
    :type output_fasta: str
    :rtype: str
    :return: Script completion message
    
    '''
    if input_gbk.find('.gbk') == 0:
        raise ValueError(f'Wrong file format in input!')
    if os.path.exists(os.path.join('.', 'Analyzed_data')) == False:
        os.mkdir(os.path.join('.', 'Analyzed_data'))
    genes_for_search = genes
    genes_gbk = []
    genes_for_search_in_gbk = []
    neighbour_genes = dict()
    if output_fasta == None:
        output_fasta = 'output_for_gbk.fasta'
    with open (input_gbk) as gbk:
        for line in gbk:
            if '/gene' in line:
                genes_gbk += [line.strip().split('=')[1]]
    for el in genes_for_search:
        genes_for_search_in_gbk += [gn for gn in genes_gbk if el in gn]
    for gene in genes_for_search_in_gbk:
        gene_index = genes_gbk.index(gene)
        if gene_index >= 0 and gene_index < (len(genes_gbk) - 1):
            for i in range(1, n_before + 1):
                neighbour_genes[(genes_gbk[gene_index - i])] = 0
            for i in range(1, n_after + 1):
                neighbour_genes[(genes_gbk[gene_index + i])] = 0
        else:
            for i in range(1, n_before + 1):
                neighbour_genes[(genes_gbk[gene_index - i])] = 0
            for i in range(1, n_after):
                neighbour_genes[(genes_gbk[0 + i])] = 0
    with open (input_gbk) as gbk:
        counter = 0
        gene_name = ''
        for line in gbk:
            for key in neighbour_genes:
                if key in line:
                    counter = 1
                    gene_name = key
                    continue
            if counter != 0 and '/translation' in line:
                counter = 0
                neighbour_genes[gene_name] = line.strip().split('=')[1]
                continue
    with open (os.path.join('.', 'Analyzed_data', output_fasta), mode = 'w') as fasta:
        for name, seq in neighbour_genes.items():
            name = '>'+name.replace('\"','')
            fasta.write(name + '\n')
            fasta.write(seq.replace('\"','') + '\n')
    return 'All sequences processed!'

