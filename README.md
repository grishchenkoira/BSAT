# BSAT - Biological Sequences Analysis Toolbox

This repository contains tools which helps you work with nucleic acid or protein sequences and with NGS-reads. It is capable to process multiple sequences, that makes analysis faster.

## Installation

To use this toolbox one need to clone repository

```shell
git clone git@github.com:grishchenkoira/BSAT.git
cd BSAT
```

### System requirements:

Key packages and programs:
- [Python](https://www.python.org/downloads/) (version >= 3.9)

## Usage

```python
# import main function
from Bio_Seq_Analysis_Tool import Bio_Seq_Analysis_Tool
```      

## Works with main functions

This section contains description of Bio_Seq_Analysis_Tool main functions.

### dna_rna_analysis(*args: str, operation: str)

This function performs a number of operations on DNA or RNA. 
Operations supported by this functions:
- transcribe - return transcribed sequence
- reverse - return reverse sequence
- complement - return complement sequence
- reverse_complement - return reverse complement sequence
- gc_calculate - return sequence GC-content in percent

**Parameters:**
- **args**: *str*

Nucleic acid sequence
- **operation**: *str*

Type of operation required

**Returns**:
- **analysis**: *str*

Analysis of nucleic acid sequence

### analyse_fastq(seqs, gc_bounds, length_bounds, quality_threshold)

Apply one of the operations described below to fastq sequences.

**Parameters:**
- **seqs**: *dict*

A dictionary consisting of fastq sequences. The structure is as follows: Key - string, sequence name. The value is a tuple of two strings: sequence and quality. The sequence is RNA or DNA.

- **gc_bounds**: *Union[int, float, Tuple [int], Tuple [float]]*

Boundary parameters for filtering sequences by GC-content. Save only reads with a GC-content between boundaries  or lower than one boundary. Lower boundary cannot be less than 0 and upper boundary cannot be greater than 100. gc_bounds default value is (0,100).

- **length_bounds** : *Union[int, Tuple [int]]*

Boundary parameters for filtering sequences by length. Works the same as gc_bounds. Lower boundary cannot be less than 0 and upper boundary cannot be greater than 2^32. length_bounds default value is (0,2^32).

- **quality_threshold** : *float*

Threshold for quality of each nucleotide in read. Quality incodes by ASCII codes. The threshold cannot be more than 40. quality_threshold default value is 0 

**Returns**:
- **analysed_seq**: *Dict[str]*

New dictionary with fastq sequence.This one consists of filtered fastq sequences and 

- **analysed_seq**: *Dict[str, str]*

New dictionary with fastq sequence. This one consists of sequences that did not pass filters.

### run_protein_analysis(*args: str)

Apply operations described below to any number of sequences with any case.

**Parameters:**
**\*args**:
- **sequences**: *str*

input coma-separated sequences in 1-letter or 3-letter code with any case (as many as you wish)
- **add_arg**: *str*

necessary parameter for certain functions (for example, specify target protein site)
- **procedure** : *str*

specify procedure you want to apply

**Returns**:
- **operation_result**: str or list

result of function work in list or str format (dependent on number of input sequences)

**Note!**
- Operation name always must be the last argument
- Additional argument must be always before operation name

## Modules

This section contains description of modules using by main functions you can find in our library.

- DNA & RNA analysis tool(#title1)
- FASTQ analysis tool(#title2)
- Amino acid sequences analysis tool(#title3)

### <a id="title1">DNA & RNA analysis</a>

This module performs a number of operations on DNA or RNA.

#### Operations

##### transcribe(seq)

Function return return transcribed sequence.

**Parameters:**
- **seq**: *str* 

DNA sequence

**Returns:**
- **rna_seq**: *str*

**Example**
```python
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
```


##### reverse(seq)

Function return return reversed sequence.

**Parameters:**
- **seq**: *str* 

DNA or RNA sequence

**Returns:**
- **reverse_seq**: *str*

**Example**
```python
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
```

##### complement(seq)

Function return return complement sequence.

**Parameters:**
- **seq**: *str* 

DNA or RNA sequence

**Returns:**
- **complement_seq**: *str*

**Example**
```python
run_dna_rna_tools('AtG', 'complement') # 'TaC'
```

##### reverse_complement(seq)

Function return return reverse complement sequence.

**Parameters:**
- **seq**: *str* 

DNA or RNA sequence

**Returns:**
- **reverse_complement_seq**: *str*

**Example**
```python
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
```

##### gc_calculate(seq)

Function return sequence GC-content in percent.

**Parameters:**
- **seq**: *str* 

DNA or RNA sequence

**Returns:**
- **gc_content**: *str*

**Example**
```python
run_dna_rna_tools ('GTAccca','gc_calculate') # '57.14'
```


### <a id="title2">FASTQ analysis tool</a>

This module contains functions for FASTQ sequnces filtration. The function allow to filter the desired reads according to three parameters: GC-content, length and reading quality.

#### Operations

##### analyse_gc(seq)

Return GC-content of DNA/RNA sequence.

**Parameters:**

- **seq**: *str*

DNA/RNA sequence

**Returns:**
- **gc_content**: *float*

##### analyse_length(seq)

Return length of DNA/RNA sequence

**Parameters:**

- **seq**: *str*

DNA/RNA sequence

**Returns:**
- **length**: *int*

##### analyse_quality(seq)

Return quality score of read, that coding by ASCII code

**Parameters:**

- **seq**: *str*

quality symbols for each nucleotide

**Returns:**
- **q_score_sum**: *float*


### <a id="title3">Amino acid sequences analysis tool</a>

This module contains functions for protein sequences analysis. You can reencode peptides sequences: 1-letter to 3-letter code and vice versa, calculate physical features, find specific sites, get predicted mRNA that coding your protein.

#### Operations

##### change_residues_encoding(seq, query='one')

Transfer amino acids from 3-letter to 1-letter code and vice versa.

**Parameters:**

- **seq**: *str*

Input protein seq in any encoding and case. If the input is a sequence of amino acids written in a three-letter code, then the amino acids must be separated by a space. If the input is a sequence of amino acids written in a single-letter code, then the amino acids may not be separated by a space.

- **encoding**: {'one', 'three'}, default: 'one'

specify target encoding

**Returns:**
- **encode_seq_registered**: *str*

same protein seq in another encoding

**Example**
```python
seq = 'AAA'
change_residues_encoding(seq, 'one', 'change_residues_encoding') # 'AAA'

seq = 'ALA ALA ALA'
change_residues_encoding(seq, 'one', 'change_residues_encoding') # 'AAA'

seq = 'AAA'
change_residues_encoding(seq, 'three', 'change_residues_encoding') # 'ALA ALA ALA'
```

##### is_protein(seq)

Check if sequence is protein or not by identify invalid seq elements, which are not presented in dicts above.

**Parameters:**
- **seq**: *str* 

input protein seq in 1-letter encoding and upper case

**Returns:**
- **verification_result**: *bool*

if seq is correct protein seq or not 

**Example**
```python
seq = 'AAA'
is_protein(seq) #True
```

##### get_seq_characteristic(seq)

Count entry of each residue type in your sequence

**Parameters:**
- **seq**: *str* 

input protein seq in 1-letter encoding and upper case

**Returns:**
- **res_count**: *dict*

each residue type in seq in 3-letter code and its amount in current seq

**Example**
```python
seq = 'AAA'
get_seq_characteristic(seq) #{'ALA': 3}
```

##### find_res(seq, res_of_interest)

Find all positions of certain residue in your seq

**Parameters:**
- **seq**: *str* 

input protein seq in 1-letter encoding and upper case
- **res_of_interest**: *str*

residue of interest in 1-letter encoding and upper case

**Returns:**
- **res_positions**: *str*

positions of specified residue in your seq

**Example**
```python
seq = 'AAA'
res = 'A'
find_res(seq, res) # 'A positions: [1, 2, 3]'
```

##### find_site(seq, site)

Find if seq contains certain site and get positions of its site

**Parameters:**
- **seq**: *str* 

input protein seq in 1-letter encoding and upper case

- **site**: *str*

specify site of interest

**Returns:**
- **site_positions**: *str*

the range of values for amino acid positions of specified site in your seq in which the last number is excluded

**Example**
```python
seq = 'AAADDDF'
site = 'AAA'
find_site(seq, site) # "Site entry in sequence = 1. Site residues can be found at positions: ['1:4']"
```

##### calculate_protein_mass(seq)

Get sum of residues masses in your seq in Da

**Parameters:**
- **seq**: *str* 

input protein seq in 1-letter encoding and upper case

**Returns:**
- **total_mass**: *float*

mass of all residues in seq in Da

**Example**
```python
seq = 'AAA'
calculate_protein_mass(seq) #267
```

##### calculate_average_hydrophobicity(seq)

Get average hydrophobicity index for protein seq as sum of index for each residue in your seq divided by its length

**Parameters:**
- **seq**: *str* 

input protein seq in 1-letter encoding and upper case

**Returns:**
- **average_hydrophobicity_idx**: *float*

average hydrophobicity index for your seq

**Example**
```python
seq = 'AAA'
calculate_average_hydrophobicity(seq) #1.8
```

##### get_mrna(seq)

Get encoding mRNA nucleotides for your seq

**Parameters:**
- **seq**: *str* 

input protein seq in 1-letter encoding and upper case

**Returns:**
- **mrna_seq**: *str*

potential encoding mRNA sequences with multiple choice for some positions

**Example**
```python
seq = 'AAA'
get_mrna(seq) # ['GCN', 'GCN', 'GCN']
```

##### calculate_isoelectric_point(seq)

Find isoelectrinc point as sum of known pI for residues in your seq

**Parameters:**
- **seq**: *str* 

input protein seq in 1-letter encoding and upper case

**Returns:**
- **pi**: *float*

isoelectric point for your seq

**Example**
```python
seq = 'AAA'
calculate_isoelectric_point(seq) # 6.01
```
##### analyze_secondary_structure(seq)

Calculates the percentage of amino acids found in the three main types of protein secondary structure: beta-turn, beta-sheet and alpha-helix in your seq

**Parameters:**
- **seq**: *str* 

input protein seq in 1-letter encoding and upper case

**Returns:**
- **result**: *list*

percentage of amino acids belonging to three types of secondary structure for seq

**Example**
```python
seq = 'AAA'
analyze_secondary_structure(seq) # [0.0, 0.0, 100.0]
``` 

## Contact

*This is the repo for the 5th homework of the BI Python 2023 course*

Author:
- *Grishenko Irina*