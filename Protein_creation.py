def dna_to_rna(dna_sequence):
    # Create a mapping of DNA bases to their complementary RNA bases
    base_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # Translate DNA sequence to RNA sequence
    rna_sequence = ''.join([base_complement[base] for base in dna_sequence])
    
    return rna_sequence

def rna_to_mrna(rna_sequence):
    # Translate RNA sequence to mRNA sequence
    mrna_sequence = rna_sequence.replace('T', 'U')
    return mrna_sequence

def translate_mrna_to_amino_acids(mrna_sequence):
    # Create a dictionary to map mRNA codons to amino acids
    codon_to_amino_acid = {
        'UUU': 'Phe', 'UUC': 'Phe',
        'UUA': 'Leu', 'UUG': 'Leu',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        'AUG': 'Met',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'UAU': 'Tyr', 'UAC': 'Tyr',
        'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
        'CAU': 'His', 'CAC': 'His',
        'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn',
        'AAA': 'Lys', 'AAG': 'Lys',
        'GAU': 'Asp', 'GAC': 'Asp',
        'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys',
        'UGG': 'Trp',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AGU': 'Ser', 'AGC': 'Ser',
        'AGA': 'Arg', 'AGG': 'Arg',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }
    
    # Translate mRNA sequence to amino acid sequence
    amino_acid_sequence = ''
    i = 0
    while i < len(mrna_sequence):
        codon = mrna_sequence[i:i+3]
        amino_acid = codon_to_amino_acid.get(codon, 'Unknown')
        amino_acid_sequence += amino_acid + ' (' + codon + ')'
        if i + 3 < len(mrna_sequence):
            amino_acid_sequence += ' – '
        i += 3
    
    return amino_acid_sequence

def count_codons(mrna_sequence):
    # Initialize a dictionary to count codons
    codon_counts = {}

    i = 0
    while i < len(mrna_sequence):
        codon = mrna_sequence[i:i+3]
        codon_counts[codon] = codon_counts.get(codon, 0) + 1
        i += 3

    return codon_counts

# Example DNA sequence
input_dna_sequence = 'TTACGATTA'

# Translate DNA to RNA
rna_sequence = dna_to_rna(input_dna_sequence)
print('Input DNA =', input_dna_sequence)
print('Complement (RNA) =', rna_sequence)

# Translate RNA to mRNA
mrna_sequence = rna_to_mrna(rna_sequence)
print('mRNA =', mrna_sequence)

# Count codons in the mRNA sequence
codon_counts = count_codons(mrna_sequence)

# Display codon frequencies
for codon, count in codon_counts.items():
    print(f'{codon} = {count}')

# Translate mRNA to amino acids
amino_acid_sequence = translate_mrna_to_amino_acids(mrna_sequence)
print('Aminoacid =', amino_acid_sequence)
