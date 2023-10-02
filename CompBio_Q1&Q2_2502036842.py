#Q1

dna = {
    "A":"U",
    "T":"A",
    "G":"C",
    "C":"G"
}

complement = {
    "A":"T",
    "T":"A",
    "G":"C",
    "C":"G"
}

aminoacid = {
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "TGT":"C", "TGC":"C",
    "GAT":"D", "GAC":"D",
    "GAA":"E", "GAG":"E",
    "TTT": "F", "TTC":"F",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "CAT":"H", "CAC":"H",
    "ATA":"I", "ATT":"I", "ATC": "I",
    "AAA":"K", "AAG":"K",
    "TTA":"L", "TTG":"L", "CTT":"L","CTC":"L", "CTA":"L", "CTG":"L",
    "ATG":"M", "AAT":"N", "AAC":"N",
    "CCT":"P", "CCC":"P", "CCA":"P","CCG":"P",
    "CAA":"Q","CAG":"Q",
    "CGT":"R", "CGC":"R","CGA":"R","CGG":"R", "AGA":"R", "AGG":"R",
    "TCT":"S", "TCC":"S","TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "TGG":"W"
}

def convert_dnacomp(dna_s):
    comp_s = "".join([complement.get(base.upper(), "none") for base in dna_s])
    return comp_s

def convert_dna_mrna(dna_s):
    mrna_s = "".join([dna.get(base.upper(), "none") for base in dna_s])
    return mrna_s

def convert_amino(mrna):
    aminoacids = []
    start_codon = "AUG"
    end_codons = ["UAA", "UAG", "UGA"]
    i = 0

    while i < len(mrna):
        codon = mrna[i:i+3]

        if codon == start_codon:
            aminoacids.append("Met")
            i += 3
        elif codon in end_codons:
            break
        else:
            aminoacid1 = aminoacid.get(codon, "Unknown")
            aminoacids.append(aminoacid1)
            i += 3

    return aminoacids

#Q2
def codon(dna_s, amino_acid = None):
    mrna_s = "".join([dna.get(base.upper(), "none") for base in dna_s])
    amino_s = convert_amino(mrna_s)

    codon_number = {}

    for i, aminoacid in enumerate(amino_s):
        if amino_acid is None or aminoacid == amino_acid:
            codon = mrna_s[i*3:i*3+3]
            if codon in codon_number:
                codon_number[codon] += 1
            else:
                codon_number[codon] = 1

#printing
dna_sequence = "AGTTCCGATA"
#Q1
print("input DNA = ", dna_sequence)

print("Complement = ", convert_dnacomp(dna_sequence))
print("mRNA = ", convert_dna_mrna(dna_sequence))
print("Aminoacid = ", convert_amino(dna_sequence))

#Q2
aminoacidresult = None
number = codon(dna_sequence, aminoacidresult)
inputacid = "GAA"
print("\namino acid = ", inputacid)