from Bio import SeqIO
import numpy as np
import pickle

aa_to_idx = {
    'A': 0,  'C': 1,  'D': 2,  'E': 3,  'F': 4,
    'G': 5,  'H': 6,  'I': 7,  'K': 8,  'L': 9,
    'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14,
    'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19,
    'B': 20, 'Z': 21, '-':22  # unknown or non-canonical
}

NUM_AAS = 22
MAXLEN = 2000

def one_hot_encode_22(seq):
    matrix = np.zeros((MAXLEN, NUM_AAS), dtype=np.float32)
    for i, aa in enumerate(seq[:MAXLEN]):
        idx = aa_to_idx.get(aa.upper(), aa_to_idx['-'])  # use '-' if unknown
        matrix[i][idx] = 1.0
    return matrix

def fasta_to_onehot_pickle_22(fasta_file, output_file):
    encodings = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        onehot = one_hot_encode_22(seq)
        encodings.append(onehot)
    with open(output_file, "wb") as f:
        pickle.dump(encodings, f)

    print(f"Saved {len(encodings)} one-hot encoded proteins to {output_file}")

# Example usage
fasta_file = "/Users/eelias13/gpcr_phospho/benchmarkers/DeepGOA_bm/human_GPCRs_all-protein_sequences-January_26_2024.fasta"
output_pkl = "your_onehot_22dim_sequences.pkl"
fasta_to_onehot_pickle_22(fasta_file, output_pkl)