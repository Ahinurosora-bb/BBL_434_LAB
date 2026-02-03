import sys
import os
import random

class PlasmidAssembler:
    def __init__(self):
        # Default BHR genes necessary for replication 
        self.BHR_DEFAULT_BACKBONE = "ATGC_REPA_REPB_REPC_SEQUENCE" #following the research paper given

        self.RESTRICTION_SITES = {
            "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
            "PstI": "CTGCAG", "SphI": "GCATGC", "SalI": "GTCGAC",
            "XbaI": "TCTAGA", "KpnI": "GGTACC", "SacI": "GAGCTC", "SmaI": "CCCGGG"
        }

        self.GENE_SEQUENCES = {
            "AmpR": "ATGAAAGCGTTGCTGATGCTGCTGCTAA", 
            "KanR (nptII/aphA)": "ATGAGCCATATTCAACGGGAAACGCTAA", 
            "TetR (tetA/tetR)": "ATGTTGACCTGCTGCTGACGATGACTAA", 
            "lacZ_alpha": "ATGACCATGATTACGCCAAGCTGCTAA"
        }

    def at_fraction(self, seq):
        return (seq.count('A') + seq.count('T')) / len(seq)

    def gc_skew(self, seq):
        g, c = seq.count('G'), seq.count('C')
        return (g - c) / (g + c) if (g + c) != 0 else 0

    def get_ori_sequence(self, dna_seq):
        """Your logic to identify the high-score 500bp segment[cite: 37]."""
        WINDOW = 500
        best_score, best_pos = -1, 0
        for i in range(0, len(dna_seq) - 1000, 100):
            window = dna_seq[i:i+1000]
            score = self.at_fraction(window) - abs(self.gc_skew(window))
            if score > best_score:
                best_score, best_pos = score, i
        return dna_seq[best_pos : best_pos + WINDOW + 1]

    def generate(self, fastafile, designfile):
        
        # 2. Extract DNA from pUC19
        with open(fastafile, 'r') as f:
            lines = f.readlines()
            dna_seq = ''.join([line.strip() for line in lines[1:]])

        # 3. Find the segment from the sequence
        ori_seq = self.get_ori_sequence(dna_seq)

        # 4. Parse design and join sequences
        with open(designfile, 'r') as f:
            for line in f:
                line = line.strip().split(', ')
                feature = line[0].replace('_site', '').replace('_gene', '')
                
                if feature in self.RESTRICTION_SITES:
                    ori_seq += self.RESTRICTION_SITES[feature]
                    # Spacer needed
                    ori_seq += ''.join(random.choice('ATGC') for _ in range(6))
                elif feature in self.GENE_SEQUENCES:
                    ori_seq += self.GENE_SEQUENCES[feature]

        final_output = self.BHR_DEFAULT_BACKBONE + ori_seq

        # Output.fa
        os.makedirs("Output", exist_ok=True)
        with open("Output/Output.fa", 'w') as f:
            f.write(">designed_plasmid\n")
            f.write(final_output + "\n")
        
        print("Success: FASTA file generated in Output/Output.fa")

if __name__ == "__main__":
    if len(sys.argv) == 3:
        app = PlasmidAssembler()
        app.generate(sys.argv[1], sys.argv[2])