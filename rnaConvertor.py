from Bio.Seq import Seq
import re
from feamodule import ORF

class RnaConvertor:

    def orf_translate(self,cds):
        return Seq(cds).translate()

    def rna_longest_orf(self,rna):
        strinfoAmbiguous = re.compile("X|B|Z|J|U", re.I)
        ptU = re.compile("U", re.I)
        rna_seq = ptU.sub("T", str(rna).strip())
        rna_seq = rna_seq.upper()
        orf_length, orf_integrity, orf_seq = ORF.ExtractORF(rna_seq).longest_ORF(start=['ATG'], stop=['TAA', 'TAG', 'TGA'])

        prot_seq = self.orf_translate(orf_seq)
        # l = len(prot_seq.strip("*"))
        prot_seq = strinfoAmbiguous.sub("", str(prot_seq))

        return prot_seq,orf_seq,orf_length,orf_integrity

    def rna2prot(self,rna):

        prot_seq, orf_seq, _, _ = self.rna_longest_orf(rna)
        return str(prot_seq.strip("*")), orf_seq

    def rnas2prots(self,rna_ls):
        prot_seq_ls = []
        orf_seq_ls = []
        for i in range(len(rna_ls)):
            rna = rna_ls[i]
            prot_seq,orf_seq = self.rna2prot(rna)
            prot_seq_ls.append(prot_seq)
            orf_seq_ls.append(orf_seq)
        return  prot_seq_ls, orf_seq_ls
