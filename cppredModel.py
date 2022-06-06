from feamodule import CTD
from feamodule import ProtParam as PP
from feamodule import ORF_length as ole
from feamodule import fickett
from feamodule  import FrameKmer
import numpy as np


class CPPredModel:

    def __init__(self, hex_file, species):

        self.hex_file = hex_file
        self.species = species

    def coding_nocoding_potential(self):
        coding = {}
        noncoding = {}
        for line in open(self.hex_file).readlines():
            fields = line.split()
            if fields[0] == 'hexamer': continue
            coding[fields[0]] = float(fields[1])
            noncoding[fields[0]] = float(fields[2])
        return coding, noncoding

    def calc_cppred_feas_all(self,rna_seqs,orf_seqs,prot_seqs):
        feas_ls = []
        for i in range(len(rna_seqs)):
            rna_seq = rna_seqs[i]
            orf_seq = orf_seqs[i]
            prot_seq = prot_seqs[i]
            feas = self.calc_cppred_feas(rna_seq,orf_seq,prot_seq)
            feas_ls.append(feas)
        return np.array(feas_ls)

    def calc_cppred_feas(self,rna_seq,orf_seq,prot_seq):

        coding, noncoding = self.coding_nocoding_potential()

        A, T, G, C, AT, AG, AC, TG, TC, GC, A0, A1, A2, A3, A4, T0, T1, T2, T3, T4, G0, G1, G2, G3, G4, C0, C1, C2, C3, C4 = CTD.CTD(rna_seq)
        insta_fe, PI_fe, gra_fe = PP.param(prot_seq)
        fickett_fe = fickett.fickett_value(rna_seq)
        hexamer = FrameKmer.kmer_ratio(orf_seq, 6, 3, coding, noncoding)
        Len, Cov, inte_fe = ole.len_cov(rna_seq)

        if self.species == "Human":
            # return [inte_fe, Cov, insta_fe, PI_fe, Len, gra_fe, hexamer, fickett_fe]
            return [inte_fe, Cov, insta_fe, T2, C0, PI_fe, Len, AC, T0, G0, C2, A4, G2, TG, A0, TC, G1, C3, T3, A1, GC, T1, G4, C1, G3, A3, gra_fe, hexamer, C4, AG, fickett_fe, A2, T4, C, G, A, T, AT]

        else:
            # return [Cov,inte_fe,insta_fe,Len,fickett_fe,PI_fe,hexamer, gra_fe]
            return [Cov,inte_fe,GC,insta_fe,Len,T0,fickett_fe,G2,C3,PI_fe,A3,C1,G3,hexamer,TG,G1,TC,A0,A1,AC,C2,G0,T4,C0,A4,G,A2,T,T3,G4,C4,gra_fe,T2,AG,AT,T1,A,C]

