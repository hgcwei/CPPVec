from Bio import SeqIO

class DataIO:


    def __init__(self,fasta_dir,fasta_ratio,label=-1,is_cut=False,cut_ratio=0.2):
        self.fasta_dir = fasta_dir
        self.label = label
        self.is_cut = is_cut
        self.cut_ratio = cut_ratio
        self.fasta_ratio = fasta_ratio

    def read_fasta(self):
        seqs = []
        labels = []
        for record in SeqIO.parse(self.fasta_dir, "fasta"):
            # x= ''.join()
            seqs.append(str(record.seq))
            labels.append(self.label)
            if len(labels) > self.fasta_ratio:
                break
        if not self.is_cut:
            return seqs,labels
        else:
            c = int(len(labels)*(1-self.cut_ratio))
            return seqs[:c],labels[:c],seqs[c:],labels[c:]









