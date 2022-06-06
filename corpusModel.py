from gensim.models.doc2vec import Doc2Vec
from gensim.models import doc2vec
import multiprocessing

import os

class CorpusModel(object):

    def __init__(self,corpus_k,corpus_s,corpus_dir):

        self.corpus_k = corpus_k
        self.corpus_s = corpus_s

        self.corpus_save_file = corpus_dir + '/prot2vec_k' + str(self.corpus_k) + 's' + str(self.corpus_s) + '.file'
        self.corpus_save_model = corpus_dir + '/prot2vec_k' + str(self.corpus_k) + 's' + str(self.corpus_s) + '.model'

    def corpus_from_prot(self, prot):
        l = len(prot)
        inds = list(range(0, l - self.corpus_k + 1, self.corpus_s))
        x = []
        for i in range(len(inds)):
            x.append(str(prot[inds[i]:inds[i] + self.corpus_k] + ' '))

        x = ''.join(x)
        return x.rstrip()

    def corpus_from_prots(self, prots):

        f = open(self.corpus_save_file, 'w')
        cps_ls = []
        for i in range(len(prots)):
            cps = self.corpus_from_prot(prots[i])
            cps_ls.append(cps)
            f.write(cps)
            f.write('\n')
        f.close()
        return cps_ls

    def train_model(self):
        docslist = doc2vec.TaggedLineDocument(self.corpus_save_file)

        model = Doc2Vec(documents=docslist, vector_size=100,
                        alpha=0.025,
                        min_alpha=0.00025,
                        window=5,
                        min_count=1,
                        workers=multiprocessing.cpu_count(),
                        dm=1)

        model.save(self.corpus_save_model)

    def test_model(self):

        model = Doc2Vec.load(self.corpus_save_model)
        feas_ls = []
        for i in range(len(model.docvecs)):
            feas_ls.append(model.docvecs[i])

        return feas_ls








