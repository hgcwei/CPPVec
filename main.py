import dataIO
import corpusModel
import featureUniform
import numpy as np
import cppredModel
import rnaConvertor
import svModel
import os

# ------------ only updata this section -------------------
#  1. train_fname = 'Integrated' test_fname = 'Integrated'
#  2. train_fname = 'Human'      test_fname = 'Human/Mouse/Fruit_fly/S.cerevisiae/Zebrafish

train_fname = 'Integrated'
test_fname = 'Integrated'



# method = 'cppvec', 'ovec', 'non-ovec'
method = 'non-ovec'
c = 300
g = 0.4
# ------------- only update this section --------------


svm_dir = 'svm/' + train_fname + '.' + test_fname

corpus_dir = 'corpus/' + train_fname + '.' + test_fname

if not os.path.isdir('svm'):
    os.makedirs('svm')

if not os.path.isdir(svm_dir):
    os.makedirs(svm_dir)

if not os.path.isdir('corpus'):
    os.makedirs('corpus')

if not os.path.isdir(corpus_dir):
    os.makedirs(corpus_dir)


fa1 = 'fasta/'+ train_fname +'.coding.train.fa'
fa2 = 'fasta/'+ train_fname +'.ncrna.train.fa'

fa3 = 'fasta/'+ test_fname +'.coding.test.fa'
fa4 = 'fasta/'+ test_fname +'.ncrna.test.fa'

di = dataIO.DataIO(fa1, fasta_ratio=60000, label=1, is_cut=False, cut_ratio=0.2)
train_seqs_pos, train_labels_pos = di.read_fasta()

di = dataIO.DataIO(fa2, fasta_ratio=60000, label=0, is_cut=False, cut_ratio=0.2)
train_seqs_neg, train_labels_neg = di.read_fasta()

di = dataIO.DataIO(fa3, fasta_ratio=60000, label=1, is_cut=False, cut_ratio=0.2)
test_seqs_pos, test_labels_pos = di.read_fasta()

di = dataIO.DataIO(fa4, fasta_ratio=60000, label=0, is_cut=False, cut_ratio=0.2)
test_seqs_neg, test_labels_neg = di.read_fasta()


train_seqs = train_seqs_pos + train_seqs_neg
train_labels = train_labels_pos + train_labels_neg

test_seqs = test_seqs_pos + test_seqs_neg
test_labels = test_labels_pos + test_labels_neg

seqs_all = train_seqs + test_seqs
labels_all = train_labels + test_labels

train_num = len(train_labels)

rc = rnaConvertor.RnaConvertor()
train_prots, train_orf_seqs = rc.rnas2prots(train_seqs)
test_prots, test_orf_seqs = rc.rnas2prots(test_seqs)

orf_seqs = train_orf_seqs + test_orf_seqs
prot_seqs = train_prots + test_prots

hex_fname = 'Hexamer/' +train_fname+ '_Hexamer.tsv'

cpm = cppredModel.CPPredModel(hex_fname,train_fname)
cppred_feas = cpm.calc_cppred_feas_all(seqs_all,orf_seqs,prot_seqs)


cm = corpusModel.CorpusModel(3,3,corpus_dir)

cm.corpus_from_prots(prot_seqs)

cm.train_model()

feas_vec = cm.test_model()

feas_vec = np.matrix(feas_vec)


label_svm = np.array(labels_all)

if method == 'cppvec':

    data_svm = np.hstack((feas_vec,cppred_feas))
elif method == 'ovec':
    data_svm = feas_vec
else:
    data_svm = cppred_feas


train_svm_label = label_svm[:train_num]
train_svm_data = data_svm[:train_num,:]

test_svm_label = label_svm[train_num:]
test_svm_data = data_svm[train_num:,:]


fu = featureUniform.FeatureUniform(train_svm_data,train_svm_label,svm_dir + '/' + method + '.train')

fu.libsvm_form()

fu = featureUniform.FeatureUniform(test_svm_data,test_svm_label,svm_dir + '/' + method + '.test')

fu.libsvm_form()


sm = svModel.SVModel(svm_dir,c=c,g=g,method= method)
sm.svm_scale()
sm.svm_train()
sm.svm_predict()
sm.svm_evaluate()