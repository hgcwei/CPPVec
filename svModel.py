import os
from sklearn.metrics import auc,roc_curve
from sklearn.metrics import f1_score,matthews_corrcoef,precision_recall_curve
from sklearn.metrics import precision_score
from sklearn.metrics import accuracy_score
import numpy as np


class SVModel:

    def __init__(self, save_dir, c, g, method):
        os.environ["PATH"] += os.pathsep + os.getcwd() + '/libsvm/windows/'

        self.range = save_dir + '/range.file'
        self.train_file = save_dir + '/' + method + '.train.data.nonscaled'
        self.scaled_train_file = save_dir + '/' + method + '.train.data.scaled'
        self.test_file = save_dir + '/' + method + '.test.data.nonscaled'
        self.scaled_test_file = save_dir + '/' + method + 'test.data.scaled'
        self.model_save_file = save_dir + '/' + method + '.model.file'
        self.predict_save_file = save_dir + '/' + method + '.predict.file'

        self.results_save_file = save_dir + '/' + method + '.results.file'
        self.c = c
        self.g = g

    def get_test_y(self):
        test_y = []
        f = open(self.test_file,'r')
        for line in f.readlines():
            test_y.append(int(line[0]))
        return test_y

    def svm_scale(self):
        cmd0 = 'svm-scale -l 0 -u 1 -s '+ self.range + ' '+ self.train_file + ' > ' + self.scaled_train_file
        cmd1 = 'svm-scale -r '+ self.range + ' '+ self.test_file + ' > ' + self.scaled_test_file
        os.system(cmd0)
        os.system(cmd1)

    def svm_train(self):
        cmd = 'svm-train -c '+ str(self.c) +' -g '+ str(self.g) +' -b 1 '+ self.scaled_train_file + ' ' + self.model_save_file
        os.system(cmd)


    def svm_predict(self):
        cmd = 'svm-predict -b 1 '+ self.scaled_test_file + ' '+ self.model_save_file + ' '+ self.predict_save_file
        os.system(cmd)
        with open(self.predict_save_file,mode='r') as f:
            line = f.readlines()
            try:
                line = line[1:]
                f = open(self.predict_save_file,mode='w')
                f.writelines(line)
                f.close()
            except:
                pass


    def eval_perf(self,label, scores):
        fpr, tpr, thresholds = roc_curve(label, scores)
        pr, re, thresholds2 = precision_recall_curve(label, scores)
        for i in range(len(thresholds)):
            if thresholds[i] < 0.5:
                return  tpr[i], 1-fpr[i], auc(fpr, tpr), auc(re, pr)

    def svm_evaluate(self):

        test_y_predict = np.loadtxt(self.predict_save_file,delimiter=' ')
        test_y = self.get_test_y()
        sn,sp,auROC,auPRC = self.eval_perf(test_y,test_y_predict[:,1])

        pre = precision_score(test_y,test_y_predict[:,0])
        acc = accuracy_score(test_y,test_y_predict[:,0])
        f1 = f1_score(test_y, test_y_predict[:,0])
        mcc = matthews_corrcoef(test_y, test_y_predict[:,0])


        print(sp,sn, pre,acc,f1,auROC,mcc)
        result = str(sp) + ' ' + str(sn) + ' ' + str(pre) + ' ' + str(acc) + ' ' + str(f1) + ' ' + str(auROC) + ' ' + str(mcc) + ' -c ' + str(self.c) + ' -g ' + str(self.g)
        with open(self.results_save_file,'a') as f:
            f.write(result)
            f.write('\n')
