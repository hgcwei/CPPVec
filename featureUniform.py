import os

class FeatureUniform:

    def __init__(self,features,labels,name):

        self.features = features
        self.labels = labels
        self.svm_save_file = name + '.data.nonscaled'


    def libsvm_form(self):

        m,n = self.features.shape

        f = open(self.svm_save_file,'w')
        for i in range(m):
            f.write(str(self.labels[i])+' ')
            for j in range(n):
                if j == n-1:
                    f.write(str(j + 1) + ':' + str(self.features[i, j]))
                else:
                    f.write(str(j+1)+':'+str(self.features[i,j])+' ')
            f.write('\n')




