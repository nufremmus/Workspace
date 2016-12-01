import numpy as np
import utilities as utils
import time
# http://aimotion.blogspot.com/2011/11/machine-learning-with-python-logistic.html
class Classifier:
    """
    Generic classifier interface; returns random classification
    """
    
    def __init__( self ):
        """ Params can contain any useful parameters for the algorithm """
        
    def learn(self, Xtrain, ytrain):
        """ Learns using the traindata """ # the random output method didn't really learn so it's empty
        
    def predict(self, Xtest):
        probs = np.random.rand(Xtest.shape[0])  #X_test is a matrix
        ytest = utils.threshold_probs(probs)
        return ytest
            
class NaiveBayes(Classifier):

    def __init__( self ):
        """ Params can contain any useful parameters for the algorithm """
        self.summaries = {}
        
    def learn(self, Xtrain, ytrain):
        # Separate by class
        separated = {}
        for tt in range(Xtrain.shape[0]):
            inputv = Xtrain[tt]
            outputy = ytrain[tt]
            if (outputy not in separated):
                separated[outputy] = []
            separated[outputy].append(inputv)        
        for classValue, instances in separated.iteritems():
            summ = [(utils.mean(attribute), utils.stdev(attribute)) for attribute in zip(*instances)]
            del summ[-1]
            self.summaries[classValue] = summ      
    
    def predict(self, Xtest):
        predictions = []
        for tt in range(Xtest.shape[0]):
            inputVector = Xtest[tt]
            probabilities = {}
            for classValue, classSummaries in self.summaries.iteritems():
                probabilities[classValue] = 1
                for i in range(len(classSummaries)):
                    mean, stdev = classSummaries[i]
                    x = inputVector[i]
                    probabilities[classValue] *= utils.calculateprob(x, mean, stdev)
            
            bestLabel, bestProb = None, -1
            for classValue, probability in probabilities.iteritems():
                if bestLabel is None or probability > bestProb:
                    bestProb = probability
                    bestLabel = classValue
            predictions.append(bestLabel)
        
        return predictions
    
    
class LogitReg(Classifier):

    def __init__( self ):
        self.weights = None
        
    def learn(self, Xtrain, ytrain):
        """ Learns using the traindata """
        self.weights = np.zeros((Xtrain.shape[1],1))
        lossfcn = lambda w: self.logit_cost(w, Xtrain,ytrain)['loss']
        direction = lambda w: self.logit_cost(w, Xtrain,ytrain)['gradient']
        self.weights = utils.fmin_simple(lossfcn,direction, self.weights)
        self.weights = np.reshape(self.weights,(np.shape(self.weights)[0],))

    def predict(self, Xtest):
        probs = utils.sigmoid(Xtest.dot(self.weights))
        ytest = utils.threshold_probs(probs)  
        return ytest
 
    def logit_cost(self, theta,X,y): 
        #print np.shape(X)
        tt = X.shape[0]
        theta = np.reshape(theta,(len(theta),1))
        y = np.reshape(y,(len(y),1))
        J = (1./tt) * (-np.transpose(y).dot(np.log(utils.sigmoid(X.dot(theta)))) - np.transpose(1-y).dot(np.log(1-utils.sigmoid(X.dot(theta)))))
        grad = (1./tt) * (-np.transpose(X).dot(y * np.divide(utils.sigmoid_der(X.dot(theta)),utils.sigmoid(X.dot(theta)))) + np.transpose(X).dot((1-y) * np.divide(utils.sigmoid_der(X.dot(theta)),1-utils.sigmoid(X.dot(theta)))))
     
        #print 'X^T theta %f' %X.dot(theta)[0]
        return {'loss':J, 'gradient':grad }