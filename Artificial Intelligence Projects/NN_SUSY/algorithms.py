import numpy as np
import utilities as utils

class Classifier:
    """
    Generic classifier interface; returns random classification
    """
    
    def __init__( self, params=None ):
        """ Params can contain any useful parameters for the algorithm """
        
    def learn(self, Xtrain, ytrain):
        """ Learns using the traindata """
        
    def predict(self, Xtest):
        probs = np.random.rand(Xtest.shape[0])
        ytest = utils.threshold_probs(probs)
        return ytest
            
class NaiveBayes(Classifier):

    def __init__( self, params=None ):
        """ Params can contain any useful parameters for the algorithm """
        # From http://aimotion.blogspot.com/2011/11/machine-learning-with-python-logistic.html
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

    def __init__( self, params=None ):
        self.weights = None
        
    def learn(self, Xtrain, ytrain):
        """ Learns using the traindata """
        self.weights = np.zeros(Xtrain.shape[1],)
        lossfcn = lambda w: self.logit_cost(w, Xtrain,ytrain)
        grad = lambda w: self.logit_cost_grad(w, Xtrain,ytrain)
        self.weights = utils.fmin_simple(lossfcn, grad, self.weights)
        print self.weights
        
    def predict(self, Xtest):
        probs = utils.sigmoid(np.dot(Xtest, self.weights))
        ytest = utils.threshold_probs(probs)  
        return ytest
 
    def logit_cost(self, theta,X,y): 
        tt = X.shape[0] # number of training examples
        theta = np.reshape(theta,(len(theta),1))

        lgsum = utils.logsumexp(X.dot(theta))
        lgsumneg = utils.logsumexp(-X.dot(theta))
        J = (1./tt) * (np.transpose(y).dot(lgsumneg) + np.transpose(1-y).dot(lgsum))

        return J

    def logit_cost_grad(self, theta,X,y): 
        tt = X.shape[0] # number of training examples
        theta = np.reshape(theta,(len(theta),1))

        sig = utils.sigmoid(X.dot(theta))
        y = np.reshape(y,(len(sig),1))

        grad = np.transpose((1./tt)*np.transpose(sig - y).dot(X))
        grad = np.ndarray.flatten(grad)
        return grad    

class NeuralNet(Classifier):
    def __init__(self, params=None):
        # Number of input, hidden, and output nodes
        # Hard-coding sigmoid transfer for this test
        self.ni = params['ni']
        self.nh = params['nh']
        self.no = params['no']
        self.transfer = utils.sigmoid
        self.dtransfer = utils.dsigmoid

        # Create random {0,1} weights to define features
        self.wi = np.random.randint(2, size=(self.nh, self.ni)) #equivalent to W(1)
        self.wo = np.random.randint(2, size=(self.no, self.nh)) #equivalent to W(2)

    def learn(self, Xtrain, ytrain):
        """ Your implementation for learning the weights """
        lossfcn = lambda a,b : self.NN_cost(a,b,Xtrain,ytrain)
        derivative1 = lambda a,b,X,y : self.NN_direction(a,b,X,y)['direction1']
        derivative2 = lambda a,b,X,y : self.NN_direction(a,b,X,y)['direction2']
        optimum = utils.f_min(lossfcn,derivative1,derivative2,self.wi,self.wo,Xtrain,ytrain)
        self.wi= optimum['theta1']
        self.wo= optimum['theta2']
        #print self.wo

    def predict(self,Xtest):
        probs = np.zeros(Xtest.shape[0],)
        for samp in range(Xtest.shape[0]):
            activations = self.evaluate(Xtest[samp,:])
            probs[samp] = activations[1][0]
        ytest = utils.threshold_probs(probs)  
        return ytest
    
    def evaluate(self, inputs):
        if inputs.shape[0] != self.ni:
            raise ValueError('NeuralNet:evaluate -> Wrong number of inputs')
        
        # hidden activations
        ah = np.ones(self.nh)
        ah = self.transfer(np.dot(self.wi,inputs))  

        # output activations
        ao = np.ones(self.no)
        ao = self.transfer(np.dot(self.wo,ah))
        
        return (ah, ao)

    def NN_cost(self,theta1,theta2,X,y): #each row in X is one sample, theta1 is W1 theta2 is W2
        D = np.dot(theta1,np.transpose(X))
        #print('shape of D is %d %d')%(D.shape[0],D.shape[1])
        #print('dtype of D is %s')%D.dtype
        if D.dtype != 'float64':
            print D
            D = np.asarray(D)
            print('dtype of D is %s')%D.dtype
            print D
        C = self.transfer(D)
        B = self.transfer(np.dot(theta2,C))
        #print B
        #B = self.transfer(np.dot(theta2,self.transfer(np.dot(theta1,X.T))))
        # B[B<0.00001] = 0.00001 # to avoid overflow
        # B[B>0.99999] = 0.99999
        loss = -( np.dot(np.log(B),y)+np.dot(np.log(1-B),1-y))
        return loss

    def NN_direction(self,theta1,theta2,XX,yy): #(XX,yy) is just one sample; theta1 is W1, theta2 is W2
        XX = XX.reshape((max(XX.shape),1)) #turn XX into a column vector
        A = (theta2.T)*self.dtransfer(np.dot(theta1,XX))  # as defined in the report
        S = self.transfer(np.dot(theta2,self.transfer(np.dot(theta1,XX))))
        direction1 = -(yy-S)*np.dot(A,XX.T)
        #print direction1.shape

        direction2 = -(yy-S)*np.transpose(self.transfer(np.dot(theta1,XX)))
        
        return {'direction1':direction1, 'direction2':direction2 }
