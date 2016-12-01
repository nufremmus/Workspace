import math
import numpy as np
from scipy.optimize import fmin
import time

def mean(numbers):
    return sum(numbers)/float(len(numbers))
 
def stdev(numbers):
    avg = mean(numbers)
    variance = sum([pow(x-avg,2) for x in numbers])/float(len(numbers)-1)
    return math.sqrt(variance)

def sigmoid(X):
    """ Compute the sigmoid function """
         
    den = 1.0 + np.exp(-1.0 * X)
 
    d = 1.0 / den
 
    return d
    
def threshold_probs(probs):
    """ Converts probabilities to hard classification """
    classes = np.ones(len(probs),)
    classes[probs < 0.5] = 0
    return classes
    
def sigmoid_der(X):
    
    nn = sigmoid(X)*sigmoid(1-X)
    
    return nn

def threshold_probs(probs):
    """ Converts probabilities to hard classification """
    #print 'probs size'
    #print np.shape(probs)
    time.sleep(2)
    classes = np.ones(len(probs),)
    #print 'new probs size'
    #print np.shape(classes)
    classes[probs < 0.5] = 0
    return classes
            
def fmin_simple(loss, direction, initparams):
    #this only works for convex functions
    #initialization
    loss0 = loss(initparams)
    direction0 = direction(initparams)
    
    previous_params=initparams
    new_params = previous_params
    previous_loss = loss0
    new_loss=previous_loss
    direction1 = direction0
    i=0 #count the number of steps
    c=1
    
    step_size = math.sqrt(sum([pow(x,2) for x in direction1]))#step size is L-2 norm of direction1
    #adjust the step_size if it's too big
    if step_size > 0.01:
        direction1 = direction1/(100*step_size)
        step_size = math.sqrt(sum([pow(x,2) for x in direction1]))
    #print 'step_size'
    #print step_size
    # print np.shape(direction1)
    # print np.shape(previous_params)
    
    #when the gradient is nonzero
    #give an upper bound for the number of total runs
    while (i< 2000) & (step_size > math.pow(10,-23)): 
        new_params = previous_params - c*direction1; #gradient descend
        new_loss = loss(new_params)
        i+=1
        #print i 
        #print 'previous_loss %f' % previous_loss 
        #print 'new_loss %f' % new_loss
        
        while  new_loss >  previous_loss : #if over shoooting occurs
            new_params+=direction1 # return to previous parameters
            direction1 = direction1*1/2 #decrease the step size
            new_params = previous_params - c*direction1
            new_loss = loss(new_params)
            
            #the while loop stops when new_loss < previous_loss
        #print i 
        #print 'previous_loss %f' % previous_loss 
        #print 'new_loss %f' % new_loss
        if 0 <= (previous_loss - new_loss)<math.pow(10,-23) : #there is no more improvement 
            break
        
        else:
            previous_loss=new_loss
            previous_params=new_params
            direction1 = direction(previous_params)
            #adjust step_size if it's too big
            step_size = math.sqrt(sum([pow(x,2) for x in direction1]))#step size is L-2 norm of direction1
            if step_size > 0.01:
                direction1 = direction1/(100*step_size)
                step_size = math.sqrt(sum([pow(x,2) for x in direction1]))
            #print 'step_size'
            #print step_size
    return new_params
    """ Temporarily simply calls fmin in scipy, should be replaced by your own personally written method """
    #return fmin(loss, initparams)
