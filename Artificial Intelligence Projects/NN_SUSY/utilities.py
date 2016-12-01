import math
import numpy as np
from scipy.optimize import fmin,fmin_bfgs

def create_susy_dataset(filenamein,filenameout,maxsamples=100000):
    dataset = np.genfromtxt(filenamein, delimiter=',')
    y = dataset[0:maxsamples,0]
    X = dataset[0:maxsamples,1:9]
    data = np.column_stack((X,y))
    
    np.savetxt(filenameout, data, delimiter=",")

def mean(numbers):
    return sum(numbers)/float(len(numbers))
 
def stdev(numbers):
    avg = mean(numbers)
    variance = sum([pow(x-avg,2) for x in numbers])/float(len(numbers)-1)
    return math.sqrt(variance)

def calculateprob(x, mean, stdev):
    exponent = math.exp(-(math.pow(x-mean,2)/(2*math.pow(stdev,2))))
    return (1 / (math.sqrt(2*math.pi) * stdev)) * exponent
    
def sigmoid(xvec): 
    """ Compute the sigmoid function """
    # Cap -xvec, to avoid overflow
    # Undeflow is okay, since it get set to zero
    xvec[xvec < -100] = -100
    
    vecsig = 1.0 / (1.0 + np.exp(np.negative(xvec)))
    
    return vecsig #dimension of vecsig is the same as xvec

def dsigmoid(xvec):
    """ Gradient of standard sigmoid 1/(1+e^-x) """
    vecsig = sigmoid(xvec)
    return vecsig * (1 - vecsig)#the output is of the same dimension as xvec

def threshold_probs(probs):
    """ Converts probabilities to hard classification """
    classes = np.ones(len(probs),)
    classes[probs < 0.5] = 0
    return classes
            
def fmin_simple(loss, grad, initparams):
    """ Lets just call fmin_bfgs, which is better """
    return fmin_bfgs(loss,initparams,fprime=grad)                

def logsumexp(a):
    """
    Compute the log of the sum of exponentials of input elements.
    Modified scipys logsumpexp implemenation for this specific situation
    """

    awithzero = np.hstack((a, np.zeros((len(a),1))))
    maxvals = np.amax(awithzero, axis=1)
    aminusmax = np.exp((awithzero.transpose() - maxvals).transpose())

    # suppress warnings about log of zero
    with np.errstate(divide='ignore'):
        out = np.log(np.sum(aminusmax, axis=1))

    out = np.add(out,maxvals)

    return out

def f_min(loss, direction1,direction2, theta1,theta2,Xtrain,ytrain):
	pre_loss = loss(theta1,theta2)#stores the previous loss
	m = Xtrain.shape[0]
	
	grad1=direction1(theta1,theta2,Xtrain[0,:],ytrain[0])
	grad2=direction2(theta1,theta2,Xtrain[0,:],ytrain[0])
	old_grad1 = grad1
	old_grad2 = grad2
	
	theta1_opt = theta1
	theta2_opt = theta2
	min = math.pow(10,10)#min records the min seen so far
	iteration =1
	num=1 #num records the number of local minimum
	c=0.1 #learning rate
	d= 0.0#d is the coefficient for momentum
	tol = 0.01; #tolerance
	max_round = 200
	
	while iteration <= max_round:
		#print('iteration =%d')%iteration
		flag=1.0 #change the step size
		round =1 #records the rounds in adapting overshooting
		# if num == 500:
			# return {'theta1':theta1_opt,'theta2':theta2_opt}
			# break
		
		#update parameters using stochastic gradient descent
		for i in range(m):
			#print('i=%d')%i 
			grad1=direction1(theta1,theta2,Xtrain[i,:],ytrain[i])
			grad2=direction2(theta1,theta2,Xtrain[i,:],ytrain[i])
			theta1 =theta1-(c*grad1+d*old_grad1)
			theta2 =theta2-(c*grad2+d*old_grad2)
			old_grad1=grad1
			old_grad2=grad2
			
		cur_loss = loss(theta1,theta2)
		
		#print pre_loss
		#print cur_loss
		#reaching local minimum
		if np.fabs(pre_loss - cur_loss)< tol:
			return {'theta1':theta1,'theta2':theta2}
			#randomly restart
			# if loss1 < min:
				# min = loss1
				# theta1_opt = theta1
				# theta2_opt = theta2
			# num+=1
			
			# print('num = %d')%num
			# if num == 10:
				# return {'theta1':theta1_opt,'theta2':theta2_opt}
			
			# theta1 = np.random.randint(2,size=theta1.shape)
			# theta2 = np.random.randint(2,size=theta2.shape)
			# cur_loss = pre_loss(theta1,theta2)
		
		pre_loss=cur_loss
		iteration +=1
		
		if (iteration % 10)==0:
			print('iteration = %d')%iteration
		
	
	theta1_opt = theta1
	theta2_opt = theta2
	return {'theta1':theta1_opt,'theta2':theta2_opt}