import csv
import random
import math
import numpy as np
import time
import algorithms as algs
 
def loadcsv(filename):
	dataset = np.genfromtxt(filename, delimiter=',')
	return dataset
 
def splitdataset(dataset, trainsize=1000, testsize=100):
	randindices = np.random.randint(0,dataset.shape[0],trainsize+testsize)
	numinputs = dataset.shape[1]-1
	Xtrain = dataset[randindices[0:trainsize],0:numinputs]
	ytrain = dataset[randindices[0:trainsize],numinputs]
	Xtest = dataset[randindices[trainsize:trainsize+testsize],0:numinputs]
	ytest = dataset[randindices[trainsize:trainsize+testsize],numinputs]
					   
	return ((Xtrain,ytrain), (Xtest,ytest))
 
 
def getaccuracy(ytest, predictions):
	correct = 0
	for i in range(len(ytest)):
		if ytest[i] == predictions[i]:
			correct += 1
	return (correct/float(len(ytest))) * 100.0


if __name__ == '__main__':
	filename = 'susysubset.csv'
	f = open('myfile.txt','a')

	dataset = loadcsv(filename)
	# trainset, testset = splitdataset(dataset)
	# print('Split {0} rows into train={1} and test={2} rows').format(len(dataset), trainset[0].shape[0], testset[0].shape[0])
	
	classalgs = {'Random': algs.Classifier(),
				 'Naive Bayes': algs.NaiveBayes(),
				 'Logistic Regression': algs.LogitReg()
				 }
				 
	rh = [3,4,5,6,7] # nh = 2**i for i in rh
	round = 100;
				 
	#store the learner names in f
	for learnername, learner in classalgs.iteritems():
		print>>f,  learnername,
	
	for i in rh:
		nh = 2**i
		print>>f, ' NN ' + str(nh),

	print >>f,'\n'
	f.close()
	
	#run this for round =100 many times
	for j in range(0,round):
		tic = time.time()
		f = open('myfile.txt','a')
		print 'round = ' + str(j)
		trainset, testset = splitdataset(dataset)
		print('Split {0} rows into train={1} and test={2} rows').format(len(dataset), trainset[0].shape[0], testset[0].shape[0])
		for learnername, learner in classalgs.iteritems():
			print 'Running learner = ' + learnername
			# Train model
			learner.learn(trainset[0], trainset[1])
			# Test model
			predictions = learner.predict(testset[0])
			accuracy = getaccuracy(testset[1], predictions)
			print 'Accuracy for ' + learnername + ': ' + str(accuracy)
			#f.write( 'Accuracy for ' + learnername + ': ' + str(accuracy))
			print>>f , str(accuracy ) + ' ',
		
		# run NeuralNet for different number of hidden nodes
		for i in rh:
			nh = 2**i
			nnparams = {'ni': trainset[0].shape[1], 'nh': nh, 'no': 1}
			classalgs_NN = {'Neural Network': algs.NeuralNet(nnparams)}
			for learnername, learner in classalgs_NN.iteritems():
				print 'Running learner = ' + learnername 
				# Train model
				learner.learn(trainset[0], trainset[1])
				# Test model
				predictions = learner.predict(testset[0])
				accuracy = getaccuracy(testset[1], predictions)
				print 'Accuracy for ' + learnername + str(nh) + ': ' + str(accuracy)
				#f.write( 'Accuracy for ' + learnername + str(nh)+ ': ' + str(accuracy) )
				print>>f,  str(accuracy) + ' ',
		
		toc = time.time()
		print 'running time for one round =' + str(toc - tic)
		print >>f, '\n'
		f.close()

	
