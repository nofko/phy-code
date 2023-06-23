import matplotlib.pyplot as plt
import numpy as np
import idx2numpy
import os
import random

import _pickle as cPickle
import gzip

## ============================================================================
## MACHINE LEARNING EXAMPLES                                                 ##
##                                                                           ##
##                         MNIST RECOGNITION FROM SCRATCH                    ##
##                                                                           ##
## ============================================================================




### LOAD DATASET ###


train = idx2numpy.convert_from_file("mnist/train-images.idx3-ubyte")

train_labels = idx2numpy.convert_from_file("mnist/train-labels.idx1-ubyte")

test = idx2numpy.convert_from_file("mnist/t10k-images.idx3-ubyte")

test_labels = idx2numpy.convert_from_file("mnist/t10k-labels.idx1-ubyte")


def vectorize(x):
    e = np.zeros(10)
    e[x] = 1.0
    return e

training_images = [np.reshape(i, (784))/255 for i in train]
training_labels = [vectorize(i) for i in train_labels]

training_set = list(zip(training_images,training_labels))

test_images = [np.reshape(i, (784))/255 for i in test]
#test_labels = [vectorize(i) for i in test_labels]
test_set = list(zip(test_images,test_labels))


### DEFINITIONS ###


class myNet():

    def __init__ (self , sizes ):
    
        self.sizes = sizes
    
        self.N = len(sizes)

        self.w = [np.random.randn(y, x) for x, y in zip( sizes [: -1] , sizes [1:]) ]
        
        self.b = [np.random.randn(i) for i in sizes[1:]]

    
    def relu(self, z):
        return np.maximum(0,z)
    
    def dxrelu(self, z):
        return np.where(z >= 0, 1, 0)

    def softmax(self,x):

        return np.exp(x)/np.sum(np.exp(x))

    def sigmoid (self,z):
        
        return 1.0/(1.0+ np.exp(-z))

    def sigmoid_prime (self,z):

        return self.sigmoid (z)*(1 - self.sigmoid (z))

    def cost_derivative (self,output_activations , y):
    
        return ( output_activations - y)

    def feedforward (self , a):

        for bb, ww in zip(self.b , self.w ):
            a = self.sigmoid (np.dot(ww, a)+bb)
            
        return a

    # def backprop(self,x,y):

    #     nabla_b = [np. zeros (bb.shape) for bb in self.b]
    #     nabla_w = [np. zeros (ww.shape) for ww in self.w]

    #     activation = x
    #     activations = [x]

    #     zs = []

    #     for bb,ww in zip(self.b,self.w):

    #         z = np.dot(ww,activation)+bb
    #         zs.append(z)

    #         activation = self.sigmoid(z)
    #         activations.append(activation)

    #     delta = self.cost_derivative(activations[-1].transpose(),y).transpose()*self.sigmoid_prime(zs[-1])
        
    #     nabla_b [-1] = delta

    #     nabla_w [-1] = np.dot(delta,activations[-2].transpose())
        
    #     for l in range (2,self.N):

    #         delta = np.dot(self.w[-l+1].transpose(),delta)*self.sigmoid_prime(zs[-l])

    #         nabla_b[-l] = delta
           
    #         nabla_w[-l] = np.dot(delta,activations[-l-1].transpose())

    #     return (nabla_b, nabla_w)


    def backprop (self , x, y):

        nabla_b = [np. zeros (bb. shape ) for bb in self.b ]
        nabla_w = [np. zeros (ww. shape ) for ww in self.w ]
    
        activation = x

        activations = [x] # list to store all the activations , layer by layer
        zs = [] # list to store all the z vectors , layer by layer
        
        for bb, ww in zip(self.b , self. w ):
            z = np.dot(ww, activation )+bb

            zs. append (z)
            activation = self.sigmoid (z)
            activations . append ( activation )
        # backward pass
        
        delta = self. cost_derivative ( activations [-1], y) * self.sigmoid_prime (zs [ -1])

        nabla_b [-1] = delta
        nabla_w [-1] = np.outer(delta , activations [ -2])
        # Note that the variable l in the loop below is used a little

        for l in range (2, self.N ):

            z = zs[-l]
            sp = self.sigmoid_prime (z)
            delta = np.dot(self.w [-l+1]. transpose () , delta ) * sp
            nabla_b [-l] = delta
            nabla_w [-l] = np.outer(delta , activations [-l -1])
            
        return (nabla_b , nabla_w )

    def update(self,mini_batch,eta):

        nabla_b = [np.zeros (bb.shape ) for bb in self.b ]
        nabla_w = [np.zeros (ww.shape ) for ww in self.w ]


        for x, y in mini_batch :
            
            delta_nabla_b , delta_nabla_w = self. backprop (x, y)
            
            nabla_b = [nb+dnb for nb , dnb in zip(nabla_b , delta_nabla_b )]
            nabla_w = [nw+dnw for nw , dnw in zip(nabla_w , delta_nabla_w )]
            
        self.w = [ww +( eta/len( mini_batch ))*nw
                         for ww, nw in zip(self.w , nabla_w )]
        self.b = [bb +( eta/len( mini_batch ))*nb
                        for bb, nb in zip(self.b, nabla_b )]
        
        # nabla_b , nabla_w = self.backprop(x, y)

        # self.b = [bb-(eta)*nb/len(training_set) for bb, nb in zip(self.b, nabla_b)]
        # self.w = [ww-(eta)*nw/len(training_set) for ww, nw in zip(self.w, nabla_w )]     

        return

    
    def gradient_descent(self,training_data,epochs,mini_batch_size,eta,test_data):

        i = 0

        n = len( training_data )

        for j in range (epochs):
            
            random.shuffle (training_data)
            
            mini_batches = [
                training_data [k:k+ mini_batch_size ]
                for k in range (0, n, mini_batch_size )]

            for mini_batch in mini_batches :
                self.update( mini_batch , eta)

            print("Epoch {0}: {1}". format (
                j, self.evaluate(test_data)))

        # for x,y in training_data:

        #     self.update(x,y,eta)

        #     if (i % 10 == 0):
                
        #         #os.system('clear') 
        #         print("Iteration: ", i)
        #         print("Accuracy: ", self.evaluate(test_data))
        #         print(y)
        #         print(self.feedforward(x))
        #     i += 1

        return

    def evaluate (self, test_data):
        
        test_results = [( np.argmax (self.feedforward (x)), y)
                        for (x, y) in test_data ]

        return sum(int(x == y) for (x, y) in test_results )
    

sizes =[28*28, 30, 10]

net = myNet(sizes)
    
net.gradient_descent(training_set,30,10,3.0,test_set)

