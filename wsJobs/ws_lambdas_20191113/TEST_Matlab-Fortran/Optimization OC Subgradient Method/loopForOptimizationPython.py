#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy
import time

from nsopy.methods.subgradient import SubgradientMethod
from nsopy.loggers import GenericMethodLogger

# Author: Renzo Caballero
# KAUST: King Abdullah University of Science and Technology
# email 1: renzo.caballerorosas@kaust.edu.sa
# email 2: CaballeroRenzo@hotmail.com
# email 3: CaballeroRen@gmail.com
# Website: None
# November 2019; Last revision: 17/11/2019


# In[ ]:


def oracle(x):
    
    sft = 0;
    comm = numpy.asarray([x[0]-sft,x[1]-sft,x[2]-sft,x[3]-sft,0,0,0,0,0,0]);
    numpy.savetxt("loopForOpt.csv", comm, delimiter=",");
    done = 0;

    while done == 0:

        time.sleep(20);
        comm = numpy.loadtxt("loopForOpt.csv", delimiter =', ');

        if float(comm[9]) == 1:
            print([comm[0],comm[1],comm[2],comm[3],comm[4],comm[5],comm[6],comm[7],comm[8]]);
            done = 1;
        else:
            print('Waiting for MATLAB...');
            
    f_x_k = comm[4];
    diff_f_xk = [comm[5],comm[6],comm[7],comm[8]];
    
    return 0, f_x_k, numpy.array(diff_f_xk)

def projection_function(x_k):
    if x_k is 0:
        return numpy.array([0,])
    else:
        return numpy.maximum(x_k, 0)


# In[ ]:


method = SubgradientMethod(oracle, projection_function, dimension=4, stepsize_0=2, stepsize_rule="1/k", sense='min')
logger = GenericMethodLogger(method)

for iteration in range(100):
    method.step()

