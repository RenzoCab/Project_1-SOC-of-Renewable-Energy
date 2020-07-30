#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy
import time
from scipy.optimize import minimize

# Author: Renzo Caballero
# KAUST: King Abdullah University of Science and Technology
# email 1: renzo.caballerorosas@kaust.edu.sa
# email 2: CaballeroRenzo@hotmail.com
# email 3: CaballeroRen@gmail.com
# Website: None
# November 2019; Last revision: 15/11/2019


# In[14]:


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
    
    return f_x_k


# In[16]:


minimize(oracle, [0,0,0,0], args=(), method='Nelder-Mead', tol=None, callback=None);


# In[ ]:




