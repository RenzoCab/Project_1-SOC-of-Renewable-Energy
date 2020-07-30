#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy
import time
from datetime import datetime, timedelta
import scipy.io
import csv
import os.path

from nsopy.methods.subgradient import SubgradientMethod
from nsopy.loggers import GenericMethodLogger

# Author: Renzo Caballero
# KAUST: King Abdullah University of Science and Technology
# email 1: renzo.caballerorosas@kaust.edu.sa
# email 2: CaballeroRenzo@hotmail.com
# email 3: CaballeroRen@gmail.com
# Website: None
# November 2019; Last revision: 07/12/2019


# In[ ]:


def oracle(x):
    
    # We load data from the previous optimization:
    path = ['Simulations/table_',actualDate,'.mat']
    path = "".join(path)
    mat = scipy.io.loadmat(path)
    array = mat.get('saveTable')
    lastRow = array[-1]
    auxSft = lastRow[0:4]
    sft = numpy.array([val for val in auxSft for _ in (0, 1)])
     
    if len(x) == 1:
        x = x[0]
    
    comm = numpy.asarray([x[0]+sft[0],x[1]+sft[1],x[2]+sft[2],x[3]+sft[3], # 4
                          x[4]+sft[4],x[5]+sft[5],x[6]+sft[6],x[7]+sft[7], # 4
                          0,0,0,0,0,0,0,0,0,0]);
    numpy.savetxt("loopForOpt_2.csv", comm, delimiter=",");
    done = 0;

    while done == 0:

        time.sleep(20);
        comm = numpy.loadtxt("loopForOpt_2.csv", delimiter =', ');

        if float(comm[-1]) == 1:
            print(comm[0:8]);
            print(comm[9:17]);
            print(comm[8]);
            done = 1;
        else:
            print('Waiting for MATLAB...');
            
    f_x_k = comm[8];
    diff_f_xk = [comm[9:17]];
    
    return 0, f_x_k, numpy.array(diff_f_xk)

def projection_function(x_k):
    if x_k is 0:
        return numpy.array([0,])
    else:
        return numpy.maximum(x_k, 0)


# In[ ]:


#initialDate = "20180609"
initialDate = "20180405"
finalDate = "20191130"
dateStringFormat = '%Y%m%d'

actualDate = initialDate

while actualDate != finalDate:
    
    path = ['Simulations/table_',actualDate,'.mat']
    path = "".join(path)
    path2 = ['Simulations/table_2_',actualDate,'.mat']
    path2 = "".join(path2)
    
    # If table 1 exists and table 2 does not, then we compute table 2:
    if os.path.isfile(path) and not(os.path.isfile(path2)):
    
        print("".join(['We will compute day: ',actualDate]))
    
        # Create .csv file with the date that we will compute:
        myFile = open('date_2.csv', 'w')
        with myFile:
            writer = csv.writer(myFile)
            writer.writerows([[actualDate]])

        # Create path and load .mat file:
        path = ['Historical/dailyData/Day_',actualDate,'.mat']
        path = "".join(path)
        mat = scipy.io.loadmat(path)
        # We use only the dams costs and create the step:
        array = mat.get("Matrix")
        array = array[0]
        array = array[2]
        div = 0.05
        step = numpy.mean(array) / div
        # If all dams have value zero, we set the step to 1:
        if step == 0:
            step = 1

        method = SubgradientMethod(oracle, projection_function, dimension=8, stepsize_0=step, stepsize_rule="1/k", sense='min')
        logger = GenericMethodLogger(method)

        for iteration in range(5):
            method.step()
        
    # We increment the actualDate in one day:
    actualDate = datetime.strptime(actualDate, dateStringFormat)
    actualDate = actualDate + timedelta(days=5)
    actualDate = actualDate.strftime(dateStringFormat)

