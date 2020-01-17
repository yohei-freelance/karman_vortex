#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 15:47:14 2019

@author: yohei
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 12:25:37 2019

@author: yohei
"""
import matplotlib.pyplot as plt
import numpy as np
import csv
from mpl_toolkits.mplot3d import Axes3D

p = np.zeros((401, 201))
C = []
with open('data9.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row != []:
            c = row[0].split(' ')[2]
            C.append(float(c))
    for i in range(1, 401):
        for j in range(1, 201):
            p[i-1][j-1] += C[201*(i-1)+j-2]

X = np.linspace(-10.0, 30.0, 401)
Y = np.linspace(-10.0, 10.0, 201)
print(p)
#print(p)
#print(len(X)*len(Y))
        
plt.pcolormesh(X, Y, p.T, cmap='jet')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
cont = plt.contour(X, Y, p.T, levels = [-0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20])
plt.title('Tstep=10000, time=200[s]')
plt.savefig('object_9.png')