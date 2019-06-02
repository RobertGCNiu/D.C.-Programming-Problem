__author__ = 'Niu Guanchong'
import cvxpy as cvx
import numpy as np
from dccp.problem import is_dccp
import matplotlib.pyplot as plt
import dccp
import math

##follow the numerical design in Fast global optimal power allocation in wiereless networks by local D.C. Programming
H = np.matrix([[0.431,0.0002,0.2605,0.0039],[0.0002,0.3018,0.0008,0.0054],[0.0129,0.0005,0.4266,0.1007],[0.0011,0.0031,0.0099,0.0634]])

M = 4
ri = 0.0
w = [1.0/6,1.0/6,1.0/3,1.0/3]

#p_max = [1.0,1.0,1.0,1.0]
p_max = [0.7,0.8,0.9,1.0]
sigmal2 = 0.0001
p = cvx.Variable(M)
t = cvx.Variable(1)
constr = []
# The constraints: (3c) and (5)
for i in range(M):
    laterValue = 0
    for j in range(M):
        if j != i:
            laterValue = laterValue + H[j,i]*p[j]
    constr.append(H[i,i]*p[i]+(1-pow(2,ri))*(laterValue+sigmal2)>=0) ## (3c)
    constr.append(p[i]<=p_max[i]) ## (5)
    constr.append(0<=p[i])
    #constr.append(p[[1,2]]==0)

## The objective function is reformulated as f(x)-t: t = g(x);
f_p = 0
g_p = 0
for i in range(M):
    laterValue_all = 0
    laterValue = 0
    for j in range(M):
        laterValue_all = laterValue_all + H[j,i]*p[j]
        if j != i:
            laterValue = laterValue + H[j,i]*p[j]
    f_p = f_p + w[i]*cvx.log(sigmal2 + laterValue_all)/cvx.log(2).value
    g_p = g_p + w[i]*cvx.log(sigmal2 + laterValue)/cvx.log(2).value


constr.append(t==g_p)

prob_value = cvx.Maximize(f_p-t)
prob = cvx.Problem(prob_value, constr)
#print is_dccp(prob)

prob.solve(method='dccp',tau=1.0)


print p.value
print prob_value.value


