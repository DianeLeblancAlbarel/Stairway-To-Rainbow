import numpy as np
from scipy import optimize
import math as ma
import os
import bisect


#Values to choose for generation

N = 2**24
t = 1000
steps = [750,800,850,900,950] + [t]
alpha = 0.95333
ell = 2
nhashingnodes =5
jobsize = 100


#Values to choose for attack

nbAttack =1000

#Value to choose if the best position of filters in wanted

H = 1/11547002. # Hash speed in computing nodes
F = 1/24456220. # Filtration speed
K = 12.7e-9 # Overhead constand
overheadFactor=0.88

########


column_file = "sarchedOrder" 
steps_file = "list_steps"
filters_file = "filters"


##### FUNCTIONS USED TO ESTIMATE THE BEST PLACEMENT FOR FILTERS ############

def mi (i,alpha,t):
    k = ((t+2)/alpha)-t
    return 2*N/(k+i)

def precomp_cost(filter_columns,alpha,t):
    P = 0
    has = 0
    cs = sorted(np.append(filter_columns, [t-1]))
    prevc = 0
    for c in cs:
        comp_cost = ((c-prevc) * mi(prevc,alpha,t))*H/nhashingnodes
        has += (c-prevc) * mi(prevc,alpha,t)
        filter_and_com_cost = mi(prevc,alpha,t)*F
        P += max(comp_cost, filter_and_com_cost)
        P += mi(prevc,alpha,t)*K
        prevc = c
    return P

def filtration_p_cost (t,alpha):
    mtmax = (2*N)/(t+2)
    if alpha==1.:
        m0 = N
    else:
        r = alpha/(1-alpha)
        m0 = r * mtmax
    mincost = (m0*t*H)/nhashingnodes
    for nfilters in range(1,10):
        x_init = [t*float(i)/nfilters for i in range(nfilters)]
        bnds = [(0, t-1)]*nfilters
        res = optimize.minimize(precomp_cost, x_init,(alpha,t), method='TNC', options={'disp':False, 'maxiter':10000, 'eps':0.1, 'ftol':0.0001}, bounds=bnds)
        if res.fun < mincost: 
            mincost = res.fun
            pos = (list(res.x))
            for i in range (0,len(pos)):
                pos[i]=int(round(pos[i]))
    return mincost,pos


def insert_steps (steps_list,filters_list):
    for i in range (0,len(steps_list)):
        if steps_list[i] not in filters_list:
            bisect.insort(filters_list,int(steps_list[i]))
    return filters_list

########## TABLE GENERATION ##############

def generation_table (steps_list,filters_list,):
    for i in range (0,len(steps_list)):
        if i==0:
            cmd = f'echo "{steps_list[i]}">{steps_file}'
            os.system(cmd)
        else:
            cmd = f'echo "{steps_list[i]}">>{steps_file}'
            os.system(cmd)
    for i in range (0,len(filters_list)):
        if i==0:
            cmd = f'echo "{filters_list[i]}">{filters_file}'
            os.system(cmd)
        else:
            cmd = f'echo "{filters_list[i]}">>{filters_file}'
            os.system(cmd)
    for j in range (0,ell):
        cmd = f'mpiexec -n {nhashingnodes} -hostfile ./host_file ./rainbow -a {round(alpha,3)} -N {int(ma.log2(N))} -t {t} -c {jobsize} -f {filters_file} -s {steps_file} -n {len(filters_list)} -l {j}'
        os.system(cmd)


########### ATTACK USING TABLES GENERATED ################

def attack (steps_list,filters_list,nbattack):
    for i in range (0,len(steps_list)):
        if i==0:
            cmd = f'echo "{steps_list[i]}">{steps_file}'
            os.system(cmd)
        else:
            cmd = f'echo "{steps_list[i]}">>{steps_file}'
            os.system(cmd)
    for i in range (0,len(filters_list)):
        if i==0:
            cmd = f'echo "{filters_list[i]}">{filters_file}'
            os.system(cmd)
        else:
            cmd = f'echo "{filters_list[i]}">>{filters_file}'
            os.system(cmd)
    cmd =  f'./attack -N {int(ma.log2(N))} -e {ell} -t {t} -a {round(alpha,3)}  -s {steps_file} -f {column_file} -n {nbattack}'
    os.system(cmd)


filter_list = insert_steps(steps,filtration_p_cost(t,alpha)[1])
generation_table(steps,filter_list)
attack(steps,filter_list,nbAttack)
