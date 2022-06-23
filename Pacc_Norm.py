import numpy as np
import matplotlib.pyplot as plt
import random
from math import exp

num_sims = 100
mut = 0.05
dose = .7

g = 0.02
K = 100
lamn = 1
b = 1
r = 0.6
cn = 0.4
a = 0.2
rho = 0.7
d = 0.01
bred = 0.05

tend = 900

norm_macro = []
pacc_macro = []
v1_macro = []
v2_macro = []
t_macro = []
extinct = []

time1 = 100
time2 = 400
time3 = 500 #500 & 800
time4 = 800

def drug1(v1,m1): #NOTE that for now, we're running this with lamn
        return m1/(lamn+b*v1)

def drug2(v2,m2):
        return m2/(lamn+b*v2)

def run_this():
        normal = [50]
        pacc = [0]
        strat1 = [0]
        strat2 = [0]
        t = [0]

        while t[-1] < tend:
            
            if t[-1]>time4:
                m1=m2=0
            elif t[-1]>time3:
                m1 = dose
                m2 = 0
            elif t[-1]>time2:
                m1=m2=0
            elif t[-1]>time1:
                m1=dose
                m2=0
            else:
                m1=m2=0
                
            current_n = normal[-1]
            current_pacc = pacc[-1]
            current_strat1 = strat1[-1]
            current_strat2 = strat2[-1]

            #broken into birth, death, switching


            rates = [current_n*r, current_n*(m1/(lamn+b*current_strat1)+m2/(lamn+b*current_strat2)+r*((current_pacc+current_n)/K)),current_n*(cn*m1/(lamn+b*current_strat1)+cn*m2/(lamn+b*current_strat2)+g),
                     current_pacc*d,current_pacc*a]

            rate_sum = sum(rates)

            if rate_sum == 0:
                    extinct.append(t[-1])
                    break

            tau = np.random.exponential(scale=1/rate_sum)

            t.append(t[-1] + tau)

            rand = random.uniform(0,1)

            #Normal cell division event
            if rand * rate_sum <= rates[0]:
                    pacc.append(pacc[-1])
                    if random.uniform(0,1)>mut:
                            strat1.append(strat1[-1])
                            strat2.append(strat2[-1])
                            normal.append(normal[-1] + 1)
                    else:
                        if m1==dose: #change if giving concurrently
                            strat2.append(strat2[-1])
                            new_strat = strat1[-1]+np.random.normal(0,0.01)
                            if drug1(new_strat,m1)<drug1(strat1[-1],m1):
                                    strat1.append(new_strat)
                                    normal.append(normal[-1] + 1)
                            elif drug1(new_strat,m1)==drug1(strat1[-1],m1):
                                    strat1.append(strat1[-1])
                                    normal.append(normal[-1] + 1)
                            else:
                                    strat1.append(strat1[-1])
                                    normal.append(normal[-1])
                        elif m2==dose:
                            strat1.append(strat1[-1])
                            new_strat = strat2[-1]+np.random.normal(0,0.01)
                            if drug2(new_strat,m2)<drug2(strat2[-1],m2):
                                    strat2.append(new_strat)
                                    normal.append(normal[-1] + 1)
                            elif drug2(new_strat,m2)==drug2(strat2[-1],m2):
                                    strat2.append(strat2[-1])
                                    normal.append(normal[-1] + 1)
                            else:
                                    strat2.append(strat2[-1])
                                    normal.append(normal[-1])
                        else:
                            strat1.append(strat1[-1])
                            strat2.append(strat2[-1])
                            normal.append(normal[-1]+1)

            #Normal cell death event
            elif rand * rate_sum > rates[0] and rand * rate_sum <= sum(rates[:2]):
                    normal.append(normal[-1] - 1)
                    pacc.append(pacc[-1])
                    strat1.append(strat1[-1])
                    strat2.append(strat2[-1])

            #Normal switch event
            elif rand * rate_sum > sum(rates[:2]) and rand * rate_sum <= sum(rates[:3]):
                    normal.append(normal[-1]-1)
                    if random.uniform(0,1)<rho:
                            pacc.append(pacc[-1] + 1)
                    else:
                            pacc.append(pacc[-1])

                    strat1.append(strat1[-1])
                    strat2.append(strat2[-1])
                    
            #PACC death event
            elif rand * rate_sum > sum(rates[:3]) and rand * rate_sum <= sum(rates[:4]):
                    normal.append(normal[-1])
                    pacc.append(pacc[-1]-1)
                    strat1.append(strat1[-1])
                    strat2.append(strat2[-1])

            #PACC switch event
            elif rand * rate_sum > sum(rates[:4]) and rand * rate_sum <= sum(rates[:5]):
                    pacc.append(pacc[-1]-1)
                    if random.uniform(0,1)>mut: #using same mutation rate, but larger breadth
                            strat1.append(strat1[-1])
                            strat2.append(strat2[-1])        
                            normal.append(normal[-1]+2)
                    else:
                        if m1==dose: #change if giving concurrently
                            strat2.append(strat2[-1])
                            new_strat = strat1[-1]+np.random.normal(0,bred)
                            if drug1(new_strat,m1)<drug1(strat1[-1],m1):
                                    strat1.append(new_strat)
                                    normal.append(normal[-1] + 2)
                            elif drug1(new_strat,m1)==drug1(strat1[-1],m1):
                                    strat1.append(strat1[-1])
                                    normal.append(normal[-1] + 2)
                            else:
                                    strat1.append(strat1[-1])
                                    normal.append(normal[-1])
                        elif m2==dose:
                            strat1.append(strat1[-1])
                            new_strat = strat2[-1]+np.random.normal(0,bred)
                            if drug2(new_strat,m2)<drug2(strat2[-1],m2):
                                    strat2.append(new_strat)
                                    normal.append(normal[-1] + 2)
                            elif drug2(new_strat,m2)==drug2(strat2[-1],m2):
                                    strat2.append(strat2[-1])
                                    normal.append(normal[-1] + 2)
                            else:
                                    strat2.append(strat2[-1])
                                    normal.append(normal[-1])
                        else:
                            strat1.append(strat1[-1])
                            strat2.append(strat2[-1])
                            normal.append(normal[-1]+2)
        
            
        norm_macro.append(normal)
        pacc_macro.append(pacc)
        v1_macro.append(strat1)
        v2_macro.append(strat2)
        t_macro.append(t)

for i in range(num_sims):
       run_this()
        
'''for i in range(10):
    for i in range(num_sims):
        run_this()
    print(len(extinct))
    extinct=[]
    norm_macro = []
    pacc_macro = []
    v1_macro = []
    v2_macro = []
    t_macro = []'''

lwi = 0.3

extinct_events = len(extinct)
#avg_extinct = np.mean(extinct)
#std_extinct = np.std(extinct)

print(extinct_events)
#print(avg_extinct)
#print(std_extinct)

plt.figure()
plt.subplot(211)
plt.title('Low Mutational Breadth: Same Tx')
plt.plot(t_macro[0],norm_macro[0], label='2N+', c='r',lw=lwi)
plt.plot(t_macro[0],pacc_macro[0], label='PACC', c='b',lw=lwi)

for i in range(1,num_sims):
    plt.plot(t_macro[i],norm_macro[i], c='r',lw=lwi)
    plt.plot(t_macro[i],pacc_macro[i], c='b',lw=lwi)
plt.xlim([0,tend])
plt.legend(loc='best')
ax = plt.gca()
ax.axvspan(time1,time2, facecolor='gold') #change this back
ax.axvspan(time3, time4, facecolor='gold') #pink
plt.ylabel('Population Size')
plt.subplot(212)
plt.plot(t_macro[0],v1_macro[0], c='k',lw=lwi,label='Drug 1 Resist')
plt.plot(t_macro[0],v2_macro[0], c='m',lw=lwi,label='Drug 2 Resist')

for i in range(1, num_sims):
    plt.plot(t_macro[i],v1_macro[i], c='k',lw=lwi)
    plt.plot(t_macro[i],v2_macro[i], c='m',lw=lwi)
plt.xlim([0,tend])
plt.legend()
ax = plt.gca()
ax.axvspan(time1,time2, facecolor='gold')
ax.axvspan(time3, time4, facecolor='gold')
plt.ylabel('Drug Resistance')
plt.tight_layout()
plt.show()
