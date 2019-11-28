#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Q1.1
import numpy as np
import matplotlib.pyplot as plt
dna1=str(np.load('dna1.npy'))
dna2=str(np.load('dna2.npy'))
def countbase1(base):
    count=0
    for i in np.arange(len(dna1)):
        if dna1[i]==base:
            count+=1
    return count
print(countbase1('A'),countbase1('T'),countbase1('C'),countbase1('G'))
print('base percentage in dna1','A',countbase1('A')/len(dna1),'T',countbase1('T')/len(dna1),'C',countbase1('C')/len(dna1),'G',countbase1('G')/len(dna1))
print(dna1.count('A'),dna1.count('T'),dna1.count('C'),dna1.count('G'))
def countbase2(base):
    count=0
    for i in np.arange(len(dna2)):
        if dna2[i]==base:
            count+=1
    return count
print(countbase2('A'),countbase2('T'),countbase2('C'),countbase2('G'))
print('base percentage in dna2','A:',countbase2('A')/len(dna2),'T:',countbase2('T')/len(dna2),'C:',countbase2('C')/len(dna2),'G:',countbase2('G')/len(dna2))
print('base counts in dna2',dna2.count('A'),dna2.count('T'),dna2.count('C'),dna2.count('G'))
#Q1.2
def countTTTT(seq):
    count=0
    for i in np.arange(len(seq)):
        if seq[i:i+4]=='TTTT':
            count+=1
           # if seq[i-1]!='T':
           #     count+=1
           # elif seq[i-4:i]=='TTTT' and seq[i-5]!='T':
           #     count+=1
           # elif seq[i-8:i]=='TTTTTTTT' and seq[i-9]!='T':
            #    count+=1
           # elif seq[i-12:i]=='TTTTTTTTTTTT' and seq[i-13]!='T':
           #     count+=1   
    print(count/len(seq),count)
countTTTT(dna1)
countTTTT(dna2)
print('Prob of TTTT in dna1:',0.29925**4,'Prob of TTTT in dna2:',0.3105**4)


# the probability of TTTT in dna1 is more similar to the predicted value given each base is independent of its precursor, while dna2 has a much higher TTTT count than expected. 
# Hense, dna2 is real and dna1 is synthetic.

# In[2]:


#2.1
def entropy(i):
    event=d[i]
    H=0
    for j in np.arange(len(d[i])):
        if event[j]==0:
            H+=0
        else:
            H-=event[j]*np.log2(event[j])
    return H
#2.2
d = np.load('annualprobs.npy')
d.shape
def entropy(i):
    event=d[i]
    H=0
    for j in np.arange(len(d[i])):
        if event[j]==0:
            H+=0
        else:
            H-=event[j]*np.log2(event[j])
    return H
for i in np.arange(6):
    ent=['0, 3.5193', '1, 3.5829', '2, 0.9602', '3, 1.9253', '4, 0.0','5, 2.6847']
    month= ('J','F','M','A','M','J','J','A','S','O','N','D')
    y_pos= np.arange(len(month))
    plt.figure()
    plt.bar(y_pos, d[i],align='center')
    plt.xticks(y_pos,month)
    plt.ylabel('prob')
    plt.title(ent[i])
    plt.show()


# In chemistry or physics, entropy denotes the probable position that a group of particles can take. In this case, it denotes the possible ways to allocate given number of occurences into events. probability distributions more evenly spread over all the events have higher entropy, vice versa. In the extreme case,if one event has prob =1 and others all have prob =0, the entropy =0 since there's only one way to allocate all the occurances. On the other hand, a uniform distribution has the highest entropy, which is simplified to -log2(1/n),(n>=0),where n is the number of events.

# In[3]:


#Q3.1
dxy = np.load ('XandY.npy')
def coeff(x,y):
    a=dxy[x]
    b=dxy[y]
    Ex=np.mean(dxy[x])
    Ey=np.mean(dxy[y])
    top=0
    de1=0
    de2=0
    for i in np.arange(len(dxy[x])):
        top+=(a[i]-Ex)*(b[i]-Ey)
        de1+=(a[i]-Ex)**2
        de2+=(b[i]-Ey)**2
    correlation=top/np.sqrt(de1*de2)
    return correlation
print('Correlation coefficient for the first set:',coeff(0,1),'for the second set:',coeff(2,3),'for the third set:', coeff(4,5),'for the forth set:', coeff(6,7))


# In[7]:


#Q3.2
xedges =np.linspace(np.min(dxy[0]),np.max(dxy[0]),11)
yedges =np.linspace(np.min(dxy[1]),np.max(dxy[1]),11)
x = dxy[0]
y = dxy[1]
H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges))
H = np.rot90(H) 
H = np.flipud(H)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(title='data dxy[0] and dxy[1]')
plt.imshow(H, interpolation='nearest', origin='low', aspect='auto',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
#dot plot
plt.figure(figsize=(5,5))
plt.plot(dxy[0],dxy[1], 'o')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
#manipulating the counts per bin
Pxy=np.flipud(H/10000)
Px=np.sum(Pxy,axis=0)
Py=np.sum(Pxy,axis=1)
Hy=-np.sum(Py*np.log2(Py))
Hx=-np.sum(Px*np.log2(Px))
Hyx=np.sum(Pxy*np.log2(np.divide(Px,Pxy,where=Pxy!=0),where=Pxy!=0)) 
MI=(Hy-Hyx)/Hy
MI, coeff(0,1)


# In[6]:


#Q3.3
def dot_plot(c,d):
    xedges =np.linspace(np.min(dxy[c]),np.max(dxy[c]),11)
    yedges =np.linspace(np.min(dxy[d]),np.max(dxy[d]),11)
    x = dxy[c]
    y = dxy[d]
    H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges))
    H = np.rot90(H) 
    H = np.flipud(H)
    #normalize Pxy
    Pxy=H/10000
    Px=np.sum(Pxy,axis=0)
    Py=np.sum(Pxy,axis=1)
    Hy=-np.sum(Py*np.log2(Py))
    Hx=-np.sum(Px*np.log2(Px))
    Hyx=np.sum(Pxy*np.log2(np.divide(Px,Pxy,where=Pxy!=0),where=Pxy!=0)) 
    #Calculated the conditional entropy
    MI=(Hy-Hyx)/Hy
    a=dxy[c]
    b=dxy[d]
    Ex=np.mean(dxy[c])
    Ey=np.mean(dxy[d])
    top=0
    de1=0
    de2=0
    for i in np.arange(len(dxy[c])):
        top+=(a[i]-Ex)*(b[i]-Ey)
        de1+=(a[i]-Ex)**2
        de2+=(b[i]-Ey)**2
    correlation=top/np.sqrt(de1*de2)
    ti=str(correlation),str(MI)
    plt.figure(figsize=(5,5))
    plt.plot(dxy[c],dxy[d], 'o')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(ti)
    plt.show()
dot_plot(0,1)
dot_plot(2,3)
dot_plot(4,5)
dot_plot(6,7)


# Correlation coefficient measures the linear correlation between two sets of data, whereas mutual information is more inclusive and measures the dependence of two sets of data, which is not necessarily linear. Correlation coefficient with larger absolute value indicates stronger linear correlation, while a value close to zero indicates little linear correlation. Hence, we observe low correlation coefficient for set 2 and 3, while set 4 has higher coefficient than Set 1 since the linear correlation is stronger in set 4. Set 3 exhibits a non-linear correlation, which is overlooked by the algorithm of correlation coefficient.
# Set 2 scores the lowest mutual information since its x and y display no dependence at all (given the x value fixed, its corresponding y values distribute almost evenly cross the y scale). Set 1 has higher mutual information as the graph indicates a linear correlation, though not as strong as that of Set 4. Set 3 and Set 4 have the highest mutual information since both exhibit obviously strong correlation, although Set 3 has a non-linear correlation yet Set 4 has a linear correlation.
