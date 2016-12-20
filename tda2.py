def c(s,p):
    if len(s)==1:
        return[]
    else:
        l=[]
        for i in range(p,len(s)):
            n=len(s)-i-1
            ss=s.copy()
            ss.remove(s[n])
            l.append(ss)
            l=l+c(ss,i)
        return l

def isLess(s,t):
    n=min(len(s),len(t))
    for i in range(0,n):
        if s[i]<t[i]:
            return True
        elif s[i]>t[i]:
            return False
    if len(s)<len(t):
        return True
    return False

import random

def MakeLex(l):
    if len(l)<=1:
        return l
    p=random.randrange(0,len(l))
    piv=l[p]
    a=[]
    b=[]
    for i in range(0,len(l)):
        y = [len(l[i])]+l[i]
        z =[len(piv)]+piv
        if isLess(y,z) and i!=p:
            a.append(l[i])
        if isLess(z,y) and i!=p:
            b.append(l[i])
    l=MakeLex(a)+[piv]+MakeLex(b)
    return l

#import sympy

from sympy import *

def choose(n,k):
    if 0<=k<=n:
        num=1
        den=1
        for t in range(1,min(k,n-k)+1):
            num*=n
            den*=t
            n-=1
        return num//den
    else:
        return 0

def j(gen,n):
    sum=0
    for i in range(0,len(gen)):
        sum=sum+gen[i][0]*choose(gen[i][1]+1,n+1)
    return sum

def par(gen):
    P=[]
    N=[]
    for i in range(0,j(gen,1)):
        N.append(i)
    for t in range(0,len(gen)):
        n=gen[t]
        for k in range(0,n[0]):
            P.append(N[0:n[1]+1])
            N=N[n[1]+1:len(N)]
    return P


def M(gen):
    P=par(gen)
    l=[]
    for i in range(0,len(P)):
        l=l+c(P[i],0)+[P[i]]
    return MakeLex(l)


def h(s): #converts a simplice to string and returns its hash value
    str1 = ','.join(str(e) for e in s)
    return hash(str1)

def location(s,x = 'r'):
    k = h(s)%len(H)
    if type(H[k][0]) == list:
        for i in range(0,len(H[k])):
            if H[k][i][2] == s:
                if x == 'r':
                    return H[k][i][0]
                else:
                    return H[k][i][1]
    else:
        if x == 'r':
            return H[k][0]
        else:
            return H[k][1]



gen=[[2,10]]
P = par(gen)
d=gen[len(gen)-1][1]
global D,H
D=[zeros(1,j(gen,0))]
for i in range(0,d):
    D.append(zeros(j(gen,i),j(gen,i+1)))




master = M(gen)
H = []
for i in range(0,2*len(master)):
    H.append(0)
c = 0
for i in range(0,len(master)):
    n = h(master[i])%len(H)
    if len(master[i]) > len(master[i-1]):
        c = 0
    if H[n] == 0:
        H[n] = [c,i,master[i]]
    elif type(H[n][0]) == list:
        H[n].append([c,i,master[i]])
    else:
        H[n] = [H[n],[c,i,master[i]]]
    c = c+1

def boundarydecomp(s,p):
    global D
    if type(s) == int:
        s = [s]
    dim = len(s) - 1
    m = D[dim]
    for x in range(0,p):
        if len(s) <= 1:
            break
        ss = s.copy()
        ss.remove(s[dim-x])
        pos = location(ss,'r')
        if dim-x%2 ==0:
            m[pos,location(s,'r')] = 1
        else:
            m[pos,location(s,'r')] = -1
    for y in range(p,len(s)):
        if len(s) <= 1:
            break
        ss = s.copy()
        ss.remove(s[dim-y])
        pos = location(ss,'r')
        if dim-y%2 ==0:
            m[pos,location(s,'r')] = 1
        else:
            m[pos,location(s,'r')] = -1
        boundarydecomp(ss,y)

for i in range(0,len(P)):
    boundarydecomp(P[i],0)
#print(D)

#print(d)
#print(j(gen,0))
#print(j(gen,1))
#print(D[1].shape)








#s=[0,1,2,3]
#t=[4,5,6]
#l=c(s,0)+[s]+c(t,0)+[t]
#print(l)
#print(MakeLex(l))
