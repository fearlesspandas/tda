#cleaned up/optimized collection of functions from Betti and TDA. This collection can currently generate the boundary matrix of 
#simplices with relations (gluings). TDA 2016-17 Landon Renzullo/Marian Anton CCSU
def complexify(s,p = 0): #generates the complex of simplex s when p is left as 0. Otherwise a subcomplex of faces with vertices < p never removed, is generated
    if len(s)==1:
        return[]
    else:
        l=[]
        for i in range(p,len(s)):
            n=len(s)-i-1
            ss=s.copy()
            ss.remove(s[n])
            l.append(ss)
            l=l+complexify(ss,i)
        return l

def isLess(s,t): #lexicographical ordering boolean function for simplices
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

def MakeLex(l): #quicksort implimented with lexicographic ordering that also preferences size of simplices
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

def choose(n,k): #standard combination function
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

def j(gen,n): #returns the number of n dimensional simplices that will be generated from the generating set gen
    sum=0
    for i in range(0,len(gen)):
        sum=sum+gen[i][0]*choose(gen[i][1]+1,n+1)
    return sum

def par(gen): #returns partition generated from generating set gen
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


def M(gen): #generates the full master list of simplices from gen
    P=par(gen)
    l=[]
    for i in range(0,len(P)):
        l=l+complexify(P[i],0)+[P[i]]
    return MakeLex(l)


def h(s): #converts a simplice to string and returns its hash value
    str1 = ','.join(str(e) for e in s)
    return hash(str1)

def location(s,x = 'r'): #uses hash table H to find quickly find the relative ('r') or absolute ordering of a simplex
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

def setGenData(s,r = []):
    global master,D,H
    master.clear()
    D.clear()
    master = M(s)
    D=[zeros(1,j(s,0))]
    d = s[len(s)-1][1]
    for i in range(0,d): #initializes boundary matrix D
        D.append(zeros(j(s,i),j(s,i+1)))
    H.clear()
    for i in range(0,2*len(master)): #initializes ordinal hash table H
        H.append(0)
    c = 0
    for i in range(0,len(master)): #fills H with simplex positions
        n = h(master[i])%len(H)
        if len(master[i]) > len(master[i-1]):
            c = 0
        if H[n] == 0:
            H[n] = [c,i,master[i],master[i]]
        elif type(H[n][0]) == list:
            H[n].append([c,i,master[i],master[i]])
def complexify(s,p = 0): #generates the complex of simplex s when p is left as 0. Otherwise a subcomplex of faces with vertices < p never removed, is generated
    if len(s)==1:
        return[]
    else:
        l=[]
        for i in range(p,len(s)):
            n=len(s)-i-1
            ss=s.copy()
            ss.remove(s[n])
            l.append(ss)
            l=l+complexify(ss,i)
        return l + [s]

def isLess(s,t): #lexicographical ordering boolean function for simplices. if s < t, returns true
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

def MakeLex(l): #quicksort implimented with lexicographic ordering that also preferences size of simplices
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

def choose(n,k): #standard combination function
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

def j(gen,n): #returns the number of n dimensional simplices that will be generated from the generating set gen
    sum=0
    for i in range(0,len(gen)):
        sum=sum+gen[i][0]*choose(gen[i][1]+1,n+1)
    return sum

def par(gen): #returns partition generated from generating set gen
    P=[]
    N=[]
    for i in range(0,j(gen,0)):
        N.append(i)
    for t in range(0,len(gen)):
        n=gen[t]
        for k in range(0,n[0]):
            P.append(N[0:n[1]+1])
            N=N[n[1]+1:len(N)]
    return P


def M(gen): #generates the full master list of simplices from gen
    P=par(gen)
    l=[]
    for i in range(0,len(P)):
        l=l+complexify(P[i],0)
    return MakeLex(l)


def h(s): #converts a simplice to string and returns its hash value
    str1 = ','.join(str(e) for e in s)
    return hash(str1)

def location(s,x = 'r'): #uses hash table H to find quickly find the relative ('r') or absolute ordering of a simplex
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

def setGenData(s,r = []):
    global master,D,H
    master.clear()
    D.clear()
    master = M(s)
    D=[zeros(1,j(s,0))]
    d = s[len(s)-1][1]
    for i in range(0,d): #initializes boundary matrix D
        D.append(zeros(j(s,i),j(s,i+1)))
    H.clear()
    for i in range(0,2*len(master)): #initializes ordinal hash table H
        H.append(0)
    c = 0
    for i in range(0,len(master)): #fills H with simplex positions
        n = h(master[i])%len(H)
        if len(master[i]) > len(master[i-1]):
            c = 0
        if H[n] == 0:
            H[n] = [c,i,master[i],master[i]]
        elif type(H[n][0]) == list:
            H[n].append([c,i,master[i],master[i]])
        else:
            H[n] = [H[n],[c,i,master[i],master[i]]]
        c = c+1

    initializeRelations(r)
    makeBoundary(s)


global D,H,master,tor,kln
master = []
D = []
H = []
tor = [[[0,2],[3,5]],[[0,1],[4,5]],[[1,2],[3,4]]]
kln = [[[0,2],[3,5]],[[0,1],[4,5]],[[1,2],[3,4]]]

def boundarydecomp(s,p = 0): #takes simplex s, and recursively calls a rightbound. Will update all boundaries of the COMPLEX of s (s and all faces)
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
        pos = location(lowestOrderRelation(ss),'r')
        if (dim-x)%2 ==0:
            m[pos,location(s,'r')] = m[pos,location(s,'r')] + 1
        else:
            m[pos,location(s,'r')] = m[pos,location(s,'r')] - 1
    for y in range(p,len(s)):
        if len(s) <= 1:
            break
        ss = s.copy()
        ss.remove(s[dim-y])
        pos = location(lowestOrderRelation(ss),'r')
        if (dim-y)%2 == 0:
            m[pos,location(s,'r')] = m[pos,location(s,'r')] + 1
        else:
            m[pos,location(s,'r')] = m[pos,location(s,'r')] - 1
        boundarydecomp(lowestOrderRelation(ss),y)

def makeBoundary(s):
    P = par(s)
    for i in range(0,len(P)):
        boundarydecomp(P[i],0)

def initializeRelations(Rel):#gets the full set of relations, and updates H with relation information for each simplex
    global H
    Rel = getAllRelations(Rel)
    for i in range(0,len(Rel)):
        s = MakeLex(Rel[i])
        lor = s[0]
        for u in range(0,len(s)):
            k = h(s[u])%len(H)
            if type(H[k][0]) == list:
                for x in range(0,len(H[k])):
                    if H[k][x][2] == s[u]:
                        lor2 = lowestOrderRelation(s[u])
                        if isLess(lor2,lor):
                            H[k][x][3] = lor2
                            l = h(s[0])%len(H)
                            b = H[l]
                            if type(b[0]) == list:
                                for y in range(0,len(b)):
                                    if b[y][2] == s[0]:
                                        b[y][3] = lor2
                            else:
                                b[3] = lor2
                            lor = lor2

                        else:
                            H[k][x][3] = lor
                        break
            else:
                lor2 = lowestOrderRelation(s[u])
                if isLess(lor2,lor):
                    H[k][3] = lor2
                    l = h(s[0])%len(H)
                    b = H[l]
                    if type(b[0]) == list:
                        for y in range(0,len(b)):
                            if b[y][2] == s[0]:
                                b[y][3] = lor2
                    else:
                        b[3] = lor2
                    lor = lor2
                else:
                    H[k][3] = lor
def lowestOrderRelation(s): #finds the lowest order relation of simplex s, using H
    k = h(s)%len(H)
    if type(H[k][0]) == list:
        for i in range(0,len(H[k])):
            if H[k][i][2] == s:
                if H[k][i][3] == s:
                    return s
                else:
                    return lowestOrderRelation(H[k][i][3])
    else:
        if H[k][3] == s:
            return s
        else:
            return lowestOrderRelation(H[k][3])
def getAllRelations(Rel):#generates all lower order relations from a relation set Rel
    fullRel = []
    for i in range(0,len(Rel)):
        lowerRels = []
        for u in range(0,len(Rel[i])):
            lowerRels.append(complexify(Rel[i][u],0))
        for k in range(0,len(lowerRels[0])):
            newRel = []
            for x in range(0,len(lowerRels)):
                newRel.append(lowerRels[x][k])
            fullRel.append(newRel)
    return Rel+fullRel

