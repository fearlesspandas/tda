def complexify(s,p = 0): #generates the complex of simplex s when p is left as 0. Otherwise a subcomplex of faces with vertices < p never removed, is generated
    if len(s)==1:
        return [s]
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

import random, BettiCalc

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

def setGenData(s,r = [],boundaryrel =None,steps = 0):
    global master,D,H,deadsimps
    master.clear()
    D.clear()
    deadsimps.clear()
    master = M(s)
    D=[zeros(1,j(s,0))]
    d = s[len(s)-1][1]
    #initializes boundary matrix D
    D = D + [zeros(j(s,i),j(s,i+1)) for i in range(0,d)]
    H.clear()
    H = [0 for i in range(0,2*len(master))]
    c = 0
    for i in range(0,len(master)): #fills H with simplex positions
        n = h(master[i])%len(H)
        flag = 0
        if len(master[i]) > len(master[i-1]):
            c = 0
        if H[n] == 0:
            flag  = 1
            H[n] = [c,i,master[i],master[i]]
        else:
            if type(H[n][0]) == list:
                flag = 2
                H[n].append([c,i,master[i],master[i]])
            else:
                flag = 3
                H[n] = [H[n],[c,i,master[i],master[i]]]
        #if master[i] == [18,19,20,21,22,23,24,25]:
        #    print("H(n)",H[n])
         #   input(flag)
        c = c+1
    if not (boundaryrel == None):
        a = getBoundary(boundaryrel[0])
        b = getBoundary(boundaryrel[1])
        if type(steps) == list:
            x = coupleSimps(a,b,steps[0])
        else:
            x = coupleSimps(a,b,steps)
        stepind = 1
        for simp in boundaryrel[2:len(boundaryrel)]:
            bound = getBoundary(simp)
            if type(steps) == list:
                x = coupleSimps(x,bound,steps[stepind])
                stepind+=1
            else:
                x = coupleSimps(x,bound,steps)
        r = r + x
    initializeRelations(r)
    makeBoundary(s)
    pruneSimps()


global D,H,master,tor,kln,deadsimps
master = []
D = []
H = []
deadsimps = []
tor = [[[0,2],[3,5]],[[0,1],[4,5]],[[1,2],[3,4]]]
kln = [[[0,2],[3,5]],[[0,1],[4,5]],[[1,2],[3,4]]]

def boundarydecomp(s,p = 0): #takes simplex s, and recursively calls a rightbound. Will update all boundaries of the COMPLEX of s (s and all faces)
    global D,deadsimps
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
        if lowestOrderRelation(ss) != ss:
            deadsimps.append(ss)
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
        if lowestOrderRelation(ss) != ss:
            deadsimps.append(ss)
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
    global H,deadsimps
    Rel = getAllRelations(Rel)
    for a in Rel:
        for j in range(1,len(a)):
            deadsimps+=complexify(a[j])
    for i in range(0,len(Rel)):
        s = MakeLex(Rel[i])
        lor = s[0]
        for u in range(0,len(s)):
            k = h(s[u])%len(H)
            #print("Simp",s[u],"hash",k,"table element",H[k],"\n")
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
def pruneSimps():
    global D
    d = MakeLex(deadsimps)
    currentdim = 0
    i = 0
    for x in d:
        if len(x)-1>currentdim:
            currentdim = len(x)-1
            i = 0
        p = location(x)
        D[currentdim].col_del(p - i)
        if currentdim < len(D) - 1:
            D[currentdim + 1].row_del(p - i)
        i = i + 1
def Betti(k,m = None): #returns codimension of m[k+1] in the the kernel of m[k]
    if m == None:
        m = D
    if k>=len(m):
        return 0
    mm = ker(k,m)
    if k +1>= len(m):
        imd = 0
    else:
        imd = len(m[k+1].rref()[1])
    return mm.shape[1] - imd

def ker(k,mat):
    return BettiCalc.ker(k,mat)

def getBoundary(s):
    b = []
    for a in s:
        ss = s.copy()
        ss.remove(a)
        b.append(ss)
    b.reverse()
    return b

def coupleSimps(s1,s2,parity = 0): #takes two lists of simplices and parity as either an integer or list, and will combine simps based on steps by parity
    start = s1[0]
    a = 0
    if type(parity) == int:
        b = parity
    else:
        b = parity[0]
    rels = []
    for i in range(0,len(s1)):
        if type(s1[a][0]) == int:
            r = [s1[a],s2[b%len(s2)]]
            rels.append(r)
        elif type(s1[a][0]) == list:
            rels.append(s1[a])
            rels[len(rels)-1].append(s2[b%len(s2)])
        a += 1
        if parity == 0:
            b+=  1
        elif type(parity) == list:
            b+= parity[(i+1)%len(parity)]
        else:
            b += parity
    return rels


def openData(n = 23,m= 13,filename = "ExampleData.txt"):
    file = open(filename,"r")
    txt = file.read()
    file.close()
    dat = []
    dat = [[txt[2*i + j] for i in range(0,m) if (txt[2*i+j] != '\t' and txt[2*i+j] != '\n') ] for j in range(0,2*m*n,2*m)]
    return Matrix(dat)

def intersect(M,a,b): #returns list as [[i,pos1,pos2]]
    intersection = []
    ca = 0
    cb = 0
    for i in range(0,M.shape[0]):
        if M[i,a] == 1:
            ca = ca + 1
        if M[i,b] == 1:
            cb = cb + 1
        if M[i,a] == 1 and M[i,b] == 1:
            intersection.append((i,ca,cb))
    return intersection


def intersectAll(M): #returns the indexes of all intersections in the complexes from incidence matrix M
    return [[intersect(M,i,j),i,j] for i in range(0,M.shape[1]) for j in range(0,i)]

def findLargeSimp(M):
    m = 0
    for i in range(0,M.shape[1]):
        size = 0
        for j in range(0,M.shape[0]):
            if M[j,i] == 1:
                size = size + 1
        if size>m:
            m = size
    return m

def rowcount(M,i):
    c = 0
    for j in range(0,M.shape[0]):
        if M[j,i] == 1:
            c = c + 1
    return c
def genSimpDataFromMat(M):
    I = intersectAll(M)
    maxdim = findLargeSimp(M)
    CompData = [[0,i] for i in range(0,maxdim)]
    for i in range(0,M.shape[1]):
        dim = rowcount(M,i) - 1
        if dim>=0:
            CompData[dim][0] = CompData[dim][0] + 1
    return CompData


def calcPosfromMat(M,i): #finds the position among simplex M[i] among homogenous simplices
    count = 0
    for j in range(0,i):
        if rowcount(M,j) == rowcount(M,i):
            count = count + 1
    return count

def simpStart(M,G,i,dim = None):#takes incidence matrix, genData and index i, and finds the unversal label of the first vertex of the simplex at column i in M
    if dim == None:
        dim = rowcount(M,i)
    #print(i,j(G[0:dim-1],0) + calcPosfromMat(M,i))
    return j(G[0:dim-1],0) + calcPosfromMat(M,i)*dim
def getRelData(M,G,I): #gets relational data from incidence matrix
    rels = []
    for i in range(0,len(I)):
        if len(I[i][0])>0:
            matches = I[i][0]
            #dim_subsimp = len(I[i][0]) - 1
            dim_supersimp_i = rowcount(M,I[i][1])
            dim_supersimp_j = rowcount(M,I[i][2])
            v_i = simpStart(M,G,I[i][1],dim_supersimp_i)
            v_j = simpStart(M,G,I[i][2],dim_supersimp_j)
            simp1 = [v_i + matches[j][1]-1 for j in range(0,len(matches))]
            simp2 = [v_j + matches[j][2]-1 for j in range(0,len(matches))]
            #print("Simmps",simp1,simp2,"\n")
            rels.append([simp1,simp2])
    return rels

def BettiAll():
    return [Betti(k) for k in range(0,len(D))]
def incToBoundary(k,filename = "ExampleData.txt"): #reads incidence matrix off of text file, then converts to generating Data, then to Boundary matrix
    m = openData(filename)                         #k bounds the maximum dimension a simplex can have
    mm = m[0:k,0:m.shape[1]]                       #returns Betti numbers of the generated complex
    i = intersectAll(mm)
    g = genSimpDataFromMat(mm)
    r = getRelData(mm,g,i)
    setGenData(g,r)
    return BettiAll()
