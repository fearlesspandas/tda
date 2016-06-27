import random



def setsimpdata(s):
    global simpdata
    simpdata = s

def countvertices(s):
    count = 0
    for i in range(0,len(s)):
        count = count + (s[i][0]*s[i][1])
    return count


def generatenaturals(n):
    x = []
    for i in range(0,n):
        x.append(i)
    return x


def getsimplicialdata_basic(s): #takes simplicial data, and the list of vertices
    vertices_n = countvertices(s)
    v = generatenaturals(vertices_n)
    simps = []
    startpos = 0
    for i in range(0,len(s)):
        if i == 0 and s[i][1] != 0:
            a = s[i][1] - 1
            for j in range(0,a):
                simps.append([])
        homogsimps = []
        u = 0
        while u < s[i][0]:
            homogsimps.append(v[startpos + u*s[i][1]:startpos + (u+1)*s[i][1]])
            u = u + 1
        startpos = startpos + u*s[i][1]
        simps.append(homogsimps)

        if i != len(s) - 1 and s[i+1][1] - 1 > s[i][1]:
            a = s[i+1][1] - s[i][1] - 1
            for j in range(0,a):
                simps.append([])

    return simps

def getmaxdim(s): #returns the maximum dimension of the generators
    maxdim= 0
    for i in range(0,len(s)):
        if s[i][1] > maxdim:
            maxdim = s[i][1]
    return maxdim - 1

def getdecomposition(s): #takes a list of homogenous (same dimension) faces
    lowersimps = []
    for i in range(0,len(s)):
        face = s[i]
        if len(s[i]) == 1 or len(s[i]) == 0:
            return []
        for u in range(i,len(face)):
            lowface = face.copy()
            lowface.remove(face[len(face) - 1 - u - i])
            lowersimps.append(lowface)
            lowersimps = lowersimps + getdecomposition([lowface])
    return lowersimps
def decomp(s,p):
    if len(s) == 1:
        return []
    else:
        lowsimps = []
        for i in range(p,len(s)):
            s_low = s.copy()
            s_low.remove(s[len(s)-1-i])
            lowsimps.append(s_low)
            lowsimps = lowsimps + decomp(s_low,i)
        return lowsimps

def getFulldecomposition(s):
    simps = []
    for i in range(0,len(s)):
        simps = simps + getdecomposition(s[i])
    return simps

def getFulldecomp(s):
    simps = []
    for i in range(0,len(s)):
        for u in range(0,len(s[i])):
            simps = simps + decomp(s[i][u],0)
    return simps
def group_by_dim(s,generators):
    a = []
    for x in range(0,getmaxdim(generators)):
        a.append([])
    for i in range(0,len(s)):
        a[len(s[i])-1].append(s[i])
    return a

def mergelists(a,b):
    a_ = len(a)
    b_ = len(b)
    x = min(a_,b_)
    for i in range(0,x):
        if a_>b_:
            a[i] = a[i] + b[i]
        else:
            b[i] = b[i] + a[i]
    if len(a)>len(b):
        return a
    else:
        return b
def islessthan(a,b):
    if type(a) == int:
        if a<b:
            return True
        else:
            return False
    if len(a) != len(b) or len(a) == 0 or len(b) == 0: #this is new, may need to be changed
        return False
    for i in range(0,len(a)):
        if len(a) == 1 and len(b) == 1:
            a = a[0]
            b = b[0]
            return islessthan(a,b)
        elif a[i] < b[i]:
            return True
        elif a[i] > b[i]:
            return False
    return False
    
def quicksort(a): #returns lexicographically sorted list
    if len(a) <= 1:
        return a
    p = random.randrange(0,len(a)) #pointer to the piv ####turn random.randrange to pyb.rng()%len(a)
    piv = a[p]
    c = []
    d = []
    for i in range(0,len(a)):
        if islessthan(a[i], piv) and i!=p:
            c.append(a[i])
        if (not islessthan(a[i],piv)) and i!=p:
            d.append(a[i])

    a = quicksort(c) + [piv] + quicksort(d)
    return a #returns sorted list

def quicksort_dim(a):
    if len(a) <= 1:
        return a
    p = random.randrange(0,len(a))

    piv = len(a[p][0])
    c = []
    d = []
    for i in range(0,len(a)):

        dim = len(a[i][0])
        if dim < piv and i != p:
            c.append(a[i])
        elif dim > piv and i !=p:
            d.append(a[i])
    a = quicksort_dim(c) + [a[p]] + quicksort(d)
    return a

def sortsimplices(s):
    for i in range(0,len(s)):
        s[i] = quicksort(s[i])
    s = quicksort_dim(s)
    return s

def getsimplicialdata(s):
    sd = getsimplicialdata_basic(s)
    sdc = getFulldecomp(sd)
    dimsdc = group_by_dim(sdc,s)
    merged = mergelists(dimsdc,sd)
    sortedsimps = sortsimplices(merged)
    return sortedsimps


