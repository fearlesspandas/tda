from sympy import *
import random,Delta

global matrix,Big_matrix,rel_matrix,mth,style,gen_data,simp_li,rel_set,simplicial_data

rel_set = []

simplicial_data = []
simp_li = []
gen_data = []

matrix = []
Big_matrix = []
rel_matrix = []
mth = 'm'
style = 'eq'

def set_gen_data(s):
    global gen_data,simp_li,simplicial_data
    gen_data = s
    simp_li = Delta.getsimplicialdata_basic(s)
    simplicial_data = Delta.getsimplicialdata(s)
    rel_matrix_init()
    boundary_init()

#m = zeros(len(simplicial_data[0]),len(simplicial_data[1]))
def init_data(title):
    rawsimps = []
    simps = []
    try:
        file = open(title,'r')
        raw = file.read()
        i = 0
        while i < len(raw):
            if raw[i] == '{':
                u = i
                prev = 0
                textsimps = [] #homogenous simplices
                while u < len(raw) and raw[u] != '}':
                    if raw[u] == '[':
                        prev = u
                    elif raw[u] == ']':
                        textsimps.append(raw[prev+1:u])
                    u = u + 1
                rawsimps.append(textsimps)
                i = u
            else:
                i = i + 1
        for x in range(0,len(rawsimps)):
            homogsimps = []
            simpstext = rawsimps[x]
            for y in range(0,len(simpstext)):
                onesimp = simpstext[y]
                goodsimp = []
                commapoints = []
                if len(onesimp) == 1:
                    goodsimp.append(int(onesimp))
                else:
                    for z in range(0,len(onesimp)):
                        if onesimp[z] == ',':
                            commapoints.append(z)
                        if z == len(onesimp) - 1:
                            commapoints.append(len(onesimp))
                    for z in range(0,len(commapoints)):
                        if z == 0:
                            goodsimp.append(int(onesimp[0:commapoints[0]]))
                        else:
                            goodsimp.append(int(onesimp[commapoints[z-1]+1:commapoints[z]]))
                homogsimps.append(goodsimp)
            simps.append(homogsimps)
    except:
        file = open(title,'w')
    file.close()
    return simps


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


def rel_matrix_init(): #creates identity matrix among simplices
    global rel_matrix
    rel_matrix = []
    for i in range(0,len(simplicial_data)):
        n = len(simplicial_data[i])
        a = []
        for x in range(0,n):
            a.append([])
            for y in range(0,n):
                if y == x:
                    a[x].append(1)
                else:
                    a[x].append(0)
        m = Matrix(a)
        rel_matrix.append(m)

def findsimp(simp,dim): #this could definitely be improved for speed
    if type(simp) == int:
        simp = [simp]
    for i in range(0,len(simplicial_data[dim])):
        if simp == simplicial_data[dim][i]:
            return i
    return False

def update_rel_matrix(r):#total equivalence among relations, takes a pair of relations only
    if type(r[0]) == int:
        r[0] = [r[0]]
        r[1] = [r[1]]
    dim = len(r[0]) - 1
    if dim < 0:
        return 0
    a = findsimp(findrel(r[0]),dim)
    b = findsimp(findrel(r[1]),dim)
    m = rel_matrix[dim]
    m[a,b] = 1
    m[b,a] = 1
    for i in range(0,m.shape[0]):
        if m[a,i] == 1:
            m[i,a] = 1
            m[i,b] = 1
            m[b,i] = 1
    for u in range(0,m.shape[0]):
        if m[u,b] == 1:
            m[a,u] = 1
            m[b,u] = 1
            m[u,a] = 1
    ##cleanup
    for x in range(0,m.shape[0]):
        for y in range(0,m.shape[0]):
            if m[x,y] == 0:
                s1 = simplicial_data[dim][x]
                s2 = simplicial_data[dim][y]
                if findrel(s1) == findrel(s2):
                    m[x,y] = 1
    rel_matrix[dim] = m

def find_local_rel(m,col):
    for i in range(0,m.shape[0]):
        if m[i,col] == 1:
            return i

def update_rel_matrix_ns(r):#update using an matrix that is not symmetric, partially transitive, takes a pair of relations only
    r1 = convert_to_equiv([r])
    r = r1[0]
    r = quicksort(r)
    if type(r[0]) == int:
        r[0] = [r[0]]
        r[1] = [r[1]]
    dim = len(r[0]) - 1
    a = findsimp(findrel(r[0]),dim)
    b = findsimp(findrel(r[1]),dim)
    m = rel_matrix[dim]
    if a == b:
        return m
    m[a,b] = 1
    for i in range(0,m.shape[0]):
        if m[i,b] == 1 and a < i:
            m[a,i] = 1
    for h in range(0,m.shape[0]):
        if m[b,h] == 1:
            m[a,h] == 1
    for u in range(0,m.shape[0]):
        if m[a,u] == 1:
            if u<b:
                m[u,b] = 1
            else:
                m[b,u] = 1
    #cleanup
    for l in range(0,m.shape[0]):
        for n in range(0,l):
            colrel = find_local_rel(m,l)
            rowrel = find_local_rel(m,n)
            if m[n,l] == 1 and colrel != rowrel :
                if colrel < rowrel:
                    m[colrel,n] = 1
                elif colrel > rowrel:
                    m[rowrel,l] = 1

    for x in range(0,m.shape[0]):
        for y in range(0,m.shape[0]):
            if m[x,y] == 0 and x < y :
                s1 = simplicial_data[dim][x]
                s2 = simplicial_data[dim][y]
                if findrel(s1) == findrel(s2):
                    m[x,y] = 1

    rel_matrix[dim] = m

def update_rel_matrix_nst(r): #not symmetric or transitive, takes a pair of relations only
    r1 = convert_to_equiv([r])
    r = r1[0]
    r = quicksort(r)
    if type(r[0]) == int:
        r[0] = [r[0]]
        r[1] = [r[1]]
    dim = len(r[0]) - 1
    a = findsimp(findrel(r[0]),dim)
    b = findsimp(findrel(r[1]),dim)
    m = rel_matrix[dim]
    if a == b:
        return m
    m[a,b] = 1
    return m
def update_rel_matrix_full(r,state):
    top = findrel(r[0])
    for i in range(1,len(r)):
        if state == 'nst':
            update_rel_matrix_nst([top,findrel(r[i])])
        elif state == 'ns':
            update_rel_matrix_ns([top,findrel(r[i])])
        else:
            update_rel_matrix([top,findrel(r[i])])

def killdeadchains():
    hitlist = [] #ordered in subgroups by dimension
    for i in range(0,len(rel_matrix)):
        m = rel_matrix[i]
        x = []
        for u in range(0,m.shape[0]):
            if find_local_rel(m,u) < u:
                x.append(u) #appends the column/row order necessary to be deleted
        hitlist.append(x)
    simpdata = simplicial_data.copy()
    for i in range(0,len(simpdata)):
        if len(hitlist) == 0:
            break
        simps = simpdata[i]
        del_set= hitlist[i]
        for h in range(0,len(del_set)):
            del_set[h] = simps[del_set[h]]
        for u in range(0,len(del_set)):
            simps.remove(del_set[u])
        simpdata[i] = simps
    return simpdata

def killdeadboundaries(matrix):
    rows = []
    col = []
    choppedmatrix = []
    for i in range(0,len(simplicial_data)):
        s_n = simplicial_data[i]
        k = []
        for u in range(0,len(s_n)):
            if findrel(s_n[u]) != s_n[u]:
                k.append(s_n[u])
        if i != 0 and i != len(simplicial_data) - 1:
            col.append(k)
            rows.append(k)
        elif i == 0:
            rows.append(k)
        elif i == len(simplicial_data) - 1:
            col.append(k)
    for j in range(1,len(matrix)):
        n = 0
        m = matrix[j]
        m_copy = m.copy()
        m_row = m_copy.shape[0]
        while n < m_row and n < len(rows[j-1]):
            placement = findsimp(rows[j-1][n],len(rows[j-1][n])-1)
            m_copy.row_del(placement - n)
            m_row = m_copy.shape[0]
            n = n + 1
        n = 0
        m_col = m_copy.shape[1]
        while n < m_col and n < len(col[j-1]):
            placement = findsimp(col[j-1][n],len(col[j-1][n])-1)
            m_copy.col_del(placement - n)
            m_col= m_copy.shape[1]
            n = n + 1
        choppedmatrix.append(m_copy)
    return choppedmatrix

def findrel(simp): #works
    if type(simp) == int:
        simp = [simp]
    dim = len(simp) - 1
    n = findsimp(simp,dim)
    m = rel_matrix[dim]
    x = n
    for i in range(0,m.shape[0]):
        if m[i,n] == 1:
            x = i
            break
    return simplicial_data[dim][x]


def boundary_init(): #updates matrix to boundary data with no relations observed
    global matrix,simplicial_data
    matrix = []
    for i in range(0,len(simplicial_data)):
        m = zeros(len(simplicial_data[i-1]),1)
        for u in range(0,len(simplicial_data[i])): #C_n loop
            high_simp = simplicial_data[i][u]
            a=[]
            for x in range(0,len(simplicial_data[i-1])):#constructing loop
                a.append(0)
            for j in range(0,len(high_simp)):
                low_simp = high_simp.copy()
                low_simp.remove(high_simp[j])
                k = findsimp(low_simp,len(low_simp)-1)
                if j%2 == 0:
                    a[k] = 1
                else:
                    a[k] = -1
            m = m.col_insert(m.shape[1]-1,Matrix(a))
        m.col_del(m.shape[1]-1)
        matrix.append(m)
    return matrix

def getgenerator(dim): #returns the generator of a given n-dimension simplicy
    for i in range(0,len(gen_data)):
        if gen_data[i][1] == dim + 1:
            return gen_data[i]
    return [0,dim + 1]

def boundarydecomp(s,p,t):#will take a linearly independent simplicie s and update its boundary
    global matrix          # as well as the boundaries of its decomposition. t is the order of
    if type(s) == int:      #s among the linearly independent simplicial set
        s = [s]
    dim = len(s) - 1
    m = matrix[dim]
    for i in range(p,len(s)):
        if len(s) <= 1:
            break
        g = getgenerator(dim - 1)
        r = i + t*(dim+1) + g[0]
        s_low = s.copy()
        s_low.remove(s[len(s)-1-i])
        print('n',g[0],"t",t,'pos',r,'i',i)
        if (dim-i)%2 == 0:
            print("+1",s_low)
            m[r,t] = 1
        else:
            print("-1",s_low)
            m[r,t] = -1
        boundarydecomp(s_low,0,r) #0 to i

def boundaryinit():
    global matrix
    matrix = []
    for x in range(0,len(simplicial_data)): #clears and prepares appropriate sized matrices
        cols = len(simplicial_data[x])
        rows = len(simplicial_data[x-1])
        matrix.append(zeros(rows,cols))
    for i in range(0,len(simp_li)):
        n_simps = simp_li[i]
        for u in range(0,len(n_simps)):
            boundarydecomp(n_simps[u],0,u)


def convert_to_equiv(r): #converts every element in a list of relations to
    for i in range(0,len(r)): #their equivalent relations
        if type(r[i]) == int:
            r[i] = [r[i]]
        x = len(r[i])
        for h in range(0,x):
            r[i][h] = findrel(r[i][h])
    return r
def prune_relation(r): #makes a list of relations with no repeats
    i = 0
    x = len(r)
    while i < x:
        c = r.count(r[i])
        while(c>1):
            r.remove(r[i])
            c = r.count(r[i])
        x = len(r)
        i = i + 1
    return r

def mtxmap(m,r):
    deadlist = []
    if type(r[0]) == int:
        topdim = 0
    else:
        topdim = len(r[0]) - 1
    rlow = compute_lower_rel(r)
    rlow = convert_to_equiv(rlow)
    rlow = prune_relation(rlow)
    for stuff in range(0,len(rlow)):
        rlow[stuff]= quicksort(rlow[stuff])
    for i in range(0,len(rlow)):
        rel = rlow[i]
        top = findrel(rel[0])
        if type(top) == int:
            dim = 0
        else:
            dim = len(top)-1
        x = findsimp(top,dim)
        a = m.row(x)
        for u in range(1,len(rel)):
            add = findrel(rel[u])
            if type(add) == int:
                lowdim = 0
            else:
                lowdim = len(add) - 1
            y = findsimp(add,lowdim)
            if x!=y:
                b = m.row(y)
                a = a + b
                deadlist.append(y)
        m.row_del(x)
        m = m.row_insert(x,a)
    for h in range(0,len(rlow)):
        update_rel_matrix(rlow[h])

    return m


def compute_lower_rel(r):
    if len(r) == 1 or type(r[0]) == int:
        return 0
    a = []
    for x in range(0,len(r[0])):
        a.append([])
    for i in range(0,len(r[0])):
        for u in range(0,len(r)):
            s = r[u].copy()
            s.remove(r[u][i])
            a[i].append(s)
    return a

def delta_rowaddition():
    global Big_matrix,matrix,rel_set
    for x in range(0,len(rel_set)):
        update_rel_matrix(rel_set[x])
    i = 0
    d = len(rel_set)
    while i < d:
        m = mtxmap(matrix[len(rel_set[i][0])-1],rel_set[i])
        matrix[len(rel_set[i][0])-1]= m
        nextrelz = compute_lower_rel(rel_set[i])
        if nextrelz !=0:
            rel_set = rel_set[0:i+1] + nextrelz + rel_set[i+1:len(rel_set)]
        Big_matrix.append(matrix)
        d = len(rel_set)
        print("Crrnt Relz", rel_set[i])
        print("___DelOperator___",matrix)
        print("Lwr Relz",nextrelz)
        print("---------------------------")
        i = i + 1
    return matrix

def delta_m(s):
    global Big_matrix,matrix,rel_set
    x = len(rel_set)
    if x == 0:
        return matrix
    i = 0
    while i < x:
        update_rel_matrix_full(rel_set[i],s) #change to nst for nst update
        nextrelz = compute_lower_rel(rel_set[i])
        #print("Crrnt Relz", rel_set[i])
        #print("___relM___",rel_matrix)
        #print("Lwr Relz",nextrelz)
        #print("---------------------------")
        rel_set = rel_set[0:i+1] + nextrelz + rel_set[i+1:len(rel_set)]
        x = len(rel_set)
        i = i + 1

    boundaries = []
    for u in range(1,len(matrix)):
        boundaries.append(rel_matrix[u-1]*matrix[u])
    return boundaries



def set_simplicial_data(s): #constructor for simplicial_data, also updates rel_matrix and boundary_init
    global simplicial_data
    simplicial_data = s
    print("completed variable change")
    rel_matrix_init()
    print("completed rel")
    boundary_init()
def set_rel_set(s):
    global rel_set
    rel_set = s

def set_preferences_to(mthd,styl): #mthd is the method for matrix multiplication and styl refers to the type of equivalence matrix used
        global mth,style
        success = False
        valid_m = ['m','ra','r']
        valid_s = ['nst','ns','eq']
        for i in range(0,len(valid_m)):
            if mthd == valid_m[i]:
                mth = mthd
                success = True
                if mthd == 'm':
                    success = False
                    for u in range(0,len(valid_s)):
                        if styl == valid_s[u]:
                            style = styl
                            success = True
        return success

def startDelta():
    m = input('Manipulation Method-->')
    s = input('Relational Style-->')
    if set_preferences_to(m,s):
        if m == 'm':
            print("___finalDel___", delta_m(s))
        else:
            print("___finalDel__", delta_rowaddition())
        print("___relM___",rel_matrix)
        print("Livingchains-->", killdeadchains())
        #print("____Matrix____",matrix)
        #print("____Big____",Big_matrix)


