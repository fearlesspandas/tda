from sympy import *

def ker(k,mat): #takes matrix boundary list, and dimension k argument
    m1 = mat[k]
    b = m1.shape[0]
    mi = eye(m1.shape[1])
    m1 = m1.col_join(mi)
    m = m1.T.rref()
    ma = m[0].T
    kermat = []
    for i in range(0,min(b,len(m[1]))):
        ma.row_del(0)
        if m[1][0] < b:
            m[1].remove(m[1][0])
    for j in range(0,len(m[1])):
        kermat.append(ma.col(ma.shape[1] -1 - j).T)
    return Matrix(kermat).T
def getdegree(mat,col): #matrix argument, not list of matrices
	deg = 0
	for i in range(0,mat.shape[0]):
		if abs(mat[i,col]) == 1:
			deg = deg + 1
	return deg
def ori(v_1,v_2):
    if v_1.shape != v_2.shape:
        return 0
    for i in range(0,v_1.shape[0]):
        if v_1[i,0] != 0 and v_2[i,0] != 0 and abs(v_1[i,0]) == abs(v_2[i,0]):
            if v_1[i,0] == v_2[i,0]:
                return -1
            else:
                return 1
    return 0
def minsearch(mat,k,b=0,s=[],m = 0,m_deg = 0): # k is the index of the column of cycle D, b bounds the list of remaining columns. Returns the m cycle hiding in basis vector mat_k
    if m == 0:
        m = mat.col(k)
        m_deg = getdegree(mat,k)
        s.append(m)
    for i in range(b,mat.shape[1]):
        if i !=k:
            x = mat.col(i)
            if i == b:
                s.append(x)
            else:
                s[len(s) - 1] = x
            y = mat.col(k) + ori(x,mat.col(k))*x
            if y != mat.col(k):
                y_deg = getdegree(y,0)
                if y_deg < m_deg and y_deg != 0:
                    m = y
                    m_deg = y_deg
                mat2 = mat.col_insert(k,y)
                mat2.col_del(k+1)
                m = minsearch(mat2,k,i+1,s,m,m_deg)
    return m
def allmins(mat):
    minmat = []
    for i in range(0,mat.shape[1]):
        x = minsearch(mat,i).T
        found = False
        for j in range(len(minmat)):
            if x == minmat[j]:
                found = True
                break
        if not found:
            minmat.append(minsearch(mat,i).T)
    minmat = Matrix(minmat).T
    return minmat

def Bettislow(dim,D):
    mat = D.copy()
    if dim>= len(mat):
        return 0
    m = allmins(ker(dim,mat))
    if dim+1<len(mat):
        im = allmins(mat[dim+1].copy())
        print(im)
    else:
        return m.shape[1]
    count = 0
    for i in range(0,im.shape[1]):
        if getdegree(im,i) > 0:
            count+=1

    return m.shape[1] - im.shape[1]
