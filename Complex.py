from sympy import *
import random
class Complex:
    def __init__(self):
        self.master = []
        self.D = []
        self.H = []
        self.deadsimps = []
        self.tor = [[[0,2],[3,5]],[[0,1],[4,5]],[[1,2],[3,4]]]
    def complexify(self,s,p = 0): #generates the complex of simplex s when p is left as 0. Otherwise a subcomplex of faces with vertices < p never removed, is generated
        if len(s)==1:
            return [s]
        else:
            l=[]
            for i in range(p,len(s)):
                n=len(s)-i-1
                ss=s.copy()
                ss.remove(s[n])
                l.append(ss)
                l=l+self.complexify(ss,i)
            return l + [s]
    def isLess(self,s,t): #lexicographical ordering boolean function for simplices. if s < t, returns true
        n=min(len(s),len(t))
        for i in range(0,n):
            if s[i]<t[i]:
                return True
            elif s[i]>t[i]:
                return False
        if len(s)<len(t):
            return True
        return False
    def MakeLex(self,l): #quicksort implimented with lexicographic ordering that also preferences size of simplices
        if len(l)<=1:
            return l
        p=random.randrange(0,len(l))
        piv=l[p]
        a=[]
        b=[]
        for i in range(0,len(l)):
            y = [len(l[i])]+l[i]
            z =[len(piv)]+piv
            if self.isLess(y,z) and i!=p:
                a.append(l[i])
            if self.isLess(z,y) and i!=p:
                b.append(l[i])
        l=self.MakeLex(a)+[piv]+self.MakeLex(b)
        return l
    def choose(self,n,k): #standard combination function
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
    def j(self,gen,n): #returns the number of n dimensional simplices that will be generated from the generating set gen
        sum=0
        for i in range(0,len(gen)):
            sum=sum+gen[i][0]*self.choose(gen[i][1]+1,n+1)
        return sum
    def par(self,gen): #returns partition generated from generating set gen
        P=[]
        N=[]
        for i in range(0,self.j(gen,0)):
            N.append(i)
        for t in range(0,len(gen)):
            n=gen[t]
            for k in range(0,n[0]):
                P.append(N[0:n[1]+1])
                N=N[n[1]+1:len(N)]
        return P
    def M(self,gen): #generates the full master list of simplices from gen
        P=self.par(gen)
        l=[]
        for i in range(0,len(P)):
            l=l+self.complexify(P[i],0)
        return self.MakeLex(l)
    def h(self,s): #converts a simplice to string and returns its hash value
        str1 = ','.join(str(e) for e in s)
        return hash(str1)
    def location(self,s,x = 'r'): #uses hash table H to find quickly find the relative ('r') or absolute ordering of a simplex
        k = self.h(s)%len(self.H)
        if type(self.H[k][0]) == list:
            for i in range(0,len(self.H[k])):
                if self.H[k][i][2] == s:
                    if x == 'r':
                        return self.H[k][i][0]
                    else:
                        return self.H[k][i][1]
        else:
            if x == 'r':
                return self.H[k][0]
            else:
                return self.H[k][1]
    def setGenData(self,s,r = [],boundaryrel =None,steps = 0):
        self.master.clear()
        self.D.clear()
        self.deadsimps.clear()
        self.master = self.M(s)
        self.D=[zeros(1,self.j(s,0))]
        d = s[len(s)-1][1]
        #initializes boundary matrix D
        self.D = self.D + [zeros(self.j(s,i),self.j(s,i+1)) for i in range(0,d)]
        self.H.clear()
        self.H = [0 for i in range(0,2*len(self.master))]
        c = 0
        for i in range(0,len(self.master)): #fills H with simplex positions
            n = self.h(self.master[i])%len(self.H)
            flag = 0
            if len(self.master[i]) > len(self.master[i-1]):
                c = 0
            if self.H[n] == 0:
                flag  = 1
                self.H[n] = [c,i,self.master[i],self.master[i]]
            else:
                if type(self.H[n][0]) == list:
                    flag = 2
                    self.H[n].append([c,i,self.master[i],self.master[i]])
                else:
                    flag = 3
                    self.H[n] = [self.H[n],[c,i,self.master[i],self.master[i]]]
            #if master[i] == [18,19,20,21,22,23,24,25]:
            #    print("H(n)",H[n])
             #   input(flag)
            c = c+1
        if not (boundaryrel == None):
            a = self.getBoundary(boundaryrel[0])
            b = self.getBoundary(boundaryrel[1])
            if type(steps) == list:
                x = self.coupleSimps(a,b,steps[0])
            else:
                x = self.coupleSimps(a,b,steps)
            stepind = 1
            for simp in boundaryrel[2:len(boundaryrel)]:
                bound = self.getBoundary(simp)
                if type(steps) == list:
                    x = self.coupleSimps(x,bound,steps[stepind])
                    stepind+=1
                else:
                    x = self.coupleSimps(x,bound,steps)
            r = r + x
        self.initializeRelations(r)
        self.makeBoundary(s)
        self.pruneSimps()

    def boundarydecomp(self,s,p = 0): #takes simplex s, and recursively calls a rightbound. Will update all boundaries of the COMPLEX of s (s and all faces)
        if type(s) == int:
            s = [s]
        dim = len(s) - 1
        m = self.D[dim]
        for x in range(0,p):
            if len(s) <= 1:
                break
            ss = s.copy()
            ss.remove(s[dim-x])
            pos = self.location(self.lowestOrderRelation(ss),'r')
            if self.lowestOrderRelation(ss) != ss:
                self.deadsimps.append(ss)
            if (dim-x)%2 ==0:
                m[pos,self.location(s,'r')] = m[pos,self.location(s,'r')] + 1
            else:
                m[pos,self.location(s,'r')] = m[pos,self.location(s,'r')] - 1
        for y in range(p,len(s)):
            if len(s) <= 1:
                break
            ss = s.copy()
            ss.remove(s[dim-y])
            pos = self.location(self.lowestOrderRelation(ss),'r')
            if self.lowestOrderRelation(ss) != ss:
                self.deadsimps.append(ss)
            if (dim-y)%2 == 0:
                m[pos,self.location(s,'r')] = m[pos,self.location(s,'r')] + 1
            else:
                m[pos,self.location(s,'r')] = m[pos,self.location(s,'r')] - 1
            self.boundarydecomp(self.lowestOrderRelation(ss),y)

    def makeBoundary(self,s):
        P = self.par(s)
        for i in range(0,len(P)):
            self.boundarydecomp(P[i],0)
    def initializeRelations(self,Rel):#gets the full set of relations, and updates H with relation information for each simplex
        global H,deadsimps
        Rel = self.getAllRelations(Rel)
        for a in Rel:
            for j in range(1,len(a)):
                self.deadsimps+=self.complexify(a[j])
        for i in range(0,len(Rel)):
            s = self.MakeLex(Rel[i])
            lor = s[0]
            for u in range(0,len(s)):
                k = self.h(s[u])%len(self.H)
                #print("Simp",s[u],"hash",k,"table element",H[k],"\n")
                if type(self.H[k][0]) == list:
                    for x in range(0,len(self.H[k])):
                        if self.H[k][x][2] == s[u]:
                            lor2 = self.lowestOrderRelation(s[u])
                            if self.isLess(lor2,lor):
                                self.H[k][x][3] = lor2
                                l = self.h(s[0])%len(self.H)
                                b = self.H[l]
                                if type(b[0]) == list:
                                    for y in range(0,len(b)):
                                        if b[y][2] == s[0]:
                                            b[y][3] = lor2
                                else:
                                    b[3] = lor2
                                lor = lor2

                            else:
                                self.H[k][x][3] = lor
                            break
                else:
                    lor2 = self.lowestOrderRelation(s[u])
                    if self.isLess(lor2,lor):
                        self.H[k][3] = lor2
                        l = self.h(s[0])%len(self.H)
                        b = self.H[l]
                        if type(b[0]) == list:
                            for y in range(0,len(b)):
                                if b[y][2] == s[0]:
                                    b[y][3] = lor2
                        else:
                            b[3] = lor2
                        lor = lor2
                    else:
                        self.H[k][3] = lor
    def lowestOrderRelation(self,s): #finds the lowest order relation of simplex s, using H
        k = self.h(s)%len(self.H)
        if type(self.H[k][0]) == list:
            for i in range(0,len(self.H[k])):
                if self.H[k][i][2] == s:
                    if self.H[k][i][3] == s:
                        return s
                    else:
                        return self.lowestOrderRelation(self.H[k][i][3])
        else:
            if self.H[k][3] == s:
                return s
            else:
                return self.lowestOrderRelation(self.H[k][3])
    def getAllRelations(self,Rel):#generates all lower order relations from a relation set Rel
        fullRel = []
        for i in range(0,len(Rel)):
            lowerRels = []
            for u in range(0,len(Rel[i])):
                lowerRels.append(self.complexify(Rel[i][u],0))
            for k in range(0,len(lowerRels[0])):
                newRel = []
                for x in range(0,len(lowerRels)):
                    newRel.append(lowerRels[x][k])
                fullRel.append(newRel)
        return Rel+fullRel
    def pruneSimps(self):
        global D
        d = self.MakeLex(self.deadsimps)
        currentdim = 0
        i = 0
        for x in d:
            if len(x)-1>currentdim:
                currentdim = len(x)-1
                i = 0
            p = self.location(x)
            self.D[currentdim].col_del(p - i)
            if currentdim < len(self.D) - 1:
                self.D[currentdim + 1].row_del(p - i)
            i = i + 1
    def coupleSimps(self,s1,s2,parity = 0): #takes two lists of simplices and parity as either an integer or list, and will combine simps based on steps by parity
        start = s1[0]                       #This function is useful for easily specifying relations over a set of ordered simplices (see makeHole fuction as an example)
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
    def getBoundary(self,s): #returns a list representing the boundary of simplex s
        b = []
        for a in s:
            ss = s.copy()
            ss.remove(a)
            b.append(ss)
        b.reverse()
        return b
    def makeHole(self,k): #takes dimension k, and set generation data to construct a complex with a k-dimensional hole
        S = [[2,k]]
        P = self.par(S)
        R = self.coupleSimps(self.getBoundary(P[0]),self.getBoundary(P[1]))
        self.setGenData(S,R)
    def ker(self,k,mat): #takes matrix boundary list, and dimension k argument
        return Matrix([a.T for a in mat[k].nullspace()]).T
    def Betti(self,k,m = None): #returns codimension of m[k+1] in the the kernel of m[k]
        if m == None:
            m = self.D
        if k>=len(m):
            return 0
        mm = self.ker(k,m)
        if k +1>= len(m):
            imd = 0
        else:
            imd = len(m[k+1].rref()[1])
        return mm.shape[1] - imd
    def BettiAll(self):
        return [self.Betti(k,self.D) for k in range(0,len(self.D))]
