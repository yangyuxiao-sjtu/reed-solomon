
from GF import  GF256 as GF
from Poly import Polynomial
import numpy as np
class rscoder:
    def __init__(self,n,k):
        """
        k is the length of message ,n is the lenth of codeword while m=n-k is the lenth of extra bit
        """
        self.n= n
        self.k =k
        self.m=n-k
        alpha = GF(3)
        """g(x)=(x-a^0)(x-a^1)..."""
        g=Polynomial([GF(1)])
        for i in range(0,n-k):
            p=Polynomial([GF(1),alpha**i])
            g=g*p
        self.g = g
        # print("g is",g.coef)
    def encode(self,message):
        Message = list(GF(x) for x in message)
        #print(list(self.g.coef))
        M=Polynomial(list(reversed(Message)))
        lss=[GF(1)]+[GF(0)]*self.m
        xm = Polynomial(lss)
        # print("XM")
        # print(xm)
      
        s= M*xm
        
       # print(list(s.coef))
        b= s%self.g
        s=s-b
        # print('s is ',s.coef)
        coef = s.coef
        return coef
        # return reversed(coef)
    def syndromes(self,r):
        n= self.n
        k= self.k
        s=[]
        for i in range(0,n-k+1):
            s.append(r.evaluate(GF(3)**i))
            # if(s[i]!=GF(0)):print("ERROR!",i)
        # print("r is",r.coef)
        tmp=GF(0)
        
        for i in r.coef:
            tmp=tmp+i
        # print("tmp is",tmp)
        # print("s is",s)
        return s
    def multiple(self,A,B):
        n=len(A[0])
        C=np.zeros((n,n),dtype=GF)
        for i in range(0,n):
            for j in range(0,n):
                for k in range(0,n):
                    C[i][j]+=A[i][k]*B[k][j]
        return C

    def gf_invert(self,A):
        n = len(A[0])
        I=np.zeros((n,n),dtype=GF)
        # print(A[0][0]*A[1][1]-A[0][1]*A[1][0])
        for i in range(n):
            I[i][i]=GF(1)
        for i in range(n):
            for j in range(n):
                if i != j:
                    ratio = A[j][i]*A[i][i].inverse()

                    for k in range(n):
                        A[j][k] = A[j][k] - ratio * A[i][k]
                        I[j][k] = I[j][k] - ratio * I[i][k]

        for i in range(n):
            divisor = A[i][i]
            for j in range(n):
                A[i][j] = A[i][j]*divisor.inverse()
                I[i][j] = I[i][j]*divisor.inverse()
        # print(A)
        return I
    def gauss_jordan(self,A, b):
        # 拼接增广矩阵
        M = np.hstack((A, b.reshape(-1, 1)))

        # 针对每一列进行消元
        for i in range(min(M.shape)):
            # 选取主元素
            max_row = i + np.argmax(np.abs(M[i:, i]))
            M[[i, max_row]] = M[[max_row, i]]

            # 将主元素所在的行除以主元素
            M[i] *= M[i, i].inverse()

            # 将主元素所在的列上下的元素消为0
            for j in range(M.shape[0]):
                if j != i:
                    for k in range(len(M[j])):
                        M[j][k] -= M[j, i] * M[i][k]

        # 返回解向量
        return M[:, -1]
    def decode(self,received):
        # print(received)
        n=self.n
        k=self.k
        m=n-k
        Received = list(GF(x) for x in received)
        r = Polynomial(Received)
        v = (n-k)//2
        S = self.syndromes(r)
        M=np.zeros((v,v),dtype=GF)
        for i in range(0,v):
            for j in range(0,v):
                    M[i][j]=S[i+j]
        L = np.zeros((v,1),dtype=GF)
        for i in range(0,v):
            L[i]=-(S[v+i])
        # print(M)
        # M_inv=self.gf_invert(M)
        # gamma = np.dot(M_inv,L)
        gamma = self.gauss_jordan(M,L)
        # print("gamma:",gamma)
        gamma_new = np.concatenate((gamma,np.array([GF(1)],dtype=GF) ),axis=0)
        # print(gamma_new)
        gamma = list()
        # for i in gamma_new:
        #     gamma.append(GF(int(gamma_new)))
        G = Polynomial(list(gamma_new))
        ret=list()
        fake_error=list()
        rev=reversed(received)
        for i in range(0,n):
            num = G.evaluate(GF(3)**(-i))
            if(num==GF(0)):
                print("error locate:",n-i-1)
                ret.append(n-i-1)
        error_length = len(ret)
        X_matrix=np.zeros((error_length,error_length),dtype=GF)
        alpha = GF(3)
        for i in range(0,error_length):
            for j in range(0,error_length):
                X_matrix[i][j]=(alpha**ret[j])**i
        S_vec = np.zeros((error_length,1),dtype=GF)
        for i in range(0,error_length):
            S_vec[i]=S[i]
        Y_vec=list()
      
        if(error_length==1):
            Y_vec.append( S_vec[0]*X_matrix[0][0].inverse())  
        else :Y_vec=self.gauss_jordan(X_matrix,S_vec)
        # Y_vec = np.linalg.lstsq(X_matrix, S_vec, rcond=None)[0]
        # print("ret",ret)
        r_true = list(received)
        for i in range(0,error_length):
            Ii=ret[i]
            r_true[Ii]-=Y_vec[i]
        return r_true



        




    






