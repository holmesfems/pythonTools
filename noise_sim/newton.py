def arraycopy(ar):
    ret=[]
    for item in ar:
        ret.append(item)
    return ret

def subtAr(a1,a2):
    ret=[]
    for i in range(0,len(a1)):
        ret.append(a1[i]-a2[i])
    return ret

def sumAr(a1,a2):
    ret=[]
    for i in range(0,len(a1)):
        ret.append(a1[i]+a2[i])
    return ret

def abs2Ar(ar):
    ret=0
    for item in ar:
        ret+=item**2
    return ret

"""
partial_derivation of "A_n" for 
scalar function of vector:F(A)
where "A" is vector of [A_i],i=0,1,2...N-1
\frac{\partial F}{\partial A_n}|_A=
\frac{F(A^+)-F(A^-)}{d}
A^+=A+[A:An<-An+0.5d]
A^-=A+[A:An<-An-0.5d]

"""
def partial_deriv(F,A,n,d):
    A1=arraycopy(A)
    A2=arraycopy(A)
    A1[n]=A1[n]+0.5*d
    A2[n]=A2[n]-0.5*d
    F1=F(A1)
    F2=F(A2)
    return (F1-F2)/d

def solve_minimize(F,A0,partial_d,err,count=0):
    MAXCOUNT=1000
    if count>=MAXCOUNT:
        return None
    N=len(A0)
    p1=[]
    p2=[]
    for i in range(0,N):
        p1.append(partial_deriv(F,A0,i,partial_d))
    for i in range(0,N):
        """
        def partial(A):
            return partial_deriv(F,A,i,partial_d)
        """
        p2.append(partial_deriv(lambda x:partial_deriv(F,x,i,partial_d),A0,i,partial_d))
    d=[]
    for i in range(0,N):
        d.append(-1.0*p1[i]/p2[i])
    A1=sumAr(d,A0)
    if abs2Ar(d)<=err:
        return A1
    else:
        return solve_minimize(F,A1,partial_d,err,count=count+1)

def demofunc(A):
    return A[0]**2+(A[1]-3)**2+(A[2]-1)**2

