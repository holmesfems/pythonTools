def matmul(a,b):
    ret=[]
    for i in range(0,len(a)):
        line=[]
        for j in range(0,len(b[0])):#len(a[0])=len(b)
            s=0
            for k in range(0,len(b)):
                s+=a[i][k]*b[k][j]
            line.append(s)
        ret.append(line)
    return ret

def nummul(c,a):
    ret=[]
    for line in a:
        add=[]
        for ele in line:
            add.append(c*ele)
        ret.append(add)
    return ret

def matadd(a,b):
    ret=[]
    for i in range(0,len(a)):
        line=[]
        for j in range(0,len(a[0])):
            line.append(a[i][j]+b[i][j])
        ret.append(line)
    return ret

def unitmat(rank):
    ret=[]
    for i in range(0,rank):
        line=[]
        for j in range(0,rank):
            if i==j:
                line.append(1.0)
            else:
                line.append(0.0)
        ret.append(line)
    return ret

def cofmat(a,i,j):
    ret=[]
    for i1 in range(0,len(a)):
        if i1!=i:
            line=[]
            for j1 in range(0,len(a[0])):
                if j1!=j:
                    line.append(a[i1][j1])
            ret.append(line)
    return ret

def det(a):
    if len(a)!=len(a[0]):
        return
    if len(a)==1:
        return a[0][0]
    s=0
    sig=1
    for i in range(0,len(a)):
        if a[0][i]!=0.0:
            s+=sig*det(cofmat(a,0,i))*a[0][i]
        sig*=-1
    return s*1.0

def matinv(a):
    d=det(a)
    ret=[]
    sig0=1
    for i in range(0,len(a)):
        line=[]
        sig=sig0
        for j in range(0,len(a[0])):
            line.append(det(cofmat(a,i,j))*sig/d)
            sig*=-1
        ret.append(line)
        sig0*=-1
    return turn(ret)

def turn(a):
    ret=[]
    for i in range(0,len(a[0])):
        line=[]
        for j in range(0,len(a)):
            line.append(a[j][i])
        ret.append(line)
    return ret

def diag(a):
    ret=[]
    for i in range(0,len(a)):
        line=[]
        for j in range(0,len(a[0])):
            if i==j:
                line.append(a[i][j])
            else:
                line.append(0.0)
        ret.append(line)
    return ret

def printmat(a):
    for i in range(0,len(a)):
        print a[i]

def demo():
    a=[
        [0,1,1],
        [1,1,0],
        [1,2,0]
        ]
    m=matinv(a)
    printmat(m)
    printmat(matmul(m,a))

if __name__=='__main__':
    demo()
