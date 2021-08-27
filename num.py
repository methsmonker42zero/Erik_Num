#code section

# n,m=col_count,row_count
# x,y=col_index,row_index

def s(a,n,x,y,new): a[x%n + y*n] = new
def ss(anm,x,y,new): s(*anm[:2],x,y,new)
g = lambda a,n,x,y: a[x%n + y*n]
gg = lambda anm,x,y: g(*anm[:2],x,y)
a_ = lambda a,b: [a[i]+b[i]for i in range(len(a))]
def a (anm,bnm):
    anm, bnm = s3(anm), s3(bnm)
    assert anm[1:]==bnm[1:]
    return (a_(anm[0],bnm[0]),*anm[1:])
xy = lambda n,i: (i%n,i//n)
dot = lambda a,b: sum(a[i]*b[i]for i in range(len(a)))
sc = lambda s, anm: ([a*s for a in anm[0]],*anm[1:])
sc_ = lambda s,a: [a*s for a in a]
col = lambda a,n,m,i: [a[k*n+i] for k in range(m)]
row = lambda a,n,i: a[i*n:i*n+n]
def srow(a,n,i,new): a[i*n:i*n+len(new)] = new
eq = lambda anm, bnm: s3(anm)==s3(bnm)
s3 = lambda a: a if len(a)==3 else (a[0],)+(a[1],)*2
def s2s(a): #s2strict, since it checks validity
    if len(a)==3:
        assert a[1]==a[2]
        return a[:2]
    elif len(a)==1:
        sqrtlen = len(a[0])**.5
        assert sqrtlen.is_integer()
        return (*a, sqrtlen)
    elif len(a)==2:
        return a
    raise Exception
vecv = lambda a: (a,1,len(a))if type(a)is list else a
vech = lambda a: (a,len(a),1)if type(a)is list else a
copy = lambda anm: (anm[0].copy(),*anm[1:])
rev = lambda anm: (list(reversed(anm[0])),*anm[1:])
roundm = lambda anm,precision=3: ([round(x,precision)for x in anm[0]],*anm[1:])

#syntactical flowing honey nectar
class Matrix(tuple):
    def _new__(s, values, n, m=None):
        if m is None:
            return tuple.__init__(Matrix, (values, n,))
        else:
            return tuple.__init__(Matrix, (values, n, m,))

def matmul(anm,bon):
    a,n,m = s3(anm)
    b,o,n = s3(bon)
    c = [None]*(m*o)
    for i in range(o):
        for k in range(m):
            s(c,o,i,k, dot(row(a,n,k),col(b,o,n,i)))
    return c,o,m

def empty_mat(n, default=0, m=None):
    if m is None: m = n
    return [default]*(n*m),n,m

iden = lambda n: diag((1,)*n)
# creates matrix with given diagonal_entries
def diag(diag_entries,default=0,base=None):
    n = len(diag_entries)
    r = [default]*n**2if base is None else base
    assert len(r)==n**2
    for i in range(n): r[i*(n+1)] = diag_entries[i]
    return r, n

# <(=): lower triangle
# >(=): upper triangle
# masks a matrix entries with (column_index condition row_index)
def mask(anm,cond,default=0,base=None):
    a,n,m = s3(anm)
    b = [default]*len(a) if base is None else base
    assert len(a)==len(b)
    for i in range(len(a)):
        if eval(('%s'+cond+'%s')%xy(n,i)):
            b[i] = a[i]
    return b,n,m

#transpose matrix
def transp(anm,base=None):
    a,n,m = s3(anm)
    b = [None]*len(a) if base is None else base
    assert len(a)==len(b)
    for i in range(n):
        for k in range(m):
            s(b,m,k,i,g(a,n,i,k))
    return b,m,n

#only works for real matrices
def strm(*ans, to=3):
    from math import log, ceil
    if not ans: return
    ans = [s3(an)for an in ans]
    la = ans[0][2]
    cw = max(4, ceil(log(max(abs(v) for an in ans for v in an[0]), 10))) # cell width (in characters) = cw + 2
    r = ''
    i = 0
    while i < la:
        for a,n,m in ans:
            for v in a[i*n:i*n+n]:
                if v==int(v): v = int(v)
                a = ('{:.%sg}'%cw).format(v)
                a = a[1:]if a[0]=='0'else a
                r += a.rjust(cw+2)
            r += ' | '
        r += '\n'
        i += 1
    return r

def printm(*ans, to=3):
    print(strm(*ans, to=to))

def swap_row(an,first,second):
    an = s2s(an)
    b = row(*an,first)
    srow(*an, first, row(*an, second))
    srow(*an, second, b)

def add_row(an,src,tar,factor):
    an = s2s(an)
    srow(*an,tar,a_(row(*an,tar),scrow(factor,row(*an,src))))

# a)
def highest_abv(an, layer):
    a,n = s2s(an)
    highest_abv = highest_abv_i = -1
    for i in range(layer, n):
        abv = abs(g(a,n,layer,i))
        if abv > highest_abv:
            highest_abv = abv
            highest_abv_i = i
    return highest_abv_i

# b)
def necessary(an, layer):
    while g(*s2s(an),layer,layer)==0: layer += 1
    return layer

def plr(an, permute=highest_abv, explicit=0):
    an = copy(s2s(an))
    a,n = an
    A = (a.copy(),n)
    perm_v = ([i for i in range(n)],1,n)
    if explicit:
        printm(perm_v,an)
    for layer in range(n-1):
        i = permute(an, layer)
        swap_row(an, layer, i)
        swap_row(perm_v, layer, i)
        if explicit:
            printm(perm_v,an)
        for j in range(layer+1,n):
            factor = - g(a,n,layer,j) / g(a,n,layer,layer)
            add_row(an,layer,j,factor)
            assert g(a,n,layer,j) == 0
            s(a,n,layer,j,-factor)
            if explicit:
                print('row %s <- row %s + %s * row %s' % (j, j, factor, layer))
                print('additive inverse of factor %s stored in cell %s,%s (row,column)' % (factor,j,layer))
                printm(an)

    L = mask(an,'<',base=iden(n)[0])
    if explicit:
        print('L =')
        printm(L)

    R = mask(an,'>=')
    if explicit:
        print('R =')
        printm(R)

    P = [[], n]
    for permv in perm_v[0]: P[0] += [0if i!=permv else 1for i in range(n)]
    if explicit:
        print('P =')
        printm(P)

    if explicit:
        pa = matmul(P,A)
        lr = matmul(L,R)
        assert pa==lr
        print('PA = LR =')
        printm(lr)

    return P, L, R

def chol(an, explicit=0):
    an = copy(s2s(an))
    a, n = an
    L = s2s(empty_mat(n))
    s(*L,0,0,g(a,n,0,0)**.5)
    if explicit:
        print('cholesky start')
        printm(L)
    for i in range(1,n):
        y = []
        for k in range(i):
            b = g(a,n,i,k)
            ssum = 0
            for j in range(k): ssum += g(*L,j,k) * y[j]
            y.append((b-ssum) / g(*L,k,k))
        srow(*L,i,y)
        s(*L,i,i,(g(a,n,i,i)-dot(y,y))**.5)
        if explicit:
            print('cholesky step')
            printm(L)
    LT = transp(L)
    LLT = matmul(L, LT)
    #assert eq(an,LLT) rounding kekw
    if explicit:
        print('L =')
        printm(L)
        print('L^T =')
        printm(LT)
        print('A = LL^T =')
        printm(LLT)
    return L

def solve_chol(an, b):
    L = chol(an)
    return solve(transp(L), solve(L, b))

def solve(LorR, b):
    b = vecv(b)
    n = LorR[1]
    if any(gg(LorR, i, j) for i in range(n) for j in range(i)):
        solve = backward_sub
    else:
        solve = forward_sub
    return solve(LorR, b)

def backward_sub(R, b):
    return rev(forward_sub(rev(R), rev(b)))

def forward_sub(L, b):
    b = vecv(b)
    n = L[1]
    x = empty_mat(n=1,m=n)
    for i in range(n):
        subway = sum(gg(L, j, i)*gg(x, 1, j) for j in range(i))
        ss(x, 1, i, (gg(b, 1, i) - subway) / gg(L, i, i))
    return x
