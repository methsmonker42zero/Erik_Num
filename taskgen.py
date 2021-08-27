from num import *
import random
import json

def genchol(n, value_range=32):
    vals = [random.random() for i in range (n**2)]
    A = (vals, n)
    A = sc(value_range, a(sc(.5, a(A, transp(A))), diag((n,)*n)))
    A = roundm(A,0)
    x = vecv([random.randint(-value_range, value_range) for i in range (n)])
    b = matmul(A, x)
    L = chol(A)
    print('A=\n%s'%strm(A))
    print('b=\n%s'%strm(b))
    print('L=\n%s'%strm(L))
    print('x=\n%s'%strm(x))
