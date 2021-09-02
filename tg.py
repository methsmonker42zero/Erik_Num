from utils.meth.num import *
import random
import json

value_range = 4**2

def random_int_ltri_matrix(n, value_range=value_range, no_zeros=1, reg=1):
    from utils.meth.rng import randint, randintnz
    return ([e for i in range(n)for e in
    [
        randintnz(-value_range,value_range)if no_zeros or j==i and reg else randint(-value_range,value_range)
        for j in range(i+1)
    ]+[0]*(n-i-1)], n)

def random_ltri_matrix(n, value_range=value_range, pos_only=0):
    from utils.meth.rng import random
    return ([e for i in range(n)for e in
    [
        (random()if pos_only else random()*2-1)* value_range
        for j in range(i+1)
    ]+[0]*(n-i-1)], n)

def random_int_posdef_matrix_using_int_tri(n, value_range=value_range, no_zeros=1):
    L = random_int_ltri_matrix(n, value_range=int(value_range**.5), no_zeros=no_zeros)
    return matmul(L, transp(L))

def random_int_posdef_matrix_using_rounding(n, value_range=value_range):
    L = random_tri_matrix(n, value_range**.5, pos_only=1)
    A = matmul(L, transp(L))
    A_rounded = a(intm(A), iden(n))
    return A, A_rounded

def random_int_posdef_matrix_using_chol(n, value_range=value_range, no_zeros=1, explicit=1):
    from utils.meth.rng import randint, randintnz
    A = Symm.__to__(empty_mat(n))
    L = Symm.__to__(empty_mat(n))
    punish = 0
    for i in range(n):
        overceil_punish = int(punish)+1
        a_i_i = randint(overceil_punish, max(value_range, overceil_punish))
        A[i,i] = a_i_i
        l_i_i = (a_i_i - punish)**.5
        L[i,i] = l_i_i
        if explicit:
            print("A's %s. diagonal entry set to %s" % (i,a_i_i))
            print("L's %s. diagonal entry set to %s" % (i,l_i_i))
        if i == n-1: break
        a_next_side = [randint(no_zeros,value_range)for j in range(i+1)]
        A[i+1,:i+1] = a_next_side
        l_next_side = forward_sub(L[:i+1,:i+1], a_next_side)
        L[i+1,:i+1] = l_next_side
        punish = dot(l_next_side, l_next_side)
        if explicit:
            print("A's %s. row set to %s" % (i+1,a_next_side))
            print("L's %s. row forced to %s" % (i+1,l_next_side))
            print("A's %s. punisher value is %s" % (i+1,punish))

    return A

def genchol(n, value_range=value_range, no_zeros=1, explicit=0, debug=1):
    from utils.meth.rng import randint
    global A,x,b,L,y
    A,x,b,L,y = (None,)*5
    try:
        fuck = 0
        A = random_int_posdef_matrix_using_chol(n, value_range=value_range)
        a = [randint(-value_range, value_range) for i in range (n)]
        x = vecv([randint(-value_range, value_range) for i in range (n)])
        b = matmul(A, x)
        L = chol(A)
        y = forward_sub(L, b)
    except:
        fuck = debug
        raise
    finally:
        if fuck:
            print("genchol fail")
            print(A,x,b,L,y)
        if fuck or explicit:
            print('A=\n%s'%strm(A))
            print('b=\n%s'%strm(b))
            print('='*13)
            print('L=\n%s'%strm(L))
            print('x=\n%s'%strm(x))
            print('y=\n%s'%strm(y))
        if fuck: print('-'*13)

if __name__=='__main__':
    for i in range(100000):
        A,AR = random_int_posdef_matrix_using_rounding(3)
        if not chol(AR,verify_posdef=1):
            printm(A,AR,names=['A','AR'])
            printm(*fund_diagon(A),names=['A diag form','A diagonalizer'])
            printm(*fund_diagon(AR),names=['AR diag form','AR diagonalizer'])
