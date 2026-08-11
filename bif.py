#Viktoriia Volkova
from __future__ import division
import sys
from math import sqrt

def read_input():
    data = sys.stdin.read().split()
    all_nums = [float(x) for x in data]
    lam = int(all_nums[-1])
    d = int(all_nums[-2])
    rest = all_nums[:-2]
    n = 1
    while True:
        k = (n+3)*(n+2)//2
        if (n+1) + n*k == len(rest):
            break
        n += 1
        if n > 20:
            raise Exception("n too large")
    p0 = list(rest[:n+1])
    coeffs = []
    off = n+1
    k = (n+3)*(n+2)//2
    for i in range(n):
        coeffs.append(list(rest[off:off+k]))
        off += k
    return n, p0, coeffs, d, lam

def basis_idx(m):
    b = [[]]
    for i in range(m): b.append([i])
    for i in range(m):
        for j in range(i,m): b.append([i,j])
    return b

def eval_poly(coeffs, p, n):
    m = n+1
    basis = basis_idx(m)
    mono = [1.0 if not bk else (p[bk[0]] if len(bk)==1 else p[bk[0]]*p[bk[1]]) for bk in basis]
    return [sum(c*mv for c,mv in zip(row,mono)) for row in coeffs]

def eval_F(coeffs, p, n):
    fp = eval_poly(coeffs, p, n)
    return [fp[i]-p[i] for i in range(n)]

def jac_f(coeffs, p, n):
    m = n+1
    basis = basis_idx(m)
    J = [[0.0]*m for _ in range(n)]
    for i, row in enumerate(coeffs):
        for c, b in zip(row, basis):
            if c == 0.0: continue
            for j in range(m):
                cnt = b.count(j)
                if cnt == 0: continue
                val = c * cnt
                rem = b[:]
                rem.remove(j)
                for idx in rem: val *= p[idx]
                J[i][j] += val
    return J

def jac_F(coeffs, p, n):
    J = jac_f(coeffs, p, n)
    for i in range(n): J[i][i] -= 1.0
    return J

def gauss(A, b):
    n = len(b)
    M = [A[i][:]+[b[i]] for i in range(n)]
    for col in range(n):
        best = max(range(col,n), key=lambda r: abs(M[r][col]))
        M[col],M[best] = M[best],M[col]
        piv = M[col][col]
        if abs(piv) < 1e-14: continue
        for row in range(n):
            if row==col: continue
            f = M[row][col]/piv
            for k in range(col,n+1): M[row][k] -= f*M[col][k]
        for k in range(col,n+1): M[col][k] /= piv
    return [M[i][n] for i in range(n)]

def nrm(v): return sqrt(sum(x*x for x in v))
def nmlz(v): s=nrm(v); return [x/s for x in v]

def tangent(J, n):
    best = None; best_n = -1.0
    for j in range(n+1):
        sub = [[J[r][c] for c in range(n+1) if c!=j] for r in range(n)]
        rhs = [-J[r][j] for r in range(n)]
        x = gauss(sub, rhs)
        u = x[:j]+[1.0]+x[j:]
        nu = nrm(u)
        if nu > best_n: best_n=nu; best=u
    return nmlz(best)

def newton(coeffs, q, fi, n):
    p = q[:]
    for _ in range(60):
        Fv = eval_F(coeffs, p, n)
        if nrm(Fv) < 1e-15: break
        J = jac_F(coeffs, p, n)
        Jr = [[J[r][c] for c in range(n+1) if c!=fi] for r in range(n)]
        dx = gauss(Jr, [-f for f in Fv])
        dx = dx[:fi]+[0.0]+dx[fi:]
        p = [p[k]+dx[k] for k in range(n+1)]
    return p

def eigs_1(M): return [M[0][0]]

def eigs_2(M):
    a,b,c,d = M[0][0],M[0][1],M[1][0],M[1][1]
    tr=a+d; disc=max(0.0,(a-d)**2+4*b*c)
    sq=sqrt(disc)
    return sorted([(tr-sq)/2,(tr+sq)/2])

def eigs_qr(M0, n):
    A = [r[:] for r in M0]
    for _ in range(3000):
        cols = [[A[r][c] for r in range(n)] for c in range(n)]
        orth = []; R = [[0.0]*n for _ in range(n)]
        for j in range(n):
            v = cols[j][:]
            for k,q in enumerate(orth):
                dt = sum(v[r]*q[r] for r in range(n))
                R[k][j]=dt; v=[v[r]-dt*q[r] for r in range(n)]
            nv=nrm(v); R[j][j]=nv
            if nv>1e-15: v=[x/nv for x in v]
            else: v=[0.0]*n; v[j]=1.0
            orth.append(v)
        newA=[[sum(R[i][k]*orth[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        A=newA
    return sorted([A[i][i] for i in range(n)])

def df_eigs(coeffs, p, n, lam):
    Jf = jac_f(coeffs, p, n)
    Jsq = [r[:n] for r in Jf]
    if n==1: evs=eigs_1(Jsq)
    elif n==2: evs=eigs_2(Jsq)
    else: evs=eigs_qr(Jsq,n)
    return evs

def closest_ev(evs, lam):
    return min(evs, key=lambda e: abs(e-lam))

def main():
    n, p0, coeffs, d0, lam = read_input()
    delta = 0.005
    p = p0[:]
    J = jac_F(coeffs, p, n)
    u = tangent(J, n)
    if u[-1]*d0 < 0: u=[-x for x in u]
    
    evs = df_eigs(coeffs, p, n, lam)
    prev_sign = 1 if closest_ev(evs,lam)-lam > 0 else -1
    prev_p = p[:]
    
    for _ in range(200000):
        q = [p[k]+delta*u[k] for k in range(n+1)]
        fi = max(range(n+1), key=lambda k: abs(u[k]))
        pe = newton(coeffs, q, fi, n)
        J2 = jac_F(coeffs, pe, n)
        v = tangent(J2, n)
        if sum(a*b for a,b in zip(u,v)) < 0: v=[-x for x in v]
        
        evs = df_eigs(coeffs, pe, n, lam)
        sv = 1 if closest_ev(evs,lam)-lam > 0 else -1
        
        if sv != prev_sign:
            pa, pb = prev_p[:], pe[:]
            for _ in range(100):
                pm = [(pa[k]+pb[k])/2.0 for k in range(n+1)]
                Jm = jac_F(coeffs, pm, n)
                vm = tangent(Jm, n)
                fim = max(range(n+1), key=lambda k: abs(vm[k]))
                pm = newton(coeffs, pm, fim, n)
                em = df_eigs(coeffs, pm, n, lam)
                sm = 1 if closest_ev(em,lam)-lam > 0 else -1
                if sm == prev_sign: pa=pm[:]
                else: pb=pm[:]
            evs_out = sorted(df_eigs(coeffs, pb, n, lam))
            print(' '.join('%.17e'%x for x in pb))
            print(' '.join('%.17e'%x for x in evs_out))
            return
        
        prev_sign=sv; prev_p=pe[:]; p=pe[:]; u=v[:]

main()
