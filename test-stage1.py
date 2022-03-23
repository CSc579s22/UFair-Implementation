from sympy import symbols, Eq, nsolve
import numpy as np
M=[[-3.035,-0.5061,1.022],[-4.85,-0.647,1.011],[-17.53,-1.048,0.9912]]
#M=[[-1,-0.5,1],[-1,-0.6,1],[-10,-1,1]]
r_max=[20,33,45,22,33]
BW=40
N=5
res=[720,360,1080,720,360]
u=[]

r = list(symbols('r0:%d' % N))
print("vector of r= ",r)

for i in range(N):
    if res[i]==1080:
        a=M[0][0]
        b=M[0][1]
        c=M[0][2]
        frac= a * (r_max[i]**b) + c
        u.append((a * (r[i]**b) + c) / frac)
    elif res[i]==720:
        a = M[1][0]
        b = M[1][1]
        c = M[1][2]
        frac = a * (r_max[i] ** b) + c
        u.append((a * (r[i]**b) + c) / frac)
    elif res[i]==360:
        a = M[2][0]
        b = M[2][1]
        c = M[2][2]
        frac = a * (r_max[i] ** b) + c
        u.append((a * (r[i]**b) + c) / frac)
#print(u)

eq2=Eq(sum(r),BW)
eq=[]
eq.append(eq2)
for i in range(N-1):
    e1=u[i]-u[i+1]
    eq.append(Eq(e1,0))
print("vector of equations= ",eq)
ov= np.ones(N)
sol= nsolve(eq,r,ov)
print("optimal r= ",sol)


