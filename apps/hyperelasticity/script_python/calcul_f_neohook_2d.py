import sympy
from sympy import *
init_printing(use_unicode=True)

mu = sympy.Symbol("mu", positive=True)
lamba = sympy.Symbol("lamba", positive = True)
alpha = sympy.Symbol("alpha", positive = True)

x = sympy.Symbol("x")
y = sympy.Symbol("y")

# deplacement a renseigner
ux = (1/lamba + alpha) * x
uy = (1/lamba - alpha/(alpha+1)) * y + alpha * sin(pi*x) * (exp(x) - 1)

#calcul F
F = sympy.Matrix([[diff(ux,x) + 1, diff(ux,y)], [diff(uy,x), diff(uy,y)+1]])
print("F =")
pprint(F)
J = simplify(F.det())
print("J =")
pprint(J)
invF = simplify(F**-1)

#calcul PK1 (depend de la loi de comportement)
PK1 = mu * F + (lamba * ln(J) - mu )* invF.transpose()

fx = - (diff(PK1[0,0],x)+ diff(PK1[0,1],y))
fy = - (diff(PK1[1,0],x)+ diff(PK1[1,1],y))

#calcul force volumique
f = sympy.Matrix([[fx], [fy]])

print("force f =")
pprint(f)
