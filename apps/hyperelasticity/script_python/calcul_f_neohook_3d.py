import sympy
from sympy import *
init_printing(use_unicode=True)

mu = sympy.Symbol("mu", positive=True)
lamba = sympy.Symbol("lamba", positive=True)
beta = sympy.Symbol("beta", positive=True)
alpha = sympy.Symbol("alpha", positive=True)



x = sympy.Symbol("x")
y = sympy.Symbol("y")
z = sympy.Symbol("z")

# deplacement a renseigner
#indentation 3D
ux = alpha * (1/lamba + 1)*x + alpha*sin(pi*y)
uy =  (alpha+beta+alpha*beta)*(1/lamba - 1/(1+alpha+beta+alpha*beta))*y
uz = beta * (1/lamba+ 1)*z + beta*sin(pi*x) ;

#calcul F
F = sympy.Matrix([[simplify(diff(ux,x) + 1.0), simplify(diff(ux,y)), simplify(diff(ux,z))], [simplify(diff(uy,x)), simplify(diff(uy,y)+1.0), simplify(diff(uy,z))], [simplify(diff(uz,x)), simplify(diff(uz,y)), simplify(diff(uz,z)+1.0)]])
print("F =")
pprint(F)
J = simplify(F.det())
pprint(J)
invF = simplify(F**-1)

#calcul PK1 (depend de la loi de comportement)
PK1 = mu * F + (lamba * ln(J) - mu )* invF.transpose()

fx = - (diff(PK1[0,0],x)+ diff(PK1[0,1],y) + diff(PK1[0,2],z))
fy = - (diff(PK1[1,0],x)+ diff(PK1[1,1],y) + diff(PK1[1,2],z))
fz = - (diff(PK1[2,0],x)+ diff(PK1[2,1],y) + diff(PK1[2,2],z))

#calcul force volumique
f = sympy.Matrix([[simplify(fx)], [simplify(fy)], [simplify(fz)]])

print("force f =")
pprint(f)
