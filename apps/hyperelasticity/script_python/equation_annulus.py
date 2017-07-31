
import sympy as sy
from sympy import *


R = sy.symbols("R", real=True, positive=True)
f = sy.symbols("f", function=True)
f(R)
f_R = f(R).diff(R)

mu = sy.Symbol("mu", positive=True)
lamda = sy.Symbol("lamba", positive=True)

Prr = sy.symbols("Prr", function=True)

#Prr = (lamda * log(f*f_R/R,10))/f_R + mu *(-1/f_R + f_R)
Prr = (lamda * log(f(R)*f_R/R))/f_R + mu *(-1/f_R + f_R)

Prr_R = factor(simplify(diff(Prr,R)))

#Poo = (lamda * log(f*f_R/R,10))/(f/R) + mu *(-R/f + f/R)

Poo = (lamda * log(f(R)*f_R/R))/(f(R)/R) + mu *(-R/f(R) + f(R)/R)


#pprint(diff(log(f*diff(f,R))

print("Prr")
pprint(Prr)

print("Poo")
pprint(Poo)

print("Prr_t")
pprint(Prr_R)


rhs = simplify((Prr -Poo)/R)

eqdiff = simplify(Prr_R + rhs)



print(" ")
print("rhs")
pprint(rhs)

print(" ")
print("Eq diff")
pprint(eqdiff)
