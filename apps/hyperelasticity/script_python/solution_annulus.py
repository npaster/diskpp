import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#solve annulus under imposed displacement
#y1 = phi
#y2 = dphi/dR

#R0 = inner radius
#R1 = outer radius

#r0 = imposed displacement

#mu, lamba = material constant

## for initial value mu =0.333 and r0 = 1.5
# lamba = 1.66644, y2_0 = 0.56974065999
# lamba = 16.6644, y2_0 =
# lamba = 166.644, y2_0 =
# lamba = 1666.44, y2_0 =
# lamba = 16664.4, y2_0 =
# lamba = 166644,  y2_0 =



def f(y, R, params):
    y1, y2 = y      # unpack current values of y
    mu, lamba = params  # unpack parameters
    f1 = y2
    t1 = mu*(1+1/np.power(y2, 2)) + (lamba*(1- np.log(y1*y2/R)))/np.power(y2, 2)
    f2 =  (lamba * np.log(y1*y2/R)*(1/y1 - 1/(R*y2)) + mu * (y1/(R*R) - 1/y1 + 1/(R*y2) - y2/R) -lamba*(1/y1 - 1/(R*y2)) )/t1
    derivs = [f1, f2]      # list of dy/dt=f functions
    return derivs

def Poo(y, R, params):
    y1, y2 = y      # unpack current values of y
    mu, lamba = params  # unpack parameters
    P = lamba * np.log(y1*y2/R)/(y1/R) + mu * (y1/R - R/y1)
    return P


def Prr(y, R, params):
    y1, y2 = y      # unpack current values of y
    mu, lamba = params  # unpack parameters
    P = lamba * np.log(y1*y2/R)/y2 + mu * (y2 - 1/y2)
    return P

# Parameters
mu = 0.333       # first materail coefficient
lamba =  1.66644      # second material coefficient
R0 = 0.5         # inner radius
R1 = 1           # outer radius


r0 = 0.8         # imposed displacement

# Initial values
y10 = r0        # initial deformation
y20 = 1.1      # initial derivate (hazard value)

# Bundle parameters for ODE solver
params = [mu, lamba]


# Make time array for solution
nb_point = 100
R = np.linspace(R0, R1, num=nb_point)
RInc = (R1 - R0)/nb_point

tole = 1E-12
a = 0.1
b = 1.1 # b = 0.2 pour r0 = 1.5,  b = 0.7 for r0 = 0.6

# Bundle initial conditions for ODE solver
y0 = [r0, y20]
# Call the ODE solver
psoln = odeint(f, y0, R, args=(params,), rtol = (1E-13), atol= 1E-13, hmin =1E-16, hmax = 1E-11, ixpr=(True),)

print("y1(R1)")
print(psoln[-1,0])
print("y2(R1)")
print(psoln[-1,1])

m = 0
#dichotomie to find an appropriate value such that Prr(R1) = 0
while np.abs(Prr(psoln[-1,:], R1, params)) >= tole:
    m = (a+b)/2
    # Bundle initial conditions for ODE solver
    y0 = [r0, m]
    # Call the ODE solver
    psoln = odeint(f, y0, R, args=(params,))
    #compute
    P = Prr(psoln[-1,:], R1, params)
    if(P > 0):
        b = m
    else:
        a = m
    print(a,b,P)




print("y1(R1)")
print(psoln[-1,0])
print("y2(R1)")
print(psoln[-1,1])


print("u(R0)")
print(psoln[0,0]-R0)

print("u(R1)")
print(psoln[-1,0]-R1)



##value Prr(R1)
Prr_R1 = Prr(psoln[-1,:], R1, params)
print("PRR(R1)")
print(Prr_R1)

Prr_sol = []
Poo_sol = []

fichier = open("sol_anneau.dat", "w")
fichier.write("#R\tPhi\tdPhi\tPrr\tPoo")
fichier.write("\n" + str(len(psoln)))

i = 0
for y in psoln:
    Rp = R[i]
    Pr = Prr(y, Rp, params)
    Po = Poo(y, Rp, params)
    Prr_sol.append(Pr)
    Poo_sol.append(Po)
    fichier.write("\n" +  str(Rp) + "\t" + str(y[0]) + "\t" + str(y[1]) + "\t" +  str(Pr) + "\t" +  str(Po) )
    i = i + 1

fichier.close()


#Plot results
fig = plt.figure(1, figsize=(8,8))
fig2 = plt.figure(2, figsize=(8,8))

# Plot phi as a function of R
ax1 = fig.add_subplot(311)
ax1.plot(R, psoln[:,0])
ax1.set_xlabel('R')
ax1.set_ylabel('phi')

## Plot dphi as a function of R
ax2 = fig.add_subplot(312)
ax2.plot(R, psoln[:,1])
ax2.set_xlabel('R')
ax2.set_ylabel('dphi')

# Plot PRR as a function of time
ax3 = fig2.add_subplot(311)
ax3.plot(R, Prr_sol)
ax3.set_xlabel('R')
ax3.set_ylabel('Prr')
#
# Plot Poo as a function of time
ax4 = fig2.add_subplot(312)
ax4.plot(R, Poo_sol)
ax4.set_xlabel('R')
ax4.set_ylabel('Poo')


plt.tight_layout()
plt.show()
