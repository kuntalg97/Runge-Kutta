# Kuntal Ghosh
# June 2021

import math

t0 = 0.0
x0 = 0.0
v0 = -1.0
tmax = 500.0
dt = 0.1
k = 1.0
m = 1.0

n = int((tmax-t0)/dt)

def slope (aa,xv):
    if aa == 1:
        deriv = (-k/m)*xv
    else:
        deriv = xv
    return deriv

def analytical (t,x,v):
    file1 = open ('analytical.dat','w')

    w = math.sqrt(k/m)

    for i in range(1,n+1):        
        x = (-1.0/w)*math.sin(w*t)
        v = -math.cos(w*t)
        t += dt

        file1.write ("%16.6f%16.6f%16.6f\n"%(t,x,v))    

def rk2 (t,x,v):
    file2 = open ('rk2_spring.dat','w')
    for i in range(1,n+1):
        ke = 0.5*m*v**2
        pe = 0.5*k*x**2
        te = ke + pe
        file2.write ("%16.6f%16.6f%16.6f%16.6f%16.6f%16.6f\n"%(t,x,v,ke,pe,te))

        k1v = dt*slope(1,x)
        k1r = dt*slope(0,v)
        k2v = dt*slope(1,x+k1r/2.0)
        k2r = dt*slope(0,v+k1v/2.0)

        v += k2v
        x += k2r
        t += dt

def rk4 (t,x,v):
    file3 = open ('rk4_spring.dat','w')
    for i in range(1,n+1):
        ke = 0.5*m*v**2
        pe = 0.5*k*x**2
        te = ke + pe
        file3.write ("%16.6f%16.6f%16.6f%16.6f%16.6f%16.6f\n"%(t,x,v,ke,pe,te))

        k1v = dt*slope(1,x)
        k1r = dt*slope(0,v)
        k2v = dt*slope(1,x+k1r/2.0)
        k2r = dt*slope(0,v+k1v/2.0)
        k3v = dt*slope(1,x+k2r/2.0)
        k3r = dt*slope(0,v+k2v/2.0)
        k4v = dt*slope(1,x+k3r)
        k4r = dt*slope(0,v+k3v)

        v += (k1v+k4v)/6.0 + (k2v+k3v)/3.0
        x += (k1r+k4r)/6.0 + (k2r+k3r)/3.0
        t += dt

def main():
    order = int(input ("Enter the order of the Runge Kutta scheme = "))

    t = t0
    x = x0
    v = v0
    analytical (t,x,v)

    t = t0
    x = x0
    v = v0    
    if order == 2:        
        rk2 (t,x,v)
    elif order == 4:
        rk4 (t,x,v)

main()
