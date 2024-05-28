#!/usr/bin/python3 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def eqinfi(x, n, a, t):
    m = 1
    h = 1        
    E = ((np.pi*n*h)**2)/(2*m*a**2)
    return np.sqrt(0.5)*np.sqrt(2/a) * np.sin(n * x * np.pi / a) * np.exp(-t*E/h)

def edelta(x,n,a,t):
    m = 1
    h = 1
    #a = 1
    E = -(m*a**2)/(2*h**2)
    return (np.sqrt(m*a)/h)*np.exp(-m*a*abs(x)/(h**2))*np.exp(-t*E/h)

def gausfree(x,n,a,t):
    m = 1
    h = 1
    norm = (2*a/np.pi)**(0.25)
    ekth = (-a*(x**2))/(1+(2*h*a*t/m))
    paronom = np.sqrt(1+(2*h*a*t/m)) 
    return norm*np.exp(ekth)/paronom

def fpfree(p,n,a,t):
    h =1 
    m = 1
    k = p/h
    ekth = -k**2/4*a
    paronom = (2*np.pi*a)*0.25
    return np.exp(ekth)/paronom

def animwf(f,xin,xf,time,nmax,fast):
    for n in range(1,nmax):
        # Set up the figure, axis, and plot element
        fig, ax = plt.subplots()
        xi = np.linspace(xin,xf, 100)
        line, = ax.plot(xi, f(xi, n, 10, 0))

        def init():
            # Initial plot setup
            line.set_ydata(np.ma.array(xi, mask=True))
            return line,

        def update(t):
            # Update the plot with new data
            y = f(xi, n, 10, t)
            line.set_ydata(y)
            ax.set_title(f"t = {t:.2f}")
            return line,

        # Creating the animation
        plt. grid()
        plt.title("infinite square dynamic")
        plt.xlabel("x")
        plt.ylabel("Ψ(x,t)")
        plt.axvline(x=0,color = "black")
        plt.axhline(y=0,color = 'black')
        ani = FuncAnimation(fig, update, frames=np.linspace(0, time, fast), init_func=init, blit=True, repeat=True)
        plt.show()

animwf(eqinfi,0,10,100,5,50)
animwf(edelta,-1,1,10,2,10**6)
animwf(gausfree,-10,10,1000,2,500)
animwf(fpfree,-10,10,1,2,50**6)
