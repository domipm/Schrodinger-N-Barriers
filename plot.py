import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.rcParams["animation.convert_path"] = "/usr/bin/convert"

pot = np.loadtxt("pot.dat") #POTENCIAL (N,V)
fonda = np.loadtxt("fonda.dat") #FUNCION DE ONDA (N, Re F, Im F) Cada N datos tiempo T+1

N = int(len(pot)) #DIVISIONES DE X
iter = int(len(fonda)/N) #ITERACIONES REALIZADAS

skip = 10

x = np.arange(0,N,1) #EJE X

fig, ax = plt.subplots()
data = []
ln, = plt.plot([], [], linestyle='-', color="Red", label="$Wavefunction \\phi$")

def init():

    ax.set_title("Schrödinger N-Barrier Potential")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_xlim(-1, N+1)
    ax.set_ylim(0, 1)

    return ln,

def update(k):

    #k = numero iteracion

    data.clear()
    
    for i in range(N):
        data.append(fonda[k*N+i])

    ln.set_data(x, data)

    return ln,

plt.plot(x, pot, label="Potential $V$")
plt.legend(loc="upper right")

anim = animation.FuncAnimation(fig, update, frames=np.arange(0,iter,skip),init_func=init,blit=False)
anim.save("schrodinger.gif", writer=animation.PillowWriter(fps=30), dpi=300)
#plt.show()