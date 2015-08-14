import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import animation


def kinematics(x, t, m, g):
    n = int(len(x)/4)
    #m = np.ones(n)

    out = np.zeros(x.shape)

    for i in range(n):
        out[2*i:2*(i+1)] = x[2*(n + i):2*(n + i + 1)]/m[i]

        xi = x[2*i:2*(i+1)]
        for j in range(n):
            xj = x[2*j:2*(j+1)]
            if i == j:
                out[2*(n + i):2*(n + i + 1)] += np.zeros(2)
            else:
                r = xj - xi
                out[2*(n + i):2*(n + i + 1)] += m[i]*m[j]*g* r/np.linalg.norm(r)**2
    return out


class Constelation:
    def __init__(self, n, g, dt, factor):
        self.n = n
        self.g = g
        self.dt = dt
        self.m = np.random.exponential(scale=0.3, size=n)
        self.x = np.random.rand(4*n)
        self.x[2*n:] = np.zeros(2*n)
        self.factor = factor
        self.patches = [plt.Circle((self.x[2*j],self.x[2*j + 1]), np.sqrt(self.factor*self.m[j]/np.pi), facecolor=np.random.rand(3,1), visible=False) for j in range(n)]

    def move(self):
        dx = kinematics(self.x, 0.0, self.m, self.g)
        self.x += dx*self.dt

    def init_animation(self):
        return self.patches

    def animate(self, i):
        if i == 1:
            for p in self.patches:
                p.set_visible(True)

        self.move()
        for j in range(n):
            self.patches[j].center = (self.x[2*j], self.x[2*j + 1])

        return self.patches

if __name__ == '__main__':
    n = 10
    g = 1e-2
    dt = 0.04
    factor = 1e-2
    time = 100
    interval = 40

    system = Constelation(n, g, dt, factor)
    print(system.m)
    print(system.x)

    fig = plt.figure()
    limit = 2
    plt.axis([-limit,limit,-limit,limit])
    ax = plt.gca()
    ax.set_aspect(1)

    for p in system.patches:
        ax.add_patch(p)

    anim = animation.FuncAnimation(fig, system.animate, init_func=system.init_animation, frames=int(time/dt), interval=interval, blit=True, repeat=True)

    #plt.plot(x[:2*n:2], x[1:2*n:2], '*')
    #for i in range(n):
    #    plt.plot(trace[:,2*i], trace[:,2*i+1], '.-', markersize=np.sqrt(m[i]/np.pi)*20)
    figManager = plt.get_current_fig_manager()
    figManager.resize(*figManager.window.maxsize())
    plt.show()
