import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk
import matplotlib.animation as animation

# Parametros de la teoria que hicimos sobre los pendulos
m1 = 1               # Masa de los pendulos
g = 9.81             # Aceleracion de la gravedad
T = 60               # Tiempo hasta que se vuelven a poner en fase
K_0 = 24             # Cantidad de oscilaciones del primer pendulo hasta que vuelven a estar en fase
                     # (tomo 21 para que el largo maximo sean 1.5 metros)
f_0 = 2*np.pi*K_0/T  # Frecuencia del primer pendulo
s = 2*np.pi/T        # Delta frecuencia
N = 12               # Numero de pendulos
np.random.seed(1000) # Seteo la semilla
error = 0.02
errorEnL = np.random.rand(N)*error*2 - error

L = []
def L_n (n):         # Largo del n-esimo pendulo
	return g/(f_0 + n*s)**2
for n in range(N):
	L.append(L_n(n))
print (L, errorEnL)
L = np.array(L) + errorEnL
#L = L * (1 + errorEnL)

def f (L):           # Frecuencia en funcion del largo
	return ((g/L)**(1/2))

# Calculo theta para todo t
theta0 = np.radians(20)
dt = 1
Tiempo = np.arange(0, (T+dt)/1, dt)
theta = []
for n in range(N):
	theta.append(theta0 * np.cos(f(L[n]) * Tiempo))

# Preparo cosas para la simulacion
nsteps = len(Tiempo)

def get_coords(th, L):
    """Return the (x, y) coordinates of the bob at angle th."""
    return L * np.sin(th), -L * np.cos(th)

fig = plt.figure()
bob_radius = max(L)/50
ax = fig.add_subplot(aspect='equal')
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)  # Label del tiempo


# Creo la lineas y circulos de cada pendulo
line = []
circle = []
colors = ['k','b','g','r','y','w']
for i in range(N):
	x0, y0 = get_coords(theta0, L[i])
	line.append(ax.plot([0, x0], [0, y0], lw=1, c='k'))
	circle.append(ax.add_patch(plt.Circle(get_coords(theta0, L[i]), bob_radius, fc=colors[i%6], zorder=3)))

# Set the plot limits so that the pendulum has room to swing!
ax.set_xlim(-max(L)*1.2, max(L)*1.2)
ax.set_ylim(-max(L)*1.2, max(L)*1.2)

# Esto es para escribir para cada tiempo que fraccion del periodo total

def fraccionMasCarcana (t):
	if 0.14 < t and t < 0.18:
		return '1/6'
	elif 0.18 < t and t < 0.225:
		return '1/5'
	elif 0.225 < t and t < 0.29:
		return '1/4'
	elif 0.29 < t and t < 0.36:
		return '1/3'
	elif 0.36 < t and t < 0.43:
		return '2/5'
	elif 0.43 < t and t < 0.57:
		return '1/2'
	elif 0.57 < t and t < 0.63:
		return '3/5'
	elif 0.63 < t and t < 0.69:
		return '2/3'
	elif 0.69 < t and t < 0.77:
		return '3/4'
	elif 0.77 < t and t < 0.82:
		return '4/5'
	elif 0.82 < t and t < 0.86:
		return '5/6'
	else:
		return '#'


def animate_i(i):
	def animate(j):
		x, y = get_coords(theta[i][j], L[i])
		line[i][0].set_data([0, x], [0, y])
		circle[i].set_center((x, y))
		time_text.set_text('time = ' + str(round(j*T/len(Tiempo), 1)) + '    ' + str(fraccionMasCarcana(j/len(Tiempo))) + ' T')
	return animate

interval = 1 # Delay entre fotogramas (en microsegundos)

# Creo las animaciones de cada pendulo
animaciones = []
for i in range(N):
	decorado = animate_i(i)
	nframe = nsteps
	ani = animation.FuncAnimation(fig, decorado, frames=nframe, repeat=False, interval=interval)
	animaciones.append(ani)
plt.show()
