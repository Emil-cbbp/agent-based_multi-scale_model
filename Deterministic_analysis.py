import numpy as np
import matplotlib.pyplot as plt
import OdeSolver as ODE
from matplotlib import cm
import matplotlib
from matplotlib.colors import ListedColormap

SMALLER_SIZE = 10
SMALL_SIZE = 14
MEDIUM_SIZE = 15
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALLER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def dydx_rates(x, params):
	R, T, G, P, X, N = x
	''' p-parameters '''
	p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10, p_11, p_12, p_13, p_14, p_15, p_16 = params[:16]
	
	
	''' gamma-parameters (/h) '''
	g_R, g_T, g_G, g_P, g_X = params[16:-1]
	
	''' knockdown '''
	kd = params[-1]

	# x: [R, T, G, P, X, N]]
	dx = np.array([	[-g_R * R	+ kd['Runx1'] * ((p_1 * R + p_2 * N) / (1 + p_1 * R + p_2 * N))],
					[-g_T * T	+ kd['Tcf7'] * ((p_3 * T + p_4 * G + (p_5 * N) / (p_6 + P)) /
								(1 + p_3 * T + p_4 * G+ p_5 * N + p_7 * P))],
					[-g_G * G	+ kd['Gata3'] * ((p_8 * T + (p_9 * N) / (p_10 + P)) /
								(1 + p_8 * T + p_9 * N + p_11 * P))],
					[-g_P * P	+ kd['PU.1'] * ((p_12 * P) / (1 + p_12 * P + p_13 * R * G * T))],
					[-g_X * X	+ 1 / (1 + p_14 * T + p_15 * G)],
					[0			+ p_16 / (1 + N)]])
	return dx

### Knockdown ###
knockdown = {'Runx1':1.0, 'Tcf7': 1.0, 'Gata3':1.0, 'PU.1':1.0}

### Parameters ###
params = np.array([	0.1, 1.0, 5.0, 1.0, 1.5, 0.01, 
					0.5, 0.7, 0.5, 1.0, 0.2, 
					2.5, 2.6, 2.0, 1.0, 0.01, 
					0.15, 0.15, 0.23, 0.06, 0.02, knockdown])

# Initial values ###
x_0_org = [round(1.3), round(1.6), round(1.3), round(4.6), round(8.0),  round(7.0)]	# [R, T, G, P, X, N]

shift = 28.056
N_CpG = 500
meth_frac = 0.75
day_length = 75				# 1 day = 75 [a.u.] in created time series
h_to_au = day_length/24		# 1 hour = 75/24 = 3.125 [a.u] in ts
T_max = (120 - shift) * h_to_au


T = 500
h = 0.1
steps = int(T / h)

ode = ODE.DetOdeSolver(	f=dydx_rates,
						x0=x_0_org,
						a=params,
						h=h,
						steps=steps)
ode.runge_kutta_4()
x_wt = ode.x
x_kd = {}
t = {}

for g, gene in enumerate(knockdown):
	knockdown = {'Runx1':1.0, 'Tcf7': 1.0, 'Gata3':1.0, 'PU.1':1.0}
	kd = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
	x_kd[gene] = []

	for k in kd:
		knockdown[gene] = k
		x_0 = x_0_org.copy()
		x_0[g] = x_0_org[g] * k
		params = np.array([	0.1, 1.0, 5.0, 1.0, 1.5, 0.01, 
							0.5, 0.7, 0.5, 1.0, 0.2, 
							2.5, 2.6, 2.0, 1.0, 0.01, 
							0.15, 0.15, 0.23, 0.06, 0.02, knockdown])
		ode = ODE.DetOdeSolver(	f=dydx_rates,
								x0=x_0,
								a=params,
								h=h,
								steps=steps)
		ode.runge_kutta_4()
		x_kd[gene].append(ode.x)



ranges = {'Runx1':(0,10), 'Tcf7':(0,50), 'Gata3':(0,5), 'PU.1':(0,20), 'X':(0,50)}
titles = ['Runx1', 'Tcf7', 'Gata3', 'PU.1', 'X']

t = np.linspace(0,T, steps)
time_days = np.array([24, 48, 72, 96, 120]) * h_to_au
time_labels = ['Day1', 'Day 2', 'Day 3', 'Day 4', 'Day 5']



nrows = len(titles)
ncols = len(knockdown)


colours = {	'Runx1': 'blue', 'Tcf7': 'red', 'Gata3': 'yellow', 'PU.1': 'purple', 
			'X': 'green', 'Notch':'grey', 'Bcl11b': 'orange'}

# create yellow colormaps
N_yel = 256
yellow = np.ones((N_yel, 4))
yellow[:, 0] = np.linspace(255/256, 1, N_yel)[::-1] # R = 255
yellow[:, 1] = np.linspace(232/256, 1, N_yel)[::-1] # G = 232
yellow[:, 2] = np.linspace(11/256, 1, N_yel)[::-1]  # B = 11
Yellows = ListedColormap(yellow)




cmaps = {	'Runx1':plt.cm.get_cmap('Blues'), 'Tcf7':plt.cm.get_cmap('Reds'), 
			'Gata3':Yellows, 'PU.1':plt.cm.get_cmap('Purples'),
			'X':plt.cm.get_cmap('Greens')}
norm = matplotlib.colors.Normalize(vmin=-2, vmax=len(kd))

lines = [	'dotted', 'dashed', 'dashdot', 'solid', 
			(0, (1, 1)),
			(0, (5, 1)), 
			(0, (3, 1, 1, 1, 1, 1))]


fig, ax = plt.subplots(figsize=(10, 14), nrows=nrows, ncols=ncols, sharex=True, sharey='row')
fig.set_tight_layout(True)
for j, gene_exp in enumerate(titles):
	for i, gene_kd in enumerate(knockdown):
		for l, k in enumerate(kd):
			ax[j][i].plot(t, x_kd[gene_kd][l][:,j], label=str(k), color=cmaps[gene_exp](norm(l)), linestyle=lines[l])
		ax[j][i].plot(t, x_wt[:,j], label='WT', color='grey')
		ax[j][i].set_ylim(ranges[gene_exp])
		ax[j][i].tick_params(axis='both', direction='in', top=True, right=True)
		if(j == 0):
			ax[j][i].set_title('{:s} KD'.format(gene_kd))
		if(i == 0):
			ax[j][i].set_ylabel('{:s} level [counts]'.format(gene_exp))
		if(i == ncols - 1):
			ax[j][i].legend(loc='center right')
		if(j == nrows - 1):	
			ax[j][i].set_xlabel('Time [days]')
			ax[j][i].set_xticks(time_days)
			ax[j][i].set_xticklabels(time_labels, rotation=45)
	




plt.show()



























#plt.show()

