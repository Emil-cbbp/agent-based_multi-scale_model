import numpy as np
import matplotlib.pyplot as plt
import Colony as col
import Cell
import Gillespie_methyl_collab as GILcoll
import timer
import os


shift = 28.056
def division_scheme(gen):
	# cell cycle length for generations 0 to 4, given as mu and sigma from normal
	# distribution
	cellcycle = np.array([	[34.-shift, 13.],
							[15., 5.],
							[13., 5.],
							[12., 4.],
							[12., 3.]])
	if(gen > 4):
		gen = 4
	div = -1
	while(div < 0):
		div = np.random.normal(cellcycle[gen][0], cellcycle[gen][1])
	return div


''' Transcriptional Level '''
### Transcription rates ###
def dydx_rates(x, params):
	R, T, G, P, X, N = x
	''' p-parameters '''
	p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10, p_11, p_12, p_13, p_14, p_15, p_16 = params[:16]
	
	
	''' gamma-parameters (/h) '''
	g_R, g_T, g_G, g_P, g_X = params[16:-1]
	
	''' knockdown '''
	kd = params[-1]

	dx = np.array([	[g_R * R, 	kd['Runx1'] * ((p_1 * R + p_2 * N) / (1 + p_1 * R + p_2 * N))],
					[g_T * T, 	kd['Tcf7'] * ((p_3 * T + p_4 * G + (p_5 * N) / (p_6 + P)) /
								(1 + p_3 * T + p_4 * G+ p_5 * N + p_7 * P))],
					[g_G * G, 	kd['Gata3'] * ((p_8 * T + (p_9 * N) / (p_10 + P)) /
								(1 + p_8 * T + p_9 * N + p_11 * P))],
					[g_P * P, 	kd['PU.1'] * ((p_12 * P) / (1 + p_12 * P + p_13 * R * G * T))],
					[g_X * X, 	1 / (1 + p_14 * T + p_15 * G)],
					[0, 		p_16 / (1 + N)]])
	return dx

### Knockdown ###
knockdown = {'Runx1':1.0, 'Tcf7': 1.0, 'Gata3':1.0, 'PU.1':1.0, 'time_on':float('inf')}		# No knockdown

### Parameters ###
params = np.array([	0.1, 1.0, 5.0, 1.0, 1.5, 0.01, 
					0.5, 0.7, 0.5, 1.0, 0.2, 
					2.5, 2.6, 2.0, 1.0, 0.01, 
					0.15, 0.15, 0.23, 0.06, 0.02, knockdown])

# Initial values ###
x_0 = [round(1.3), round(1.6), round(1.3), round(4.6), round(8.0),  round(7.0)]	# [R, T, G, P, X, N]



''' Epigenetic Level '''
def collaborative(y, x, t, params):
	k_1, k_2, k_3, alpha, beta, gamma, delta, epsilon = params
	R, X, N = x

	rates =	[	k_1 * X,			# u -> h
				k_1 * X,			# h -> m
				k_2 * N + k_3 * R,	# m -> h
				k_2 * N + k_3 * R,	# h -> u
				alpha,				# u -> h, if med=m
				beta,				# h -> m, if med=m
				gamma,				# h -> u, if med=u
				delta,				# m -> h, if med=u
				epsilon]			# h -> m, if med=h
	return rates


### Parameters ###
params_m = np.array([0.28, 0.20, 0.20, 0.002, 0.002, 0.0005, 0.0005, 0.002])
N_CpG = 500				# Total number of promotor CpG sites
mediator_range = N_CpG	# Range of mediator. If same as CpG sites, spatial effects eliminated

### Initial values
y_0_m = [N_CpG, 0, 0]		# Completely methylated



day_length = 75				# 1 day = 75 [a.u.] in created time series
h_to_au = day_length/24		# 1 hour = 75/24 = 3.125 [a.u] in ts

measure_points = np.array([59, 134, 209, 284]) #* h_to_au	# already converted
stuff_to_measure = {'R': 0, 'T': 1, 'G': 2, 'P': 3, 'X': 4, 'N': 5}

T_max = (120 - shift) * h_to_au

names = ['Runx1', 'Tcf7', 'Gata3', 'PU.1', 'X', 'Notch']
gene_dic = {g: i for i, g in enumerate(names)}




tot_time_start = timer.start_time()

wt_simulations = ['wt']		# wt simulations
no_kd = [-1]	#no knockdown initialisation time

kd_simulations = ['Runx1', 'Tcf7' , 'Gata3', 'PU.1']	# genes to knockdown
kd_init = [50, 100, 150, 200, 250, 0]	# knockdown times

# pick simualtion type and kd times
simulations = wt_simulations
kd_times = no_kd
#simulations = kd_simulations
#kd_times = kd_init



N_colonies = 3		# Number of colonies to simulate per simulation type
for KD in simulations:
	for t_kd in kd_times:
		print('\n\n\nPerforming knockdown experiments with {:s} at time {:d}.'.format(KD, t_kd))

		PATH = os.getcwd() + '/Results/'
		list_subfolders = [f.name for f in os.scandir(PATH) if f.is_dir()]
		
		
		knockdown = {k:1.0 for k in knockdown}
		
		if(KD == 'wt'):
			kd_frac = 1.0
			name = 'wt/'
			t_kd = float('inf')
		else:
			kd_frac = 0.2
			knockdown[KD] = kd_frac
			name = 'KD_{:s}_{:s}_timeOn_{:d}/'.format(KD, str(kd_frac).replace('.', ''), t_kd)
		knockdown['time_on'] = t_kd
		
		PATH += name
		if(name[:-1] not in list_subfolders):	#[:-1] to exclude the slash
			os.mkdir(PATH)
		
		
		for n in range(N_colonies):
			print('\n\nIteration {:d}\n'.format(n))
			colony = col.Colony(x_0,division_scheme, dydx_rates, params, gene_dic, PATH=PATH, methyl='advanced', day_length=day_length, 
							rates_epi=collaborative, params_epi=params_m,
							y_0=y_0_m, N_CpG=N_CpG, mediator_range=mediator_range, 
							gillespie=GILcoll.Gillespie, knockdown=knockdown)
					
			time_s = timer.start_time()
			colony.run_max_time(T_max)
			timer.time_elapsed(time_s)

		#### Measure
			print('Measurements')
			colony.measure(	measure_points=measure_points, measure_levels=stuff_to_measure, 
						measure_methylation=True, meth_frac=0.75)

			
			save = True

			if(save):
				colony.save_colony(sparse=False)

			t, ts = colony.create_tree(binary=True, save=save)
			t, ts = colony.create_tree(binary=True, save=save, mark_X=True)
			t, ts = colony.create_tree(binary=False, save=save)

			timer.time_elapsed(tot_time_start, current_time=True)


timer.time_elapsed(tot_time_start)










