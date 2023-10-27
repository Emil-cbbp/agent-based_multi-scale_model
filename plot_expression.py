import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import os
import timer

SMALLER_SIZE = 10
SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALLER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title






N_CpG = 500
meth_frac = 0.75


shift = 28.056
N_CpG = 500
meth_frac = 0.75
day_length = 75				# 1 day = 75 [a.u.] in created time series
h_to_au = day_length/24		# 1 hour = 75/24 = 3.125 [a.u] in ts
T_max = (120 - shift) * h_to_au

time_days = np.array([24, 48, 72, 96, 120]) * h_to_au
time_labels = ['Day1', 'Day 2', 'Day 3', 'Day 4', 'Day 5']



	
	

""" Plot selected lineages of selected colony """
plot_order = [1,2,3,4,5,6]
genes = ['Runx1', 'Tcf7', 'Gata3', 'PU.1', 'X', 'Notch']
gene_map = {genes[x]:x  for x in range(len(genes))}
colours = {'Runx1': 'blue', 'Tcf7': 'red', 'Gata3': 'yellow', 'PU.1': 'purple', 'X': 'green', 'Notch':'grey', 'Bcl11b': 'orange'}
ls = ['-', ':', '--', '-.'] * 5



KD = 'wt'

PATH_cwd = os.getcwd()
PATH = '/export/scratch/emila/T-Cells/Results/Initial_pseudo/'
PATH = PATH_cwd + '/Results/'
save_extention = ''
title_appendum = ''

if(KD == 'wt'):
	kd_frac = 1.0	#k_f
	PATH += 'wt/'
	save_extention += 'wt'
	title_appendum += 'Wild Type'
else:
	kd_frac = 0.2		#####
	PATH += ''	# chose a colony
	save_extention += 'KD_{:s}_{:s}'.format(KD, str(kd_frac).replace('.', ''))
	title_appendum += '{:s} knockdown {:.1f}'.format(KD, kd_frac)



fol = PATH + 'gen_8_run_1'


if('Plots' not in [f.name for f in os.scandir(os.getcwd()) if f.is_dir()]):
	os.mkdir('Plots')
if('Expression_plots' not in [f.name for f in os.scandir(os.getcwd() + '/Plots/') if f.is_dir()]):
	os.mkdir('Plots/Expression_plots')


tot_time_start = timer.start_time()




timer.time_elapsed(tot_time_start)
colony = pd.read_csv(fol + '/colony_labels.tsv', delimiter='\t')
colony.rename(columns={'Bcl11b_frac': 'Bcl11b_final'}, inplace=True)

final_cells = [187, 205, 241]		# Pick some interesting lineages



cells_id = final_cells.copy()
cells = colony.loc[cells_id]
i = 0
while(True):
	c = cells.loc[cells_id[i]]
	if(c.mother == -1):
		break
	if(c.mother not in cells_id):
		cells_id.append(c.mother)
		cells = cells.append(colony.loc[c.mother])
	i += 1


lineages = {c: [c] for c in final_cells}
for c in lineages:
	i = 0
	while(True):
		cell = cells.loc[lineages[c][i]]
		if(cell.mother == -1):
			break
		else:
			lineages[c].append(cell.mother)
		i += 1
	lineages[c].reverse()


gens = cells['gen']
t_div = cells['t_div']

trails = [[] for _ in range(6)]
cell_development = {}

bcl11b = []

for c in cells.index:
	cell = cells.loc[c]
	ef = cell['expression_file']
	expr = np.loadtxt(fol + '/' + ef + '.tsv', delimiter='\t', skiprows=1)
	t = expr[:,0]
	for g in genes:
		trails[gene_map[g]].append([t,expr[:,plot_order[gene_map[g]]]])
	cell_development[c] = {g: [t,expr[:,plot_order[gene_map[g]]]] for g in genes}
	meth = np.loadtxt((fol + '/' + ef + '.tsv').replace('Expression', 'Methylation'), delimiter='\t', skiprows=1)
	bcl11b.append([meth[:,0],meth[:,-1] / N_CpG])
	cell_development[c]['Bcl11b'] = [meth[:,0],meth[:,-1] / N_CpG]
	


#'''
### One plot containing all chosen lineages
fig, ax1 = plt.subplots(figsize=(32,16))
for i,g in enumerate(plot_order):
	for c in range(len(trails[i])):
		if(c == 0):
			ax1.plot(trails[i][c][0], trails[i][c][1], color=colours[genes[g-1]], label=genes[g-1])
		else:
			ax1.plot(trails[i][c][0], trails[i][c][1], color=colours[genes[g-1]])
ax2 = ax1.twinx()
for c in range(len(bcl11b)):
	if(c == 0):
		ax2.plot(bcl11b[c][0], bcl11b[c][1], color='orange', label='Bcl11b')
	else:
		ax2.plot(bcl11b[c][0], bcl11b[c][1], color='orange')

for i,g in enumerate(gens):
	ax1.vlines(t_div[cells_id[i]], 0, 80, colors='black', linestyles=ls[g])

ax1.set_ylim(0,80)
ax2.set_ylim(0,1.18)
ax1.set_ylabel('Expression level [counts]')
ax2.set_ylabel('Fraction of open Bcl11b regulatory sites')
ax1.set_xlabel('Time [a.u.]')
ax1.grid()
fig.legend(loc='upper left', ncol=2)
fig.tight_layout()


plt.savefig('Plots/Expression_plots/All_chosen_lineages-{:s}-{:s}.png'.format(save_extention, '_'.join([str(s) for s in final_cells])), format='png')
#'''




#'''
### Splitting lineages in subplots
fig, ax1 = plt.subplots(figsize=(14.8,5), nrows=1, ncols=3, sharey='row')
ax1 = ax1.flatten()
fig.set_tight_layout(True)

max_expr = 0
for j, c in enumerate(lineages):
	ax = ax1[j]
	for g in genes:
		expr = []
		time = []
		for lin in lineages[c]:
			time += list(cell_development[lin][g][0])
			expr += list(cell_development[lin][g][1])
		ax.plot(time, expr, color=colours[g], label=g)
		max_expr = max(max_expr, max(expr))
	ax2 = ax.twinx()
	g = 'Bcl11b'
	expr = []
	time = []
	for lin in lineages[c]:
		time += list(cell_development[lin][g][0])
		expr += list(cell_development[lin][g][1])
	ax2.plot(time, expr, color=colours[g], label=g)
	max_expr = max(max_expr, max(expr))


	ax2.set_ylim(0,1.05)
	if(j == 0):
		ax.set_ylabel('Expression level [counts]')
	if(j == len(lineages) - 1):
		ax2.set_ylabel('Fraction of open\nBcl11b regulatory sites')
	else:
		ax2.set_yticklabels([])
	ax.set_xlabel('Time [days]')
	ax.tick_params(axis='both', direction='in', top=True)
	ax2.tick_params(axis='y', direction='in')
	
	ax.set_xticks(time_days)
	ax.set_xticklabels(time_labels, rotation=45)
	ax.set_title('Final cell {:d}'.format(c))
	
	lines, labels = ax.get_legend_handles_labels()
	lines2, labels2 = ax2.get_legend_handles_labels()
	ax.legend(lines + lines2, labels + labels2, loc='center left', ncol=2)

for j, c in enumerate(lineages):
	ax = ax1[j]

	gens = cells[cells['cell_ID'].isin(lineages[c])]['gen']
	t_div = cells[cells['cell_ID'].isin(lineages[c])]['t_div']
	for i,gen in enumerate(gens):
		ax.vlines(t_div.loc[t_div.index[i]], 0, max_expr * 1.1, colors='black', linestyles=ls[gen])
		lab = cells.loc[t_div.index[i]]['Labels']
		if('closed post' in lab):
			lab = lab.replace('closed post', 'closed post \n')
		elif('closed pre' in lab):
			lab = lab.replace('closed pre', 'closed \npre')
		elif('closed off' in lab):
			lab = lab.replace('closed off', 'closed \noff')
		ax.text(t_div.loc[t_div.index[i]], max_expr * 1.1, lab + ' ',
					color='black', fontsize='xx-small', rotation='vertical',
					horizontalalignment='right', verticalalignment='top')

	ax.set_ylim(0,max_expr * 1.1)


plt.savefig('Plots/Expression_plots/Panels_chosen_lineages-{:s}-{:s}.png'.format(save_extention, '_'.join([str(s) for s in final_cells])), format='png')


#'''


#'''
### Subplot for each gene
fig, ax1 = plt.subplots(figsize=(14, 5), nrows=1, ncols=4, sharey=False)#True
ax1 = ax1.flatten()
fig.set_tight_layout(True)

skip = 0
skip_genes = ['Notch', 'Runx1', 'Gata3']
for j, g in enumerate(colours):
	if(g in skip_genes):
		skip += 1
		continue
	max_expr = 0
	ax = ax1[j - skip]
	for c in lineages:
		expr = []
		time = []
		for lin in lineages[c]:
			time += list(cell_development[lin][g][0])
			expr += list(cell_development[lin][g][1])
		ax.plot(time, expr, label=c)
		max_expr = max(max_expr, max(expr))

	gens = cells[cells['cell_ID'].isin(lineages[c])]['gen']
	t_div = cells[cells['cell_ID'].isin(lineages[c])]['t_div']

	ax.set_ylim(0,max_expr * 1.1)
	if(j == skip):
		ax.set_ylabel('Expression level [counts]')
	ax.set_xlabel('Time [days]')
	ax.tick_params(axis='both', direction='in', top=True, right=True)
	ax.set_xticks(time_days)
	ax.set_xticklabels(time_labels, rotation=45)
	ax.legend(loc='upper left', ncol=1)
	ax.set_title(g)

skip = 0
for j, g in enumerate(colours):	
	if(g in skip_genes):
		skip += 1
		continue
	max_expr = 0
	ax = ax1[j - skip]
	for k, c in enumerate(lineages):
		for l in lineages[c]:
			ax.scatter([cells.loc[l]['t_div']], [cells.loc[l][g+'_final']], 
						color=ax.lines[k].get_color(), marker='o', s=100, 
						edgecolors='black', zorder=2)

plt.savefig('Plots/Expression_plots/Panels_genes_chosen_lineages-{:s}-{:s}.png'.format(save_extention, '_'.join([str(s) for s in final_cells])), format='png')


#'''








plt.show()









	
	
	
	
	
