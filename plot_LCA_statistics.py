import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import os
import timer

SMALL_SIZE = 25
MEDIUM_SIZE = 20
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title




def plot_expression(data, labs, nobs, measure=None, names=None, save='', title='', fmt='png'):
	if(measure == None):
#		measure = ['Runx1_final', 'Tcf7_final', 'Gata3_final', 'Pu1_final', 'X_final', 'Bcl11b_frac', 'M_final', 'H_final', 'U_final']
		measure = ['Runx1_final', 'Tcf7_final', 'Gata3_final', 'PU.1_final', 'X_final', 'Bcl11b_frac', 'C_final', 'I_final', 'O_final']
	if(names == None):
		names = ['Runx1', 'Tcf7', 'Gata3', 'PU.1', 'X', 'Bcl11b', 'Closed regulatory sites', 'Intermediate regulatory sites', 'Open regulatory sites']

	colours = {'Runx1': 'blue', 'Tcf7': 'red', 'Gata3': 'yellow', 'PU.1': 'purple', 'X': 'green', 'Notch':'grey', 'Bcl11b': 'orange', 'Closed regulatory sites':'pink', 'Intermediate regulatory sites':'brown', 'Open regulatory sites':'cyan'}


	n_cols = int(np.floor(np.sqrt(len(names))))
	n_rows = int(np.ceil(np.sqrt(len(names))))

	size_x = 32 if n_cols >= 3 else 15.5
	size_y = 16 if n_rows >= 3 else 11.5
	size = (size_x, size_y)



	fig, ax1 = plt.subplots(figsize=size, nrows=n_rows, ncols=n_cols, sharex=True)
	ax1 = ax1.flatten()
	fig.set_tight_layout(True)
	for i, g in enumerate(names):
		ax = ax1[i]
		sns.barplot(x='Labels', y=measure[i], data=data, order=labs, ax=ax, color=colours[g])	#, errorbar='sd'

		pos = range(len(nobs)) 
		
		ax.set_xlabel('')
		if(i % n_rows == 0):
			ax.set_ylabel('Expression level [counts]')
		else:
			ax.set_ylabel('')
		ax.set_title(g)
		if('M_' in measure[i] or 'H_' in measure[i] or 'U_' in measure[i]):
			ax.set_ylim(0,510)
		elif('Bcl11b_frac' == measure[i]):
			ax.set_ylim(0, 1.05)
		
		for tick in pos:
			ax.text(pos[tick], 0.05 * ax.get_ylim()[1], nobs[tick], horizontalalignment='center', verticalalignment='bottom', size='medium', color='gray', weight='semibold', rotation=90)

		if(i // n_rows == n_rows - 1):
			ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
		ax.tick_params(axis='both', direction='in', top=True, right=True)
	
		if(title):
			plt.suptitle(title)
		if(save):
			plt.savefig(save + '.' + fmt, format=fmt)



N_CpG = 500
meth_frac = 0.75

n_rows = 2
n_cols = 2
		
kd_exp = ['wt', 'Runx1', 'Tcf7', 'Gata3', 'PU.1']
for KD in kd_exp:	#
	print('\n\n\nAnalysing knockdown experiments with {:s}.'.format(KD))
	PATH = 'Results/'
	save_extention = ''
	title = ''
	if(KD == 'wt'):
		kd_frac = 1.0	#k_f
		save_extention += '_wt'
		title += 'Wild Type'
		variants = ['wt']
	else:
		kd_frac = 0.2
		save_extention += '_KD_{:s}_{:s}'.format(KD, str(kd_frac).replace('.', ''))
		save_copy = save_extention
		title += '{:s} knockdown {:.1f}'.format(KD, kd_frac)
		variants = [f.name for f in os.scandir(PATH) if (f.is_dir() and (KD in f.name))]

	for var in variants:
		PATH = 'Results/'
		PATH += var + '/'
		folder = [f.name for f in os.scandir(PATH) if (f.is_dir() and f.name.startswith('gen'))]
		
		if('timeOn' in var):
			t_kd = int(var.split('_')[-1])
			title += 'at time {:d}'.format(t_kd)
			save_extention = save_copy + '_{:d}'.format(t_kd)
			print('Knockdown at time {:d}.'.format(t_kd))

		N_colonies = len(folder)
		title += ' -- {:d} Colonies'.format(N_colonies)

		statistics = pd.DataFrame({'Generation': [], 'State': [], 'Open': [], 'Closed': [], 'Same': [], 'Opposite': []})

		tot_time_start = timer.start_time()
		for k, fol in enumerate(folder):
			print('Colony {:d}'.format(k))
			print(fol)
			colony = pd.read_csv(PATH + fol + '/colony_labels.tsv', delimiter='\t')
			
			kd_type = [KD] * len(colony)
			colony['Experiment'] = kd_type
			fracs = np.array(colony['Bcl11b_frac'])
			states = ['open' if i >= 0.75 else 'closed' for i in fracs]
			colony['State'] = states
			
			if(k == 0):		#KD == 'wt'):	# exp for LCA_plot, KD for KD_plot
				total = colony.copy()
			else:
				total = pd.concat([total, colony], ignore_index=True)

		labels_org = np.unique(total['Labels']).astype(str)
		labels_org = labels_org[np.where(np.char.find(labels_org, 'pre') == -1)[0]]
		labels_org = labels_org[np.where(np.char.find(labels_org, 'open') == -1)[0]]
		labels_org = labels_org[np.where(np.char.find(labels_org, 'offspring') == -1)[0]]
		post_lca = sorted(labels_org[np.where(np.char.find(labels_org, 'post') > -1)[0]], key = lambda s: int(s.split(' ')[-1].split()[0]))
		lca = reversed(sorted(np.setdiff1d(	labels_org[np.where(np.char.find(labels_org, 'LCA') > -1)[0]], 
						labels_org[np.where(np.char.find(labels_org, 'post') > -1)[0]]), key = lambda s: int(s.split(' ')[-1].split()[0])))
		labels = ['closed pre LCA'] + list(lca) + list(post_lca) + ['open', 'closed offspring']

		obs = total['Labels'].value_counts(sort=False) 
		nobs = obs.values 
		lobs = list(obs.keys()) 
		nobs = [str(nobs[lobs.index(x)]) if x in lobs else '0' for x in labels]
		#nobs = ["n: " + i for i in nobs] 


	LCA_lab = [s for s in labels if s.startswith('LCA')]
	LCA_lab.insert(0, 'closed pre LCA')
	LCA_lab.append('open')
	LCA = total[total['Labels'].isin(LCA_lab)]
	
	post_lab = [s for s in labels if s.startswith('closed post')]
	post_lab.insert(0, 'closed pre LCA')
	post_lab.append('open')
	post = total[total['Labels'].isin(post_lab)]
	
	obs = total['Labels'].value_counts(sort=False) 
	nobs = obs.values 
	lobs = list(obs.keys()) 
	nobs_lca = [str(nobs[lobs.index(x)]) if x in lobs else '0' for x in LCA_lab]
	nobs_post = [str(nobs[lobs.index(x)]) if x in lobs else '0' for x in post_lab]
	
	measure = ['Tcf7_final', 'PU.1_final', 'X_final', 'Bcl11b_frac']
	names = ['Tcf7', 'PU.1', 'X', 'Bcl11b']
	
	if('Plots' not in [f.name for f in os.scandir(os.getcwd()) if f.is_dir()]):
		os.mkdir('Plots')
	if('LCA_statistics' not in [f.name for f in os.scandir(os.getcwd() + '/Plots/') if f.is_dir()]):
		os.mkdir('Plots/LCA_statistics')
	save_path = os.getcwd() + '/Plots/LCA_statistics/' + 'Statistics' + save_extention

	plot_expression(LCA, LCA_lab, nobs_lca, title=title, fmt='png', save=save_path+'-LCA-all')
	plot_expression(post, post_lab, nobs_post, title=title, fmt='png', save=save_path+'-closed_post_LCA-all')
	
	plot_expression(LCA, LCA_lab, nobs_lca, measure=measure, names=names, title=title, fmt='png', save=save_path+'-LCA-'+'_'.join(names))
	plot_expression(post, post_lab, nobs_post, measure=measure, names=names, title=title, fmt='png', save=save_path+'-closed_post_LCA-'+'_'.join(names))
	










