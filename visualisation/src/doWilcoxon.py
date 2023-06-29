import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as smm
from itertools import repeat
import os

# Load data
dataPathBase = 'Projects/LaSalle/FEMVoQ/Corpus/OPENGLOT/specialIssue'
experimentDir = 'OpenGlotIExtended'
dataDir = 'optimizationResults'
dataPath = os.path.join(os.path.expanduser('~'), dataPathBase, experimentDir, dataDir)
resultsDir = 'statisticsResults'
resultsFileName = 'WilcoxonSpecialIssue2023' + experimentDir
resultsPath = os.path.join(os.path.expanduser('~'), dataPathBase, experimentDir, resultsDir)
if not os.path.isdir(resultsPath):
	os.makedirs(resultsPath)
dataFileNames = [
    'repository_I_extended_Original-IAIF_grid-search_pulse.csv',
    'repository_I_extended_IOP-IAIF_grid-search_pulse.csv',
    'repository_I_extended_GFM-IAIF_grid-search_pulse.csv',
    'repository_I_extended_QCP_grid-search_pulse.csv',
    'repository_I_extended_ST-QCP_grid-search_pulse.csv'
]

f0bins = [0, 190, 280, 600]
f0Labels = ['low', 'mid', 'high']

df = pd.DataFrame()
for dataFileName in dataFileNames:
	dataFilePath = os.path.join(dataPath, dataFileName)
	dfTemp = pd.read_csv(dataFilePath)
	df = pd.concat([df, dfTemp])

df.loc[(df['GIFmethod'] == 'QCP') & (df['STcompensation'] == 1), 'GIFmethod'] = 'ST-QCP'

df.loc[(df['f0'] <= f0bins[1]), 'f0r'] = f0Labels[0]
df.loc[(df['f0'] > f0bins[1]) & (df['f0'] <= f0bins[2]), 'f0r'] = f0Labels[1]
df.loc[(df['f0'] > f0bins[2]) & (df['f0'] <= f0bins[3]), 'f0r'] = f0Labels[2]

# Define cfg order
cfgOrder1=["Original-IAIF", "Original-IAIF", "Original-IAIF", "Original-IAIF", "IOP-IAIF", "IOP-IAIF", "IOP-IAIF", "GFM-IAIF", "GFM-IAIF", "QCP"]
cfgOrder2=["IOP-IAIF", "GFM-IAIF", "QCP", "ST-QCP", "GFM-IAIF", "QCP", "ST-QCP", "QCP", "ST-QCP", "ST-QCP"]
tuples=list(zip(*[cfgOrder1,cfgOrder2]))

gridRows = ['rmseNorm',
            'naq',
            'h1h2',
            'hrfdB',
            'errorSpectralTilt'
            ]

gridCols = ['GIFmethod', 'vocalEffort', 'f0r', 'vowel']
colOrder = [[], ['creaky','normal','breathy','whispery'], ['low','mid','high'] ,['I','E','AE','A','U','O']]
genders = ['Female', 'Male']

nRow = len(gridRows);
nCol = len(gridCols);

## WILCOXON
ncfg=len(cfgOrder1)

v=[]
p=[]
pc=[]
groupingVar=[]
group=[]
variable=[]

# global comparison
ngroups=1

#for each measure
for rowId in range(nRow):
	measure = gridRows[rowId]
	v_=[]
	p_=[]
	# perform the comparisons between GIF methods
	for cfg1, cfg2 in zip(cfgOrder1,cfgOrder2):
		dataA = df[df.GIFmethod==cfg1][measure]
		dataB = df[df.GIFmethod==cfg2][measure]
		dataA = dataA.replace([np.inf, -np.inf], np.nan)
		dataB = dataB.replace([np.inf, -np.inf], np.nan)
		dataA_filt = dataA[dataA.notna() & dataB.notna()]
		dataB_filt = dataB[dataA.notna() & dataB.notna()]
		v__,p__=stats.wilcoxon(dataA_filt, dataB_filt)
		v_.append(v__)
		p_.append(p__)
	# Apply p correction Holm-Bonferroni
	rej,pc_=smm.multipletests(p_,method='holm')[:2]
	v.extend(v_)
	p.extend(p_)
	pc.extend(pc_)
	groupingVar.extend( ["None"] *(len(p_)) )
	group.extend( ["None"] *(len(p_)) )
	variable.extend( [measure] *(len(p_)) )

# for each measure
for rowId in range(nRow):
	measure = gridRows[rowId]
	# for each grouping variable
	for colId in range(1,nCol):
		groups = colOrder[colId]
		if rowId==0:
			ngroups = ngroups + len(groups)
		# for each group
		for group_ in groups:
			v_=[]
			p_=[]
			df_ = eval("df[df." + gridCols[colId] + "=='" + group_ + "']")
			for cfg1, cfg2 in zip(cfgOrder1,cfgOrder2):
				dataA = df[df.GIFmethod==cfg1][measure]
				dataB = df[df.GIFmethod==cfg2][measure]
				dataA = dataA.replace([np.inf, -np.inf], np.nan)
				dataB = dataB.replace([np.inf, -np.inf], np.nan)
				dataA_filt = dataA[dataA.notna() & dataB.notna()]
				dataB_filt = dataB[dataA.notna() & dataB.notna()]
				v__,p__=stats.wilcoxon(dataA_filt, dataB_filt)
				v_.append(v__)
				p_.append(p__)
			# Apply p correction Holm-Bonferroni
			rej,pc_=smm.multipletests(p_,method='holm')[:2]
			v.extend(v_)
			p.extend(p_)
			pc.extend(pc_)
			groupingVar.extend( [gridCols[colId]] *(len(p_)) )
			group.extend( [group_] *(len(p_)) )
			variable.extend( [measure] *(len(p_)) )
d={'measure':variable,'groupingVar':groupingVar,'group':group,'cfg1':cfgOrder1*(ngroups*nRow),'cfg2':cfgOrder2*(ngroups*nRow),'p':p, 'v':v, 'pc':pc}

resultsFilePath = os.path.join(resultsPath, resultsFileName + '.csv')
wcx = pd.DataFrame(data=d)
wcx.to_csv(resultsFilePath)
