import pandas as pd
from pandas.api.types import CategoricalDtype
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as smm
import itertools
import os

# Load data
#dataPathBase = 'Projects/LaSalle/FEMVoQ/Corpus/OPENGLOT/specialIssue'
dataPathBase = '/data/Dropbox/Salle/2021FEMVoQ/Corpus/OPENGLOT/specialIssue' #MF
experimentDir = 'OpenGlotIExtended'
dataDir = 'optimizationResults'
#dataPath = os.path.join(os.path.expanduser('~'), dataPathBase, experimentDir, dataDir)
dataPath = os.path.join(dataPathBase, experimentDir, dataDir) #MF
resultsDir = 'statisticsResults'
resultsFileName = 'WilcoxonSpecialIssue2023' + experimentDir
#resultsPath = os.path.join(os.path.expanduser('~'), dataPathBase, experimentDir, resultsDir)
resultsPath = os.path.join(dataPathBase, experimentDir, resultsDir) #MF
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

# Delete 80 Hz
df = df[df['f0']!=80]

df.loc[(df['GIFmethod'] == 'QCP') & (df['STcompensation'] == 1), 'GIFmethod'] = 'ST-QCP'
#df.loc[df['GIFmethod'] == 'Original-IAIF','GIFmethod'] = 'IAIF'

# discretize f0 into 3 categories
df["f0r"]=pd.cut(df["f0"],bins=f0bins,labels=f0Labels,right=False)

# Define cfg comparisons
GIFMethods = ['Original-IAIF','IOP-IAIF','GFM-IAIF','QCP','ST-QCP']
GIFMethodsPairs=list(itertools.combinations(GIFMethods,2))

gridRows = ['rmseNorm',
            'naq',
            'h1h2',
            'hrfdB',
            'errorSpectralTilt'
            ]
tableDict={'rmseNorm':'RMSE','naq':'NAQ','h1h2':'H1H2','hrfdB':'HRF','errorSpectralTilt':'Spectral Tilt','Original-IAIF':'IAIF'}
gridCols = ['GIFmethod', 'vocalEffort', 'f0r', 'vowel']
colOrder = [[], ['creaky','normal','breathy','whispery'], ['low','mid','high'], ['I','E','AE','A','U','O']]
groups = ['None','creaky','normal','breathy','whispery','low','mid','high','I','E','AE','A','U','O']
groupAbr = {'None':'','vocalEffort':'c n b w','f0r':'l m h','vowel':'i e æ a u o'}
groupingVars = ['None','vocalEffort','f0r','vowel']
genders = ['Male','Female']

# Define categorical variables to preserve order
cat_groupVar_order = CategoricalDtype(groupingVars,ordered=True)
cat_group_order = CategoricalDtype(groups,ordered=True)
cat_measure_order = CategoricalDtype(gridRows,ordered=True)
cat_cfg_order = CategoricalDtype(GIFMethods,ordered=True)
cat_gender_order = CategoricalDtype(genders,ordered=True)

nRow = len(gridRows);
nCol = len(gridCols);

## WILCOXON
# ncfg=len(cfgOrder1)
for idx, gender in enumerate(genders):
	df_gender = df[df['gender'] == gender]
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
		for cfg1, cfg2 in GIFMethodsPairs:
			dataA = df_gender[df_gender.GIFmethod==cfg1][measure]
			dataB = df_gender[df_gender.GIFmethod==cfg2][measure]
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
				df_gender_ = eval("df_gender[df_gender." + gridCols[colId] + "=='" + group_ + "']")
				for cfg1, cfg2 in GIFMethodsPairs:
					dataA = df_gender_[df_gender_.GIFmethod==cfg1][measure]
					dataB = df_gender_[df_gender_.GIFmethod==cfg2][measure]
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

	cfgOrder1,cfgOrder2=list(zip(*GIFMethodsPairs))
	d={'measure':variable,'groupingVar':groupingVar,'group':group,'cfg1':cfgOrder1*(ngroups*nRow),'cfg2':cfgOrder2*(ngroups*nRow),'p':p, 'v':v, 'pc':pc}

	resultsFilePath = os.path.join(resultsPath, resultsFileName + '_' + gender + '.csv')
	wcx = pd.DataFrame(data=d)
	wcx.to_csv(resultsFilePath)
	
	wcx['gender']=gender
	if idx==0:
		wcxAll=wcx.copy()
	else:
		wcxAll=pd.concat([wcxAll,wcx])

# Apply categorical variables to preserve order
wcxAll['groupingVar']=wcxAll['groupingVar'].astype(cat_groupVar_order)
wcxAll['group']=wcxAll['group'].astype(cat_group_order)
wcxAll['measure']=wcxAll['measure'].astype(cat_measure_order)
wcxAll['cfg1']=wcxAll['cfg1'].astype(cat_cfg_order)
wcxAll['cfg2']=wcxAll['cfg2'].astype(cat_cfg_order)
wcxAll['gender']=wcxAll['gender'].astype(cat_gender_order)
wcxAll=wcxAll.sort_values(['groupingVar','measure','cfg1','gender','cfg2','group'])

wcxAst=wcxAll.copy()
wcxAst.loc[wcxAll['pc']>=0.05,'pc']='n'
wcxAst.loc[(wcxAll['pc']>=0.01) & (wcxAll['pc']<0.05),'pc']='*'
wcxAst.loc[wcxAll['pc']<0.01,'pc']='*'

# for each grouping variable
for groupingVar in groupingVars:
	wcxTable=pd.pivot_table(wcxAst[wcxAst.groupingVar==groupingVar],values='pc',index=['measure','cfg1'],columns=['gender','cfg2'],aggfunc=[lambda x: " ".join(v for v in x)])
	# drop lambda level
	wcxTable.columns=wcxTable.columns.droplevel()
	if groupingVar != 'None':
		wcxTable.loc[('',groupingVar),:]=[groupAbr[groupingVar]]*len(wcxTable.columns)
		idx=wcxTable.index.tolist()
		idtmp=idx.pop(len(wcxTable)-1)
		wcxTable=wcxTable.reindex([idtmp] + idx)
	# renaming before saving the table
	wcxTable=wcxTable.rename(index=tableDict)
	# save the table
	resultsFilePath = os.path.join(resultsPath, resultsFileName + '_' + groupingVar + '.xlsx')
	wcxTable.to_excel(resultsFilePath)
	resultsFilePath = os.path.join(resultsPath, resultsFileName + '_' + groupingVar + '.tex')
	wcxTable.to_latex(resultsFilePath, na_rep='')
