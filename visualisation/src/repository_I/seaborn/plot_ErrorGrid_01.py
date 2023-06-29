import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


sns.set_theme(style="whitegrid")


# gridRows = ['Error_RMSE_pulse', 'Error_NAQ', 'Error_H1H2', 'Error_HRFdB', 'Error_ST']
# errLabels = ["RMS distance\n(%)", "NAQ error\n(%)", "H1H2 error\n(dB)", "HRF error\n(dB)", "ST error\n(dB/decade)"]
gridRows = ['rmseNorm',
            'naq',
            'h1h2',
            'hrfdB',
            'errorSpectralTilt'
            ]
errLabels = ["RMS distance\n(%)",
             "NAQ error\n(%)",
             "H1H2 error\n(dB)",
             "HRF error\n(dB)",
             "ST error\n(dB/decade)"
             ]

gridCols =  ['GIFmethod', 'vocalEffort', 'f0', 'vowel']
colLabels = ['GIF method', 'phonation type', 'f0', 'vowel']
colOrder = [[], ['creaky', 'normal', 'breathy', 'whispery'], ['low', 'mid', 'high'], ['I', 'E', 'AE', 'A', 'U', 'O']]
genres = ['Female', 'Male']
f0bins = [0, 190, 280, 600]
f0Labels = ['low', 'mid', 'high']

# Paths
# fileIn = os.path.join(os.path.expanduser('~'), 'Projects/FEMVoQ/Corpus/OPENGLOT/repository_I/results', 'Study_GIF_ResultsP_Iberspeech2022.csv')
# pathToFigure = os.path.join(os.path.expanduser('~'), 'Projects/FEMVoQ/Corpus/OPENGLOT/repository_I/results/figures', 'GridErrors_RMSE_NAQ_H1H2_HRF_ST.svg')

# df = pd.read_csv(fileIn,delimiter=';')
# #df = df[df.f0<250]

dataPathBase = 'Projects/LaSalle/FEMVoQ/Corpus/OPENGLOT/specialIssue'
experimentDir = 'OpenGlotIExtended'
dataDir = 'optimizationResults'
dataPath = os.path.join(os.path.expanduser('~'), dataPathBase, experimentDir, dataDir)
dataFileNames = [
    'repository_I_extended_Original-IAIF_grid-search_pulse.csv',
    'repository_I_extended_IOP-IAIF_grid-search_pulse.csv',
    'repository_I_extended_GFM-IAIF_grid-search_pulse.csv',
    'repository_I_extended_QCP_grid-search_pulse.csv',
    'repository_I_extended_ST-QCP_grid-search_pulse.csv'
]

# Read input files into DataFrame
df = pd.DataFrame()
for dataFileName in dataFileNames:
    dataFilePath = os.path.join(dataPath, dataFileName)
    dfTemp = pd.read_csv(dataFilePath)
    df = pd.concat([df, dfTemp])

for index, row in df.iterrows():
	if int(row['f0']) < f0bins[1]:
		df.at[index, 'f0'] = f0Labels[0]
	elif int(row['f0']) < f0bins[2]:
		df.at[index, 'f0'] = f0Labels[1]
	else:
		df.at[index, 'f0'] = f0Labels[2]

# df.loc[(df['f0'] <= f0bins[1]), 'f0r'] = f0Labels[0]
# df.loc[(df['f0'] > f0bins[1]) & (df['f0'] <= f0bins[2]), 'f0r'] = f0Labels[1]
# df.loc[(df['f0'] > f0bins[2]) & (df['f0'] <= f0bins[3]), 'f0r'] = f0Labels[2]


figuresDir = 'figures'
figuresPath = os.path.join(dataPath, figuresDir)
# figureNameBase = 'GridErrors_RMSE_NAQ_H1H2_HRF_ST_noYlim'
# figureNameBase = 'GridErrors_RMSE_NAQ_H1H2_HRF_ST'
figureNameBase = 'GridErrors_RMSE_NAQ_H1H2_HRF_ST_noYlim_standWhis'
# figureNameBase = 'GridErrors_RMSE_NAQ_H1H2_HRF_ST_standWhis'
if not os.path.isdir(figuresPath):
    os.makedirs(figuresPath)


# Convert error NAQ to %
# df['Error_NAQ']=100*df['Error_NAQ']
df['naq'] = 100 * df['naq']
df.loc[(df['GIFmethod'] == 'QCP') & (df['STcompensation'] == 1), 'GIFmethod'] = 'QCP-ST'


nRow = len(gridRows);
nCol = len(gridCols);
axes = []
axId = 0

yAxesLimArray = [[0, 80], [0, 100], [0,30], [0, 15], [0, 40]]

for idx, genre in enumerate(genres):
	fig, axes = plt.subplots(nRow, nCol, gridspec_kw={'width_ratios': [5,8,6,10]},sharey='row')
	df_genre = df[df['gender'] == genres[idx]]
	for rowId in range(nRow):
		for colId in range(nCol):
			if len(colOrder[colId])>0:
				g = sns.boxplot(data=df_genre,
								x=gridCols[colId],
								y=gridRows[rowId],
								hue='GIFmethod',
								order=colOrder[colId],
								# whis=[5, 95],
								showfliers=False,
								ax=axes[rowId][colId]
								)
				# axes[rowId][colId].set_ylim(yAxesLimArray[rowId])
			else:
				g = sns.boxplot(data=df_genre,
								x=gridCols[colId],
								y=gridRows[rowId],
								hue='GIFmethod',
								# whis=[5, 95],
								showfliers=False,
								ax=axes[rowId][colId],
								dodge=False
								)
				# axes[rowId][colId].set_ylim(yAxesLimArray[rowId])
			
			if colId==0:
				axes[rowId][colId].set_ylabel(errLabels[rowId])
			else:
				g.set(ylabel=None)
				axes[rowId][colId].tick_params(labelleft=False)
			if rowId==nRow-1:
				axes[rowId][colId].tick_params(axis='x',labelrotation=45,pad=0)
				axes[rowId][colId].set_xlabel(colLabels[colId])
			else:
				g.set(xlabel=None)
				axes[rowId][colId].tick_params(labelbottom=False)

			axes[rowId][colId].legend_.remove()
			axId = axId+1

	plt.suptitle(figureNameBase + '_' + genre, y=0.03)
	fig.set_size_inches([15, 12])
	fig.align_labels() 
	plt.subplots_adjust(left=0.1, bottom=0.13, right=0.99, top=0.98, hspace=0.12, wspace=0.03)
	#plt.show()
	pathToFigure = os.path.join(figuresPath, figureNameBase + '_' + genre + '.svg')
	plt.savefig(pathToFigure)
