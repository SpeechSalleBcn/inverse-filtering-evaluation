import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


sns.set_theme(style="whitegrid")


# gridCols = ['Error_RMSE_pulse'   , 'Error_NAQ'      ,'Error_H1H2'       ,'Error_HRFdB'    ,'Error_ST']
# colLabels = ["RMS distance (%)" , "NAQ error (%)" ,"H1H2 error (dB)" ,"HRF error (dB)","ST error (dB/decade)"]

gridCols = ['rmseNorm',
            'naq',
            'h1h2',
            'hrfdB',
            'errorSpectralTilt'
            ]
colLabels = ["RMS distance\n(%)",
             "NAQ error\n(%)",
             "H1H2 error\n(dB)",
             "HRF error\n(dB)",
             "ST error\n(dB/decade)"
             ]

gridRows =  ['GIFmethod','GIFmethod']
rowLabels = ['GIF method','GIFmethod']

boxOrder = [[], ['creaky','normal','breathy','whispery'], ['low','mid','high'] ,['I','E','AE','A','U','O']]

# Paths
# fileIn = os.path.join(os.path.expanduser('~'), 'Projects/LaSalle/FEMVoQ/Corpus/OPENGLOT/repository_I/2023_05_11_LJ/optimizationResults', 'Study_GIF_ResultsP_Iberspeech2022.csv')
# pathToFigure = os.path.join(os.path.expanduser('~'), 'Projects/LaSalle/FEMVoQ/Corpus/OPENGLOT/repository_I/results/figures', 'GridErrors_RMSE_NAQ_H1H2_HRF_ST_global.svg')

dataPathBase = 'Projects/LaSalle/FEMVoQ/Corpus/OPENGLOT/repository_I'
experimentDir = '2023_05_11_LJ'
dataDir = 'optimizationResults'
dataPath = os.path.join(os.path.expanduser('~'), dataPathBase, experimentDir, dataDir)
dataFileNames = [
    'repository_I_Original-IAIF_grid-search_pulse.csv',
    'repository_I_IOP-IAIF_grid-search_pulse.csv',
    'repository_I_GFM-IAIF_grid-search_pulse.csv',
    'repository_I_QCP_grid-search_pulse.csv',
    'repository_I_ST-QCP_grid-search_pulse.csv'
]

# Read input files into DataFrame
df = pd.DataFrame()
for dataFileName in dataFileNames:
    dataFilePath = os.path.join(dataPath, dataFileName)
    dfTemp = pd.read_csv(dataFilePath)
    df = pd.concat([df, dfTemp])

figuresDir = 'figures'
figuresPath = os.path.join(dataPath, figuresDir)
figureNameBase = 'GridErrors_RMSE_NAQ_H1H2_HRF_ST'
if not os.path.isdir(figuresPath):
    os.makedirs(figuresPath)

# df = pd.read_csv(fileIn,delimiter=';')
#df = df[df.f0<250]
# Convert error NAQ to %
df['naq'] = 100 * df['naq']
df.loc[(df['GIFmethod'] == 'QCP') & (df['STcompensation'] == 1), 'GIFmethod'] = 'QCP-ST'

nRow = len(gridRows);
nCol = len(gridCols);
#axes = []
axId = 0

fig, axes = plt.subplots(nRow,nCol)

for rowId in range(nRow):
	for colId in range(nCol):
		g = sns.boxplot(data=df, y=gridCols[colId], x=gridRows[rowId], hue='GIFmethod', whis=[5, 95], showfliers=False,ax=axes[rowId][colId], dodge=False)
		axes[rowId][colId].set_ylabel(colLabels[colId])
		axes[rowId][colId].tick_params(labelbottom=False)
		g.set(xlabel=None)
		if rowId==0:
			axes[rowId][colId].legend_.remove()

		axId = axId+1

fig.set_size_inches([12,6])
fig.align_labels() 
plt.subplots_adjust(left=0.05, bottom=0.024, right=0.99, top=0.98,hspace=0.12,wspace=0.35)
#plt.show()
pathToFigure = os.path.join(figuresPath, figureNameBase + '.svg')
plt.savefig(pathToFigure)