import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


sns.set_theme(style="whitegrid")


gridCols = ['Error_RMSE_pulse'   , 'Error_NAQ'      ,'Error_H1H2'       ,'Error_HRFdB'    ,'Error_ST']
colLabels = ["RMS distance (%)" , "NAQ error (%)" ,"H1H2 error (dB)" ,"HRF error (dB)","ST error (dB/decade)"]

gridRows =  ['vocalEffort'   ,'f0r','vowel']
rowLabels = ['phonation type','f0' ,'vowel']

boxOrder = [['creaky','normal','breathy','whispery'], ['low','mid','high'] ,['I','E','AE','A','U','O']]

# Paths
# fileIn = '../Study_GIF_ResultsP_Iberspeech2022.csv'
# currPath = os.path.dirname(__file__)
fileIn = os.path.join(os.path.expanduser('~'), 'Projects/FEMVoQ/Corpus/OPENGLOT/repository_I/results', 'Study_GIF_ResultsP_Iberspeech2022.csv')
pathToFigure = os.path.join(os.path.expanduser('~'), 'Projects/FEMVoQ/Corpus/OPENGLOT/repository_I/results/figures', 'GridErrors_RMSE_NAQ_H1H2_HRF_ST_perGroups.svg')



df = pd.read_csv(fileIn,delimiter=';')
#df = df[df.f0<250]
# Convert error NAQ to %
df['Error_NAQ']=100*df['Error_NAQ']

nRow = len(gridRows);
nCol = len(gridCols);
#axes = []
axId = 0

fig, axes = plt.subplots(nRow,nCol)

for rowId in range(nRow):
	for colId in range(nCol):
		if len(boxOrder[rowId])>0:
			g = sns.boxplot(data=df, y=gridCols[colId], x=gridRows[rowId], hue='GIFmethod', whis=[5, 95], showfliers=False,ax=axes[rowId][colId], order=boxOrder[rowId],linewidth=0.7)
		else:
			g = sns.boxplot(data=df, y=gridCols[colId], x=gridRows[rowId], hue='GIFmethod', whis=[5, 95], showfliers=False,ax=axes[rowId][colId], dodge=False)
	
		axes[rowId][colId].set_ylabel(colLabels[colId])
		axes[rowId][colId].tick_params(rotation=45,pad=-5)
		g.set(xlabel=None)
		if rowId<8:
			axes[rowId][colId].legend_.remove()

		axId = axId+1

fig.set_size_inches([12,9])
fig.align_labels() 
plt.subplots_adjust(left=0.05, bottom=0.03, right=0.99, top=0.98,hspace=0.22,wspace=0.29)
#plt.show()
plt.savefig(pathToFigure)
