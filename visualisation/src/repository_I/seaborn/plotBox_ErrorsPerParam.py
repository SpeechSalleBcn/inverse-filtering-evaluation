import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import os

sns.set_theme(style="whitegrid")
paletteColors = ["#E69F00", "#0072B2", '#C0392B', "#F0E442", "#009E73", "#D55E00"]
palette = sns.color_palette(paletteColors)
sns.set_palette(palette)

filterF0Gender = False
standWhis = False
fixY = False
errorSign = False
hrfdB = True

titlePlot = False


hrf = 'hrfdB' if hrfdB else 'hrf'
hrdSign = 'hrfdBSign' if hrfdB else 'hrfSign'

gridRowsArray = [['rmseNorm', 'naq', 'h1h2', hrf, 'errorSpectralTilt'],
			['rmseNorm', 'naqSign', 'h1h2Sign', hrdSign, 'errorSpectralTiltSign']]

errLabels = ["RMS dist.\n(%)",
             "NAQ error\n(%)",
             "H1H2 error\n(dB)",
             "HRF error\n(dB)",
             "ST error\n(dB/decade)"
             ]

plots =  ['GIFmethod', 'vocalEffort', 'f0r', 'vowel']
plotsTitles =  ['GIF Methods', 'Vocal Effort', 'F0 Range', 'Vowels']
colLabels = ['GIF method', 'phonation type', 'f0', 'vowel']
colOrder = [[], ['creaky', 'normal', 'breathy', 'whispery'], ['low', 'mid', 'high'], ['I', 'E', 'AE', 'A', 'U', 'O']]
vowelLabels = ['i', 'e', 'æ', 'a', 'u', 'o']
f0bins = [0, 190, 280, 600]
f0Labels = ['low', 'mid', 'high']
genders = ['Male', 'Female']
# filter by Genders F0
F0Range = {"Male":[100, 240], "Female":[220, 360]}

# Paths
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

# discretize f0 into 3 categories
df["f0r"]=pd.cut(df["f0"],bins=f0bins,labels=f0Labels,right=False)

figuresDir = 'figures/paperVersion/BoxPlots'
# figuresDir = 'figures/paperVersion/BoxPlots/'
figuresPath = os.path.join(dataPath, figuresDir)
if not os.path.isdir(figuresPath):
    os.makedirs(figuresPath)

figureNameBase = 'BoxPlots_RMSE_NAQ_H1H2_HRF_ST'
if fixY:
	figureNameBase = figureNameBase + '_YFixed'
if filterF0Gender:
	figureNameBase = figureNameBase + '_F0Range'
if errorSign:
	figureNameBase = figureNameBase + '_Sign'
if hrfdB:
	figureNameBase = figureNameBase + '_hrfdB'
if standWhis:
	figureNameBase = figureNameBase + '_standWhis'
else:
	figureNameBase = figureNameBase + '_5-95'


# Convert error NAQ to %
df['naq'] = 100 * df['naq']
df.loc[(df['GIFmethod'] == 'QCP') & (df['STcompensation'] == 1), 'GIFmethod'] = 'ST-QCP'
df.loc[(df['GIFmethod'] == 'Original-IAIF'), 'GIFmethod'] = 'IAIF'
df = df[df['f0'] != 80]
# print('remove NaNs in naq')
df = df[~df['naq'].isin([np.inf, -np.inf])]
df = df[~df['naq'].isnull()]

gridRows = gridRowsArray[1] if errorSign else gridRowsArray[0]
nRow = len(gridRows);
nPlots = len(plots);
nGenders = len(genders)
axes = []
yAxesLimArray = [[0, 80], [0, 100], [0,30], [0, 15], [0, 40]]
wisker = 1.5 if standWhis else [5, 95]

for plotId, plot in enumerate(plots):
	fig, axes = plt.subplots(nRow, nGenders, sharey='row')
	for rowId in range(nRow):
		for colId in range(nGenders):
			df_gender = df[df['gender'] == genders[colId]]
			if(filterF0Gender):
				df_gender = df_gender[df_gender['f0'] >= F0Range[genders[colId]][0]]
				df_gender = df_gender[df_gender['f0'] <= F0Range[genders[colId]][1]]
			if len(colOrder[plotId])>0:
				g = sns.boxplot(data=df_gender,
								x=plots[plotId],
								y=gridRows[rowId],
								hue='GIFmethod',
								order=colOrder[plotId],
								whis=wisker,
								showfliers=False,
								ax=axes[rowId][colId]
								)
				if fixY:
					axes[rowId][colId].set_ylim(yAxesLimArray[rowId])
			else:
				g = sns.boxplot(data=df_gender,
								x=plots[plotId],
								y=gridRows[rowId],
								hue='GIFmethod',
								whis=wisker,
								showfliers=False,
								ax=axes[rowId][colId],
								dodge=False
								)
				if fixY:
					axes[rowId][colId].set_ylim(yAxesLimArray[rowId])
			
			if colId == 0:
				axes[rowId][colId].set_ylabel(errLabels[rowId])
			else:
				g.set(ylabel=None)
				axes[rowId][colId].tick_params(labelleft=False)

			if rowId == nRow-1:
				if (plot == 'vowel'):
					axes[rowId][colId].set_xticks(range(len(colOrder[plotId])))
					axes[rowId][colId].set_xticklabels(vowelLabels);
				if (plot == 'GIFmethod'):
					axes[rowId][colId].tick_params(labelbottom=False)
				axes[rowId][colId].set_xlabel(genders[colId])
			else:
				g.set(xlabel=None)
				axes[rowId][colId].tick_params(labelbottom=False)

			axes[rowId][colId].legend_.remove()
			axes[rowId][colId].yaxis.set_major_locator(MultipleLocator(20))
			axes[rowId][colId].yaxis.set_minor_locator(AutoMinorLocator(4))
			if rowId > 1:
				axes[rowId][colId].yaxis.set_major_locator(MultipleLocator(10))
			if rowId == 3:
				axes[rowId][colId].yaxis.set_major_locator(MultipleLocator(2))

			axes[rowId][colId].grid(visible=True, which='major', axis='y', linestyle='-')
			axes[rowId][colId].grid(visible=True, which='minor', axis='y', linestyle=':')

	if titlePlot:
		plt.suptitle(plotsTitles[plotId], y=0.98)
	fig.set_size_inches([10, 9])
	plt.subplots_adjust(left=0.09, bottom=0.07, right=0.99, top=0.95, hspace=0.12, wspace=0.03)
	plt.legend(loc='upper center', bbox_to_anchor=(0.116, 5.8), ncol=6, fancybox=True, shadow=True)
	if plot == 'GIFmethod':
		plt.subplots_adjust(bottom=0.06)
	fig.align_labels() 
	# plt.show()
	pathToFigureSVG = os.path.join(figuresPath, plot + '_' + figureNameBase + '.svg')
	pathToFigurePNG = os.path.join(figuresPath, plot + '_' + figureNameBase + '.png')
	pathToFigurePDF = os.path.join(figuresPath, plot + '_' + figureNameBase + '.pdf')
	plt.savefig(pathToFigureSVG)
	plt.savefig(pathToFigurePNG)
	plt.savefig(pathToFigurePDF)
