import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import os
sns.set_theme(style="whitegrid")
sns.set(style="whitegrid", rc={"lines.linewidth": 1})
paletteColors = ["#E69F00", "#0072B2", '#C0392B', "#F0E442", "#009E73", "#D55E00"]
palette = sns.color_palette(paletteColors)
sns.set_palette(palette)

filterF0Gender = False
fixY = False
hrfdB = True
errorSign = False
piStand = False

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
gridCols =  ['f0']
colLabels = ['GIF method', 'phonation type', 'f0', 'vowel']
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

figuresDir = 'figures/paperVersion/PointPlots'
figuresPath = os.path.join(dataPath, figuresDir)
if not os.path.isdir(figuresPath):
    os.makedirs(figuresPath)

figureNameBase = 'F0ErrorsBar_RMSE_NAQ_H1H2_HRF_ST' # pi percentile intercal
if fixY:
	figureNameBase = figureNameBase + '_YFixed'
if filterF0Gender:
	figureNameBase = figureNameBase + '_F0Range'
if errorSign:
	figureNameBase = figureNameBase + '_Sign'
if hrfdB:
	figureNameBase = figureNameBase + '_hrfdB'
if piStand:
	figureNameBase = figureNameBase + '_piStand'
else:
	figureNameBase = figureNameBase + '_pi50'

# Read input files into DataFrame
df = pd.DataFrame()
for dataFileName in dataFileNames:
    dataFilePath = os.path.join(dataPath, dataFileName)
    dfTemp = pd.read_csv(dataFilePath)
    df = pd.concat([df, dfTemp])

# Convert error NAQ to %
df['naq'] = 100 * df['naq']
df.loc[(df['GIFmethod'] == 'QCP') & (df['STcompensation'] == 1), 'GIFmethod'] = 'ST-QCP'
df.loc[(df['GIFmethod'] == 'Original-IAIF'), 'GIFmethod'] = 'IAIF'
df = df[~df['naq'].isin([np.inf, -np.inf])];
df = df[~df['naq'].isnull()]
df = df[df['f0'] != 80]
colOrder = df.f0.unique().tolist()
colOrder = sorted(colOrder)


gridRows = gridRowsArray[1] if errorSign else gridRowsArray[0]

nRow = len(gridRows)
nCol = len(gridCols)
nGenders = len(genders)
axes = []
yAxesLimArray = [[0, 60], [0, 100], [0,30], [0, 10], [0, 40]]
errorbar = "pi" if piStand else ('pi', 50)

fig, axes = plt.subplots(nRow, nGenders, sharey='row')
for colId, gender in enumerate(genders):
	df_gender = df[df['gender'] == genders[colId]]
	if(filterF0Gender):
		df_gender = df_gender[df_gender['f0'] >= F0Range[gender][0]]
		df_gender = df_gender[df_gender['f0'] <= F0Range[gender][1]]
	colOrder = df_gender.f0.unique().tolist()
	colOrder = sorted(colOrder)
	for rowId in range(nRow):
		g = sns.pointplot(
				data = df_gender,
				x = gridCols[0],
				y = gridRows[rowId],
				hue = 'GIFmethod',
				estimator='median',
				order = colOrder,
				errwidth = 2,
				errorbar = errorbar,
				dodge = 0.5,
				ax = axes[rowId][colId]
				)
		if fixY:
			axes[rowId][colId].set_ylim(yAxesLimArray[rowId])

		if colId==0:
			axes[rowId][colId].set_ylabel(errLabels[rowId])
			# axes[rowId][colId].axvspan(-0.5, 8, color='y', alpha=0.5, lw=0)
		else:
			g.set(ylabel=None)
			axes[rowId][colId].tick_params(labelleft=False)
			# axes[rowId][colId].axvspan(8, 13, color='y', alpha=0.5, lw=0)

		axes[rowId][colId].set_xlabel(genders[colId])
		axes[rowId][colId].legend_.remove()
		axes[rowId][colId].yaxis.set_major_locator(MultipleLocator(20))
		axes[rowId][colId].yaxis.set_minor_locator(AutoMinorLocator(4))
		if rowId > 1:
			axes[rowId][colId].yaxis.set_major_locator(MultipleLocator(10))
		if rowId == 3:
			axes[rowId][colId].yaxis.set_major_locator(MultipleLocator(2))

		axes[rowId][colId].grid(visible=True, which='major', axis='y', linestyle='-')
		axes[rowId][colId].grid(visible=True, which='minor', axis='y', linestyle=':')

# plt.legend(loc='lower center', bbox_to_anchor=(0, -1), ncol=3, fancybox=True, shadow=True)
plt.legend(loc='upper center', bbox_to_anchor=(0.284, 5.7), ncol=6, fancybox=True, shadow=True)

if titlePlot:
	plt.suptitle(figureNameBase)
fig.set_size_inches([12, 15])
fig.align_labels()
plt.subplots_adjust(left=0.075, bottom=0.05, right=0.99, top=0.96, hspace=0.12, wspace=0.03)
# plt.show()
pathToFigureSVG = os.path.join(figuresPath, figureNameBase + '.svg')
pathToFigurePNG = os.path.join(figuresPath, figureNameBase + '.png')
pathToFigurePDF = os.path.join(figuresPath, figureNameBase + '.pdf')
plt.savefig(pathToFigureSVG)
plt.savefig(pathToFigurePNG)
plt.savefig(pathToFigurePDF)
