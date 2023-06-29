import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

sns.set_theme(style="whitegrid")

gridRows = ['rmseNorm',
            # 'naq',
            # 'h1h2',
            # 'hrfdB',
            'errorSpectralTilt'
            ]
errLabels = ["RMS distance\n(%)",
             # "NAQ error\n(%)",
             # "H1H2 error\n(dB)",
             # "HRF error\n(dB)",
             "ST error\n(dB/decade)"
             ]

gridCols =  ['GIFmethod']
colLabels = ['a', 'b', 'c', 'd', 'e', 'f']
genres = ['Female', 'Male']

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
figuresDir = 'figures/ViolinPlots'
figuresPath = os.path.join(dataPath, figuresDir)
# figureNameBase = 'GridErrorsViolin_RMSE_NAQ_H1H2_HRF_ST'
figureNameBase = 'GridErrorsViolin_RMSE_NAQ_H1H2_HRF_ST_YFixed'
if not os.path.isdir(figuresPath):
    os.makedirs(figuresPath)

# Read input files into DataFrame
df = pd.DataFrame()
for dataFileName in dataFileNames:
    dataFilePath = os.path.join(dataPath, dataFileName)
    dfTemp = pd.read_csv(dataFilePath)
    df = pd.concat([df, dfTemp])

df.loc[(df['GIFmethod'] == 'QCP') & (df['STcompensation'] == 1), 'GIFmethod'] = 'QCP-ST'
# remove nans and infs from naq
df['naq'] = 100 * df['naq']
df = df[~df['naq'].isin([np.inf, -np.inf])];
df = df[~df['naq'].isnull()]

df = df[~(df['rmseNorm'] > 50)]

nRow = len(gridRows);
nCol = len(gridCols);
axes = []
axId = 0

yAxesLimArray = [[0, 80], [0, 100], [0,30], [0, 15], [0, 40]]

for idx, genre in enumerate(genres):
	fig, axes = plt.subplots(nRow, nCol)
	df_genre = df[df['gender'] == genres[idx]]
	for rowId in range(nRow):
		sns.violinplot(data = df_genre,
					x = gridCols[0],
					y = gridRows[rowId],
					hue = 'GIFmethod',
					width=0.3,
					# whis = [5, 95],
					# showfliers = False,
					ax = axes[rowId]
					)
		# axes[rowId].set_ylim(yAxesLimArray[rowId])
		axes[rowId].legend_.remove()
		# axes[rowId].axis('off')
		axes[rowId].get_xaxis().set_visible(False)
		# if rowId != 1:
		# 	axes[rowId].set_xlabel('')

		# plt.axis('off')
	plt.suptitle(figureNameBase + '_' + genre, y=0.03)
	# fig.set_size_inches([8, 10])
	fig.align_labels() 
	# fig.update_layout(
	#     violinmode='group',
	#     violingap=0, 
	#     violingroupgap=0,
	#     width=600,
	#     height=500,
	#     margin=dict(
	#         l=0,
	#         r=0,
	#         t=0,
	#         b=0
	#     ),
	    # plot_bgcolor='white'
	# )

	plt.subplots_adjust(left=0.1, bottom=0.13, right=0.99, top=0.98, hspace=0.12, wspace=0.03)
	plt.show()
	# pathToFigure = os.path.join(figuresPath, figureNameBase + '_' + genre + '.svg')
	# plt.savefig(pathToFigure)
