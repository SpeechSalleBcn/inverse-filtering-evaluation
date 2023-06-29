import numpy as np
import pandas as pd
# import seaborn as sns
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.colors import n_colors
import os

# sns.set_theme(style="whitegrid")

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

GIFMethods = pd.unique(df[gridCols[0]]);
print(GIFMethods)
# print(df[gridRows[1]])
# exit()
	
colors = n_colors('rgb(5, 200, 200)', 'rgb(200, 10, 10)', len(GIFMethods), colortype='rgb')

for idx, genre in enumerate(genres):
	fig = make_subplots(rows = nRow,
					    cols = nCol,
					    specs = [[{"type": "violin"}],
					             [{"type": "violin"}],
					             [{"type": "violin"}],
					             [{"type": "violin"}],
					             [{"type": "violin"}]])
	fig2 = make_subplots(rows = nRow,
					    cols = nCol,
					    specs = [[{"type": "box"}],
					             [{"type": "box"}],
					             [{"type": "box"}],
					             [{"type": "box"}],
					             [{"type": "box"}]])
	df_genre = df[df['gender'] == genres[idx]]
	for rowId in range(nRow):
		for GIFMethod, color in zip(GIFMethods, colors):
			fig.add_trace(
			    go.Violin(
						x = df_genre[gridCols[0]][df_genre[gridCols[0]] == GIFMethod],
						y = df_genre[gridRows[rowId]][df_genre[gridCols[0]] == GIFMethod],
						name = GIFMethod,
						# box_visible = True,
						meanline_visible = True,
						line_color = color
						),
		    	row = rowId + 1,
		    	col = nCol
		    )
			fig2.add_trace(
			    go.Box(
						x = df_genre[gridCols[0]][df_genre[gridCols[0]] == GIFMethod],
						y = df_genre[gridRows[rowId]][df_genre[gridCols[0]] == GIFMethod],
						name = GIFMethod,
						line_color = color
						),
		    	row = rowId + 1,
		    	col = nCol
		    )

	fig.show()
	fig2.show()



