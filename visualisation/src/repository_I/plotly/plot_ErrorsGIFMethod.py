import numpy as np
import pandas as pd
# import seaborn as sns
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.colors import n_colors
import plotly.express as px
import os

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

plotCols =  ['GIFmethod']
# colLabels = ['a', 'b', 'c', 'd', 'e', 'f']
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

nRow = len(gridRows);
nCol = len(genres);
axes = []
axId = 0

GIFMethods = pd.unique(df[plotCols[0]]);
colors = [n_colors('rgb(210, 120, 60)', 'rgb(100, 40, 20)', len(GIFMethods), colortype='rgb'),
		  n_colors('rgb(50, 170, 80)', 'rgb(30, 90, 40)', len(GIFMethods), colortype='rgb')]
sides = ['positive', 'negative']

fig = make_subplots(rows = nRow,
				    cols = nCol,
				    specs = [[{"type": "violin"}, {"type": "violin"}],
				             [{"type": "violin"}, {"type": "violin"}],
				             [{"type": "violin"}, {"type": "violin"}],
				             [{"type": "violin"}, {"type": "violin"}],
				             [{"type": "violin"}, {"type": "violin"}]])
fig2 = make_subplots(rows = nRow,
				    cols = nCol,
				    specs = [[{"type": "box"}, {"type": "box"}],
				             [{"type": "box"}, {"type": "box"}],
				             [{"type": "box"}, {"type": "box"}],
				             [{"type": "box"}, {"type": "box"}],
				             [{"type": "box"}, {"type": "box"}]])
fig3 = make_subplots(rows = nRow,
				    cols = 1,
				    specs = [[{"type": "violin"}],
				             [{"type": "violin"}],
				             [{"type": "violin"}],
				             [{"type": "violin"}],
				             [{"type": "violin"}]])

yAxesLimArray = [[0, 80], [0, 100], [0,30], [0, 15], [0, 40]]

for rowId in range(nRow):
	for GIFId, GIFMethod in enumerate(GIFMethods):
	# for colId in range(nCol):
		for genreId, genre in enumerate(genres):
			df_genre = df[df['gender'] == genres[genreId]]
			fig.add_trace(
			    go.Violin(
						x = df_genre[plotCols[0]][df_genre[plotCols[0]] == GIFMethod],
						y = df_genre[gridRows[rowId]][df_genre[plotCols[0]] == GIFMethod],
						name = GIFMethod,
						# box_visible = True,
						meanline_visible = True,
						line_color = colors[genreId][GIFId]
						),
		    	row = rowId + 1,
		    	col = genreId + 1
		    )
			fig.update_yaxes(range=yAxesLimArray[rowId], row = rowId + 1, col = genreId + 1)


			fig2.add_trace(
			    go.Box(
						x = df_genre[df_genre[plotCols[0]] == GIFMethod][plotCols[0]],
						y = df_genre[df_genre[plotCols[0]] == GIFMethod][gridRows[rowId]],
						name = GIFMethod,
						# boxpoints = False,
						line_color = colors[genreId][GIFId]
						),
			    # layout_yaxis_range=[-1,80],
		    	row = rowId + 1,
		    	col = genreId + 1
			)
			fig2.update_yaxes(range=yAxesLimArray[rowId], row = rowId + 1, col = genreId + 1)

			fig3.add_trace(
		    	go.Violin(
						x = df_genre[plotCols[0]][df_genre[plotCols[0]] == GIFMethod],
						y = df_genre[gridRows[rowId]][df_genre[plotCols[0]] == GIFMethod],
						name = GIFMethod,
						# box_visible = True,
						meanline_visible = True,
						side = sides[genreId],
						line_color = colors[genreId][GIFId]
						),
		    	row = rowId + 1,
		    	col = 1
		    )
			fig3.update_yaxes(range=yAxesLimArray[rowId], row = rowId + 1, col = genreId + 1)

fig.show()
fig2.show()
fig3.show()
