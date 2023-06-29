import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

sns.set_theme(style="whitegrid")

adduction = ['small', 'medium', 'large']
genres = ['Female', 'Male']
f0_female = [175, 196, 220, 294]
f0_male = [82, 110, 156, 220]
vowel = ['AE', 'A', 'I', 'U']
GIFmethod = ['Original-IAIF', 'IOP-IAIF', 'GFM-IAIF', 'QCP', 'QCP-ST']
semitonesThreshold = 4

# filesInPathBase = 'Projects/FEMVoQ/Corpus/OPENGLOT/repository_II/results'
dataPathBase = 'Projects/LaSalle/FEMVoQ/Corpus/OPENGLOT/specialIssue'
experimentDir = 'OpenGlotII16k'
dataDir = 'optimizationResults'
dataPath = os.path.join(os.path.expanduser('~'), dataPathBase, experimentDir, dataDir)
dataFileNames = [
    'repository_II_Original-IAIF_grid-search_pulse.csv',
    'repository_II_IOP-IAIF_grid-search_pulse.csv',
    'repository_II_GFM-IAIF_grid-search_pulse.csv',
    'repository_II_QCP_grid-search_pulse.csv',
    'repository_II_ST-QCP_grid-search_pulse.csv'
]
filterFileName = 'F0andSemitones_pulse_wavs_16k.csv'
filterFilePath = os.path.join(os.path.expanduser('~'), dataPathBase, filterFileName)

figuresDir = 'figures'
figuresPath = os.path.join(dataPath, figuresDir)
figureNameBase = 'GridErrors_RMSE_NAQ_H1H2_HRF_ST'
# figureNameBase = 'GridErrors_RMSE_NAQ_H1H2_HRF_ST_standWhis'
# figureNameBase = 'GridErrors_RMSE_NAQ_H1H2_HRF_ST_noYlim'
# figureNameBase = 'GridErrors_RMSE_NAQ_H1H2_HRF_ST_noYlim_standWhis'
if not os.path.isdir(figuresPath):
    os.makedirs(figuresPath)

# Read input files into DataFrame
dfFilter = pd.read_csv(filterFilePath)
dfFilterClean = dfFilter[~dfFilter['semitonesDifferenceGT'].str.contains('i', regex=False)] # openglotII 16k
# dfFilterClean = dfFilter[~dfFilter['semitonesDifferenceGT'].isin([np.inf, -np.inf])] # openglotII 8k

df = pd.DataFrame()
for dataFileName in dataFileNames:
    dataFilePath = os.path.join(dataPath, dataFileName)
    # print(dataFilePath)
    dfTemp = pd.read_csv(dataFilePath)
    # print(dfTemp.shape)
    # print('remove complex in semitonesDifferenceGT')
    dfTemp = dfTemp[~dfFilter['semitonesDifferenceGT'].str.contains('i', regex=False)]; # openglotII 16k
    # print('remove inf in semitonesDifferenceGT')
    # dfTemp = dfTemp[~dfFilter['semitonesDifferenceGT'].isin([np.inf, -np.inf])]; # openglotII 8k
    # print(dfTemp.shape)
    # print('remove semitonesDifferenceGT higher than'  + str(semitonesThreshold))
    dfTemp = dfTemp[dfFilterClean['semitonesDifferenceGT'].astype(float) < semitonesThreshold]
    # print(dfTemp.shape)
    # print('remove NaNs in naq')
    dfTemp = dfTemp[~dfTemp['naq'].isnull()]
    # print(dfTemp.shape)
    df = pd.concat([df, dfTemp])

# Convert error NAQ to %
df['naq'] = 100 * df['naq']
df.loc[(df['GIFmethod'] == 'QCP') & (df['STcompensation'] == 1), 'GIFmethod'] = 'QCP-ST'

# Plots by genre

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
gridCols = ['GIFmethod', 'adduction', 'f0', 'vowel']
colLabels = ['GIF method', 'adduction', 'f0', 'vowel']
nRow = len(gridRows);
nCol = len(gridCols);
axes = []

yAxesLimArray = [[0, 80], [0, 100], [0,30], [0, 15], [0, 40]]

for idx, genre in enumerate(genres):
    colOrder = [[], adduction, eval('f0_' + genre.lower()), vowel]
    df_genre = df[df['genre'] == genres[idx]]
    fig, axes = plt.subplots(nRow, nCol, gridspec_kw={'width_ratios': [5, 8, 6, 10]}, sharey='row')

    for rowId in range(nRow):
        for colId in range(nCol):
            if len(colOrder[colId]) > 0:
                g = sns.boxplot(data=df_genre,
                                x=gridCols[colId],
                                y=gridRows[rowId],
                                hue='GIFmethod',
                                order=colOrder[colId],
                                whis=[5, 95],
                                showfliers=False,
                                ax=axes[rowId][colId]
                                )
                axes[rowId][colId].set_ylim(yAxesLimArray[rowId])
            else:
                g = sns.boxplot(data=df_genre,
                                x=gridCols[colId],
                                y=gridRows[rowId],
                                hue='GIFmethod',
                                whis=[5, 95],
                                showfliers=False,
                                ax=axes[rowId][colId],
                                dodge=False
                                )
                axes[rowId][colId].set_ylim(yAxesLimArray[rowId])
            if colId == 0:
                axes[rowId][colId].set_ylabel(errLabels[rowId])
            else:
                g.set(ylabel=None)
                axes[rowId][colId].tick_params(labelleft=False)
            if rowId == nRow - 1:
                axes[rowId][colId].tick_params(axis='x', labelrotation=45, pad=0)
                axes[rowId][colId].set_xlabel(colLabels[colId])
            else:
                g.set(xlabel=None)
                axes[rowId][colId].tick_params(labelbottom=False)

            axes[rowId][colId].legend_.remove()

    plt.suptitle(figureNameBase + '_' + genre, y=0.03)
    fig.set_size_inches([15, 12])
    fig.align_labels()
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.99, top=0.98, hspace=0.12, wspace=0.03)
    # plt.show()
    pathToFigure = os.path.join(figuresPath, figureNameBase + '_' + genre + '.svg')
    plt.savefig(pathToFigure)

# Plots by gridCols

# variableNames = ['GIFmethod', 'adduction', 'vowel']

# for variableName in variableNames:
#     variableOptionsList = eval(variableName)
#     variableNumOptions = len(variableOptionsList)
#     axes = []
#     fig, axes = plt.subplots(nRow, variableNumOptions, sharey='row')
#     for rowId in range(nRow):
#         for colId, option in enumerate(variableOptionsList):
#             dfVariableName = df[df[variableName] == option]
#             sub = sns.boxplot(data=dfVariableName,
#                               x='GIFmethod',
#                               y=gridRows[rowId],
#                               hue='genre',
#                               whis=[5, 95],
#                               showfliers=False,
#                               # dodge = False,
#                               ax=axes[rowId][colId]
#                               )
#             if colId == 0:
#                 axes[rowId][colId].set_ylabel(errLabels[rowId])
#             else:
#                 sub.set(ylabel=None)
#                 axes[rowId][colId].tick_params(labelleft=False)

#             if rowId == nRow - 1:
#                 # axes[rowId][colId].tick_params(axis = 'x', labelrotation = 45, pad = 0)
#                 sub.set(xlabel=option)
#             else:
#                 sub.set(xlabel=None)
#                 axes[rowId][colId].tick_params(labelbottom=False)

#             if colId == variableNumOptions - 1 and rowId == nRow - 1:
#                 sns.move_legend(axes[rowId][colId], "upper right", bbox_to_anchor=(1, -0.25))
#             else:
#                 axes[rowId][colId].legend_.remove()

#     plt.suptitle(variableName.capitalize(), y=0.05)
#     fig.set_size_inches([13, 12])
#     plt.subplots_adjust(left=0.071, bottom=0.13, right=0.99, top=0.98, hspace=0.12, wspace=0.03)
#     # plt.show()
#     pathToFigure = os.path.join(figuresPath,
#                                 figureNameBase + '_' +
#                                 variableName.capitalize() + '.svg')
#     plt.savefig(pathToFigure)
