#
# Copyright (C) 2023 HER Speech, La Salle - Universitat Ramon Llull
#
# This file is part of the project Inverse Filtering Evaluation repository
# https://github.com/SpeechSalleBcn/inverse-filtering-evaluation
#
# License
# This file is under the LGPL license,  you can
# redistribute it and/or modify it under the terms of the GNU Lesser General 
# Public License as published by the Free Software Foundation, either version 3 
# of the License, or (at your option) any later version. This file is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
# PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# https://www.gnu.org/licenses/lgpl-3.0.html
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# Authors
# Marc Freixes - marc.freixes@salle.url.edu
#

import pandas as pd
import os

gridRows = ['rmseNorm',
            'naq',
            'h1h2',
            'hrfdB',
            'errorSpectralTilt'
            ]
genders = ['Female', 'Male']

# Load data
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

resultsDir = 'tables'
resultsPath = os.path.join(os.path.expanduser('~'), dataPathBase, experimentDir, dataDir, resultsDir)
if not os.path.isdir(resultsPath):
	os.makedirs(resultsPath)
resultsFileName = 'TableSpecialIssue2023_' + experimentDir + '_medianValues'
resultsFilePath = os.path.join(resultsPath, resultsFileName + '.csv')

df = pd.DataFrame()
for dataFileName in dataFileNames:
	dataFilePath = os.path.join(dataPath, dataFileName)
	dfTemp = pd.read_csv(dataFilePath)
	df = pd.concat([df, dfTemp])

df['naq'] = 100 * df['naq']
df = df[df['f0'] != 80]

df.loc[(df['GIFmethod'] == 'QCP') & (df['STcompensation'] == 1), 'GIFmethod'] = 'QCP-ST'

df.groupby("GIFmethod")[gridRows].median().to_csv(resultsFilePath)

for idx, gender in enumerate(genders):
	df_gender = df[df['gender'] == genders[idx]]
	resultsFilePath = os.path.join(resultsPath, resultsFileName + '_' + gender + '.csv')
	df_gender.groupby("GIFmethod")[gridRows].quantile([0.25,0.5,0.75]).to_csv(resultsFilePath)
    # df_gender.groupby("GIFmethod")[gridRows].median().to_csv(resultsFilePath)
