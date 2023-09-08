# FEMVoQ Inverse Filtering (IF) evaluation

## About
This repository contains all code for the evaluation of Inverse Filtering algorithms reviewed from literature.


##  1. corpora
In this folder you can find the code necesary to generate the corpora used to evaluate the GIF methods.

#### OPENGLOT_I

This directory contains all the necesary code to generate the extended version of the OPENGLOT repository I.

This is based on the code from [OPENGLOT repository I](http://research.spa.aalto.fi/projects/openglot/)

## 2. inverse-filtering-code

This folder contains the MATLAB code developed for running the optimization and the analysis of GIF algorithms.
Algorithms currently implemented:
1. IAIF
2. GFM-IAIF
3. IOP-IAIF
4. QCP
5. QCP with spectral tilt correction.

## 3. visualisation

Python code for the visualisation and statistical analysis of the results obtained for the evaluation of the GIF methods.

## 4. experiment

In this folder you can find the csv files with the results used in our paper: [Evaluation of Glottal Inverse Filtering Techniques on OPENGLOT Synthetic Male and Female Vowels](https://doi.org/10.3390/app13158775).
Also included are the `.ini` files used for our experiments.
The structure of subfolders is explained in the section `Using the code`.

# Using the code

This projects allows to run in two modes:
Optimisation: it runs the optimisation of parameters for the implemented GIF methods.
Analysis: executes the GIF methods given a set of parameters.

In both cases a collection of error measures are computed.

An example folder with the necessary files to execute the code is shared, called `experiment`.

## Execution

To execute the code the following comand is executed in MATLAB's Command Window:

`executeInverseFilterAnalysis('path_to_corpus/experiment_parent_directory/config.ini')`

## Experiments

The code to execute the experiments relay on configuration files to read the characteristics of said experiments.
There are two different types of configuration files.
- A general config file called  `config.ini`. This files defines the corpus and the paths and directory names.
- The GIF method configuration files, one for each method to be executed.

All the configuration files have the extension `.ini`.

### Folder structure

The main `config.ini` file needs to be in a the main experiment folder (experiment_parent_directory). Subfolders will contain the rest of files.
A first subfolder (GIF_config_files_dir) will contain the necesary GIF method configuration files.

Analysis: In the of case executing an analysis another subfolder is needed (GIF_parameters_files_dir). This one contains the csv files (one per each GIF method) with the parameters per audio file in the corpus.

The results of the execution together with the generated log files will be persisted in other subfolders (results_directory_name and logs_dir).

The structure resulting is the following.

```
experiment_parent_directory
│	config.ini
└─	GIF_config_files_dir 			(contains the GIF configuration files)
└─	GIF_parameters_files_dir		(only analysis, parameters to be used for each audio file)
└─	results_dir_name
		└─ mat 				(only optimisation, automatically created)
└─	logs_dir 				(automatically created)
audio_files_dir_name
```

### config.ini

The experiment will run following a configuration file named `config.ini`, which contains the names of folders and subfolders containing the rest of files.

```
[experiment]

repository = repository_name			; used for the results files
optimization = true				; false for analysis
filesFormat = (?<a>\w+)_(?<b>\w+)		; regex of the audio files

[paths]

corpusDir = audio_files_dir_name		; e.g. wavs

invFilterConfigDir = GIF_config_files_dir 	; e.g. inverseFilterOptimization or inverseFilterAnalysis
paramFilesDir = GIF_parameters_files_dir 	; only for analysis execution, e.g. analysisParams

resultsDir = results_dir_name			; e.g. optimizationResults or analysisResults

logDir = logs_dir
logFileName = logFile.txt

```


### GIF\_config\_files\_dir

For each GIF method a second configuration file is needed.
Examples can be found in `experiment/inverseFilterOptimization`.
This examples can be used as templates for further experiments.

### GIF\_parameters\_files\_dir

This is necesary for the analysis execution.
Inside this folder a csv file is mandatory for each GIF method to be executed.
This csv folder needs to contain a row for each audio file, with the parameters information per column.

Note: the output csv files from the optimisation work as input for the analysis. You can use the provided csv examples in `experiment/optimizationResults`.

### results\_directory\_name

The results are error measures computed after the inverse filtering is executed. Results will be persisted in this directory as csv files.
These measures are computed at a pulse basis (computed from GCIs).

For each gif method two csv files will be stored, one with the medians and means of the erros measures per audio file and another one with the raw error measures per pulse.

#### mat

In the case of optimisation, a subfolder named `mat` will be created where `.mat` files with all the grid-search results will be stored.


# Reference this work

If you use any part of this code, please reference our work as:
```bib
 @article {Freixes2023,
	AUTHOR = {Freixes, Marc and Joglar-Ongay, Luis and Socoró, Joan Claudi and Alias-Pujol, Francesc},
	TITLE = {Evaluation of Glottal Inverse Filtering Techniques on OPENGLOT Synthetic Male and Female Vowels},
	JOURNAL = {Applied Sciences},
	VOLUME = {13},
	YEAR = {2023},
	NUMBER = {15},
	ARTICLE-NUMBER = {8775},
	URL = {https://www.mdpi.com/2076-3417/13/15/8775},
	ISSN = {2076-3417},
	ABSTRACT = {Current articulatory-based three-dimensional source-filter models,
		which allow the production of vowels and diphtongs, still present very limited
		expressiveness. Glottal inverse filtering (GIF) techniques can become 
		instrumental to identify specific characteristics of both the glottal source 
		signal and the vocal tract transfer function to resemble expressive speech. 
		Several GIF methods have been proposed in the literature; however, their 
		comparison becomes difficult due to the lack of common and exhaustive 
		experimental settings. In this work, first, a two-phase analysis methodology 
		for the comparison of GIF techniques based on a reference dataset is 
		introduced. Next, state-of-the-art GIF techniques based on iterative adaptive 
		inverse filtering (IAIF) and quasi closed phase (QCP) approaches are 
		thoroughly evaluated on OPENGLOT, an open database specifically designed to 
		evaluate GIF, computing well-established GIF error measures after extending 
		male vowels with their female counterparts. The results show that GIF methods 
		obtain better results on male vowels. The QCP-based techniques significantly 
		outperform IAIF-based methods for almost all error metrics and scenarios and 
		are, at the same time, more stable across sex, phonation type, F0, and vowels. 
		The IAIF variants improve the original technique for most error metrics on 
		male vowels, while QCP with spectral tilt compensation achieves a lower 
		spectral tilt error for male vowels than the original QCP.},
	DOI = {10.3390/app13158775}
}
```

# Code Authors
Marc Freixes - [marc.freixes@salle.url.edu](marc.freixes@salle.url.edu)

Luis Joglar-Ongay - [luis.joglar@salle.url.edu](luis.joglar@salle.url.edu)

Joan Claudi Socoró - [joanclaudi.socoro@salle.url.edu](joanclaudi.socoro@salle.url.edu)

# License

This repository is under license [LGPL](https://choosealicense.com/licenses/lgpl-3.0/).

# Third Party Code modified

This repository makes use of modifed versions of software developed by other authors.

## corpora

### OPENGLOT_I

The code from [OPENGLOT repository I](http://research.spa.aalto.fi/projects/openglot/)
has been the base for the extended version.
In our version we add the possibility to generate the same set of vowels but with female vocal tract formant frequencies as described by:
```
Peterson, G.E.; Barney, H.L. 
Control methods used in a study of the vowels.
The Journal of the acoustical society of America 1952,
545 24, 175–184. [DOI](https://doi.org/10.1121/1.1906875.)
```

## inverse-filtering-code
### inverse filtering methods
The GIF methods in this project are based on open source versions as follow:

#### IAIF and IOP-IAIF
Based on the methods implemented in [Covarep](http://covarep.github.io/covarep/) with license LGPL

#### GFM-IAIF
Based on the implementation by [Olivier Perrotin](olivier.perrotin@gipsa-lab.grenoble-inp.fr) with license LGPL

#### QCP
Based on the implementation in [AaltoAparat](http://research.spa.aalto.fi/projects/aparat/legacy/Aalto_Aparat_sources.zip) an open source project (without explicit license)

### utils/vendor
#### Ini Config
This projects relies on a modified version of [INI Config](https://es.mathworks.com/matlabcentral/fileexchange/24992-ini-config) MATLAB code by [Evgeny Pr](https://es.mathworks.com/matlabcentral/profile/authors/1764615)

#### get\_vq\_params
This function from [Covarep](http://covarep.github.io/covarep/) computes NAQ, QOQ, H1H2, HRF and PSP. The modifed version computes the Spectral Tilt (ST) as well.

# Dependencies

Other methods from the following libraries and toolbox are used in this project.
To be able to execute this code you will need access to these libraries (accessed on June 2023):
#### [AaltoAparat](http://research.spa.aalto.fi/projects/aparat/legacy/Aalto_Aparat_sources.zip)
#### [Covarep](http://covarep.github.io/covarep/)
#### [VOICEBOX](http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html)
#### [egifa_toolbox](https://languageandvoice.files.wordpress.com/2017/03/egifa.zip)
