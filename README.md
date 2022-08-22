# queryMetaboLights
### 18 August 2022
R code to query Kujawinski lab data repositories at MetaboLights
Krista Longnecker

This repository is a first pass at setting up a way to query untargeted metabolomics data at MetaboLights with an elemental formula. This question has come up multple times, most recently at a meeting in late July. Voila, ''queryMetaboLights'' was started. This is a work in progress, but feel free to play with it and let me know what you think.

## Quick version of what to do with this:
``queryMetabolites_latest.R`` --> this is R script that will run a search based on an elemental formula you enter. Run this in your R environment (I use R Studio), when prompted, enter the elemental formula you want, and wait. To run enter this: ```source('~/Dropbox/XX_MetaboLights_access/queryMetaboLights_latest.R')```\
One example of a compound I know you will find in some, but not all, studies this one: C24H41N3O16 (triacetylchitotriose)

## Other files that are in this repository
* ``queryMetaboLights_getFilenames.R`` --> takes a list of MTBLS numbers, and downloads the file names from MetaboLights. The output is a series of CSV files, one for each MTBLS number. These are messy, and need to be pruned manually to figure out where the untargeted data files are, and which ion mode you are searching. Automating this is harder, and for only a few datasets it was easier to do this manually. I have already done this for the Kujawinski lab datasets.

* ``MetaboLights_untargetedLookup_curated.csv`` --> is the set of MTBLS###s and filenames for data submitted by the Kujawinski laboratory. This is based on the pruning done from the data from step #1. You will need this code in the same directory as the query script.


## Some notes
This is not the fastest program in the world because it goes to MetaboLights to get the data files. For just the Kujawinski lab files it will tak <10 minutes. Expanding across more of MetaboLights will get progressively slower and slower, and at some point it's probably worth migrating this to an HPC computing environment where you can run each study in it's own little 'computing box'
