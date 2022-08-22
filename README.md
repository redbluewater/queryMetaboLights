# queryMetaboLights
### 18 August 2022
R code to query Kujawinski lab data repositories at MetaboLights

This repository is a first pass at setting up a way to query untargeted metabolomics data at MetaboLights with an elemental formula. This question has come up multple times, most recently at a meeting in late July. Voila, ''queryMetaboLights'' was started. This is a work in progress, but feel free to play with it and let me know what you think.

There are two R scripts here and a CSV file. 
1. queryMetaboLights_getFilenames.R --> takes a list of MTBLS numbers, and downloads the file names from MetaboLights. The output is a series of CSV files, one for each MTBLS number. These are messy, and need to be pruned manually to figure out where the untargeted data files are, and which ion mode you are searching. Automating this is harder, and for only a few datasets it was easier to do this manually.
2. MetaboLights_untargetedLookup_curated.csv --> is the set of MTBLS###s and filenames for data submitted by the Kujawinski laboratory. This is based on the pruning done from the data from step #1
3. queryMetabolites_latest.R --> this is the actual query.
