<<<<<<< HEAD
#Working on how to query MetaboLights by elemental formula 
#seems obvious, but not readily available (and a question that keeps coming up)
#
#input from user: 
#code will ask for an elemental formula, this should be provide like this: C24H41N3O16
#(e.g., no spaces, dots, etc.) and is the neutral elemental formula ([M])
#output
#
#the code will send out a CSV file with all the hits in the study_id, positive or negative 
#ion mode, assuming [M+H]+ in positive ion mode and [M-H]- in negative ion mode
#Next steps: 
#(1) figure out the best way to set this up with some/all of the Kuj repositories
#(2) Consider moving this to the HPC as it may get slow with too many repositories
#
#Krista Longnecker, 17 August 2022
#Woods Hole Oceanographic Institution

library(metabolighteR)
library(pracma)
library(tibble)
library(reshape2)
library(Rdisop)
library(MetaboCoreUtils) #from Jo Rainer - has all handy math (mz,formula, adducts)

#tidy up first
rm(list = ls.str())

#start with setting the study_id, later will need to loop through 
#different studies; this is a great example of changes in repositories over time!
study_id <- 'MTBLS144' #Tps4  

#other studies to add later
# #study_id <- 'MTBLS154'#Tps6 (example of targeted and untargeted)
# MTBLS366 #Lohmann cruise
# MTBLS1024 #DeepDOM DOS project (direct infusion)
# MTBLS1820 #Jardines de la Reina
# MTBLS461 #dilution experiment neg/untargeted, separate file with targeted
# MTBLS293 #NB pos/untar, neg/untar, targeted
# MTBLS155 #Syn2 pos/untar, neg/untar, targeted
# MTBLS295 #Cara's MP1
# MTBLS157 #Rpom

#set the search window in ppm
maxErr <- 10

#note that I added this to download_study_file so that it would stop sending up error about encoding
#KL August 9, 2022; used trace; doesn't seem to have this edit?
#trace("download_study_file",edit=TRUE)
#was: file_content <- httr::content(file_download, "text")
#is: file_content <- httr::content(file_download, "text", encoding = "UTF-8")

#Define the function to do the searching up here 
FindAllWithinError <- function(inputList, targetMZ, maxErr) {
  #remember this will return the index to the closest value
  #KL notes: targetMZ is the ONE mz value I am searching for
  #inputList is the full list of mz values that are possibilities
  
  #sort inputList 
  se <- sort(inputList,index.return = TRUE)
  inputSorted <- se$x
  idxSortInput <- se$ix
  
  if (length(targetMZ)==1) {
    ti = 1 
    } else {
    stop('Something is wrong - should only be one targetMZ value')
  }
  
  #will only do this for the ppm calculation
  #do the math in two steps to make sure I get this right
  numerator <- abs(inputSorted - targetMZ[ti])
  e <- numerator/targetMZ[ti]
  ii <- which(e<=maxErr/1e6)

  #return indices wrt original (unsorted) input
  idxIn = idxSortInput[ii] #idxIn works for the original, unsorted list
  
  return(idxIn)
}

#other useful functions in library(metabolighteR)
#gs <- get_studies()
#get_study_contacts('MTBLS144')
#get_study_factors(study_id)
#get_study_samples(study_id) #this spits back an error; do manually

## spits back list of all studies and the technology used, can use this to
## filter out the NMR studies
#gst <- get_study_tech()

#this will provide the names and types of the files, but not the actual data files
d <- get_study_files(study_id, raw_data = FALSE)
fileList <- d$file
rm(d)

#here - need to use the assay file to figure out the maf file name AND ion mode
#(can have no maf file listed)
#pull the a_ files so I search through them
g <- grepl("^a_",fileList)
idxA <- which(g)

#this makes the list of a_ files the following:
fileList_a <- fileList[idxA]

track <- data.frame(matrix(NA, nrow = length(idxA), ncol = 3))
colnames(track) <- c('aFile','ionMode','mafFile')

for (a in idxA)
{
  #put in the filename first...that's the easy step
  oneFile <- fileList_a[idxA[a]]
  track[a,'aFile'] <- oneFile
  #now go into that file and pull the ion mode
  getFile <- download_study_file(study_id, fileList_a[idxA[a]])
  #this will rely on single ion mode in the data file...could be a problem
  track[a,'ionMode'] <- unique(getFile$Parameter.Value.Scan.polarity.)
  
  #and see if I have a maf file
  track[a,'mafFile'] <- unique(getFile$Metabolite.Assignment.File)
}

#now I know both the name of the maf files and which a_ have a maf file
#probably easier to just delete rows with no maf files
track <- na.omit(track)

#get the elemental formula, neutral version before I start messing around with different ion modes
frml <- readline(prompt ="enter your neutral elemental formula, no periods: ")

#no sense in recreating the wheel, MetaboCoreUtils will do the math steps
#frml <- "C9H17NO5" #pantothenic acid
#frml <- "C24H41N3O16" #triacetylchitotriose #this is in Tps4, pos mode
exactMass <- calculateMass(frml)

#now go through track, should be two iterations of the loop, one per ion mode
for (idxRow in 1:nrow(track))
{
  #idxRow <- 1  #fixed for testing one iteration of the loop
  #now set up the look to go through one maf file at a time 
  one_ionMode <- track[idxRow,'ionMode']
  
  one_A <- download_study_file(study_id, track[idxRow,'aFile'])
  one_maf <- download_study_file(study_id, track[idxRow,'mafFile'])
  
  #pull sample names
  sampleNames <- one_A$Sample.Name
  
  #will need to set ion mode in the next row based on what I learn above about ion mode
  chargedMass <- mass2mz(exactMass,adduct = adductNames(polarity=c(one_ionMode)))
  cn <- colnames(chargedMass)
  
  if (one_ionMode=="positive") {
      idxA <- which(cn=="[M+H]+")
    } else {
      idxA <- which(cn=="[M-H]-")
    }
  
  targetMZ <- chargedMass[idxA] 
 
  #spit back the study_id, file, abundances in a sample
  inputList <- one_maf$mass_to_charge
  
  #now search in the maf file
  idx <- FindAllWithinError(inputList,targetMZ,maxErr)
  
  #have extra columns here - in the maf file, 
  keep <- c('mass_to_charge','retention_time',sampleNames)
  
  ##now that I have the index returned, use that to pull the row I want in
  #the original datafile
  getData <- one_maf[idx,keep]
  
  #I think the best way to export this is to use this format:
  #study / mz / rt / sample / area/height
  id.vars <- c('mass_to_charge','retention_time')
  getData <- melt(getData, id.vars)
  getData <- add_column(getData,ionMode = one_ionMode, .before = "mass_to_charge")
  getData <- add_column(getData,study_id = study_id, .before = "ionMode")
  
  #now that I have the data from one ion mode, stuff it into data.frame
  #that will be exported after both iterations of the loop are done
  if (exists("forExport")) {
    #add more data
    forExport <- rbind(forExport,getData)
  } else {
      forExport <- getData
  }
  
}##the loop for rows in track ends here

#save forExport, out of R so I can gather up the pieces later
#make a simple file name base
fbase <- paste0('export',study_id)

#saveRDS(forExport, file = paste0(fbase,'.rds'))
write.csv(forExport,file  = paste0(fbase,'.csv'))
=======
#Working on how to query MetaboLights by elemental formula 
#seems obvious, but not readily available (and a question that keeps coming up)
#input from user: 
#code will ask for an elemental formula, this should be provide like this: C24H41N3O16
#(e.g., no spaces, dots, etc.) and is the neutral elemental formula ([M])
#output
#the code will send out a CSV file with all the hits in the study_id, positive or negative 
#ion mode, assuming [M+H]+ in positive ion mode and [M-H]- in negative ion mode
#Next steps: 
#(1) figure out the best way to set this up with some/all of the Kuj repositories
#(2) Consider moving this to the HPC as it may get slow with too many repositories
#
#Krista Longnecker, 17 August 2022
#Woods Hole Oceanographic Institution

library(metabolighteR)
library(pracma)
library(tibble)
library(reshape2)
library(Rdisop)
library(MetaboCoreUtils) #from Jo Rainer - has all handy math (mz,formula, adducts)

#tidy up first
rm(list = ls.str())

#start with setting the study_id, later will need to loop through 
#different studies; this is a great example of changes in repositories over time!
study_id <- 'MTBLS144' #Tps4  

#other studies to add later
# #study_id <- 'MTBLS154'#Tps6 (example of targeted and untargeted)
# MTBLS366 #Lohmann cruise
# MTBLS1024 #DeepDOM DOS project (direct infusion)
# MTBLS1820 #Jardines de la Reina
# MTBLS461 #dilution experiment neg/untargeted, separate file with targeted
# MTBLS293 #NB pos/untar, neg/untar, targeted
# MTBLS155 #Syn2 pos/untar, neg/untar, targeted
# MTBLS295 #Cara's MP1
# MTBLS157 #Rpom

#set the search window in ppm
maxErr <- 10

#note that I added this to download_study_file so that it would stop sending up error about encoding
#KL August 9, 2022; used trace; doesn't seem to have this edit?
#trace("download_study_file",edit=TRUE)
#was: file_content <- httr::content(file_download, "text")
#is: file_content <- httr::content(file_download, "text", encoding = "UTF-8")

#Define the function to do the searching up here 
FindAllWithinError <- function(inputList, targetMZ, maxErr) {
  #remember this will return the index to the closest value
  #KL notes: targetMZ is the ONE mz value I am searching for
  #inputList is the full list of mz values that are possibilities
  
  #sort inputList 
  se <- sort(inputList,index.return = TRUE)
  inputSorted <- se$x
  idxSortInput <- se$ix
  
  if (length(targetMZ)==1) {
    ti = 1 
    } else {
    stop('Something is wrong - should only be one targetMZ value')
  }
  
  #will only do this for the ppm calculation
  #do the math in two steps to make sure I get this right
  numerator <- abs(inputSorted - targetMZ[ti])
  e <- numerator/targetMZ[ti]
  ii <- which(e<=maxErr/1e6)

  #return indices wrt original (unsorted) input
  idxIn = idxSortInput[ii] #idxIn works for the original, unsorted list
  
  return(idxIn)
}

#other useful functions in library(metabolighteR)
#gs <- get_studies()
#get_study_contacts('MTBLS144')
#get_study_factors(study_id)
#get_study_samples(study_id) #this spits back an error; do manually

## spits back list of all studies and the technology used, can use this to
## filter out the NMR studies
#gst <- get_study_tech()

#this will provide the names and types of the files, but not the actual data files
d <- get_study_files(study_id, raw_data = FALSE)
fileList <- d$file
rm(d)

#here - need to use the assay file to figure out the maf file name AND ion mode
#(can have no maf file listed)
#pull the a_ files so I search through them
g <- grepl("^a_",fileList)
idxA <- which(g)

#this makes the list of a_ files the following:
fileList_a <- fileList[idxA]

track <- data.frame(matrix(NA, nrow = length(idxA), ncol = 3))
colnames(track) <- c('aFile','ionMode','mafFile')

for (a in idxA)
{
  #put in the filename first...that's the easy step
  oneFile <- fileList_a[idxA[a]]
  track[a,'aFile'] <- oneFile
  #now go into that file and pull the ion mode
  getFile <- download_study_file(study_id, fileList_a[idxA[a]])
  #this will rely on single ion mode in the data file...could be a problem
  track[a,'ionMode'] <- unique(getFile$Parameter.Value.Scan.polarity.)
  
  #and see if I have a maf file
  track[a,'mafFile'] <- unique(getFile$Metabolite.Assignment.File)
}

#now I know both the name of the maf files and which a_ have a maf file
#probably easier to just delete rows with no maf files
track <- na.omit(track)

#get the elemental formula, neutral version before I start messing around with different ion modes
frml <- readline(prompt ="enter your neutral elemental formula, no periods: ")

#no sense in recreating the wheel, MetaboCoreUtils will do the math steps
#frml <- "C9H17NO5" #pantothenic acid
#frml <- "C24H41N3O16" #triacetylchitotriose #this is in Tps4, pos mode
exactMass <- calculateMass(frml)

#now go through track, should be two iterations of the loop, one per ion mode
for (idxRow in 1:nrow(track))
{
  #idxRow <- 1  #fixed for testing one iteration of the loop
  #now set up the look to go through one maf file at a time 
  one_ionMode <- track[idxRow,'ionMode']
  
  one_A <- download_study_file(study_id, track[idxRow,'aFile'])
  one_maf <- download_study_file(study_id, track[idxRow,'mafFile'])
  
  #pull sample names
  sampleNames <- one_A$Sample.Name
  
  #will need to set ion mode in the next row based on what I learn above about ion mode
  chargedMass <- mass2mz(exactMass,adduct = adductNames(polarity=c(one_ionMode)))
  cn <- colnames(chargedMass)
  
  if (one_ionMode=="positive") {
      idxA <- which(cn=="[M+H]+")
    } else {
      idxA <- which(cn=="[M-H]-")
    }
  
  targetMZ <- chargedMass[idxA] 
 
  #spit back the study_id, file, abundances in a sample
  inputList <- one_maf$mass_to_charge
  
  #now search in the maf file
  idx <- FindAllWithinError(inputList,targetMZ,maxErr)
  
  #have extra columns here - in the maf file, 
  keep <- c('mass_to_charge','retention_time',sampleNames)
  
  ##now that I have the index returned, use that to pull the row I want in
  #the original datafile
  getData <- one_maf[idx,keep]
  
  #I think the best way to export this is to use this format:
  #study / mz / rt / sample / area/height
  id.vars <- c('mass_to_charge','retention_time')
  getData <- melt(getData, id.vars)
  getData <- add_column(getData,ionMode = one_ionMode, .before = "mass_to_charge")
  getData <- add_column(getData,study_id = study_id, .before = "ionMode")
  
  #now that I have the data from one ion mode, stuff it into data.frame
  #that will be exported after both iterations of the loop are done
  if (exists("forExport")) {
    #add more data
    forExport <- rbind(forExport,getData)
  } else {
      forExport <- getData
  }
  
}##the loop for rows in track ends here

#save forExport, out of R so I can gather up the pieces later
#make a simple file name base
fbase <- paste0('export',study_id)

#saveRDS(forExport, file = paste0(fbase,'.rds'))
write.csv(forExport,file  = paste0(fbase,'.csv'))
>>>>>>> 29340b9526206c0ae9491f0fe19a78532ca629d0
