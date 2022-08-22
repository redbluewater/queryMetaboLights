#Working on how to query MetaboLights by elemental formula 
#seems obvious, but not readily available (and a question that keeps coming up)
#input from user: 
#code will ask for an elemental formula, this should be provide like this: C24H41N3O16
#(e.g., no spaces, dots, etc.) and is the neutral elemental formula ([M])
#output
#"C24H41N3O16" #triacetylchitotriose will work in MTBLS144
#the code will send out a CSV file with all the hits in the study_id, positive or negative 
#ion mode, assuming [M+H]+ in positive ion mode and [M-H]- in negative ion mode
#Next steps: 
#(1) figure out the best way to set this up with some/all of the Kuj repositories
#(2) Consider moving this to the HPC as it may get slow with too many repositories
#
#This version loads in the pre-sorted Kujawinsk lab submissions from the CSV file,
#small-scale hack, but not awesome for all MTBLS repositories
#Krista Longnecker, 19 August 2022
#Woods Hole Oceanographic Institution

library(metabolighteR)
library(pracma)
library(tibble)
library(reshape2)
library(Rdisop)
library(MetaboCoreUtils) #from Jo Rainer - has all handy math (mz,formula, adducts)

#tidy up first
rm(list = ls.str())

#set the search window in ppm
maxErr <- 10

#note that I added this to download_study_file so that it would stop sending up error about encoding
#KL August 9, 2022; used trace; doesn't seem to have this edit?
#trace("download_study_file",edit=TRUE)
#was: file_content <- httr::content(file_download, "text")
#is: file_content <- httr::content(file_download, "text", encoding = "UTF-8")
#apparently trace does not persist across sessions, so define my own function with the desired edit
download_study_fileKL <- function (study_id, filename) 
{
  file_download <- httr::GET(paste0(getOption("BASE_URL"), 
                      "/studies/", study_id, "/download/public?file=", filename))
  
  file_content <- httr::content(file_download, "text", encoding = "UTF-8") %>% 
    textConnection() %>% readLines() %>% utils::read.delim(text = ., 
                      sep = "\t") %>% tibble::as_tibble()
  return(file_content)
}

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

#get the elemental formula, neutral version before I start messing around with different ion modes
frml <- readline(prompt ="enter your neutral elemental formula, no periods: ")

#no sense in recreating the wheel, MetaboCoreUtils will do the math steps
#frml <- "C9H17NO5" #pantothenic acid
#frml <- "C24H41N3O16" #triacetylchitotriose #this is in Tps4, pos mode
exactMass <- calculateMass(frml)

#read in the CSV file with the pre-sorted file information I need
mfInfo <- read.csv("MetaboLights_untargetedLookup_prunedC24H41N3O16.csv")

#only search for MTBLS numbers found in the lookup table
uM <- unique(mfInfo$MTBLS)

#loop through the unique studies here
for (a in 1:length(uM))
{
  study_id <- uM[a]
  print(paste0("Working on ",study_id," --> total studies to search: ",length(uM)))
  #figure out which row(s) I need in mfInfo
  r <- which(mfInfo$MTBLS==study_id)
  
  #go through what I find for this study, usually two iterations of the loop, one per ion mode
  #it's possible there will be a case with only one ion mode, I honestly do 
  #not remember if we have done that
  for (useR in 1:length(r))
  {
    idxRow <- r[useR]
    #now set up the look to go through one maf file at a time 
    one_ionMode <- mfInfo[idxRow,'ionMode']
    
    one_A <- download_study_fileKL(study_id, mfInfo[idxRow,'aFile'])
    one_maf <- download_study_fileKL(study_id, mfInfo[idxRow,'mFile'])
    
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
    
  }##the loop for one ion mode ends here
  
  #save forExport, out of R so I can gather up the pieces later
  #make a simple file name base
  #fbase <- paste0('export',study_id)
  #saveRDS(forExport, file = paste0(fbase,'.rds'))

} #the loop for unique MTBLS numbers (uM) ends here

#export the result
write.csv(forExport,file  = paste0('allHits','.csv'))

