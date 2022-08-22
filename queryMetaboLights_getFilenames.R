#this fill pulls the file names for each MTBLS number so it can be hacked into 
#the lookup table
#Krista Longnecker, 19 August 2022
#Woods Hole Oceanographic Institution
library(metabolighteR)

#tidy up first
rm(list = ls.str())

#apparently trace does not persist across sessions, so define my own function 
#with the desired edit so I don't get the warning about encoding
download_study_fileKL <- function (study_id, filename) 
{
  file_download <- httr::GET(paste0(getOption("BASE_URL"), 
                     "/studies/", study_id, "/download/public?file=", filename))
  
  file_content <- httr::content(file_download, "text", encoding = "UTF-8") %>% 
    textConnection() %>% readLines() %>% utils::read.delim(text = ., 
                        sep = "\t") %>% tibble::as_tibble()
  return(file_content)
}

#loop through different studies; this is a great example of changes in 
#repositories over time!

allKuj = c('MTBLS144','MTBLS154','MTBLS366','MTBLS1024','MTBLS1820',
           'MTBLS461','MTBLS293','MTBLS155','MTBLS295','MTBLS157')

#quick details on Kuj studies
#'MTBLS144' #Tps4  
#'MTBLS154'#Tps6 (example of targeted and untargeted)
#'MTBLS366' #Lohmann cruise
#'MTBLS1024' #DeepDOM DOS project (direct infusion)
#'MTBLS1820' #Jardines de la Reina
#'MTBLS461' #dilution experiment neg/untargeted, separate file with targeted
#'MTBLS293' #NB pos/untar, neg/untar, targeted
#'MTBLS155' #Syn2 pos/untar, neg/untar, targeted
#'MTBLS295' #Cara's MP1
#'MTBLS157' #Rpom

for (a in 1:length(allKuj))
{
  study_id <- allKuj[a]
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
  
  track <- data.frame(matrix(NA, nrow = length(idxA), ncol = 4))
  colnames(track) <- c('MTBLS','ionMode','aFile','mafFile')
  
  for (a in idxA)
  {
    #put in the filename first...that's the easy step
    oneFile <- fileList_a[idxA[a]]
    track[a,'MTBLS'] <- study_id
    track[a,'aFile'] <- oneFile
    #now go into that file and pull the ion mode
    getFile <- download_study_fileKL(study_id, fileList_a[idxA[a]])
    #this will rely on single ion mode in the data file...could be a problem
    track[a,'ionMode'] <- unique(getFile$Parameter.Value.Scan.polarity.)
    
    #and see if I have a maf file
    track[a,'mafFile'] <- unique(getFile$Metabolite.Assignment.File)
  }
  
  #now I know both the name of the maf files and which a_ have a maf file
  #probably easier to just delete rows with no maf files
  track <- na.omit(track)
  
  #save forExport, out of R so I can gather up the pieces later
  #make a simple file name base
  fbase <- paste0('exportTrack_',study_id)
  
  #saveRDS(forExport, file = paste0(fbase,'.rds'))
  write.csv(track,file  = paste0(fbase,'.csv'))
} #end looping through all Kuj files here

#could (but don't at the moment) read these in and concatenate into the
#look up table. However, I know there are files in here that are not
#untargeted data...leave for now