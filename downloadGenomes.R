#### downloadGenomes.R 
# Andrea Martinez Vernon   OSS Lab
#  asmvernon@gmail.com
# 
# MODIFIED: 11/11/2016
#
############################################################################################################################################
##  Download NCBI FASTA files for complete reference and repressentative genomes
#  
# INPUT
#   Paths         -  list
#       $dirPath  -  path directory 
#       $NAS      -  NCBI_Bact_Genomes/ in OSS NAS, where files are to be saved 
# 
#   Generates files:
#     * output/downloadGenomes_LOG.txt
# 
############################################################################################################################################

downloadGenomes <- function(Paths,server=T){
  setwd(Paths$dirPath)
  Start_time <- Sys.time()
  library(R.utils,quietly = T,warn.conflicts = F) # gunzip
  useDate <- Sys.Date()
  cat("***\tDownloading genomes from NCBI\t***", fill = T)
###------------------------------------------------------------------------------------------------------------------------------------###
####      INITIATE LOG FILE ----
###------------------------------------------------------------------------------------------------------------------------------------###
  LOG_FILE <- "output/downloadGenomes_LOG.txt"
  write(x = paste("downloadGenomes.R last modified on:",date(),"\n"), file = LOG_FILE,append = T)
###------------------------------------------------------------------------------------------------------------------------------------###
####      SAVE FILE TO NAS ----
###------------------------------------------------------------------------------------------------------------------------------------###
  if(server==F){
    genomesPath <- paste(Paths$NAS,"genomes_",useDate,sep="")
    if(is.element(paste("genomes_",useDate,sep=""),dir(Paths$NAS))==F){ 
      dir.create(genomesPath)
    }
  }else{
    genomesPath <- paste("data/","genomes_",useDate,sep="")
    if(is.element(paste("genomes_",useDate,sep=""),dir('data'))==F){ 
      dir.create(genomesPath)
    }
  }
###------------------------------------------------------------------------------------------------------------------------------------###
####      DOWNLOAD ASSEMBLY FILE FROM NCBI ----
###------------------------------------------------------------------------------------------------------------------------------------###
  # Download README file
  URL_fileREADME     <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt"
  
  outFile_fileREADME <- paste(genomesPath,"/README_assembly_summary_",useDate,".txt",sep="")
  download.file(url = URL_fileREADME,destfile = outFile_fileREADME,quiet = T)
  cat("\tDownloaded README file", fill = T)
  
  
  # ASSEMBLY SUMMARY
  URL <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
  # URL <- "ftp://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt" # INCLUDES EUKARYOTES AS WELL
  assembly_file <- paste(genomesPath,'/assembly_summary_NCBI_',useDate,'.txt',sep="")
  download.file(URL,assembly_file,quiet = T)
  cat("\tDownloaded assembly_summary file", fill = T)
  
  assembly_summary <- read.delim(assembly_file,stringsAsFactors = F,header = F,comment.char = "#",sep="\t")
  names(assembly_summary)   <- c("assembly_accession", "bioproject","biosample","wgs_master","refseq_category","taxid","species_taxid",
                               "organism_name","infraspecific_name","isolate","version_status","assembly_level","release_type","genome_rep",
                               "seq_rel_date","asm_name","submitter","gbrs_paired_asm","paired_asm_comp","ftp_path","excluded_from_refseq")
  
#### NUMBERS ----
  totalGenomes            <- nrow(assembly_summary)
  totalrep.refGenomes     <- length(which(assembly_summary$refseq_category!="na"))
  totalrefGenomes         <- length(which(assembly_summary$refseq_category=="reference genome"))
  totalrepGenomes         <- length(which(assembly_summary$refseq_category=="representative genome"))
  
  totalrep.refCompleteGenomes <- length(intersect(which(assembly_summary$refseq_category!="na"),  
                                                  which(assembly_summary$assembly_level=="Complete Genome")))
  totalrefCompleteGenomes     <- length(intersect(which(assembly_summary$refseq_category=="reference genome"),  
                                                  which(assembly_summary$assembly_level=="Complete Genome")))
  totalrepCompleteGenomes     <- length(intersect(which(assembly_summary$refseq_category=="representative genome"),  
                                                  which(assembly_summary$assembly_level=="Complete Genome")))
  
  ## WRITE NUMBERS TO LOG FILE 
  write(x = "", file = LOG_FILE,append = T)
  write(x = paste("Total Genomes \t\t\t\t\t\t\t\t\t\t-\t",totalGenomes), file = LOG_FILE,append = T)
  write(x = paste("Total Reference & Representative Genomes \t\t\t-\t",totalrep.refGenomes), file = LOG_FILE,append = T)
  write(x = paste("Total Reference  Genomes \t\t\t\t\t\t\t-\t",totalrefGenomes), file = LOG_FILE,append = T)
  write(x = paste("Total Representative Genomes \t\t\t\t\t\t-\t",totalrepGenomes), file = LOG_FILE,append = T)
  write(x = "", file = LOG_FILE,append = T)
  write(x = paste("Total Reference & Representative Complete Genomes \t-\t",totalrep.refCompleteGenomes), file = LOG_FILE,append = T)
  write(x = paste("Total Reference Complete Genomes \t\t\t\t\t-\t",totalrefCompleteGenomes), file = LOG_FILE,append = T)
  write(x = paste("Total Representative Complete Genomes \t\t\t\t-\t",totalrepCompleteGenomes), file = LOG_FILE,append = T)
  write(x = "", file = LOG_FILE,append = T)
  
###------------------------------------------------------------------------------------------------------------------------------------###
####      MODIFY ASSEMBLY SUMMARY FILE ----
###------------------------------------------------------------------------------------------------------------------------------------###
  ## Exclude those that aren't reference or representative
  assembly_summary <- assembly_summary[which(assembly_summary$refseq_category!="na"),]
  
  ## Exclude those that aren't reference or representative AND that aren't Complete
  assembly_summary <- assembly_summary[which(assembly_summary$assembly_level=="Complete Genome"),]
  
  
  ###------------------------------------------------------------------------------------------------------------------------------------###
  ##      DOWNLOAD ASSEMBLY FILE FROM NCBI
  ###------------------------------------------------------------------------------------------------------------------------------------###
  
  # ADD MORE INFO
  ftps  <- as.data.frame(sapply(assembly_summary$ftp_path, strsplit,split="/",simplify = T))
  assembly_summary$genomeIdentifier     <-    t(ftps[nrow(ftps),])
  assembly_summary$fastaDownload        <- NA     # Logical indicator of whether fasta file was downloaded
  
  
  refSeqcat   <- gsub(" ","_",assembly_summary$refseq_category) # replace spaces with "_"
  
  # CREATE FOLDERS
  for(r in unique(refSeqcat)){
    # REF SEQ CAT FOLDER
    refSeqCat_folder <- paste(genomesPath,"/",gsub(" ","_",r),sep="")
    if(is.element(gsub(" ","_",r),dir(genomesPath))==F) dir.create(refSeqCat_folder)
    
    
    # FOLDER TO SAVE FASTA FILES
    fastaPath <- paste(genomesPath,"/",gsub(" ","_",r),"/fasta",sep="")
    if(is.element("fasta",dir(refSeqCat_folder))==F) dir.create(fastaPath)
  }
  
  # DOWNLAOD FASTA FILES
  notifyEvery <- 50
  cat("", fill = T)
  cat("\tStarting download...... ",format(Sys.time()), fill = T)
  cat("\t\t< status update every", notifyEvery, "downloads >", fill = T)
  cat("", fill = T)
  failed <- NULL

  for(i in 1:nrow(assembly_summary)){
  # for(i in 1:10){
    if(i==1||i%%notifyEvery==0||i==nrow(assembly_summary)) {
      cat("\tDownloading\t",i, "\tof\t", nrow(assembly_summary), "genomes\t......\tRun time:",format(difftime(Sys.time(),Start_time),digits=3), fill = T)
    }
    
    URL     <- paste(assembly_summary$ftp_path[i],"/",assembly_summary$genomeIdentifier[i],"_genomic.fna.gz",sep="") # point to fasta file in ftp folder 
    outFile <- paste(genomesPath,"/",refSeqcat[i],"/fasta","/",assembly_summary$genomeIdentifier[i],".fna.gz",sep="")
    
    # Download .gz file
    try(out <- download.file(url = URL,destfile = outFile,quiet = T)) # output other than 0 indicates error downloading 
    
    try(assembly_summary$fastaDownload[i]  <- out)
    
    # Extract file
    try(gunzip(outFile))
    Sys.sleep(0.5)
  } # genome download loop

  failed <- c(which(assembly_summary$fastaDownload!=0),which(is.na(assembly_summary$fastaDownload)))  
        # Identify entries that are non-zero, including NA
  
  # Retry those that failed
  if(length(failed)>0){
    cat("*** \tRetrying failed downloads\t***", sep="",fill = T)
    for(i in failed){
      cat("\tDownloading\t",i, "\tof\t", length(failed), "genomes\t......\tRun time:",format(difftime(Sys.time(),Start_time),digits=3), fill = T)
      
      URL     <- paste(assembly_summary$ftp_path[i],"/",assembly_summary$genomeIdentifier[i],"_genomic.fna.gz",sep="")
      outFile <- paste(genomesPath,"/",refSeqcat[i],"/fasta","/",assembly_summary$genomeIdentifier[i],".fna.gz",sep="")
      
        # Download .gz file
        try(out <- download.file(url = URL,destfile = outFile,quiet = T))

        try(assembly_summary$fastaDownload[i]  <- out)
      
        # Extract file
        try(gunzip(outFile))
      Sys.sleep(0.5)
    } # Retry failed 
  } # end if failed 
  
  ## Find the genome sequence that were not downloaded
  failed2 <- c(which(assembly_summary$fastaDownload!=0),which(is.na(assembly_summary$fastaDownload)))   # works with 0 and "0"  
  
  if(length(failed2)>0){
    # write(x = paste(assembly_summary$assembly_accession[failed_index],collapse = ","), file = LOG_FILE,append = T)   
    write(x = paste(assembly_summary$assembly_accession[failed2],collapse = ","), file = LOG_FILE,append = T)   
  }
  assembly_summary$fastaDownload <- ifelse(test = assembly_summary$fastaDownload==0,yes = "YES", no = "NO") # successfull download == 0
  assembly_summary <- assembly_summary[which(assembly_summary$fastaDownload=="YES"),]
  
  # DATA CLEAN UP & COMPATILBILITY ASSURANCE FOR CSV
  assembly_summary$organism_name      <- gsub(",","",assembly_summary$organism_name)
  assembly_summary$infraspecific_name <- gsub(",","",assembly_summary$infraspecific_name)
  for(n in 1:ncol(assembly_summary)) assembly_summary[,n] <- gsub(",",";",assembly_summary[,n]) # Make it csv compatible
  
  # ORDER ALPHABEtICALLY BY REF_SEQ CATEGORY
  assembly_summary <- assembly_summary[order(assembly_summary$genomeIdentifier),]
  assembly_summary <- assembly_summary[order(assembly_summary$refseq_category),]
  
  ## WRITE CSV FILE
  write.csv(assembly_summary,file = paste(genomesPath,'/assembly_summary_',useDate,'.csv',sep=""),row.names = F)
  
  ## WRITE TAB DELIMITED FILE
  write(paste(names(assembly_summary),collapse = "\t"),file = paste(genomesPath,'/assembly_summary_',useDate,'.txt',sep=""),append = F)
  write(paste(assembly_summary,collapse = "\t"),file = paste(genomesPath,'/assembly_summary_',useDate,'.txt',sep=""),append = T)
  
  write(x = "-----------------------------------------------------------------------------------------", file = LOG_FILE,append = T)   
}


