
## #####################################################################################
##
##  File: htt_code.R
##
##  Code for analyzing target gene repeat expansions in fastq files (reads of insert)
##
##  Author: Adam Ameur, SciLifeLab/Uppsala University
##
##  Dependencies: 'seqinr' R/Bioconductor package
##
## #####################################################################################


## #######################################################################
##
## Definition of targets to be analyzed (i.e. the HTT gene in this study)
##
## #######################################################################

targets <- list()
targets$HTT <- list()
targets$HTT$repeatElements <- c("CAG","CCG") # Repeats to be analyzed (CAG:main repeat, CCG: flanking repeat)
targets$HTT$start <- "AGTCCTTCCAGCAG" # 14bp overlapping repeat start (CAGCAG)
targets$HTT$end <- "CCGCCGCCTCCTCA" # 14bp overlapping repeat end (CCGCCG)
targets$HTT$basesBefore <- c(2641,905) ## Expected read length(s) before HTT repeat start site
targets$HTT$basesAfter <- c(114) ## Expected read length(s) before HTT repeat end site


## #############################################################################################################
##
##  Function: extractRepeatsInHTTsamples
##
##  This function executes the analysis of HTT samples. For each sample, a CCS fasta file is taken as input
##
##  The following 4 files are generated as output for each of the samples:
##
##  a) ontarget_uncorrected.fasta: complete repeat sequences extracted from the on-target reads
##
##  b) ontarget_corrected.fasta: error correction of sequences in ontarget.fasta
##
##  c) uncorrected_repeats.fasta: file where repeats in 'uncorrected' file are replaced by 1 (CAG) and 2 (CCG)
##
##  d) corrected_repeats.fasta: file where repeats in 'corrected' file are replaced by 1 (CAG) and 2 (CCG)
##
## ##############################################################################################################

extractRepeatsHTTsamples <- function(){
  require(seqinr)

  analyzeRepeatsInSample(fastaFile="CCS_reads/genscript_rep1_CCS.fasta", outdir="analysis_results", sampleName="genscript_rep1",targetList=c("HTT"))
  analyzeRepeatsInSample(fastaFile="CCS_reads/genscript_rep2_CCS.fasta", outdir="analysis_results", sampleName="genscript_rep2",targetList=c("HTT"))
  analyzeRepeatsInSample(fastaFile="CCS_reads/genscript_rep3_CCS.fasta", outdir="analysis_results", sampleName="genscript_rep3",targetList=c("HTT"))
  analyzeRepeatsInSample(fastaFile="CCS_reads/genscript_rep4_CCS.fasta", outdir="analysis_results", sampleName="genscript_rep4",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample2_CCS.fasta", outdir="analysis_results", sampleName="Sample2",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample2_CCS.fasta", outdir="analysis_results", sampleName="Sample2",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample2_CCS.fasta", outdir="analysis_results", sampleName="Sample2",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample3_CCS.fasta", outdir="analysis_results", sampleName="Sample3",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample4_CCS.fasta", outdir="analysis_results", sampleName="Sample4",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample5_CCS.fasta", outdir="analysis_results", sampleName="Sample5",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample6_CCS.fasta", outdir="analysis_results", sampleName="Sample6",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample7_CCS.fasta", outdir="analysis_results", sampleName="Sample7",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample8_CCS.fasta", outdir="analysis_results", sampleName="Sample8",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample9_CCS.fasta", outdir="analysis_results", sampleName="Sample9",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample10_CCS.fasta", outdir="analysis_results", sampleName="Sample10",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample11_CCS.fasta", outdir="analysis_results", sampleName="Sample11",targetList=c("HTT"))
  #analyzeRepeatsInSample(fastaFile="CCS_reads/Sample12_CCS.fasta", outdir="analysis_results", sampleName="Sample12",targetList=c("HTT"))

}


analyzeRepeatsInSampleAndTarget <- function(reads, target, outfile, allowIndelInRecognitionSites=TRUE, maxD=0.1){

    ## Goes through a vector of sequences and returns the first one matching in the read
    getMatchingSequenceInRead <- function(read, sequences){

        for(sequence in sequences){
            if(length(grep(sequence, read)) == 1){
                return(sequence)
            }
        }

        return(NA)
    }


    ## Trim a read and remove all bases not between 'targetStart' and 'targetEnd'.
    ## Make sure repeat sequence is not occurring in flanked sequences!
    trimRead <- function(read, targetStart, targetEnd, repeatElement=NA, maxRepeatsInFlanks=5){

        repeatSequence <- NA
        longestMatchingSeq <- ""
        seqBefore <- ""
        seqAfter <- ""

        if(!is.na(repeatElement)){
            repeatSequence <- paste(rep(repeatElement,maxRepeatsInFlanks+1),collapse="")
        }

        if(length(grep(targetStart, read)) == 1){
            tmpRes <- strsplit(read,targetStart)[[1]]
            beginningOfRead <- tmpRes[1]
            tmpRead <- tmpRes[2]

            if(length(grep(targetEnd, tmpRead)) == 1){
                tmpRes <- strsplit(tmpRead,targetEnd)[[1]]
                trimmedRead <- tmpRes[1]
                endOfRead <- tmpRes[2]
                if( (length(grep(repeatSequence, beginningOfRead)) == 0) && (length(grep(repeatSequence, endOfRead)) == 0)){
                    trimmedRead <- paste(targetStart, trimmedRead, targetEnd, sep="")
                    if(nchar(trimmedRead) > nchar(longestMatchingSeq)){
                        longestMatchingSeq <- trimmedRead
                        seqBefore <- beginningOfRead
                        seqAfter <- endOfRead
                    }
                }
            }
        }

        return(list("trimmedRead"=longestMatchingSeq,"seqBefore"=seqBefore,"seqAfter"=seqAfter))
    }


    ## Check if a read has correct length before/after the repeat recognition sites
    checkSizeOfFlanks <- function(observedLenBefore, observedLenAfter, targetLengthsBefore, targetLengthsAfter, maxDist, revComp=FALSE){

      sizeBeforeOk <- FALSE
      sizeAfterOk <- FALSE

      if(revComp == TRUE){
        tmp <- observedLenBefore
        observedLenBefore <- observedLenAfter
        observedLenAfter <- tmp

      }

      ##print(maxDist)
      ##print(observedLenBefore)

      for(len in targetLengthsBefore){
        if( (findInterval(observedLenBefore,c(len-(maxDist*len),len+(maxDist*len))) == 1) ){
          sizeBeforeOk <- TRUE
        }
      }
      for(len in targetLengthsAfter){
        if( (findInterval(observedLenAfter,c(len-(maxDist*len),len+(maxDist*len))) == 1) ){
          sizeAfterOk <- TRUE
        }
      }

      return((sizeBeforeOk && sizeAfterOk))

    }


    ## Main function start

    if(file.exists(outfile)){
        file.remove(outfile)
    }

    repeatNr <- list()
    trimmedSequence <- list()

    startsFound <- 0
    endsFound <- 0
    targetsFound <- 0

    targetRepeatElements <- targets[[target]]$repeatElements

    targetStarts <- targets[[target]]$start
    targetEnds <- targets[[target]]$end

    if(allowIndelInRecognitionSites==TRUE){
        targetStarts <- introduceTwoIndels(targetStarts)
        targetEnds <- introduceTwoIndels(targetEnds)
    }

    targetStartsRC <- sapply(targetStarts, revComp)
    targetEndsRC <- sapply(targetEnds, revComp)

    totalReads <- length(names(reads))

    currentReadNr <- 0

    for(readName in names(reads)){

        read <- as.character(reads[readName])

        currentReadNr <- currentReadNr+1

        cat("\r - Target: ",target," [processing read nr ",currentReadNr," of ",totalReads,"]", sep="")

        matchStart <- FALSE
        matchEnd <- FALSE
        matchStartRC <- FALSE
        matchEndRC <- FALSE

        matchStart <- getMatchingSequenceInRead(read, targetStarts)
        matchEnd <- getMatchingSequenceInRead(read, targetEnds)
        matchStartRC <- getMatchingSequenceInRead(read, targetStartsRC)
        matchEndRC <- getMatchingSequenceInRead(read, targetEndsRC)

        if(!is.na(matchStart)){
            startsFound <- startsFound+1
        }
        if(!is.na(matchEnd)){
            endsFound <- endsFound+1
        }
        if(!is.na(matchStartRC)){
            startsFound <- startsFound+1
        }
        if(!is.na(matchEndRC)){
            endsFound <- endsFound+1
        }

        if( !is.na(matchStart) && !is.na(matchEnd) ){

            res <- trimRead(read,matchStart,matchEnd,targetRepeatElements[1])
            trimmedRead <- res[["trimmedRead"]]
            seqBefore <- res[["seqBefore"]]
            seqAfter <- res[["seqAfter"]]

            if(nchar(trimmedRead)>0){

                flankingBasesOk <- checkSizeOfFlanks(nchar(seqBefore), nchar(seqAfter), targets[[target]]$basesBefore, targets[[target]]$basesAfter, maxD, revComp=FALSE)

                if(flankingBasesOk){
                    targetsFound <- targetsFound + 1
                    ##cat(targetsFound," F - ", matchStart," ",matchEnd,"\r",sep="")
                    trimmedSequence[[readName]] <- trimmedRead
                    cat(">",target,"|",paste(targetRepeatElements,collapse=","),"|fwd|b",nchar(seqBefore),"|a",nchar(seqAfter),"|",readName,"\n",sep="", file=outfile, append=TRUE)
                    cat(trimmedSequence[[readName]],"\n",sep="",file=outfile, append=TRUE)
                }
                ##else{
                ##    cat("Filtered:",trimmedRead,"-",readName,"\n")
                ##}
            }
        }
        if( !is.na(matchStartRC) && !is.na(matchEndRC) ){

            res <- trimRead(read,matchEndRC,matchStartRC,revComp(targetRepeatElements[1]))
            trimmedRead <- res[["trimmedRead"]]
            seqBefore <- res[["seqBefore"]]
            seqAfter <- res[["seqAfter"]]

            if(nchar(trimmedRead)>0){

                flankingBasesOk <- checkSizeOfFlanks(nchar(seqBefore), nchar(seqAfter), targets[[target]]$basesBefore, targets[[target]]$basesAfter, maxD, revComp=TRUE)

                if(flankingBasesOk){
                    targetsFound <- targetsFound + 1
                    ##cat(targetsFound," R - ",matchStartRC," ",matchEndRC,"\r",sep="")
                    trimmedSequence[[readName]] <- trimmedRead
                    targetRepeatElementsRC <- sapply(targetRepeatElements, revComp)
                    cat(">",target,"|",paste(targetRepeatElementsRC,collapse=","),"|rev|b",nchar(seqBefore),"|a",nchar(seqAfter),"|",readName,"\n",sep="", file=outfile, append=TRUE)
                    cat(trimmedSequence[[readName]],"\n",sep="",file=outfile, append=TRUE)
                }
                ##else{
                ##   cat("FilteredREV:",trimmedRead,"-",readName,"\n")
                ##}
            }
        }
    }
    cat(" - Done.\n")
}



analyzeRepeatsInSample <- function(fastaFile="coco_coriell/095082_A01_NA03620_CCS.fasta", sampleName="NA03620", outdir="results", targetList = c("HTT"), allowIndelInRecognitionSites=TRUE, rerun=FALSE){

    cat("Analyzing sample: ",sampleName,"\n",sep="")

    results <- list()

    ## Create on-target fasta files (if they don't already exist)
    for(target in targetList){

        onTargetFastaFile <- paste(outdir,"/",sampleName,"_ontarget_uncorrected.fasta",sep="")

        if( !file.exists(onTargetFastaFile) || (rerun==TRUE)){
            require(seqinr)
            fastaReads <- read.fasta(fastaFile, as.string=TRUE)
            reads <- toupper(as.character(fastaReads))
            names(reads) <- names(fastaReads)
            analyzeRepeatsInSampleAndTarget(reads, target, onTargetFastaFile)
        }
    }


    results <- list()

    ## Go through fasta files and do error correction
    for(target in targetList){

        results$target <- list()

        onTargetFastaFile <- paste(outdir,"/",sampleName,"_ontarget_uncorrected.fasta",sep="")

        if(!file.exists(onTargetFastaFile)){
            file.create(onTargetFastaFile)
        }

        repeatFileUncorrected <- paste(outdir,"/",sampleName,"_ontarget_uncorrected_repeats.fasta",sep="")
        if(file.exists(repeatFileUncorrected)){
            file.remove(repeatFileUncorrected)
        }

        onTargetFastaFileCorrected <- paste(outdir,"/",sampleName,"_ontarget_corrected.fasta",sep="")
        if(file.exists(onTargetFastaFileCorrected)){
            file.remove(onTargetFastaFileCorrected)
        }

        repeatFileCorrected <- paste(outdir,"/",sampleName,"_ontarget_corrected_repeats.fasta",sep="")
        if(file.exists(repeatFileCorrected)){
            file.remove(repeatFileCorrected)
        }


        if(file.info(onTargetFastaFile)$size == 0){
            results[[target]]$nrReadsOnTarget <- 0
            results[[target]]$nrRepeatsRaw <- 0
            results[[target]]$nrRepeatsCor <- 0
        }
        else{
            data <- read.fasta(onTargetFastaFile, as.string=TRUE)
            OTreads <- toupper(as.character(data))
            names(OTreads) <- names(data)
            nrRepeats <- array(NA,0)
            nrRepeatsCorrected <- array(NA,0)

            for(i in 1:length(OTreads)){

                header <- as.character(names(OTreads)[i])
                sequence <- as.character(OTreads[header])

                headerVec <- strsplit(header,"\\|")[[1]]

                repeatSequence <- headerVec[2]
                strand <- headerVec[3]

                repeatElements <- strsplit(repeatSequence,",")[[1]]

                ## Replace repeat in uncorrected sequence
                sequenceReplaced <- sequence

                for(i in 1:length(repeatElements)){
                    sequenceReplaced <-  replaceRepeat(sequenceReplaced, repeatElements[i], i)
                }

                repeatSequenceUncorrected <- sequenceReplaced

                ## Replace repeat in corrected sequence
                sequenceCorrected <- sequence

                for(i in 1:length(repeatElements)){
                    sequenceCorrected <- doReadCorrection(sequenceCorrected, repeatElements[i], nrounds=5)
                }

                sequenceReplaced <- sequenceCorrected

                for(i in 1:length(repeatElements)){
                    sequenceReplaced <-  replaceRepeat(sequenceReplaced, repeatElements[i], i)
                }

                repeatSequenceCorrected <- sequenceReplaced

                cat(">",header,"\n",sep="", file=onTargetFastaFileCorrected, append=TRUE)
                cat(sequenceCorrected,"\n",sep="",file=onTargetFastaFileCorrected, append=TRUE)

                cat(">",header,"\n",sep="", file=repeatFileUncorrected, append=TRUE)
                cat(repeatSequenceUncorrected,"\n",sep="",file=repeatFileUncorrected, append=TRUE)

                cat(">",header,"\n",sep="", file=repeatFileCorrected, append=TRUE)
                cat(repeatSequenceCorrected,"\n",sep="",file=repeatFileCorrected, append=TRUE)
            }
        }
    }
}



## Computes the longest repeat stretch in a DNA sequence, with most 'maxGapBases' are interrupting the repeat
replaceRepeat <- function(seq, repeatElement, replaceElement="X"){

    replacedSeq <- gsub(repeatElement,replaceElement,seq)

    return(replacedSeq)

}


doReadCorrection <- function(read, repeatElement, nrounds=5){

    correctedRead <- read

    isDone <- FALSE
    iterations <- 1

    while(!isDone && iterations<=nrounds){

        iterations <- iterations+1
        indelFound <- FALSE
        indels <- generateRepeatIndels(repeatElement)

        for(i in 1:length(indels)){
            replacement <- names(indels)[i]
            indel <- indels[i]
            if(length(grep(indel, correctedRead))>0){
                indelFound <- TRUE
                correctedRead <- gsub(indel, replacement, correctedRead)
            }
        }

        if(!indelFound){
            isDone <- TRUE
        }
    }

    return(correctedRead)
}


generateRepeatIndels <- function(repeatElement, nrRepsBefore=1, nrRepsAfter=1){

   repeatVec <- strsplit(repeatElement,"")[[1]]
   repeatBefore <- paste(rep(repeatElement, nrRepsBefore), collapse="")
   repeatAfter <- paste(rep(repeatElement, nrRepsAfter), collapse="")

   deletions <- array(NA,0)
   deletionNames <- array(NA,0)
   insertions <- array(NA,0)
   insertionNames <- array(NA,0)

   for(i in 1:length(repeatVec)){
       #deletionNames[i] <- paste(repeatElement,"D",i,sep="")
       deletionNames[i] <- repeatElement
       deletions[i] <- paste(repeatVec[-i],collapse="")
       #insertionNames[i] <- paste(repeatElement,"I",i,sep="")
       insertionNames[i] <- repeatElement
       if(i < length(repeatVec)){
           insertions[i] <- paste(c(repeatVec[1:i],repeatVec[i],repeatVec[(i+1):length(repeatVec)]),collapse="")
       }
       if(i == length(repeatVec)){
           insertions[i] <- paste(c(repeatVec,repeatVec[length(repeatVec)]),collapse="")
       }
   }

   names(deletions) <- deletionNames
   names(insertions) <- insertionNames

   indels <- c(insertions,deletions)

   doubleIndels <- array(NA,0)
   doubleIndelNames <- array(NA,0)

   for(i in 1:length(indels)){
       for(j in 1:length(indels)){
           doubleIndel <- paste(indels[i],indels[j],sep="")
           doubleIndelName <- paste(names(indels)[i],names(indels)[j],sep="")

           doubleIndels <- c(doubleIndels, doubleIndel)
           doubleIndelNames <- c(doubleIndelNames, doubleIndelName)

       }
   }

   names(doubleIndels) <- doubleIndelNames

   indels <- c(indels,doubleIndels)

   indelsFinal <- paste(repeatBefore,indels,repeatAfter,sep="")
   indelNamesFinal <- paste(repeatBefore,names(indels),repeatAfter,sep="")
   names(indelsFinal) <- indelNamesFinal

   #indels <- paste(repeatBefore,indels,repeatAfter,sep="")
   #indelNames <- paste(repeatBefore,names(indels),repeatAfter,sep="")
   #names(indels) <- indelNames

   return(indelsFinal)
}



introduceIndels <- function(sequence){

    sequenceVec <- strsplit(sequence,"")[[1]]

    allSequences <- c(sequence)

    for(pos in 1:length(sequenceVec)){
        deletionSeq <- paste(sequenceVec[-pos],collapse="")
        allSequences <- c(allSequences, deletionSeq)
        if(pos < length(sequenceVec)){
            insertionSeq <- paste(c(sequenceVec[1:pos],sequenceVec[pos],sequenceVec[(pos+1):length(sequenceVec)]),collapse="")
            allSequences <- c(allSequences, insertionSeq)
        }

    }

    return(allSequences)
}


introduceTwoIndels <- function(sequence){

    sequencesRound1 <- introduceIndels(sequence)

    sequencesRound2 <- array(NA,0)

    for(seq in sequencesRound1){
        tmpSeqs <- introduceIndels(seq)

        sequencesRound2 <- c(sequencesRound2, tmpSeqs)
    }

    finalSequences <- unique(c(sequencesRound1,sequencesRound2))

    return(finalSequences)


}


## Returns the reverse complement of a IUPAC or consensus sequence
revComp <- function(seq){

    baseComp <- array(NA,0)
    baseComp[c("A","C","G","T","R","Y","K","M","S","W","B","V","D","H","N","[","]","(",")")] <- c("T","G","C","A","Y","R","M","K","S","W","V","B","H","D","N","]","[",")","(")

    return(paste(rev(sapply(strsplit(toupper(seq),"")[[1]],function(x){
        if(!is.na(baseComp[x]))
        {
            return(baseComp[x])
        }else{
            return(x)
        }
    },USE.NAMES = FALSE)),collapse = ""))
}



## Computes the longest repeat stretch in a DNA sequence, with most 'maxGapBases' are interrupting the repeat
computeNrRepeatsInSeq <- function(seq, repeatElement="1", maxGapBases=3, minRepsAtStartAndEnd=2){


    stretchStartingAtPos <- function(baseVec, pos, baseToSearchFor, maxGapBases){

        i <- pos

        subVec <- baseVec[i:length(baseVec)]

        ##print(subVec)

        mmVec <- cumsum(!(subVec == baseToSearchFor))

        if(mmVec[1]<=maxGapBases){
            ##print(subVec)
            endpos <- max(which(mmVec <= maxGapBases))
            return(baseVec[i:(i+endpos-1)])
        }
        else{
            return(array(NA,0))
        }

    }


    seqvec <- strsplit(seq,"")[[1]]

    longestStretch <- array(NA,0)

    startRepeat <- paste(rep(repeatElement,minRepsAtStartAndEnd),collapse="")

    for(i in 1:length(seqvec)){
        if(seqvec[i] == repeatElement){
            tmpvec <- stretchStartingAtPos(seqvec, i, repeatElement, maxGapBases)

            tmpstr <- paste(tmpvec,collapse="")

            ## make sure string starts and ends with 'startRepeat'
            basesToRemoveStart <- strsplit(tmpstr,startRepeat)[[1]][1]
            tmpstr <- sub(basesToRemoveStart,"",tmpstr)

            tmpstrRev <- paste(rev(strsplit(tmpstr,"")[[1]]),collapse="")
            basesToRemoveEnd <- strsplit(tmpstrRev,startRepeat)[[1]][1]
            tmpstrRev <- sub(basesToRemoveEnd,"",tmpstrRev)

            tmpvec <- rev(strsplit(tmpstrRev,"")[[1]])

            if(length(tmpvec) > length(longestStretch)){
                longestStretch <- tmpvec
            }
        }
    }

    nrRepeats <- length(which(longestStretch == repeatElement))

    return(nrRepeats)

}


computeNrRepeatsInSample <- function(fastaFile, repeatElements=c("1","2")){

    require(seqinr)

    fastaEntries <- read.fasta(fastaFile, as.string=TRUE)
    seqs <- toupper(as.character(fastaEntries))
    names(seqs) <- names(fastaEntries)

    repeatResults <- list()
    sequences <- list()


    for(repeatElement in repeatElements){
        repeatResults[[repeatElement]] <- NULL
    }


    for(seqName in names(seqs)){
        seq <- seqs[seqName]
        for(repeatElement in repeatElements){
            tmpRepeatNr <- computeNrRepeatsInSeq(seq, repeatElement)
            repeatResults[[repeatElement]] <- c(repeatResults[[repeatElement]], tmpRepeatNr)
        }

        if(length(grep("rev",seqName)) == 1){
            seq <- revComp(seq)
        }

        sequences[[seqName]] <- as.character(seq)
    }

    return(list(sequences,repeatResults))
}



plotRepeatsInSample <- function(fastaFile, repeatElements=c("1","2"), sampleName=""){

    res <- computeNrRepeatsInSample(fastaFile, repeatElements)

    sequences <- res[[1]]
    seqNames <- names(sequences)
    repeatCounts <- res[[2]]

    CAGrepeatCountVec <- repeatCounts[[1]]
    idsCluster1 <- order(CAGrepeatCountVec)[1:floor(length(res[[1]])/2)]
    idsCluster2 <- order(CAGrepeatCountVec)[ceiling(length(res[[1]])/2):length(CAGrepeatCountVec)]
    CAGcounts1 <- repeatCounts[[1]][idsCluster1]
    CAGcounts2 <- repeatCounts[[1]][idsCluster2]
    CCGcounts1 <- repeatCounts[[2]][idsCluster1]
    CCGcounts2 <- repeatCounts[[2]][idsCluster2]

    CAGcount1 <- names(which(table(CAGcounts1) == max(table(CAGcounts1))))
    CAGcount2 <- names(which(table(CAGcounts2) == max(table(CAGcounts2))))
    CCGcount1 <- names(which(table(CCGcounts1) == max(table(CCGcounts1))))
    CCGcount2 <- names(which(table(CCGcounts2) == max(table(CCGcounts2))))

    repeatString <- paste("Allele1: ",CAGcount1,"xCAG, ",CCGcount1,"xCCG; Allele2: ",CAGcount2,"xCAG, ",CCGcount2,"xCCG",sep="")

    #if(CAGcount2 < CAGcount1){
    #    repeatString <- paste("Allele1: ",CAGcount2,"xCAG, ",CCGcount2,"xCCG; Allele2: ",CAGcount1,"xCAG, ",CCGcount1,"xCCG",sep="")
    #}

    sequences <- sapply(names(sequences), function(x){gsub("1","345",sequences[[x]])})
    sequences <- sapply(names(sequences), function(x){gsub("2","678",sequences[[x]])})

    seqOrder <- order(sapply(sequences, nchar))

    dataVec <- 1:10
    names(dataVec) <- c("A","C","G","T","3","4","5","6","7","8")

    colVec <- c("white",grey(0.2),grey(0.4),grey(0.6),grey(0.8),"red1","red2","red3","blue1","blue2","blue3")

    maxLen <- max(sapply(sequences,nchar))

    dataMatrix <- matrix(0, nrow=length(sequences), ncol=maxLen)

    for(i in 1:length(seqNames)){
        seqName <- seqNames[seqOrder][i]
        seq <- strsplit(as.character(sequences[seqName]),"")[[1]]
        dataMatrix[i,1:length(seq)] <- as.numeric(dataVec[seq])
    }


    pdfFile <- sub(".fasta",".pdf",fastaFile)

    pdf(pdfFile, height=6, width=14)

    par(mfrow=c(1,2))

    plot(table(repeatCounts[[1]]), type="h", bty="n", main=paste("CAG repeat distribution, ",length(repeatCounts[[1]])," reads on target",sep=""), ylab="nr reads", xlab="repeat count", lwd=5, las=2)

    image(x=1:maxLen,y=1:length(seqNames),t(dataMatrix), col=colVec, breaks=c(0,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10), xlab="bases",ylab="reads", bty="n", axes=FALSE, main=repeatString)

    axis(1, line=NA)
    axis(2, line=NA)

    dev.off()

}


## #########################################################################################################
##
##  Function: plotRepeatsInHTTsamples
##
##  This function generates result summary plots for the HTT samples, both 'corrected' and 'uncorrected'
##
## #########################################################################################################

plotRepeatsInHTTamples <- function(){
    plotRepeatsInSample("analysis_results/genscript_rep1_ontarget_corrected_repeats.fasta")
    plotRepeatsInSample("analysis_results/genscript_rep1_ontarget_uncorrected_repeats.fasta")
    plotRepeatsInSample("analysis_results/genscript_rep2_ontarget_corrected_repeats.fasta")
    plotRepeatsInSample("analysis_results/genscript_rep2_ontarget_uncorrected_repeats.fasta")
    plotRepeatsInSample("analysis_results/genscript_rep3_ontarget_corrected_repeats.fasta")
    plotRepeatsInSample("analysis_results/genscript_rep3_ontarget_uncorrected_repeats.fasta")
    plotRepeatsInSample("analysis_results/genscript_rep4_ontarget_corrected_repeats.fasta")
    plotRepeatsInSample("analysis_results/genscript_rep4_ontarget_uncorrected_repeats.fasta")
}
