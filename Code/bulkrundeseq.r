########################################
#
# Program to run DESeq2 on the count data
#
########################################

setwd("~/jen/Files/")
require(DESeq2)
require(plyr)
require(biomaRt)

##### Part 1: Read in the data

alldirectories = c("run60msrnaseqwtsplnlungbbreg","run61msrnaseqdkospbbregwtlnbreg","run68mousernaseqtregsplnlung","run69mousernaseqtregs")
countdirectory="/home/trevor/jen/counts/"
summarytabledirectory="/home/trevor/jen/summarytables/"

# Count data
curfilesc = NULL
filenamesc = NULL
curtablefilec = vector("list",length(alldirectories))
removesummaryinfo <- function(data) {
	newdata = data[-which(substr(data$V1,1,2)=="__"),]
return(newdata)
}
for(i in 1:length(alldirectories)) {
	curdirectory = alldirectories[i]
	filenames = list.files(paste(countdirectory,curdirectory,"/",sep=""), pattern="*.txt", full.names=FALSE)
	curfiles = lapply(paste(countdirectory,curdirectory,"/",filenames,sep=""), read.table)
	curfilesc = c(curfilesc,lapply(curfiles,removesummaryinfo))
	filenamesc = c(filenamesc,paste(substr(curdirectory,1,5),filenames,sep=""))
	# Meta data
	curtablefilec[[i]] = read.table(paste(summarytabledirectory,curdirectory,"table.txt",sep=""), skip=3, sep="\t", header=TRUE)
	curtablefilec[[i]]$Sample.IDmod = paste(substr(curdirectory,1,5),curtablefilec[[i]][,3],sep="")
}
# Comparisons data
fc = file(paste(summarytabledirectory,"mousecomparisonswreps.txt",sep=""))
comparisonsfile = strsplit(readLines(fc), "\t")
close(fc)

##### Part 2: Format for DESeq analysis

countData = join_all(curfilesc,by="V1")
names(curtablefilec[[3]]) = names(curtablefilec[[1]])
allcolData = do.call(rbind,curtablefilec)
colnames(countData) = c("genes",sapply(strsplit(filenamesc,split="\\."),"[",1))
rownames(countData) = countData[,1]; countData2 = countData[,-1]
# Remove samples that were 0 in all runs
countData3 = countData2[-which(apply(countData2,1,sum)==0),]

##### Part 3: Run differential expression analysis on each comparison

allresults = vector("list",length(comparisonsfile))
for(i in 1:length(comparisonsfile)) {
	print(i)
	curcompsi = comparisonsfile[[i]]
	curcomps = NULL
	cursets = NULL
	for(j in 1:length(curcompsi)) {
		pullcomps = unlist(strsplit(curcompsi[[j]],split=" "))
		curcomps = c(curcomps,pullcomps)
		cursets = c(cursets,rep(j,length(pullcomps)))
	}
	pullcounts = which(colnames(countData3)%in%curcomps)
	pullcols = which(allcolData$Sample.IDmod%in%curcomps)
	curcolData = data.frame(condition=allcolData$Sample.Description[pullcols],group=as.factor(cursets))
	rownames(curcolData) = allcolData$Sample.IDmod[pullcols]
	curcountdata = countData3[,pullcounts]
	curcountdataorder = curcountdata[,match(colnames(curcountdata),rownames(curcolData))]
	if(length(unique(curcolData$group))>1) {
	dds = DESeqDataSetFromMatrix(countData = curcountdataorder, colData = curcolData, design = ~ group)
	}
	if(length(unique(curcolData$group))==1) { # Fall back on no replicates method if there are no replicates
	dds = DESeqDataSetFromMatrix(countData = curcountdataorder, colData = curcolData, design = ~ condition)
	}
	dds = DESeq(dds)
	res = results(dds)
	ressort = res[order(res$padj,decreasing=FALSE),]
	allresults[[i]] = ressort
}

##### Part 4: Annotate the results

for(i in 1:length(allresults)) {
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
annodata = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"),filters="ensembl_gene_id",values=allresults[[i]]@rownames,mart=ensembl)
allresults[[i]]@listData = cbind(DataFrame(allresults[[i]]),annodata[match(allresults[[i]]@rownames,annodata$ensembl_gene_id),-1])
}

##### Part 5: Save the data

save(allresults,file="mouseanalysisresuls.rdata")
