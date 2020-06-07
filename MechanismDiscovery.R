#install.packages("devtools")
#install.packages("rlang")
#install.packages("igraph")
#devtools::install_github("funnell/reactomefi")

REVEALERInput="REVEALER LIHC HCV.txt"
chromosomeLoc="mart_export.txt"
CPMinput="LIHC HCV Virus Only New.txt"
geneList="LIHC HCV Combined.txt"
cancer="HCV"

setwd("C:/Users/Harrison/Documents/PanCancerVirus/Step5-REVEALERmechanism")
library(reactomefi)
library(igraph)

REVEALERmech<-function(REVEALERInput,chromosomeLoc,CPMinput,geneList,cancer){ #direction 
  REVEALER<-read.delim(REVEALERInput, header=TRUE)
  CPM<-read.delim(CPMinput, header=FALSE)
  match<-read.delim(chromosomeLoc, header=FALSE)
  IAandSurr<-read.delim(geneList, header=FALSE)
  comb34<-paste(match[,3],match[,4],sep="") #combines 3rd and 4th column
  outMat<-matrix(,nrow=1,ncol=6) 
  for (j in 1:ncol(REVEALER)){
    comb<-matrix(ncol=3)
    REVEALERregionConv<-matrix(ncol=2)
    for (i in 1:nrow(REVEALER)){
      title<-strsplit(as.character(REVEALER[i,j]),"_")
      if(!is.na(title[[1]][2])){
        if(as.character(title[[1]][2])=="DEL"){
          eachcolumn=title[[1]][1]
          ind<-which(as.character(match[,5])==eachcolumn)
          locus<-comb34[ind[1]]
          ind2<-which(as.character(comb34)==locus)
          CytBandgenes<-unique(match[ind2,5])   #names of genes in the genomic region (cytoband)
          IAgene<-colnames(REVEALER)[j] #pull gene name in column header
          split<-strsplit(IAgene,"_")
          IAgeneName<-split[[1]][1]
          direction<-split[[1]][2]
          indIAgene<-which(as.character(CPM[,1])==IAgeneName)
          exprIAgene<-CPM[indIAgene,-1]
          mat<-matrix(,nrow=length(CytBandgenes),ncol=3)
          n=1
          for (k in CytBandgenes){
            indGene<-which(as.character(CPM[,1])==as.character(k))
            if(length(indGene)!=0){
              exprCytbandGene<-CPM[indGene[1],-1]
              cor <- cor.test(as.numeric(as.matrix(exprIAgene)),as.numeric(as.matrix(exprCytbandGene)),method="spearman",use="complete.obs")
              mat[n,1]<-k
              mat[n,2]<-cor$estimate
              mat[n,3]<-cor$p.value
              n=n+1
            }
          }
          if(direction=="NEG"){
            indMatching<-which(mat[,2]>0)
            matching<-mat[indMatching,]
            comb=rbind(comb,matching)
          } else if(direction=="POS"){
            indMatching<-which(mat[,2]<0)
            matching<-mat[indMatching,]
            comb=rbind(comb,matching)
          }
          recordRegion<-c(eachcolumn,locus)
          REVEALERregionConv<-rbind(REVEALERregionConv,recordRegion)
        } else if(title[[1]][2]=="AMP"){
          eachcolumn=title[[1]][1]
          ind<-which(as.character(match[,5])==eachcolumn)
          locus<-comb34[ind[1]]
          ind2<-which(as.character(comb34)==locus)
          CytBandgenes<-unique(match[ind2,5])   #names of genes in the genomic region (cytoband)
          IAgene<-colnames(REVEALER)[j] #pull gene name in column header
          split<-strsplit(IAgene,"_")
          IAgeneName<-split[[1]][1]
          direction<-split[[1]][2]
          indIAgene<-which(as.character(CPM[,1])==IAgeneName)
          exprIAgene<-CPM[indIAgene,-1]
          mat<-matrix(,nrow=length(CytBandgenes),ncol=3)
          n=1
          for (k in CytBandgenes){
            indGene<-which(as.character(CPM[,1])==as.character(k))
            if(length(indGene)!=0){
              exprCytbandGene<-CPM[indGene[1],-1]
              cor <- cor.test(as.numeric(as.matrix(exprIAgene)),as.numeric(as.matrix(exprCytbandGene)),method="spearman",use="complete.obs")
              mat[n,1]<-k
              mat[n,2]<-cor$estimate
              mat[n,3]<-cor$p.value
              n=n+1
            }
          }
          if(direction=="NEG"){
            indMatching<-which(mat[,2]<0)
            matching<-mat[indMatching,]
            comb=rbind(comb,matching)
          } else if(direction=="POS"){
            indMatching<-which(mat[,2]>0)
            matching<-mat[indMatching,]
            comb=rbind(comb,matching)
          }
          recordRegion<-c(eachcolumn,locus)
          REVEALERregionConv<-rbind(REVEALERregionConv,recordRegion)
        } else if(title[[1]][2]=="MUT"){
          mutGeneName=title[[1]][1]
          indMutGene<-which(as.character(CPM[,1])==mutGeneName)
          if(length(indMutGene)!=0){
            exprMutgene<-CPM[indMutGene,-1]
            IAgene<-colnames(REVEALER)[j] #pull gene name in column header
            split<-strsplit(IAgene,"_")
            IAgeneName<-split[[1]][1]
            direction<-split[[1]][2]
            indIAgene<-which(as.character(CPM[,1])==IAgeneName)
            exprIAgene<-CPM[indIAgene,-1]
            cor <- cor.test(as.numeric(as.matrix(exprIAgene)),as.numeric(as.matrix(exprMutgene)),method="spearman",use="complete.obs")
            mat<-matrix(,nrow=1,ncol=3)
            mat[1,1]<-mutGeneName
            mat[1,2]<-cor$estimate
            mat[1,3]<-cor$p.value
            comb=rbind(comb,mat)
          }
        } 
      }
    } #end for loop of i (REVEALER features)
    
    indSignificant<-which(as.numeric(comb[,3])<0.05)
    combSig<-comb[indSignificant,]
    combSigNum<-matrix(as.numeric(combSig),nrow=nrow(combSig),ncol=ncol(combSig)) #first column becomes NA due to words being converted to numbers
    combSigNum[,1]<-1:nrow(combSigNum) #fills first row with indices of combSig
    sortedComb<-combSigNum[order(combSigNum[,3]),]  #order numeric matrix by p-value ascending
    small<-min(10,nrow(sortedComb))
    top10ind<-sortedComb[1:small,1] #top 10 rows are best genes
    top10names<-combSig[top10ind,1] #take gene names
    
    #put into Cytoscape and find shortest path to IA gene
    genes<-rbind(as.matrix(IAandSurr),as.matrix(top10names)) #input all gene names including IA gene+ surrounding genes+ top 10 genes
    genesChar<-as.character(genes)
    FInetwork<-ReactomeFINetwork("2013",genesChar,use.linkers = TRUE)
    edges<-as.matrix(FInetwork@fis) #extract list of edges
    
    #loop through each of the top 10 genes and find the 1 gene with shortest path to IA gene
    igraph<-graph_from_edgelist(edges,directed=FALSE) #create igraph object from list of edges
    mat2<-matrix(,nrow=length(top10names),ncol=1) #create matrix to store distances
    m=1
    for (gene in top10names){
      try(mat2[m,1]<-distances(igraph,v=gene, to=IAgeneName))
      m=m+1
    }
    mat2[is.na(mat2)]=Inf  #set NAs in mat2 to infinity
    if(min(mat2)!=Inf){  #to filter out IA genes where no path could be found
      bestInd<-min(which(mat2==min(mat2))) #minimum distance ones are best
      bestGene<-top10names[bestInd]
      path<-shortest_paths(igraph,from=bestGene,to=IAgeneName) #find shortest path
      pathVec<-row.names(as.matrix(path$vpath[[1]])) #convert path to vector format
      pathString<-paste(pathVec,collapse = ",")
      
      #store results in matrix
      recordDat<-matrix(,nrow=1,ncol=6)
      recordDat[,1]<-IAgeneName
      recordDat[,2]<-bestGene
      ind<-which(as.character(match[,5])==bestGene)
      locus2<-comb34[ind[1]]
      recordDat[,3]<-locus2
      ind<-which(as.character(REVEALERregionConv[,2])==locus2)
      if(length(ind)!=0){
        recordDat[,4]<-REVEALERregionConv[ind[1],1]
      } else {
        recordDat[,4]<-"MUT"
      }
      recordDat[,5]<-sortedComb[bestInd,2]
      recordDat[,6]<-pathString
      outMat<-rbind(outMat,recordDat)
    }
  } #end for loop of j (IA gene)
  outMat[1,]<-c("IA Gene Name","Best Genome Alter","Genomic Region","Name in REVEALER","R^2","Path to IA Gene")
  write.table(outMat,paste("REVEALERmechanismOutput-",cancer,".txt", sep = ""),quote=FALSE,sep="\t",row.names = FALSE,col.names = FALSE)
}

REVEALERmech("REVEALER-CESC.txt","mart_export.txt","CESCVirus.txt","CESCCombined.txt","CESC")
REVEALERmech("REVEALER HNSCC.txt","mart_export.txt","HNSCC CPM Virus New.txt","HNSCC Combined.txt","HNSCC")
REVEALERmech("REVEALER LIHC HBV.txt","mart_export.txt","LIHC HBV Virus Only.txt","LIHC HBV Combined.txt","HBV")
REVEALERmech("REVEALER LIHC HCV.txt","mart_export.txt","LIHC HCV Virus Only New.txt","LIHC HCV Combined.txt","HCV")
REVEALERmech("REVEALER STAD.txt","mart_export.txt","STAD only Virus New.txt","STAD Combined.txt","STAD")
