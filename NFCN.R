##----------------------------------------------------------------------------##
## By D. H. Heiland 
## MILO Laboratory, themilolab.com
## 
## Neares functional connected neighbour (NFCN)
## This algorithm calculates from a base and a target data set (each an expression matrix with cells 
## as colnames and genes as rownames) the most connected cells with respect to a defined Receptor-Ligand couple.
##----------------------------------------------------------------------------##
##
## Input Parameter Matrix Basis (matrix), Matrix Target (matrix), 
## Receptor: A string of the receptor genes. Example: c("IL10RA", "IL10RB")
## Ligand: A string of the ligand genes. Example: c("IL10")
## GeneSets (data.frame) with: 
## first row: Genes Induction (a geneSet that is responsible to induce the ligand expression)
## second row: Genes Response (a geneSet that respond to a ligand activation)

# Get bash args
args <- commandArgs(TRUE)
#print(args)

argumente=args[grepl("^--", args)]
#print(argumente)
#print(length(argumente))
# --help

if("--help" %in% args){
  print(cat("\n
            \n ####   This is the MILO Software Neares functional connected neighbour (NFCN)  #### 
            \n ####   By D. H. Heiland
            \n
            \n --Ligand: A string of the ligand genes. Example: IL10
            \n --Receptor: A string of the receptor genes. Example: IL10RA,IL10RB
            \n --Gensets: GeneSets: 
            \n            first row: Genes Induction (a geneSet that is responsible to induce the ligand expression)
            \n            second row: Genes Response (a geneSet that respond to a ligand activation)
            \n --Matrix_basis: Gene Expression Matrix with cells as colnames and rawnames as rownames
            \n --Matrix_target: Gene Expression Matrix with cells as colnames and rawnames as rownames
            \n --DimRed: Matrix of Dimensional reduction of your cells (UMAP//TSNE...)
            \n --Cores: Number of used cores Default=8
            \n --Output: Output Folder
            \n --quantil_test: Which quantil of connected cells should be used for further anaysis Default=0.8
            \n --VisLabOutput: Return a Output of DE for VisLab Default=F (if T->  adapt --quantil_test to more than 0.9)
            "))
  q()
}else{
  
  #Call Arguments from bash
  
  get_arguments_in_list=lapply(1:length(argumente), function(i){
    if(argumente[i]=="--Ligand"){
      Ligand_inp=args[which(args=="--Ligand")+1]
      Ligand_inp=unlist(strsplit(Ligand_inp, ","))
      #print(paste0("The used Ligand is: ", Ligand_inp))
      return(Ligand_inp)
    }else{
    if(argumente[i]=="--DimRed"){
      DimRed_inp=args[which(args=="--DimRed")+1]
        return(DimRed_inp)
      }else{
    if(argumente[i]=="--Receptor" ){
      Receptor_inp=args[which(args=="--Receptor")+1]
      Receptor_inp=unlist(strsplit(Receptor_inp, ","))
      #print(paste0("The used Receptor is: ", Receptor_inp))
      return(Receptor_inp)
    }else{
    if(argumente[i]=="--Gensets" ){
      Gensets_inp=args[which(args=="--Gensets")+1]
      #print(paste0("The used Gensets is: ", Gensets_inp))
      return(Gensets_inp)
    }else{
    if(argumente[i]=="--Matrix_basis" ){
      Matrix_basis_inp=args[which(args=="--Matrix_basis")+1]
      #print(paste0("The used Matrix_basis is: ", Matrix_basis_inp))
      return(Matrix_basis_inp)
    }else{
    if(argumente[i]=="--Cores" ){
        Cores_inp=args[which(args=="--Cores")+1]
        #print(paste0("Cores: ", Cores_inp))
        return(Cores_inp)
      }else{
    if(argumente[i]=="--Output" ){
      Output_inp=args[which(args=="--Cores")+1]
          #print(paste0("Output: ", Output_inp))
          return(Output_inp)
        }else{
    if(argumente[i]=="--quantil_test" ){
      quantil_test_inp=args[which(args=="--quantil_test")+1]
            #print(paste0("Output: ", quantil_test_inp))
            return(quantil_test_inp)
          }else{
    if(argumente[i]=="--VisLabOutput" ){
     VisLabOutput_inp=args[which(args=="--VisLabOutput")+1]
              #print(paste0("Output: ", VisLabOutput_inp))
              return(VisLabOutput_inp)
            }else{
    if(argumente[i]=="--Matrix_target" ){
      Matrix_target_inp=args[which(args=="--Matrix_target")+1]
      #print(paste0("The used Matrix_target is: ", Matrix_target_inp))
      return(Matrix_target_inp)
    }else{ stop("--- Arguments incorrect listed ----")
    }}}}}}}}}}
    
    
  })
  names(get_arguments_in_list)=argumente
  #print(!"--Cores" %in% names(get_arguments_in_list))
  

  
  
  #check Argumets
  if( length(unique(c("--Ligand","--DimRed", "--Receptor", "--Gensets", "--Matrix_basis", "--Matrix_target") %in% names(get_arguments_in_list)))==2 ) stop("Arguments Missing")
    
    
}

if(is.null(get_arguments_in_list[["--Cores"]])){cores=8}else{cores=get_arguments_in_list[["--Cores"]]}
if(is.null(get_arguments_in_list[["--quantil_test"]])){quantil_test=0.8}else{quantil_test=as.numeric(get_arguments_in_list[["--quantil_test"]])}
if(is.null(get_arguments_in_list[["--VisLabOutput"]])){VisLabOutput=F}else{VisLabOutput=get_arguments_in_list[["--VisLabOutput"]]}
if(is.null(get_arguments_in_list[["--Output"]])){Output=getwd()}else{Output=get_arguments_in_list[["--Output"]]}

out=paste0("\n
           \n Input Parameters:
           \n --Ligand: ",get_arguments_in_list[["--Ligand"]],"
           \n --Receptor: ",get_arguments_in_list[["--Receptor"]],"
           \n --Gensets: ",get_arguments_in_list[["--Gensets"]],"
           \n --Matrix_basis: ",get_arguments_in_list[["--Matrix_basis"]],"
           \n --Matrix_target: ",get_arguments_in_list[["--Matrix_target"]],"
           \n --DimRed: ",get_arguments_in_list[["--DimRed"]],"
           \n 
           \n --Cores: ",cores,"
           \n --Output: ",Output,"
           \n --quantil_test: ",quantil_test,"
           \n --VisLabOutput: ",VisLabOutput,"
           ")

print(cat(out))
output_f=Output
dimReduction=read.table(get_arguments_in_list[["--DimRed"]])

## Create input

#Scores
print(cat("##### First Data Set Loaded ######"))
score=read.table(get_arguments_in_list[["--Gensets"]], stringsAsFactors = F)
Activation_score=unique(score[,2])[!unique(score[,2])==""]
Induction_score=unique(score[,1])[!unique(score[,1])==""]

#Input lists
R_Input=list(get_arguments_in_list[["--Receptor"]], Activation_score)
L_Input=list(get_arguments_in_list[["--Ligand"]], Induction_score)




#Defining classes for input data 
data_counts=setClass("data_counts", slots = c(counts="matrix", 
                                              norm_Exp = "matrix", 
                                              batch_corrected = "matrix", 
                                              used="matrix"))

R_L_Set <- setClass("R_L_Set", slots = c(  R_Input="list",
                                           L_Input = "list",
                                           data = "data_counts" ,
                                           spatial_info="data.frame",
                                           Cluster_to_check="list",
                                           UMAP="data.frame",
                                           data_fdat="data.frame"
))


setMethod("initialize",signature = "R_L_Set", definition = function(.Object, R_Input, L_Input,data){
  .Object@R_Input=R_Input
  .Object@L_Input=L_Input
  .Object@data@norm_Exp=data
  return(.Object)
})
print(cat("##### Read Basis ######"))
basis=as.matrix(read.table(get_arguments_in_list[["--Matrix_basis"]]))
cells_R=colnames(basis)

print(cat("##### Read Target ######"))
target=as.matrix(read.table(get_arguments_in_list[["--Matrix_target"]]))
cells_L=colnames(target)


if(!identical(rownames(basis), rownames(target))) stop("Genes from Input matrix do not match ")
exp_db=cbind(basis,target)

R_L_Set=R_L_Set(R_Input, L_Input, data=exp_db)

#First ealuate the z-scored enrichment of Induktion and Response Genesets 
library(GSVA)

print(cat("##### Start Analysis Step 1/2 ######"))

Activation_score=R_L_Set@R_Input[[2]]
Induction_score=R_L_Set@L_Input[[2]]
Ligant=R_L_Set@L_Input[[1]]
Receptor=R_L_Set@R_Input[[1]]

gs=list(Activation_score, Induction_score)
es.max <- gsva(exp_db, gs, mx.diff=T, method=c("zscore"), parallel.sz=cores)

#Input data used for map
ExpR=colMeans(exp_db[Receptor, cells_R])
ExpL=exp_db[Ligant, cells_L]

Exp_TC=exp_db[, cells_R]
Exp_M=exp_db[Ligant, cells_L]


### plot

plot_PL_Trajectory=function(es.max,
                            modus=c("spatial", "scRNAseq")[2],
                            ExpR,ExpL, 
                            quantil_test=.8,
                            DimRed,
                            Return_vislab=F,
                            plot_model=T,
                            plot_connectivity=F,
                            plot_barplots=F,...){
  
  #required function
  norm_s=function(x){(x-min(x))/(max(x)-min(x))}
  map2color<-function(x,pal,limits=NULL){
    if(class(x)=="numeric"){
      if(is.null(limits)) limits=range(x)
      pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
    }else{
      print(x[!duplicated(x)])
      da=data.frame(Terms=x[!duplicated(x)], Nr=seq(1,length.out = length(x[!duplicated(x)])))
      da$col=colorRampPalette(pal)(nrow(da))[da[,2]]
      daf=data.frame(x=x, col=1)
      for(i in 1:length(x)){
        daf[i, ]$col=da[da$Terms==daf[i, "x"], ]$col
        
      }
      
      return(list(daf$col, da))
    }
    
  }
  library(viridis)
  library(scales)
  if(modus=="scRNAseq"){
    
    message("Start Fitting the Model by using sc_RNAseq Data")
    #Fit model
    score1=es.max[1,cells_R]
    upper=quantile(score1, .995)
    score1[score1>upper]=jitter(rep(upper,length(score1[score1>upper])), factor=upper)
    y1=2-norm_s(score1)
    
    
    upper=quantile(ExpR, .995)
    ExpR[ExpR>upper]=jitter(rep(upper,length(ExpR[ExpR>upper])), factor=upper)
    x1=2-(norm_s(ExpR))
    
    score=es.max[2,cells_L]
    upper=quantile(score, .995)
    score[score>upper]=jitter(rep(upper,length(score[score>upper])), factor=upper)
    y2=norm_s(score)
    
    upper=quantile(ExpL, .995)
    ExpL[ExpL>upper]=jitter(rep(upper,length(ExpL[ExpL>upper])), factor=upper)
    x2=(norm_s(ExpL))
    
    
    #Start ploting
    message("Start Create Plots")
    
    #dev.off()
    
    pdf(useDingbats = F, paste0(output_f,"/",Sys.time(),"Output_Analysis_NFCN.pdf"))
    
    if(plot_model==T){
      
      scale_f=0.4
      plot(NA,xlim = c(0-scale_f,2+scale_f), ylim=c(0-scale_f,2+scale_f), axes=F, xlab="", ylab="")
      lo <- loess(c(y1,y2) ~ c(x1,x2))
      xl <- seq(0.05, 1.9, length.out = 1000)
      
      
      
      points(c(x1,x2), jitter(predict(lo,c(x1,x2)), factor=900), 
             col=c(map2color(x1, rev(viridis(20))),map2color(x2, (viridis(20))) ), pch=19,
             cex=c(y1-1,y2+0.3)*2)
      
      points(xl, predict(lo,xl), type="l", lty=2, lwd=1.5, col="black")
      
      #Add scales
      z_point=0
      up_z_point=2
      arrows(z_point,z_point, 1, z_point,  length=0.1) 
      arrows(z_point,z_point, z_point, 1,  length=0.1)
      arrows(up_z_point,up_z_point, 1, up_z_point,  length=0.1) 
      arrows(up_z_point,up_z_point, up_z_point, 1,  length=0.1)
      
      threshold=quantil_test
      
      polygon(x=c(threshold,1+(1-threshold), (1-threshold)+1, threshold ), 
              y=c(threshold,threshold,(1-threshold)+1, (1-threshold)+1), lty=2, lwd=1.5)
      
      #lables
      text(x=z_point+0.4, y=z_point-0.1, labels = "Expression Ligand")
      text(x=z_point-0.1, y=z_point+0.4, labels = "GSVA-Score Induction", srt=90)
      
      text(x=up_z_point-0.4, y=up_z_point+0.1, labels = "Expression Receptor")
      text(x=up_z_point+0.1, y=up_z_point-0.4, labels = "GSVA-Score Activation",srt=90 )
      
      #Show top cells in cluster
      
      #quantil_test=0.8
    }
    if(plot_barplots==T){
      
      #plot
      #-> Calculate percentage
      A=table(TC@meta.data[names(ExpR[ExpR>quantile(ExpR,quantil_test)]), ]$seurat_clusters)
      B=table(TC@meta.data$seurat_clusters)
      C=A/B
      C=prop.table(C)
      barplot((C), col=brewer.pal(5, "Set1"), main = "Cluster Receptor")
      
      
      A=table(ds@meta.data[names(ExpL[ExpL>quantile(ExpL,quantil_test)]), ]$seurat_clusters)
      B=table(ds@meta.data$seurat_clusters)
      C=A/B
      barplot(C, col=brewer.pal(5, "Set1"), main = "Cluster Ligant")
      
    }
    if(plot_connectivity==T){
      
      plot(DimRed, pch=19, col=alpha("black", 0.2), bty="n", axes=F, xlab="", ylab="")
      #points(ds@reductions$umap@cell.embeddings[names(ExpL[ExpL>quantile(ExpL,quantil_test)]), ],col="navy", pch=19,cex=0.3 )
      #plot(ds@reductions$umap@cell.embeddings, pch=19, col=alpha("grey", 0.2), bty="n", axes=F, xlab="", ylab="")
      #points(ds@reductions$umap@cell.embeddings[names(ExpR[ExpR>quantile(ExpR,quantil_test)]), ],col="darkgreen", pch=19, cex=0.3 )
      
      #Add connection lines
      
      #(1) Rank gene set
      ExpL_r=y1[y1>quantile(y1,quantil_test)]
      ExpR_r=y2[y2>quantile(y2,quantil_test)]
      
      ExpR_r=ExpR_r[order(ExpR_r)]
      ExpL_r=ExpL_r[order(ExpL_r)]
      
      df_r=data.frame(ExpR_r,QT=ecdf(ExpR_r)(ExpR_r), row.names = names(ExpR_r))
      df_l=data.frame(ExpL_r,QT=ecdf(ExpL_r)(ExpL_r), row.names = names(ExpL_r))
      col=alpha(map2color(df_r$QT, inferno(50)), 0.7)
      
      
      
      for(i in 1:nrow(df_r)){#Loop for connected lines defined by neares neighbour based on quantil similary
        inp=df_r[i,2]
        x=df_l[,2]
        FROM=rownames(df_r)[i]
        TO=rownames(df_l)[which.min(abs(x - inp))]
        
        FROM=gsub('\\.', '-', FROM)
        TO=gsub('\\.', '-', TO)
        
        start=DimRed[FROM, ]
        end=DimRed[TO, ]
        curveMaker <- function(x1, y1, x2, y2, ...){
          curve( plogis( x, scale = 0.9, loc = (x1 + x2) /2 ) * (y2-y1) + y1, 
                 x1, x2, add = TRUE, ...)
        }
        
        #curveMaker(start[1], start[2], end[1], end[2], col=col[i], lwd = inp*2)
        
        polygon(x=c(start[1], end[1]),y=c(start[2], end[2]),border = col[i], lwd = inp*2)
        
      }
      
      
    }
    
    dev.off()
    
    
    
    
    
    
    #Map Tcells
    
    #plot(TC@reductions$umap@cell.embeddings, pch=19, col=alpha("grey", 0.2), bty="n", axes=F, xlab="", ylab="")
    #points(TC@reductions$umap@cell.embeddings[names(ExpR[ExpR>quantile(ExpR,quantil_test)]), ],col=viridis(10)[6], pch=19, cex=0.8 )
    
    
    #compare top connected and not connected 
    
    if(Return_vislab==T){
      
      ExpL_r_up=names(y2[y2>quantile(y2,0.975)])
      ExpL_r_down=names(y2[y2<quantile(y2,0.025)])
      
      df_De=data.frame(Sample=c(ExpL_r_up,ExpL_r_down ), row.names = c(ExpL_r_up,ExpL_r_down ))
      df_De$group="up"
      df_De[ExpL_r_down,]$group="down"
      df_De$group=as.factor(df_De$group)
      df_De$treatment=as.factor(df_De$group)
      
      DeSeq_matrix=ds@assays$RNA@counts[, rownames(df_De)]
      dim(DeSeq_matrix)
      
      Data.character=list(ExpID=paste0(Ligant,"_",Receptor[1]), 
                          Researcher="DHH", 
                          Bioinformatic="DHH",
                          Institute="MILO",
                          Date_of_Analysis=Sys.Date(),
                          Datatype="scRNA-seq",
                          Sequencer="Illumina",
                          Species="Human")
      
      #source("~/Desktop/Projekt_Metabolom/Bioinformatic Tests/VIS_Lab_v1.4.1/pipeline_Create_Vizfile.R")
      pipline_Create_Vizfile=function(Modus=c("fromMatrix", "fromcOutcounts")[2], 
                                      CountMatrix=NA,
                                      f_data, 
                                      countfiles, 
                                      Data.character, 
                                      HOMER=F, 
                                      Addgene=F, 
                                      minCounts=5,
                                      norm=c("rlog", "vst")[1],
                                      Outputname=Data.character[[1]],
                                      export){
        
        message("Start Pipeline")
        ##----------------------------------------------------------------------------##
        ## Create the new class
        ##----------------------------------------------------------------------------##
        
        data_counts=setClass("data_counts", slots = c(counts="matrix", 
                                                      norm_Exp = "matrix", 
                                                      batch_corrected = "matrix", 
                                                      used="matrix"))
        
        VIZ_MILO <- setClass("VIZ_MILO", slots = c(Data.character="list",
                                                   fdata = "data.frame", 
                                                   DE="list",
                                                   condition_Analysis="character",
                                                   data =   "data_counts" ,
                                                   AutoEncoder="matrix",
                                                   UMAP="matrix",
                                                   cluster="data.frame",
                                                   Marker_genes="data.frame",
                                                   GSEA_sample_list="matrix",
                                                   GSEA="list",
                                                   Gene_Sets_used="data.frame"
        ))
        
        setMethod("initialize",signature = "VIZ_MILO", definition = function(.Object, count_matrix, f_data){
          .Object@fdata=f_data
          return(.Object)
        })
        
        ##----------------------------------------------------------------------------##
        ## Reqired Functions
        ##----------------------------------------------------------------------------##
        input_from_OutCounts=function(Out_counts,Samples_discription ){
          samples=nrow(Samples_discription)
          All_genes=unlist(lapply(1:samples,function(i){Out_counts[[i]]$HUGO}))
          All_genes=as.data.frame(as.character(All_genes[!duplicated(All_genes)]))
          names(All_genes)="genes"
          print("--------------- Start Merge Counts ------------------------")
          pb=txtProgressBar(min = 1, max = samples, initial = 1, char = "#",
                            width = 40, title=paste0("Progress Config Countfiles"))
          
          
          for(i in 1:samples){
            
            setTxtProgressBar(pb, i)
            l=Out_counts[[i]]
            l[,1]=as.character(l[,1])
            All_genes[,i+1]=0
            for(ix in 1:nrow(All_genes)){
              count=sum(l[l[,1]==as.character(All_genes[ix,1]),2])
              All_genes[ix,i+1]=count
              
            }
            
          }
          
          close(pb)
          
          input=All_genes
          return(input)
        }
        
        DESEQ_DHH=function(dds){
          
          
          print("--------------- DESEQ_DHH  ------------------------")
          dds <- dds[ rowSums(counts(dds)) >= minCounts, ]
          dim(dds)
          if(norm=="rlog"){
            rld <- rlog(dds, blind=FALSE)
          }else{rld <- varianceStabilizingTransformation(dds)}
          
          if(length(unique(sample@fdata$Batch))>1){assay(rld) <- sva::ComBat(assay(rld), batch=sample@fdata$Batch)}
          sample@data@norm_Exp=as.matrix(assay(rld))
          dds=DESeq(dds)
          return(dds)
        }
        
        
        ##----------------------------------------------------------------------------##
        ## Start input basic informations
        ##----------------------------------------------------------------------------##
        
        
        sample=VIZ_MILO(f_data=f_data)
        rownames(sample@fdata)=sample@fdata$Sample 
        sample@Data.character=Data.character
        
        
        if(Modus=="fromcOutcounts"){
          
          ##----------------------------------------------------------------------------##
          ## Nanopore Seq
          ##----------------------------------------------------------------------------##
          ##----------------------------------------------------------------------------##
          ## Input Coount files
          ##----------------------------------------------------------------------------##
          
          
          
          if(length(countfiles)>1){
            
            Out_counts=do.call(c,lapply(1:length(countfiles), function(z) readRDS(countfiles[z]) ))
            
          }else(Out_counts=readRDS(countfiles))
          
          Samples_discription=sample@fdata
          counts=input_from_OutCounts(Out_counts,Samples_discription)
          counts=data.frame(counts[2:ncol(counts)], row.names = counts$genes)
          sample@data@counts=as.matrix((na.omit(counts)))
          colnames(sample@data@counts)=sample@fdata$Sample 
          
        }
        
        if(Modus=="fromMatrix"){
          
          sample@data@counts=as.matrix((na.omit(CountMatrix)))
          colnames(sample@data@counts)=sample@fdata$Sample 
          
          counts=as.matrix(CountMatrix)
          
        }
        
        ##----------------------------------------------------------------------------##
        ## DE
        ##----------------------------------------------------------------------------##
        
        
        counts=as.matrix(counts)
        feat=sample@fdata
        
        
        library(DESeq2)
        feat$group=as.factor(sample@fdata$group)
        feat$treatment=as.factor(sample@fdata$treatment)
        
        
        
        dds_group=DESEQ_DHH(DESeqDataSetFromMatrix(countData = sample@data@counts ,DataFrame(sample@fdata), ~ group))
        dds_treat=DESEQ_DHH(DESeqDataSetFromMatrix(countData = sample@data@counts ,DataFrame(sample@fdata), ~ treatment))
        
        dds=DESeqDataSetFromMatrix(countData = sample@data@counts ,DataFrame(sample@fdata), ~ group)
        dds <- dds[ rowSums(counts(dds)) >= minCounts, ]
        dim(dds)
        if(norm=="rlog"){
          rld <- rlog(dds, blind=FALSE)
        }else{rld <- varianceStabilizingTransformation(dds)}
        
        if(length(unique(sample@fdata$Batch))>1){assay(rld) <- sva::ComBat(assay(rld), batch=sample@fdata$Batch)}
        sample@data@norm_Exp=as.matrix(assay(rld))
        
        
        sample@DE=list(dds_group, dds_treat)
        sample@fdata=feat
        names(sample@DE)=c("group", "treatment")
        sample@condition_Analysis=c("group", "treatment")
        
        diff_gene=results(dds_group)
        
        
        ##----------------------------------------------------------------------------##
        ## Data ready for the GSEA
        ##----------------------------------------------------------------------------##
        
        ### Add new GS or load the old ones
        
        
        #Add new
        
        ### Selected GS for Analysis
        library(clusterProfiler)
        library(gtools)
        library(DESeq2)
        library(clusterProfiler)
        
        
        if(Addgene==T){
          newGS=readRDS("~/Desktop/Projekt_Metabolom/Bioinformatic Tests/VIS_Lab_v1.4.1/MILO_Geneset.RDS")
          
          nGS_TCells=read.csv("~/Desktop/Projekt_Metabolom/Bioinformatic Tests/VIS_Lab_v1.4.1/New_Genesets/TC_GS.csv", sep=";")
          nGS_TCells=nGS_TCells[,c("ont", "gene")]
          
          #Stimmulation Datasets
          HOMER_DE_to_gmt=function(folder, GS_name){
            message(paste0("Number of new GS lists were found in folder:   ", length(dir(folder))))
            
            GS_fromfolder=as.data.frame(do.call(rbind,lapply(1:length(dir(folder)), function(i){
              y=read.csv(dir(folder)[i])
              
              a=gsub("group", "",dir(folder)[i] )
              a=gsub(".txt", "",a)
              ont=paste0("MILO_",GS_name,"_",a)
              
              out=data.frame(ont=ont, gene=as.character(y$V1))
              return(out)
            })))
            
            return(GS_fromfolder)
          }
          
          nGS_TCells_Stim=HOMER_DE_to_gmt(folder="~/Desktop/RNA-Seq_Analysis/Nico/Merged/HOMER/Input",
                                          GS_name="TCELL_STIM")
          
          
          names(newGS)==names(nGS_TCells_Stim)
          
          MsigDB1=rbind(MsigDB1,newGS, nGS_TCells, nGS_TCells_Stim)
          
          saveRDS(MsigDB1, "MILO_MsigDB1_14_Feb.RDS")
          
          
        }
        
        
        #load old
        MsigDB1=readRDS("~/Desktop/Projekt_Metabolom/Bioinformatic Tests/VIS_Lab_v1.4.1/MILO_MsigDB1_14_Feb.RDS")
        
        
        sample@Gene_Sets_used=MsigDB1
        #Get new genesets
        
        
        #Input explanation:  res=data.frame, GMT -> GeneSets
        GSEA_DHH=function(res, GMT){
          res=res[order(res$log2FoldChange, decreasing=T), ]
          ranks=data.frame(ID=as.character(rownames(res)),t=as.numeric(res$log2FoldChange))
          ranks <- setNames(ranks$t, ranks$ID)
          egmt1 <- GSEA(ranks, TERM2GENE=GMT, verbose = T, pvalueCutoff= 1)
        }
        
        
        
        #Matrix with all possible GSEA conditions
        to_test=levels(sample@fdata[,sample@condition_Analysis[1]])
        df_test=permutations(n=length(to_test), r=2, v=to_test, repeats.allowed=F)
        cond_mat=matrix(NA,nrow(df_test), 2)
        colnames(cond_mat)=sample@condition_Analysis
        for(z in 1:nrow(df_test)){
          cond_mat[z,sample@condition_Analysis[1]]=paste0(df_test[z,1],"_",df_test[z,2])
          
        }
        to_test=levels(sample@fdata[,sample@condition_Analysis[2]])
        df_test=permutations(n=length(to_test), r=2, v=to_test, repeats.allowed=F)
        for(z in 1:nrow(df_test)){
          cond_mat[z,sample@condition_Analysis[2]]=paste0(df_test[z,1],"_",df_test[z,2])
          
        }
        
        
        #Start list with GSEA of all conditions
        GSEA_out=lapply(1:length(sample@condition_Analysis), function(i){
          
          to_test=levels(sample@fdata[,sample@condition_Analysis[i]])
          df_test=permutations(n=length(to_test), r=2, v=to_test, repeats.allowed=F)
          
          GSEA_out2=lapply(1:nrow(df_test), function(j){
            
            DE_in=sample@DE[[sample@condition_Analysis[i]]]
            res=data.frame(results(DE_in, contrast = c(sample@condition_Analysis[i],df_test[j,1],df_test[j,2])))
            return(GSEA_DHH(res, MsigDB1))
            
          })
          
        })
        
        #Add correct name to lists
        names(GSEA_out)=sample@condition_Analysis
        names(GSEA_out[[sample@condition_Analysis[1]]])=as.character(na.omit(cond_mat[,1]))
        names(GSEA_out[[sample@condition_Analysis[2]]])=as.character(na.omit(cond_mat[,2]))
        
        #Export into class
        sample@GSEA=GSEA_out
        sample@GSEA_sample_list=cond_mat
        
        
        ##----------------------------------------------------------------------------##
        ## HOMER
        ##----------------------------------------------------------------------------##
        
        if(HOMER==T){
          message("HOMER Motif Enrnrichment Analysis")
          
          path=getwd()
          dir.create("HOMER")
          setwd(paste0(getwd(), "/HOMER"))
          dir.create("Input")
          dir.create("Output")
          message("HOMER Data will be saved in: ")
          message(getwd())
          
          #save genelist
          length_motiv=100
          setwd(paste0(getwd(), "/Input"))
          
          bash=c("#!/bin/bash")
          
          HOMER=lapply(1:length(sample@condition_Analysis), function(i){
            
            to_test=levels(sample@fdata[,sample@condition_Analysis[i]])
            df_test=permutations(n=length(to_test), r=2, v=to_test, repeats.allowed=F)
            
            HOMER2=lapply(1:nrow(df_test), function(j){
              
              setwd(paste0(path, "/HOMER/Input"))
              print(paste0(sample@condition_Analysis[i],df_test[j,1],"_UP_",df_test[j,2],length_motiv,".txt"))
              
              DE_in=sample@DE[[sample@condition_Analysis[i]]]
              res=data.frame(results(DE_in, contrast = c(sample@condition_Analysis[i],df_test[j,1],df_test[j,2])))
              
              genes_UP=head(rownames(res[order(res$log2FoldChange, decreasing = F), ]),length_motiv)
              write.table(as.matrix(genes_UP), 
                          file=paste0(sample@condition_Analysis[i],df_test[j,1],"_UP_",df_test[j,2],length_motiv,".txt"),
                          quote=FALSE, sep='\t', row.names = F)
              
              
              pathtopl="/Users/HenrikHeiland/homer/"
              
              command=paste0(paste0("findMotifs.pl"), " ",
                             paste0(path, "/HOMER/Input/", sample@condition_Analysis[i],df_test[j,1],"_UP_",df_test[j,2],length_motiv,".txt")
                             ," human ", 
                             paste0(path, "/HOMER/Output/",sample@condition_Analysis[i],df_test[j,1],"_UP_",df_test[j,2],length_motiv) ,
                             " -start -400 -end 100 -len 6,8,10,12 -p 8  ")
              
              bash<<-c(bash,command)
              
              
              
              
            })
            
            
          })
          
          write.table(data.frame(bash), "HOMER.sh", quote = F, row.names = F, col.names = F)
          
          
        }
        
        
        
        
        
        ##----------------------------------------------------------------------------##
        ## Export Data
        ##----------------------------------------------------------------------------##
        
        
        sample@data@norm_Exp
        sample@fdata=feat
        
        print(paste0(export,"/",Outputname,Sys.Date(), ".RDS"))
        
        
        saveRDS(sample,paste0(export,"/",Outputname,Sys.Date(), ".RDS"))
        
      }
      
      
      pipline_Create_Vizfile(Modus="fromMatrix", 
                             CountMatrix=DeSeq_matrix,
                             f_data=df_De, 
                             Data.character=Data.character, 
                             HOMER=F, 
                             Addgene=F, 
                             norm="vst",
                             minCounts=1,
                             Outputname=Data.character[[1]],
                             export=output_f)
      
      
    }
    
    
  }
  
  
}


print(cat("##### Start Analysis Step 2/2 ######"))

plot_PL_Trajectory(es.max,ExpR,ExpL,modus="scRNAseq", 
                   plot_model=T,
                   plot_connectivity=T,
                   quantil_test=quantil_test, 
                   DimRed=dimReduction, 
                   Return_vislab=VisLabOutput)

print(cat("##### Finish Create Output ######"))

#Rscript NFCN.R --Ligand IL10 --Receptor IL10RA,IL10RB --Gensets GS.txt --Matrix_basis Basis.txt --Matrix_target Target.txt --DimRed DimRed.txt

q()
