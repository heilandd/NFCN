### An optimized Version for the NFCN


##All packages
#packages=info$packages
#input is a df of session_info()$packages
#saveRDS(sessionInfo()$packages, file="packages.RDS")

checkP_DHH=function(df_packages){
  for(i in 1:nrow(df_packages)){
    #First check if package is installed
    package=df_packages$package[i]
    if(package %in% rownames(installed.packages()))
      {
      print(paste0("#######################--  Package: ",package, " is checked and installed--#######################  " ))
      
    }else{
      #check source
      source=df_packages$source[i]
      source=unlist(strsplit(source," "))[1]
      if(source=="CRAN"){
        ##source is CRAn
        install.packages(package)
      }else{
        if(source=="Bioconductor"){
          ##source is Bioconductor
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
            BiocManager::install(package)
        }else{
          ##Source is github
          library(devtools)
          #get githubinfo
          source=df_packages$source[i]
          source2=unlist(strsplit(source," "))[2]
          source3=gsub("\\(|\\)", "", source2)
          install_github(source3)
        }
      }
    }
  }
}
checkP_DHH(readRDS("packages.RDS"))



