
![Image](https://github.com/heilandd/NFCN/blob/master/Img.png)




This algorithm identifies the nearest connected cell from a base and a target data set with respect to a defined receptor-ligand pair. 
Data should be in a defined format: expression matrices with cells as colon names and genes as down names, of base cells and target pupulation respectively.  
Data generated by milolab.com



## How to install the package from GitHub

Install in Mac/Linux
Requirements: 

Xcode Command Line Tools or other C++ Compiler
R-Software (https://cran.r-project.org)

```
git clone https://github.com/heilandd/NFCN

cd NFCN

#Example:
Rscript NFCN.R --help

```

#Input Parameter: 

- --Ligand: A string of the ligand genes. Example: IL10
- --Receptor: A string of the receptor genes. Example: IL10RA,IL10RB
- --Gensets: GeneSets: 
            first row: Genes Induction (a geneSet that is responsible to induce the ligand expression)
            second row: Genes Response (a geneSet that respond to a ligand activation)
 - --Matrix_basis: Gene Expression Matrix with cells as colnames and rawnames as rownames
 -  --Matrix_target: Gene Expression Matrix with cells as colnames and rawnames as rownames
 - --DimRed: Matrix of Dimensional reduction of your cells (UMAP//TSNE...)
 - --Output: Output Folder
 - --quantil_test: Which quantil of connected cells should be used for further anaysis Default=0.8
 - --VisLabOutput: Return a Output of DE for VisLab Default=F (if T->  adapt --quantil_test to more than 0.9)


In any case of problems send me the Error per mail.
Dieter.henrik.heiland@uniklinik-freiburg.de

### Tutorial and Examples

In the following section we provide an example of how you can use NFCN on your MAC/Linux (Windows is not supported) First clone the github and installation requirements. In the requirements file (.RDS) you will all requirements. The script install.R contains a function to install all packages automatically.


#### Open your Terminal:

```
Rscript install.R

```

Next, we select the data we want to analyze. First, a set of cells that reflect your basic cells. The algorithm uses this set of cells and scans for connected cells of your target data set. In other words, you filter a specific subset of cells you are interested in and find connected cells in the rest of your data set. We need two gene expression matrices in which the column names reflect the cells (if the cells overlap between data sets, you will get errors) and the row names are HUGO symbols. 

Since we are trying to simulate a defined biological process, we have to define information about activations in the downstream and upstream pathways accurately. This information is sometimes difficult to obtain, but it is an essential part of your success. An example: 

For example, if we want to study the interaction based on interferon gamma, we first need to define the interaction partners: the pair IFNG and IFNGR1 & IFNGR2. In addition, we need a set of genes that describe the induction of IFNG (by the releasing cell) and the IFNG response (by the receiving cell). The data has to be in a table where the first column contains induction genes and the second column contains response genes. Sometimes GeneSets overlap (this does not lead to errors). 

| Induction  | Response    
| ---------- |:---------:| 
| IFNAR1     | ADAR      | 
| IFNAR2     | APOL6     |  
| NFKB1      | CCL2      |    
| JUN        | CD274     |
| IRF1       | CXCL9     |
|  ...       | ...       |

If you want to plot connections into your TSNE//UMAP, provide a table with Dim1 and Dim2 (HUGO genes as rownames). If you use seurat you can export the DimRed by: 

```
write.table(yourseurat@reductions$umap@cell.embeddings, "DimRed.txt")
```
From monocle:

```
write.table(reducedDims(yourcds)[["UMAP"]], "DimRed.txt")
```

Now, you have all your data, lets run the script.


#### Open your Terminal:

```
Rscript NFCN.R --Ligand IFNG --Receptor IFNGR1,IFNGR2 --Gensets pathto/GS.txt --Matrix_basis pathto/Basis.txt --Matrix_target pathto/Target.txt --DimRed pathto/DimRed.txt


```

You can export the most and less connected cells for DE analysis. We integrated an output file for the VisLab (a visualizer of gene expression data). This is helpful to further analze your data. Visit https://github.com/tumormetabolismfr/Vis_Lab1.5 for more information. 

```
 --VisLabOutput T

```


#### Pitfalls: 
- Never put space if you use more than one receptor, allways use a "," to seperate.
- Expression matrix, use the write.table() function from R to export  
- Never use more cores than total nr of cores - 2 
- If you use --VisLabOutput T also adapt the  --quantil_test parameter to 0.95 or higher. If your data set is very large (>10tsd cells) use 0.995 or 0.999 
- If you use --VisLabOutput T install VisLab first!






### Update:




### Authors

D. H. Heiland  The MILO Laboratoy, Medical-Center Freiburg, University of Freiburg, themilolab.com
