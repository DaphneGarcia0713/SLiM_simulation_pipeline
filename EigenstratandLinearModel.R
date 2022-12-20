#install.packages('tidyverse')
#install.packages('AssocTests')
library(AssocTests) # for eigenstrat and smt()
library(ggplot2)    # for plotting eigenstrat
library(LEA)        #for making the genofile
#library(reticulate) # for running python from wihin Rstudio
#detach("package:AssocTests", unload=TRUE)

args = commandArgs(trailingOnly=TRUE)
name =  args[1]
print(name)


print("hi")
#LEA
getwd()
output = vcf2geno(paste("vcf/",name,".vcf",sep = ""), paste("Genofiles/",name,".geno",sep = "")) #remember this file for real effect values!


# EIGENSTRAT 

EigenValues <- eigenstrat(genoFile = paste("Genofiles/",name,".geno",sep = ""), 
           outFile.Robj = paste("EigenstratCovariates/",name,".list",sep = ""),   ##CHANGE FILE NAME
           outFile.txt = paste("EigenstratCovariates/",name,".txt",sep = ""), rm.marker.index = NULL,
           rm.subject.index = NULL, miss.val = 9, num.splits = 10,
           topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE,
           iter.outlier = 5, sigma.thresh = 6)



#PLOTTING EIGENSTRAT RESULTS      

eigenplot=read.table(paste("EigenstratCovariates/",name,".txt",sep=""))                                                             
ggplot(eigenplot, aes(V1, V2)) +    #Only use V2 and V3 bc theyre the closest to (---) 1. create new folder for Eigenstrat plots 2. Check what happens w/ V1 and V2?
    geom_point() +
    geom_smooth(method="lm") +
    labs(x="x coordinate", y="y coordinate", title = paste("Eigenstrat coordinates for",name,"population's trait values"), caption="using logarithmic scale")
    #Make sure the plot on the right looks like how you want it to be saved!!!!!!!!
ggsave(paste("EigenstratPlots/",name,"Plot.pdf",sep=""))                                                           




## LINEAR MODEL 

eigenCovariates = read.table(paste("EigenstratCovariates/",name,".txt",sep = ""))            #covariates      
traits = read.table(paste("meanTraitFiles/",name,".csv",sep=""))                      ##CHANGE FILE
traits = c(t(t(traits)))
genofilescan = scan(paste("Genofiles/",name,".geno",sep=""),what = "character")     ##CHANGE FILE
eigenCovariatesLong <- nrow(eigenCovariates) + 600


i=1                                                      # is 1, bc genofilescan starts w/ [1] 
for(x in genofilescan) {                                    #for each line in genotable
    genotype = strsplit(genofilescan[i],split='')  #name genocharacters the list of each 0,1,2, etc.
    genotype = unlist(genotype)                   #unlisted coupled with numeric turns class() list to numeric
    g <- as.numeric(genotype)                           #    ^^
    
    mod <- lm(traits ~ eigenCovariates$V1 + eigenCovariates$V2 + g)
    #summary(mod)
    pvalue = summary(mod)$coefficients[,4]   # p value
    traitvalue = summary(mod)$coefficients[,1]   # est. trait value
    names(pvalue) <- NULL
    names(traitvalue) <- NULL
    if(length(pvalue) == 4){  #added this 11/15 bc some values have pvalue[4] as NA
        if(pvalue[4] < 0.01){
            out = c(pvalue[4],traitvalue[4], i)    #this is the result! 3 columns, the p value, the trait value, and i, which is number of mutations
            
            write(out, file = paste("EstimatedEffectValueFiles/",name,".estim",sep=""),
                  ncolumns = 3, append = TRUE, sep = " ")
        }
    }
    i=i+1
}  

#This normally takes 3 min with the Nov9sm genome size and sample size and etc.


#GO TO RealvsEstimatePlotting.py FOR NEXT STEP




##  TWO WAYS TO MAKE THE PLOTS #################

#RealvsEstPlot = read.table(paste("RealvsEstimatePlotCoordinates/",name,".coord",sep=""))

#thisplot = plot(RealvsEstPlot$V1, RealvsEstPlot$V2,
#                main=paste("Real vs Estimated Effect Values for",name),
#                xlab="Estimated Values",
#                ylab="Real Values",
#                type="p",
#                col="blue")


##### Another version of this same plot, but in ggplot!
#ggplot(RealvsEstPlot, aes(x= V1, y = V2)) +
#    geom_point() +
#    geom_smooth(method="lm") +
#    labs(x="Estimated Values", y="Real Values", title = paste("Real vs Estimated Effect Values for",name) )
    #Make sure the plot on the right looks like how you want it to be saved!!!!!!!!
ggsave(paste("RealvsEstimatePlots/",name,"RvEPlot.pdf",sep=""))  
#############his is the graph to use!!


