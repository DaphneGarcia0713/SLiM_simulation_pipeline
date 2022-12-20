
library(ggplot2)    # for plotting eigenstrat

args = commandArgs(trailingOnly=TRUE)
name =  args[1]
print(name)

RealvsEstPlot = read.table(paste("RealvsEstimatePlotCoordinates/",name,".coord",sep=""))

pdf(paste("RealvsEstimatePlots/",name,".pdf",sep=""))
EvRplot <- ggplot(RealvsEstPlot, aes(x= V1, y = V2)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Estimated Values", y="Real Values", title = paste("Real vs Estimated Effect Values for",name) )
#Make sure the plot on the right looks like how you want it to be saved!!!!!!!!
print(EvRplot)
dev.off()
