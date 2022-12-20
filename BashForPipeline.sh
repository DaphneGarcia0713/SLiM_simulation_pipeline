
echo Start of the Bash Script: input name
# Input species name without run number, same one as slim script
read FilesName




echo "Starting real effect values . py " $FilesName;

    
#There would be a function here that adds $varname and the variable for slim's num
#RealEffectValues.py makes a file in RealEffectValueFiles
python3 RealEffectValues.py $FilesName;
echo this is after RealEffectValues.py; 


#EigenstratandLinearModel.R uses LEA, ggplot2, Assoctests to make genofile, EigenstratPlots, and RealvsEstimatePlotCoordinates
 Rscript EigenstratandLinearModel.R $FilesName;
echo this is after EigenstratandLinearModel.R;

#RealvsEstimatePlotting.py makes a file in RealvsEstimatePlotCoordinates
 python3 RealvsEstimatePlotting.py $FilesName;
echo this is after RealvsEstimatePlotting.py;

#RealvsEstimatePlot.R uses ggplot to make the plot
 Rscript RealvsEstimatePlot.R $FilesName;
echo this is after RealvsEstimatePlot.R; 
