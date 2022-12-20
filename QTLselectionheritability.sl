initialize() {
 initializeMutationRate(1e-6);

 initializeMutationType("m1", 0.5, "f", 0.0); // neutral
 initializeMutationType("m2", 0.5, "n", 0.0, 1.0); // QTL
 m2.convertToSubstitution = F;  // for m2, don't ignore fixed mutations
 initializeGenomicElementType("g1", c(m1, m2), c(1, 0.05)); 
 initializeGenomicElement(g1, 0, 1e6 - 1); //g1 has 1,000,000 bp
 initializeRecombinationRate(1e-7);
 cat(num);
 
 
}
1 late() {
sim.addSubpop("p1", 500);  
sim.addSubpop("p2", 500);
sim.addSubpop("p3", 500);
sim.addSubpop("p4", 500);
sim.addSubpop("p5", 500);
p1.setMigrationRates(p2, 0.02);   
p2.setMigrationRates(p3, 0.02);
p3.setMigrationRates(p4, 0.02);
p4.setMigrationRates(p5, 0.02);
p2.setMigrationRates(p1, 0.02);                
p3.setMigrationRates(p2, 0.02);
p4.setMigrationRates(p3, 0.02);
p5.setMigrationRates(p4, 0.02);
} 


1: late() {                                     
    csvlines = readFile("SelectionHeritabilityFiles/FormattedSelectionH2.csv");
    csvlines = csvlines[substr(csvlines, 0,1) != "//"];
    
    linecount = 0;
    Y = 0;
    h2 = 0;
    for (line in csvlines)
     {
        if (linecount == 0)             ////////////////////////////////////////////////////////////////////CHANGE THIS: choose which species from csv line to use for selection and h2 
         {
            fields = strsplit(line, ",");
            Y = asFloat(fields[10]);                         // gamma as quadratic selection value
            h2 = asFloat(fields[21]);                       //h2 
            //cat("The quad var is ");
            //cat(Y);
            //cat("\n");
        }
        linecount = linecount + 1;
     }
                                                // evaluate POPULATION 1 QTLs to get phenotypes and fitness
 inds1 = p1.individuals;
 additive1 = inds1.sumOfMutationsOfType(m2);
 V_A1 = sd(additive1)^2;
 V_E1 = (V_A1 - h2 * V_A1) / h2;
 env1 = rnorm(size(inds1), 0.0, sqrt(V_E1));
 phenotypes1 = additive1 + env1;
 inds1.fitnessScaling = 1.0 + (Y)^2*(phenotypes1);
 inds1.tagF = phenotypes1;
                                                // evaluate POPULATION 2 QTLs to get phenotypes and fitness
 inds2 = p2.individuals;
 additive2 = inds2.sumOfMutationsOfType(m2);
 V_A2 = sd(additive2)^2;
 V_E2 = (V_A2 - h2 * V_A2) / h2;
 env2 = rnorm(size(inds2), 0.0, sqrt(V_E2));
 phenotypes2 = additive2 + env2;
 inds2.fitnessScaling = 1.0 + (Y)^2*(phenotypes2);
 inds2.tagF = phenotypes2;
                                                 // evaluate POPULATION 3 QTLs to get phenotypes and fitness
 inds3 = p3.individuals;
 additive3 = inds3.sumOfMutationsOfType(m2);
 V_A3 = sd(additive3)^2;
 V_E3 = (V_A3 - h2 * V_A3) / h2;
 env3 = rnorm(size(inds3), 0.0, sqrt(V_E3));
 phenotypes3 = additive3 + env3;
 inds3.fitnessScaling = 1.0 + (Y)^2*(phenotypes3);
 inds3.tagF = phenotypes3; 
                                                 // evaluate POPULATION 4 QTLs to get phenotypes and fitness
 inds4 = p4.individuals;
 additive4 = inds4.sumOfMutationsOfType(m2);
 V_A4 = sd(additive4)^2;
 V_E4 = (V_A4 - h2 * V_A4) / h2;
 env4 = rnorm(size(inds4), 0.0, sqrt(V_E4));
 phenotypes4 = additive4 + env4;
 inds4.fitnessScaling = 1.0 + (Y)^2*(phenotypes4);
 inds4.tagF = phenotypes4;
                                                 // evaluate POPULATION 5 QTLs to get phenotypes and fitness
 inds5 = p5.individuals; 
 additive5 = inds5.sumOfMutationsOfType(m2);
 V_A5 = sd(additive5)^2;
 V_E5 = (V_A5 - h2 * V_A5) / h2;
 env5 = rnorm(size(inds5), 0.0, sqrt(V_E5));
 phenotypes5 = additive5 + env5;
 inds5.fitnessScaling = 1.0 + (Y)^2*(phenotypes5);
 inds5.tagF = phenotypes5;
}


fitness(m2) {
return 1.0; // QTLs are neutral; fitness effects are handled below. make fitness of mutations not affect fitness, so only phenotype fitness affects it.
}


1:5000 late() { 
 name = "Acroceph0."+num;                 ////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 meanPhenotype1 = mean(p1.individuals.tagF);
 meanPhenotype2 = mean(p2.individuals.tagF);
 meanPhenotype3 = mean(p3.individuals.tagF); 
 meanPhenotype4 = mean(p4.individuals.tagF); 
 meanPhenotype5 = mean(p5.individuals.tagF); 
 phenos = asString( format("%.2f", meanPhenotype1)) + ", " + asString( format("%.2f", meanPhenotype2)) + ", " + asString( format("%.2f", meanPhenotype3)) + ", " + asString( format("%.2f", meanPhenotype4)) + ", " + asString( format("%.2f", meanPhenotype5));
 writeFile("meanPhenotype/"+name+".csv", phenos, append = T);   // Run until we reach the fitness peak
 } 


 

5000 late() {

// This is to get the VCF
name = "Acroceph0."+num;                  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

myInds1 = sample(p1.individuals,200); 
myInds2 = sample(p2.individuals,200);
myInds3 = sample(p3.individuals,200);
myInds4 = sample(p4.individuals,200);
myInds5 = sample(p5.individuals,200);
genos = c(myInds1.genomes, myInds2.genomes, myInds3.genomes, myInds4.genomes, myInds5.genomes);
genos.outputVCF(filePath="vcf/"+name+".vcf");

// This is list of those sample indiv's trait value for each pop
SampleTraits = asString(+ myInds1.tagF) +' ' +  asString(+ myInds2.tagF) + ' ' + asString(+ myInds3.tagF)  + ' ' + asString(+ myInds4.tagF)  + ' ' + asString(+ myInds5.tagF);
writeFile("meanTraitFiles/"+name+".csv", SampleTraits, append = T);
}






