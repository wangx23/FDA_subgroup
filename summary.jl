using DelimitedFiles
using StatsBase

resultsim1 = readdlm("../result/resultsim1.csv",',')
resultsim1v2 = readdlm("../result/resultsim1v2.csv",',')
resultsim11 = readdlm("../result/resultsim11.csv",',')
resultsim12 = readdlm("../result/resultsim12.csv",',')
resultsim12v2 = readdlm("../result/resultsim12v2.csv",',')

resultsim1lam2 = readdlm("../result/resultsim1lam2.csv",',')

mean(resultsim1lam2[:,2])



resultsim12lam0 = readdlm("../result/resultsim12lam0.csv",',')
resultsim12lam2 = readdlm("../result/resultsim12lam2.csv",',')
resultsim12lam5 = readdlm("../result/resultsim12lam5.csv",',')


mean(resultsim12lam0[:,2])
mean(resultsim12lam2[:,2])
mean(resultsim12lam5[:,2])

sum(resultsim12lam0[:,5].==1)

findall(resultsim1lam2[:,2].==0)
findall(resultsim1lam2[:,6].==1)
