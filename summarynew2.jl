using DelimitedFiles
using StatsBase

resultsim3m10 = readdlm("../resultnew/resultsim3m10.csv",',')

mean(resultsim3m10[:,3])
mean(resultsim3m10[:,7])

mean(resultsim3m10[:,1])
mean(resultsim3m10[:,5])

resultsim3m20 = readdlm("../resultnew/resultsim3m20.csv",',')

mean(resultsim3m20[:,3])
mean(resultsim3m20[:,7])

mean(resultsim3m20[:,1])
mean(resultsim3m20[:,5])



resultsim3m30 = readdlm("../resultnew/resultsim3m30.csv",',')

mean(resultsim3m30[:,3])
mean(resultsim3m10[:,7])

mean(resultsim3m30[:,1])
mean(resultsim3m30[:,5])



resultsim3s10 = readdlm("../resultnew/resultsim3s10.csv",',')

mean(resultsim3s10[:,2])
mean(resultsim3s10[:,4])

mean(resultsim3s10[:,1])
mean(resultsim3s10[:,3])

resultsim3s20 = readdlm("../resultnew/resultsim3s20.csv",',')

mean(resultsim3s20[:,2])
mean(resultsim3s20[:,4])

mean(resultsim3s20[:,1])
mean(resultsim3s20[:,3])



resultsim3s30 = readdlm("../resultnew/resultsim3s30.csv",',')

mean(resultsim3s30[:,2])
mean(resultsim3s30[:,4])

mean(resultsim3s30[:,1])
mean(resultsim3s30[:,3])
