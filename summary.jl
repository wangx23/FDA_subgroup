using DelimitedFiles
using StatsBase

resultsim1 = readdlm("../result/resultsim1.csv",',')
resultsim1v2 = readdlm("../result/resultsim1v2.csv",',')
resultsim11 = readdlm("../result/resultsim11.csv",',')
resultsim12 = readdlm("../result/resultsim12.csv",',')
resultsim12v2 = readdlm("../result/resultsim12v2.csv",',')

resultsim1lam0 = readdlm("../result/resultsim1lam0.csv",',')
resultsim1lam2 = readdlm("../result/resultsim1lam2.csv",',')
resultsim1lam5 = readdlm("../result/resultsim1lam5.csv",',')
resultsim1lam10 = readdlm("../result/resultsim1lam10.csv",',')

mean(resultsim1[:,2])
mean(resultsim1lam0[:,2])
mean(resultsim1lam2[:,2])
mean(resultsim1lam5[:,2])
mean(resultsim1lam10[:,2])

resultsim12lam0 = readdlm("../result/resultsim12lam0.csv",',')
resultsim12lam2 = readdlm("../result/resultsim12lam2.csv",',')
resultsim12lam5 = readdlm("../result/resultsim12lam5.csv",',')
resultsim12lam10 = readdlm("../result/resultsim12lam10.csv",',')

mean(resultsim12lam0[:,2])
mean(resultsim12lam2[:,2])
mean(resultsim12lam5[:,2])
mean(resultsim12lam10[:,2])

sum(resultsim12lam0[:,5].==1)

findall(resultsim1lam2[:,2].==0)
findall(resultsim1lam2[:,6].==1)


uniqtm = unique(tm)
Bmt = orthogonalBsplines(tm, knots)
Bmi = orthogonalBsplines(uniqtm,knots)


p = size(Bmi)[2]
Dm1 = zeros(p-1, p)
for j=1:(p-1)
    Dm1[j,j] = 1
    Dm1[j,j+1] = -1
end

Dm2 = Dm1[1:p-2,1:p-1] * Dm1

Dml = 10*transpose(Dm2) * Dm2

Hmat = Bmi*inv(transpose(Bmi) * Bmi + 10*Dml) * transpose(Bmi)
yi = y[1:50]
Ini = diagm(0 => ones(50))
yh = (Ini - Hmat) * yi
transpose(yh) * yh/tr(Ini - Hmat)
