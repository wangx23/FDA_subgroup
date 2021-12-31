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
resultsim1init5 = readdlm("../result/resultsim1init5.csv",',')

mean(resultsim1[:,2])
mean(resultsim1lam0[:,2])
mean(resultsim1lam2[:,2])
mean(resultsim1lam5[:,2])
mean(resultsim1lam10[:,2])
mean(resultsim1init5[:,2])

mean(resultsim1init5[:,6])
mean(resultsim1lam10[:,6])



resultsim12lam0 = readdlm("../result/resultsim12lam0.csv",',')
resultsim12lam2 = readdlm("../result/resultsim12lam2.csv",',')
resultsim12lam5 = readdlm("../result/resultsim12lam5.csv",',')
resultsim12lam10 = readdlm("../result/resultsim12lam10.csv",',')
resultsim12init5 = readdlm("../result/resultsim12init5.csv",',')

mean(resultsim12lam0[:,2])
mean(resultsim12lam2[:,2])
mean(resultsim12lam5[:,2])
mean(resultsim12lam10[:,2])
mean(resultsim12init5[:,2])



mean(resultsim12lam0[:,5])
mean(resultsim12init5[:,6])


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



resultsim3 = readdlm("../result/resultsim3.csv",',')
mean(resultsim3[:,3])
mean(resultsim3[:,7])

mean(resultsim3[:,1])
mean(resultsim3[:,8])

resultsim3m20 = readdlm("../result/resultsim3m20.csv",',')

mean(resultsim3m20[:,3])
mean(resultsim3m20[:,7])

mean(resultsim3m20[:,1])
mean(resultsim3m20[:,5])

resultsim3m20ncl100 = readdlm("../result/resultsim3m20ncl100.csv",',')

mean(resultsim3m20ncl100[:,1])
mean(resultsim3m20ncl100[:,3])

mean(resultsim3m20ncl100[:,4])
mean(resultsim3m20ncl100[:,6])

resultsim3m15 = readdlm("../result/resultsim3m15.csv",',')

mean(resultsim3m15[:,3])
mean(resultsim3m15[:,7])

mean(resultsim3m15[:,1])
mean(resultsim3m15[:,5])


resultsim3m10 = readdlm("../result/resultsim3m10.csv",',')

mean(resultsim3m10[:,3])
mean(resultsim3m10[:,7])

mean(resultsim3m10[:,1])
mean(resultsim3m10[:,5])


resultsim3m10v2 = readdlm("../result/resultsim3m10v2.csv",',')

mean(resultsim3m10v2[:,3])
mean(resultsim3m10v2[:,7])

mean(resultsim3m10v2[:,1])
mean(resultsim3m10v2[:,5])


resultsim3m20ncl100 = readdlm("../result/resultsim3m20ncl100.csv",',')


resultsim3s = readdlm("../result/resultsim3s.csv",',')
mean(resultsim3s[:,1])
mean(resultsim3s[:,2])

mean(resultsim3s[:,3])
mean(resultsim3s[:,4])



resultsim3s30 = readdlm("../result/resultsim3s30.csv",',')
mean(resultsim3s30[:,1])
mean(resultsim3s30[:,2])

mean(resultsim3s30[:,3])
mean(resultsim3s30[:,4])

resultsim3s20 = readdlm("../result/resultsim3s20.csv",',')
mean(resultsim3s20[:,1])
mean(resultsim3s20[:,2])

mean(resultsim3s20[:,3])
mean(resultsim3s20[:,4])


resultsim3s15 = readdlm("../result/resultsim3s15.csv",',')
mean(resultsim3s15[:,1])
mean(resultsim3s15[:,2])

mean(resultsim3s15[:,3])
mean(resultsim3s15[:,4])


resultsim3s10 = readdlm("../result/resultsim3s10.csv",',')
mean(resultsim3s10[:,1])
mean(resultsim3s10[:,2])

mean(resultsim3s10[:,3])
mean(resultsim3s10[:,4])
