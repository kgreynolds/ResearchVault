using Plots, DataFrames, CSV, CurveFit,Peaks
plotly()

WL_Nofilter = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/LASERS/Femtosecond/ADMIN/Femtosecond Calibration/WhiteLightCalibration_2023_07_09/WL.txt", DataFrame)
Plots.plot(WL_Nofilter[:,1], WL_Nofilter[:,2], size = (900,600))

WLpks, vals = findmaxima(WL_Nofilter[:,2])
WLpks, proms = peakproms!(WLpks, WL_Nofilter[:,2], minprom = 25)
plotpeaks!(WL_Nofilter[:,1], WL_Nofilter[:,2], peaks=WLpks)

WLpks,vals = findminima(WL_Nofilter[:,2])
WLpks, proms = peakproms!(WLpks, WL_Nofilter[:,2], minprom = 25)
plotpeaks!(WL_Nofilter[:,1], WL_Nofilter[:,2], peaks=WLpks)


WL_filter = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/LASERS/Femtosecond/ADMIN/Femtosecond Calibration/WhiteLightCalibration_2023_07_09/WL_Filter.txt", DataFrame)
Plots.plot!(WL_filter[:,1],WL_filter[:,2])

WLFilterpks, vals = findmaxima(WL_filter[:,2])
WLFilterpks, proms = peakproms!(WLFilterpks, WL_filter[:,2], minprom = 25)
plotpeaks!(WL_filter[:,1], WL_filter[:,2], peaks=WLFilterpks)

WLFilterpks, vals = findminima(WL_filter[:,2])
WLFilterpks, proms = peakproms!(WLFilterpks, WL_filter[:,2], minprom = 25)
plotpeaks!(WL_filter[:,1], WL_filter[:,2], peaks=WLFilterpks)




WL_Ocean = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/LASERS/Femtosecond/ADMIN/Femtosecond Calibration/WhiteLightCalibration_2023_07_09/WL_Ocean.txt", DataFrame)
Plots.plot(WL_Ocean[:,1],WL_Ocean[:,2], xflip = true, size = (900,600))

WLOceanPKS, vals = findmaxima(WL_Ocean[:,2])
WLOceanPKS, proms = peakproms!(WLOceanPKS, WL_Ocean[:,2], minprom = 25)
plotpeaks!(WL_Ocean[:,1],WL_Ocean[:,2], peaks=WLOceanPKS)


WL_Ocean_Filter = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/LASERS/Femtosecond/ADMIN/Femtosecond Calibration/WhiteLightCalibration_2023_07_09/WL_OceanFilter.txt", DataFrame)
Plots.plot!(WL_Ocean_Filter[:,1],WL_Ocean_Filter[:,2], xflip = true)

WLOceanFilterPKS, vals = findmaxima(WL_Ocean_Filter[:,2])
WLOceanFilterPKS, proms = peakproms(WLOceanFilterPKS,WL_Ocean_Filter[:,2], minprom = 25)
plotpeaks!(WL_Ocean_Filter[:,1],WL_Ocean_Filter[:,2], peaks=WLOceanFilterPKS)

WLOceanFilterPKS, vals = findminima(WL_Ocean_Filter[:,2])
WLOceanFilterPKS, proms = peakproms(WLOceanFilterPKS,WL_Ocean_Filter[:,2], minprom = 25)
plotpeaks!(WL_Ocean_Filter[:,1],WL_Ocean_Filter[:,2], peaks=WLOceanFilterPKS)