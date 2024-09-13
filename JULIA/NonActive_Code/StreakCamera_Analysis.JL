using CSV, DataFrames, ColorSchemes,CairoMakie; CairoMakie.activate!()
using Optim
using Plots


#Import Data to Plot
file = "/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/LASERS/Femtosecond/DATA/2023_06_24/H_Chlorin_590nmPummp_10000avg.dac"
Streak_camera_data = (CSV.read(file,header = false,DataFrame))

function StreakCameraPlot(file)

    Data = (CSV.read(file,header = false,DataFrame))
    wl = collect(Data[1,2:end])
    t = parse.(Float64,collect(Data[2:end,1]))
    Mat = Matrix(Data[2:end,2:end])

    # Set start and end time for dataset
    uppertimebound = 15.0
    t0 = 5.5

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound
    # pick wavelength range for kinetics integration
    wlbin = 625 .< wl .< 725

    # Resize Matrix to new time domain
    Mat1 = Mat[cut,:]

    KineticTrace = sum(Mat1[:,wlbin],dims = 2)[:,1]

    Spectrum = sum(Mat1,dims = 1)[1,:]


    f = Figure(resolution=(600, 600))
    ax = CairoMakie.Axis(f[1, 1], xlabel = "Time (ns)",ylabel = "Wavelength (nm)")
    hm = CairoMakie.heatmap!(f[1,1],t,wl,Mat,colormap = :matter)
    cb = Colorbar(f[1, 2], hm)
    ax.yreversed=true



    # Integrated Spectrum
    ax2 = CairoMakie.Axis(f[2,1],title = "Spectrum",xlabel = "Wavelength (nm)",ylabel = "counts")
    lines!(ax2,wl,Spectrum, color = :blue, label = "Spectrum")

    # Integrated Kinetics
    ax3 = CairoMakie.Axis(f[3,1],title = "Kinetics", xlabel = "time (ns)", ylabel = "counts")
    lines!(ax3, t[cut], KineticTrace, color = :red, label  = "Kinetics")
    # fit = fitexp(t[cut],KineticTrace, n=1)
    # lines!(ax3,fit.x,fit.y, label = "fit")

    rowsize!(f.layout,1,200)
    rowsize!(f.layout,2,100)
    rowsize!(f.layout,3,100)
    resize_to_layout!(f)

    f
end

StreakCameraPlot(file)


# Kinetics Fitting using Optim

    Data = (CSV.read(file,header = false,DataFrame))
    wl = collect(Data[1,2:end])
    t = parse.(Float64,collect(Data[2:end,1]))
    Mat = Matrix(Data[2:end,2:end])

    # Set start and end time for dataset
    uppertimebound = 15.0
    t0 = 5.0

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound
    # pick wavelength range for kinetics integration
    wlbin = 625 .< wl .< 725

    # Resize Matrix to new time domain
    Mat1 = Mat[cut,:]

    KineticTrace = sum(Mat1[:,wlbin],dims = 2)[:,1]
    time = t[cut]

    mexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3]
    p0 = [3000, 11.0, 100.0]

    # bexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] 
    # p0 = [3000, 9.0,100.0, 15.0, -500.0]

    Plots.scatter(time, KineticTrace)
    Plots.plot!(time, mexp(p0))

    peakfit0(p) = sum((KineticTrace .- mexp(p)).^2)

    res = optimize(peakfit0, p0,BFGS(),Optim.Options(time_limit = 2000))
    
    values = Optim.minimizer(res)  # Kinetics results
    mexp(values)


    Plots.plot!(time, mexp(values), linewidth = 4)


  # Fit Plots  
  f = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
  width = 2

  ax1 = CairoMakie.Axis(f[1,1], title = "Fit",
  xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
  ylabel =  "Intensity", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
  # xlims!(ax1,350,650)
  # ylims!(nothing,)
  
  lines!(f[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
  lines!(f[1,1],time, mexp(values), linewidth = 2, label = "Fit", color = :red)
  lines!(f[1,1],time,mexp(p0), linewidth = 2, label = "Guess", color = :green)
  axislegend(ax1, position = :rt, framevisible = false)

  ax2 = CairoMakie.Axis(f[2,1], title = "Residuals",
  xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
  ylabel =  "Intensity", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

  lines!(f[2,1],time, KineticTrace .- mexp(values), linewidth = width, color = :black)

  colsize!(f.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
  rowsize!(f.layout, 1, 300)
  rowsize!(f.layout,2,100)
  resize_to_layout!(f)

f

MonoFitIRF_StreakCamera(Streak_camera_data,3.0,12.0,[2128.23, 3.413, 0.961, 39.392, 10.231, 1.728],true,false)




function MonoFitIRF_StreakCamera(file,t0,uppertimebound,coeffs, residuals, ReturnFit)

  function Fitfunc(x,coeffs)
  
      Gaussian(p) = exp.(-4 .* log(2) .* (((x .- p[2]) ./ p[3]) .^2))
      mexp(p) = (p[4] .* exp.((-(x) ./ p[5])) .+ p[6])
  
      IRF = Gaussian(coeffs)
      Exp = mexp(coeffs)
  
      function DirectLinearConvolution(f,g)
          N = length(f) # for This to work both f and g need to be the same size
          M = length(g)
      
          Conv = zeros(N)
          for n=1:N
              for m=1:M
                  if n-m+1 > 0
                      Conv[n] = Conv[n] + f[m] * g[n-m+1]
                  end
                  # n+1 <= m
              end
          end
          return Conv
      end
  
      Convolution = DirectLinearConvolution(IRF,Exp) * (coeffs[1]/maximum(DirectLinearConvolution(IRF,Exp)))
  
      return Convolution
  end



  Data = file
  wl = collect(Data[1,2:end])
  t = parse.(Float64,collect(Data[2:end,1]))
  Mat = Matrix(Data[2:end,2:end])
  
  # Create Bit-vectors for desrired time range and wavelength range for integration
  wlbin = 625 .< wl .< 725
  cut = t0 .< t .< uppertimebound
  I = Mat[cut,:]

  KineticTrace = sum(I[:,wlbin],dims = 2)[:,1]
  time = t[cut]


  ## FITTING FUNCTION
  peakfit0(coeffs) = sum((KineticTrace .- Fitfunc(time,coeffs)).^2)


  res = optimize(peakfit0, coeffs,LBFGS(),Optim.Options(time_limit = 200))
  print("Fit Resuluts:", res)
  values = Optim.minimizer(res) # Kinetics results
  residual = KineticTrace .- Fitfunc(time,values)
  print("Fit Coeff:", values)

  df = DataFrame(Time = time,RawData = KineticTrace, Fit = Fitfunc(time,values),Residual = residual)

  if ReturnFit == true
      return df
  
  elseif residuals == true
      Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
      width = 2

      ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
      xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
      ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
      # xlims!(ax1,350,650)
      # ylims!(nothing,)

      lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)      
      lines!(Fit_Fig[1,1],time,Fitfunc(time,coeffs), linewidth = 2, label = "Guess", color = :green)
      lines!(Fit_Fig[1,1],time, Fitfunc(time,values), linewidth = 2, label = "Fit", color = :red)

      axislegend(ax1, position = :rt, framevisible = false)

      ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
      xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
      ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

      lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

      colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
      rowsize!(Fit_Fig.layout, 1, 300)
      rowsize!(Fit_Fig.layout,2,100)
      resize_to_layout!(Fit_Fig)

      Fit_Fig
      # save("JULIA/Outputs/MonoFit_Fig.png",Fit_Fig)
  else

      Plots.plot(time, KineticTrace, color = :blue, label = "Data")
      Plots.plot!(time, Fitfunc(time,coeffs), color = :green, label = "Guess")
      Plots.plot!(time, Fitfunc(time,values), linewidth = 3, color = :red, label = "Fit")

  end

  # return fitfunction(values)
end

function BiFit(file,column,t0,uppertimebound, residuals, ReturnFit)

  Data = file
  t = Data[:,1]
  I = Data[:,column]

  # Create Bit-vectors for desrired time range and wavelength range for integration
  cut = t0 .< t .< uppertimebound

  time = t[cut]
  KineticTrace = I[cut]

  ## FITTING FUNCTION
  fitfunction = bexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] 
  p0 = [0.01, 9.0, 0.02, 15.0, 0.0]

  peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)

  res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
  print("Fit Resuluts:", res)
  values = Optim.minimizer(res) # Kinetics results
  residual = KineticTrace .- fitfunction(values)
  print("Fit Coeff:", values)

  if ReturnFit == true
      return hcat(time,KineticTrace,fitfunction(values),residual)
  

  elseif residuals == true
      Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
      width = 2

      ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
      xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
      ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
      # xlims!(ax1,350,650)
      # ylims!(nothing,)

      lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
      lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "Fit", color = :red)
      lines!(Fit_Fig[1,1],time,fitfunction(p0), linewidth = 2, label = "Guess", color = :green)
      axislegend(ax1, position = :rt, framevisible = false)

      ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
      xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
      ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

      lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

      colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
      rowsize!(Fit_Fig.layout, 1, 300)
      rowsize!(Fit_Fig.layout,2,100)
      resize_to_layout!(Fit_Fig)

      Fit_Fig
      # save("JULIA/Outputs/BiFit_Fig.png",Fit_Fig)
  else

      Plots.plot(time, KineticTrace, color = :blue, label = "Data")
      Plots.plot!(time, fitfunction(p0), color = :green, label = "Guess")
      Plots.plot!(time, fitfunction(values), linewidth = 3, color = :red, label = "Fit")

  end

  # return fitfunction(values)
end

function BiFitIRF(file,column,t0,uppertimebound,coeffs, residuals, ReturnFit, reduce)

  function Fitfunc(x,coeffs)
  
      Gaussian(p) = exp.(-4 .* log(2) .* (((x .- p[2]) ./ p[3]) .^2))
      Biexp(p) = (p[4] .* exp.((-(x) ./ p[5])) .+ (p[6] .* exp.((-(x) ./ p[7]))) .+ p[8])
  
      IRF = Gaussian(coeffs)
      Exp = Biexp(coeffs)
  
      function DirectLinearConvolution(f,g)
          N = length(f) # for This to work both f and g need to be the same size
          M = length(g)
      
          Conv = zeros(N)
          for n=1:N
              for m=1:M
                  if n-m+1 > 0
                      Conv[n] = Conv[n] + f[m] * g[n-m+1]
                  end
                  # n+1 <= m
              end
          end
          return Conv
      end
  
      Convolution = DirectLinearConvolution(IRF,Exp) * (coeffs[1]/maximum(DirectLinearConvolution(IRF,Exp)))
  
      return Convolution
  end

  Data = file

  rebin = round(length(Data[:,1])/reduce)
  stepsize = (Data[end,1]-Data[1,1])/rebin
  InterpData = LinearInterpolation(Data[:,1],Data[:,2])
  
  t = collect(Data[1,1]:stepsize:Data[end,1])
  I = InterpData(t)

  # Create Bit-vectors for desrired time range and wavelength range for integration
  cut = t0 .< t .< uppertimebound

  time = t[cut]
  KineticTrace = I[cut]

  ## FITTING FUNCTION
  peakfit0(coeffs) = sum((KineticTrace .- Fitfunc(time,coeffs)).^2)


  res = optimize(peakfit0, coeffs,LBFGS(),Optim.Options(time_limit = 200))
  print("Fit Resuluts:", res)
  values = Optim.minimizer(res) # Kinetics results
  residual = KineticTrace .- Fitfunc(time,values)
  print("Fit Coeff:", values)

  df = DataFrame(Time = time,RawData = KineticTrace, Fit = Fitfunc(time,values),Residual = residual)

  if ReturnFit == true
      return df
  
  elseif residuals == true
      Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
      width = 2

      ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
      xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
      ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
      # xlims!(ax1,350,650)
      # ylims!(nothing,)

      lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
      lines!(Fit_Fig[1,1],time, Fitfunc(time,values), linewidth = 2, label = "Fit", color = :red)
      lines!(Fit_Fig[1,1],time,Fitfunc(time,coeffs), linewidth = 2, label = "Guess", color = :green)
      axislegend(ax1, position = :rt, framevisible = false)

      ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
      xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
      ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

      lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

      colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
      rowsize!(Fit_Fig.layout, 1, 300)
      rowsize!(Fit_Fig.layout,2,100)
      resize_to_layout!(Fit_Fig)

      Fit_Fig
      # save("JULIA/Outputs/MonoFit_Fig.png",Fit_Fig)
  else

      Plots.plot(time, KineticTrace, color = :blue, label = "Data")
      Plots.plot!(time, Fitfunc(time,coeffs), color = :green, label = "Guess")
      Plots.plot!(time, Fitfunc(time,values), linewidth = 3, color = :red, label = "Fit")

  end

  # return fitfunction(values)
end

function TriFit(file,column,t0,uppertimebound, residuals, ReturnFit)

  Data = file
  t = Data[:,1]
  I = Data[:,column]

  # Create Bit-vectors for desrired time range and wavelength range for integration
  cut = t0 .< t .< uppertimebound

  time = t[cut]
  KineticTrace = I[cut]

  ## FITTING FUNCTION
  fitfunction = triexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] .* exp.((-time) ./ p[6]).+ p[7] 
  p0 = [0.01, 2.19, 0.01, 8.7, 0.01, 19.0, 0.01]

  peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)

  res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
  print("Fit Resuluts:", res)
  values = Optim.minimizer(res) # Kinetics results
  residual = KineticTrace .- fitfunction(values)
  print("Fit Coeff:", values)

  df = DataFrame(Time = time,RawData = KineticTrace, Fit = Fitfunc(time,values),Residual = residual)

  if ReturnFit == true
      return df
  

  elseif residuals == true
      Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
      width = 2

      ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
      xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
      ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
      # xlims!(ax1,350,650)
      # ylims!(nothing,)

      lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
      lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "Fit", color = :red)
      lines!(Fit_Fig[1,1],time,fitfunction(p0), linewidth = 2, label = "Guess", color = :green)
      axislegend(ax1, position = :rt, framevisible = false)

      ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
      xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
      ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

      lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

      colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
      rowsize!(Fit_Fig.layout, 1, 300)
      rowsize!(Fit_Fig.layout,2,100)
      resize_to_layout!(Fit_Fig)

      Fit_Fig
      # save("JULIA/Outputs/TriFit_Fig.png",Fit_Fig)
  else

      Plots.plot(time, KineticTrace, color = :blue, label = "Data")
      Plots.plot!(time, fitfunction(p0), color = :green, label = "Guess")
      Plots.plot!(time, fitfunction(values), linewidth = 3, color = :red, label = "Fit")

  end

  # return fitfunction(values)
end