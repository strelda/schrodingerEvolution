print("precompiling libraries...")
using PlotlyJS #it's heavy, I just wanted to learn with it
using ProgressBars
println("done")

#0 for free particle, 1 for Gaussian
const switch=1
#if you just want to compare Gaussian without the potential, set testingMode=0
const testingMode=0

#edited thomas algorithm for matrix inversion, complexity O(n), numerically stable if diagonal elements are bigger than the rest of the matrix
function thomasUni(a::ComplexF64, b::Vector{ComplexF64}, d::Vector{ComplexF64}, n::Int64)
  x = copy(d)
  cp = Array{ComplexF64}(undef,n)
  for i in 1:n
    cp[i]=a
  end
  cp[1] /= b[1]
  x[1] /= b[1]

  for i = 2:n
    w = 1.0/(b[i]-cp[i-1]*a)
    cp[i] *= w
    x[i] = (x[i]-a*x[i-1])*w
  end
  for i = n-1:-1:1
    x[i] -= (cp[i]*x[i+1])
  end

  return x
end

function braket(P::Vector{ComplexF64},P0::Vector{ComplexF64},h::Float64) #Simpson rule for h even, O(h^4)
  n = length(P)
  I = (conj(P[1])*P0[1] + 4*conj(P[n-1])*P0[n-1] + conj(P[n])*P0[n])
  for i in 1:(round(Int,n/2)-1)
    I += 2*conj(P[2i-1])*P0[2i-1]
    I += 4*conj(P[2i])*P0[2i]
  end
  I *= h/3
  return I
end

function plotCompare(y::Vector{Float64}, yE::Vector{Float64})
  p1 = scatter(;y=y, mode="lines", name="ψ(t)")
  p2 = scatter(; y=yE, mode="lines", name="ψ(t) exact")
  p12() = plot([p1,p2])

  pErr = scatter(;y=y-yE, mode="lines", name="error = ψExact-ψ")
  we() = plot([pErr])

  return [p12(),we()]
end

function plotCompareReIm(yR::Vector{Float64}, yI::Vector{Float64})
  p1 = scatter(;y=yR, mode="lines", name="ℜ ψ(t)")
  p2 = scatter(; y=yI, mode="lines", name="ℑ ψ(t)")
  p12() = plot([p1,p2])

  return [p12()]
end


#some constants default: 2, 1.5, 0.5, 0, -5
const μ=1
const ω=1.5
const σ=0.5
const p0=0
const x0=-5

#quality, reasonable are: {10000,60,1e-2,50}, {10000,1000,1e-2,500}
const n=5000
const scale=10 #scale up to fit the wave on the screen
const dt=1e-2
const tExact=10
const h=-2*x0/n #centralizes the potential for gauss
const iter = round(Int,tExact/dt)



##exact solution 

#eigen
const k  = 0.1
const En = 5
exactψeigen(x::Float64,t::Float64) = exp(-im*En*t+im*k*x)
 
#Gauss
const X(t::Float64)  = p0*t/μ
const Σ2(t::Float64) = σ*σ + t*t/(4μ*μ*σ*σ)
const φ(x::Float64,t::Float64)= p0*(x-X(t)) + p0*p0*t/(2μ) + t*(x-X(t))^2/(8*μ*σ*σ*Σ2(t)) + angle(1/sqrt(μ+im*t/(2σ*σ)))
const exactψgauss(x::Float64,t::Float64)= (2π*Σ2(t))^(-1/4) * exp( -(x-X(t))^2/(4*Σ2(t)) +im* φ(x,t))



##functions for approximate solution
const V(x::Float64) = 0.5*μ*ω*ω*x*x
v  = Vector{Float64}(undef,n)
ψ0eigen = Vector{ComplexF64}(undef,n)
ψ0gauss = Vector{ComplexF64}(undef,n)


##initial conditions

#b is diagonal of tridiagonal matrix H, a are all offdiagonal elements
z = Vector{ComplexF64}(undef,n)
exactψgrid = Array{ComplexF64}(undef,n)


if switch==0
  #eigen
  b = Array{ComplexF64}(undef,n)
  const a = 0
  for i in 1:n
    x=scale*(x0+i*h)
    ψ0eigen[i] = exp(im*k*x)
    b[i] = En
  end
  ψ = copy(ψ0eigen)
  const ψ0= copy(ψ0eigen)
  exactψ = exactψeigen
else 
  #gauss
  for i in 1:n
    x=scale*(x0+i*h)
    if testingMode==1
      v[i]=0
    else
      v[i]=V(x)
    end 
    ψ0gauss[i] = (2π*σ*σ)^(-1/4)*exp(-x*x/(4σ*σ) + im*p0*x)
  end
  const a = -1/(2μ*h^2*scale^2)
  const b = 1/(μ*h^2)/scale^2 .+ v
  ψ = copy(ψ0gauss)
  const ψ0= copy(ψ0gauss)
  exactψ = exactψgauss
end


##iteration
const K = 0.5*im*dt
autocor = Array{ComplexF64}(undef,iter)
cnt=0
for tim in ProgressBar(1:iter)
  global ψ, z
  ψ[1]=0
  ψ[n]=0

  for j in 2:n-1
    z[j] = ψ[j] - K*(a*ψ[j-1]+b[j]*ψ[j]+a*ψ[j+1])
  end

  z[1] = ψ[1] - K*(b[1]*ψ[1] + a*ψ[2])  
  z[n] = ψ[n] - K*(a*ψ[n-1] + a*ψ[n])
  ψ = thomasUni(K*a,1 .+ K*b,z,n) # = inv(Id+0.5*im*dt*H)*z
  
  if mod(tim,round(Int,iter/100))==0
    for i in 1:n
      x=scale*(x0+i*h)
      global exactψgrid[i]=exactψ(x,tim*dt)
    end  
    if testingMode==1
      savefig(plotCompare(real.(ψ),real.(exactψgrid)), "plots/schrodinger"*lpad(cnt,4,"0")*".jpeg")
    else
      savefig(plotCompareReIm(real.(ψ),imag.(ψ)), "plots/schrodinger"*lpad(cnt,4,"0")*".jpeg")
    end
      global cnt += 1
  end

  global autocor[tim] = braket(ψ, ψ0, h)
end


#printing autoCorrelation function
pAutocor = scatter(;x=range(1,length=iter),y=real.(autocor), mode="lines", name="Autocorrelation function ℜ")
pAutocorIm = scatter(;x=range(1,length=iter),y=imag.(autocor), mode="lines", name="Autocorrelation function ℑ")

pAuto() = plot([pAutocor, pAutocorIm])
savefig(pAuto(), "autoCorrelation.jpeg")

#printing integral error, just for testing
if testingMode==1
  for i in 1:n
    x=scale*(x0+i*h)
    global exactψgrid[i]=exactψ(x,iter*dt)
  end
  plotCompare(real.(ψ),real.(exactψgrid))

  diff = (ψ-exactψgrid)[200:end-200]
  intDiff = real(braket(diff,diff,h))

  open("intDiff", "a") do io
    println(io, intDiff )
  end   
end
