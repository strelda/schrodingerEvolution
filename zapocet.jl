print("precompiling libraries...")
using PlotlyJS #it's heavy, I know...did it for the fun of learning with it
using ProgressBars
println("done")

#edited thomas algorithm for matrix inversion, complexity O(n), numerically stable for bigger than the rest of the matrix
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

function braket(P::Vector{ComplexF64},P0::Vector{ComplexF64},h::Float64) #only for h even
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


#some constants
const μ=1
const ω=1.5
const σ=0.5
const p0=0
const x0=-5

#quality
const n=10000
const scale=60 #scale up to fit the wave on the screen
const dt=1e-2
const h=-2*x0/n #centralizes the potential
const tExact=50
const iter = round(Int64,tExact/dt)

#exact solution
X(t::Float64)  = p0*t/μ
Σ2(t::Float64) = σ*σ + t*t/(4μ*μ*σ*σ)
φ(x::Float64,t::Float64)= p0*(x-X(t)) + p0*p0*t/(2μ) + t*(x-X(t))^2/(8*μ*σ*σ*Σ2(t)) + angle(1/sqrt(μ+im*t/(2σ*σ)))
exactψ(x::Float64,t::Float64)= (2π*Σ2(t))^(-1/4) * exp( -(x-X(t))^2/(4*Σ2(t)) +im* φ(x,t))
exactψgrid = Array{ComplexF64}(undef,n)


#functions for approximate solution
V(x::Float64) = 0.5*μ*ω*ω*x*x
v  = Vector{Float64}(undef,n)
ψ0 = Vector{ComplexF64}(undef,n)

for i in 1:n
  x=scale*(x0+i*h)
  v[i]=V(x)
  ψ0[i] = (2π*σ*σ)^(-1/4)*exp(-x*x/(4σ*σ) + im*p0*x)
end

#b is diagonal of tridiagonal matrix H, a are all offdiagonal elements
b = (1/(μ*h^2) .+ v)/scale^2
a = -1/(2μ*h^2*scale^2)

autocor = Array{ComplexF64}(undef,iter)

#initial conditions
ψ = ψ0
K = 0.5*im*dt
z = Vector{ComplexF64}(undef,n)


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

    savefig(plotCompare(real.(ψ),real.(exactψgrid)), "plots/schrodinger"*lpad(count,4,"0")*".jpeg")

  end
  autocor[tim] = braket(ψ, ψ0, h)
end



# pp1=plot(real.(ψ0), name="ψ(t=0)")
# pp2=plot(v, name="potential v")

# pAutocor = scatter(;x=range(1,length=iter),y=real.(autocor), mode="lines", name="Autocorrelation function ℜ")
# pAutocorIm = scatter(;x=range(1,length=iter),y=imag.(autocor), mode="lines", name="Autocorrelation function ℑ")
# pAuto() = plot([pAutocor, pAutocorIm])

# savefig(pAuto(), "autoCorrelation.jpeg")

# for i in 1:n
#   x=scale*(x0+i*h)
#   global exactψgrid[i]=exactψ(x,50.) 
# end

# y = real.(ψ)
# yE = real.(exactψgrid)

# p1 = scatter(;y=y, mode="lines", name="ψ(t)")
# p2 = scatter(; y=yE, mode="lines", name="ψ(t) exact")
# p12() = plot([p1,p2])

# pErr = scatter(;y=(y-yE)/maximum(yE), mode="lines", name="(ψExact-ψ)/max(ψExact)")
# we() = plot([pErr])


# [p12(),we(),pAuto()]


