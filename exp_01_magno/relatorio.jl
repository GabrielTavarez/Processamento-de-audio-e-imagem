### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 1434c786-af94-11ec-269f-f7938620ef5a
begin
	using DSP
	using Pkg
	using Plots
	using LinearAlgebra
	using Statistics
	plotly()
end

# ╔═╡ f71662ab-bd42-4174-b5c1-bb67939b589e


# ╔═╡ 3487e451-6344-40cc-9897-cfd4319d1ee0
begin
	M = 2
	N = 500
	n = 0:N-1
	σ² = 0.01
	s = sqrt(σ²)*rand(N)
	ϕ = 2π*rand()
	x = sin.(2π/10*n .+ π/6 .+ ϕ)
	u = 5sin.(2π/10*n .+ ϕ)
end

# ╔═╡ 3f93a713-4850-49cd-a303-c43911929f6b
md"Letra A"

# ╔═╡ d78b8006-066d-4c71-8e19-4284ffff805f
begin
	R =[25/2 25/2*cos(2π/10) ;
		25/2*cos(2π/10) 25/2]
	
	p = [5/2*cos(2π/12) ;
		5/2*cos(2π*11/60)]
end

# ╔═╡ 91382807-6609-4b0a-ae5f-d337cb3edfbf
begin
	wₒ = inv(R)*p
end

# ╔═╡ bae96fd7-dcc7-4353-8f6a-9945d9d90ae6
begin
	#Jmin
	σd² = σ² + 0.5
	Jmin = σd² - p'*wₒ
end

# ╔═╡ cd4eb269-73a7-437a-bff7-2d453558b925
function pow2db(x)
	20*log10(x)
end

# ╔═╡ aa8689c2-59ea-472e-86f8-d22067545db6
begin
	ω = range(0, π, length=500)
	filtro = PolynomialRatio(wₒ, [1])
	resp = freqz(filtro, ω)
	plot(ω./π,pow2db.(real(resp)))
end

# ╔═╡ 288bdde9-4dba-4b38-b718-3ee2bc29c5d9
md"letra b"

# ╔═╡ 1784753b-4117-485f-b3ff-b3c95c9844ec
begin
	λ = eigvals(R)
	λmax = maximum(λ)
	μmax = 2/λmax
end

# ╔═╡ 87a2deb2-179d-4ff4-ac56-fee3b07cc9b8
md"letra c"

# ╔═╡ 1986fd1b-30cd-4db2-9376-ca6afab1ba1e
function LMS(x, d, M, N, μ)
	X = zeros(M,1)
	W = zeros(N,M)
	erro = zeros(N,1)
	
	for n in 1:N-1
		X = [x[n]; X[1:M-1]]
		y = W[n,:]'*X
		erro[n] = d[n]-y
		W[n+1, :] = W[n,:] +μ*erro[n]*X
	end
	return W, erro
end

# ╔═╡ 01095d12-b298-4e86-ac7d-6ef429f44861
begin
	d = s+x
	W_LMS, erro = LMS(u, d, M, N, 0.03)
end

# ╔═╡ 28b2e8a9-60bc-4597-996e-673fbbe5d063
plot(W_LMS[1,:],W_LMS[2,:], seriestype = :scatter)

# ╔═╡ 15aebe76-9e6d-4085-bd2c-a8d5d25ea111
W_LMS

# ╔═╡ 4cddabf9-412a-4290-8210-3b7504c7d134
p

# ╔═╡ Cell order:
# ╠═1434c786-af94-11ec-269f-f7938620ef5a
# ╠═f71662ab-bd42-4174-b5c1-bb67939b589e
# ╠═3487e451-6344-40cc-9897-cfd4319d1ee0
# ╠═3f93a713-4850-49cd-a303-c43911929f6b
# ╠═d78b8006-066d-4c71-8e19-4284ffff805f
# ╠═91382807-6609-4b0a-ae5f-d337cb3edfbf
# ╠═bae96fd7-dcc7-4353-8f6a-9945d9d90ae6
# ╠═aa8689c2-59ea-472e-86f8-d22067545db6
# ╠═cd4eb269-73a7-437a-bff7-2d453558b925
# ╠═288bdde9-4dba-4b38-b718-3ee2bc29c5d9
# ╠═1784753b-4117-485f-b3ff-b3c95c9844ec
# ╠═87a2deb2-179d-4ff4-ac56-fee3b07cc9b8
# ╠═1986fd1b-30cd-4db2-9376-ca6afab1ba1e
# ╠═01095d12-b298-4e86-ac7d-6ef429f44861
# ╠═28b2e8a9-60bc-4597-996e-673fbbe5d063
# ╠═15aebe76-9e6d-4085-bd2c-a8d5d25ea111
# ╠═4cddabf9-412a-4290-8210-3b7504c7d134
