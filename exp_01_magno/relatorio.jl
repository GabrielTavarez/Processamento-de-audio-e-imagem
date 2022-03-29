### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 1434c786-af94-11ec-269f-f7938620ef5a
begin
	using DSP
	using Pkg
end

# ╔═╡ 3487e451-6344-40cc-9897-cfd4319d1ee0
begin
	M = 2
	N = 100
	n = 0:N-1
	σ² = 0.01
	s = sqrt(σ²)*rand(N)
	ϕ = 2π*rand()
	x = sin.(2π/10*n .+ π/6 .+ ϕ)
	u = 5sin.(2π/10*n .+ ϕ)
end

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


# ╔═╡ Cell order:
# ╠═1434c786-af94-11ec-269f-f7938620ef5a
# ╠═3487e451-6344-40cc-9897-cfd4319d1ee0
# ╠═d78b8006-066d-4c71-8e19-4284ffff805f
# ╠═91382807-6609-4b0a-ae5f-d337cb3edfbf
# ╠═bae96fd7-dcc7-4353-8f6a-9945d9d90ae6
