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
	using Noise
	plotly()
end

# ╔═╡ 3487e451-6344-40cc-9897-cfd4319d1ee0
begin
	M = 2
	N = 500
	n = 0:N-1
	σ² = 0.01
	s = add_gauss(zeros(N),σ²) 
	#s = sqrt(σ²)*rand(N) .- sqrt(σ²)/2
	ϕ = 2π*rand()
	f_int = 1/10
	x = sin.(2π/10*n .+ π/6 .+ ϕ)
	u = 5sin.(2π/10*n .+ ϕ)
end

# ╔═╡ 3f93a713-4850-49cd-a303-c43911929f6b
md"
#### Letra A

- o vetor de coeficientes ótimos wo
- o erro quadrático médio mínimo
- a resposta em frequência do filtro ótimo na frequência da interferência ecompare com x(n) e u(n).
"

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
	Jmin = σd² - wₒ'*p
end

# ╔═╡ aa8689c2-59ea-472e-86f8-d22067545db6
begin
	#Resposta em frequência do filtro ideal
	ω = range(0, π, length=500)
	filtro = PolynomialRatio(wₒ, [1])
	resp = freqz(filtro, ω)	
	resp_int = freqz(filtro, f_int)
end

# ╔═╡ 288bdde9-4dba-4b38-b718-3ee2bc29c5d9
md"
#### Letra b

Com a matriz R e a função eig.m do Matlab, calcule a faixa de valores do passo
de adaptação µ que garante a convergência do algoritmo Steepest Descent.


"

# ╔═╡ 1784753b-4117-485f-b3ff-b3c95c9844ec
begin
	λ = eigvals(R)
	λmax = maximum(λ)
	μmax = 2/λmax
end

# ╔═╡ 87a2deb2-179d-4ff4-ac56-fee3b07cc9b8
md"
### Letra c

Aplique o algoritmo LMS com µ = 0,03 e N = 500 iterações.


- observe inicialmente os sinais de entrada u(n), de erro e(n) e s(n) em gráficos na mesma escala;
- compare os coeficientes do filtro adaptativo com os coeficientes ótimos calculados no item a), fazendo um gráfico dos coeficientes ao longo das iterações;
- trace as curvas de nível da superfície de erro e sobre elas, a trajetória dos coeficientes;
- trace a curva do erro quadrático e²(n) em dB.

"

# ╔═╡ 4cddabf9-412a-4290-8210-3b7504c7d134
md"
#### Letra d

Determine experimentalmente o valor máximo de µ para convergência do algoritmo LMS e compare-o com o valor calculado no item b) para o algoritmo Steepest Descent.
"

# ╔═╡ 029312c5-192a-420a-ac85-83a1dbd82711
md" 
#### Letra E

 Obtenha uma aproximação para J(n) = E{e²(n)} considerando uma média de500 realizações de e2(n). Note que em cada realização, um novo valor de φ e um novo s(n) devem ser considerados. Pede-se:
"

# ╔═╡ e8bd4fac-df94-4025-922b-27f7567e75d1
md" ## Functions "

# ╔═╡ cd4eb269-73a7-437a-bff7-2d453558b925
function pow2db(x)	
	20*log10(x)
end

# ╔═╡ 3a309359-d21a-4a75-824c-86b14901b4d8
begin
	plot(ω./π,pow2db.(abs.(resp)), label = "Resp Freq")
	plot!([f_int/π], pow2db.(abs.([resp_int])), marker = :circle, label = "Freq interferência") 
	plot!(title = "Resposta em frequência do filtro")
	
end

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
	#Erro e coeficientes ao longo do tempo
	d = s+x
	W_LMS, erro = LMS(u, d, M, N, 0.03)
end

# ╔═╡ 670064be-021b-42c5-8e6f-3c06eab831cd
begin
	plot(n, u, label = "u")
	plot!(n,erro, label = "e")
	plot!(n,s, label = "s")
	
	plot!(title = "Sinais do filtro")
end

# ╔═╡ 28b2e8a9-60bc-4597-996e-673fbbe5d063
begin
	plot(W_LMS[:,1],W_LMS[:,2], marker = :circle, label = "Filtro real", markersize = 2)
	#plot!(wₒ[1],wₒ[2])
	plot!([wₒ[1],wₒ[1]],[wₒ[2],wₒ[2]], seriestype = :scatter, color = "red", markersize=5, label = "Filtro ideal")
	plot!(title="Coeficientes do filtro ao longo do tempo")
end

# ╔═╡ 15aebe76-9e6d-4085-bd2c-a8d5d25ea111
begin
	#Erro quadratico médio
	plot(erro.^2)
	plot!(title="Erro quadrático médio ao longo do tempo")
end

# ╔═╡ 67bd5e06-164f-4ca7-8a69-21bb1b2ca3d6
begin
	μ_max_experimental = 0.08
	plot(LMS(u, d, M, N, μ_max_experimental+0.01)[2].^2, label="μ = 0.09", alpha=0.6)
	plot!(LMS(u, d, M, N, μ_max_experimental)[2].^2, label="μ = 0.8", width = 2)
	plot!(ylim = [0,10], title = "Erro quadrático médio para diferentes μ")
end

# ╔═╡ efc7d7d7-dd9f-447b-832f-04ee5d1a1bc1
begin
	realizacoes = 1:500
	erro_realizacoes_001 = zeros(N)
	erro_realizacoes_003 = zeros(N)
	erro_realizacoes_005 = zeros(N)
	
	
	for realizaco in realizacoes
		σ² = 0.01
		s = sqrt(σ²)*rand(N) .- sqrt(σ²)/2
		ϕ = 2π*rand()
		f_int = 1/10
		x = sin.(2π/10*n .+ π/6 .+ ϕ)
		u = 5sin.(2π/10*n .+ ϕ)
		
		W_LMS, erro_001 = LMS(u, d, M, N, 0.01)
		W_LMS, erro_003 = LMS(u, d, M, N, 0.03)
		W_LMS, erro_005 = LMS(u, d, M, N, 0.05)
		
		
		erro_realizacoes_001 += erro_001.^2
		erro_realizacoes_003 += erro_003.^2
		erro_realizacoes_005 += erro_005.^2
		
	end
	erro_realizacoes_001 = erro_realizacoes_001./realizacoes
	erro_realizacoes_003 = erro_realizacoes_003./realizacoes
	erro_realizacoes_005 = erro_realizacoes_005./realizacoes
	
end

# ╔═╡ c5a83862-062a-4472-80ff-8a291a122fa1
begin
	plot(erro_realizacoes_001, ylim = [0,10])	
	plot!(erro_realizacoes_003, ylim = [0,10])
	plot!(erro_realizacoes_005, ylim = [0,10])	
end

# ╔═╡ 534e7124-9889-443d-aef4-0d7d82d35019


# ╔═╡ Cell order:
# ╠═1434c786-af94-11ec-269f-f7938620ef5a
# ╠═3487e451-6344-40cc-9897-cfd4319d1ee0
# ╠═3f93a713-4850-49cd-a303-c43911929f6b
# ╠═d78b8006-066d-4c71-8e19-4284ffff805f
# ╠═91382807-6609-4b0a-ae5f-d337cb3edfbf
# ╠═bae96fd7-dcc7-4353-8f6a-9945d9d90ae6
# ╠═aa8689c2-59ea-472e-86f8-d22067545db6
# ╠═3a309359-d21a-4a75-824c-86b14901b4d8
# ╟─288bdde9-4dba-4b38-b718-3ee2bc29c5d9
# ╠═1784753b-4117-485f-b3ff-b3c95c9844ec
# ╠═87a2deb2-179d-4ff4-ac56-fee3b07cc9b8
# ╠═01095d12-b298-4e86-ac7d-6ef429f44861
# ╠═670064be-021b-42c5-8e6f-3c06eab831cd
# ╠═28b2e8a9-60bc-4597-996e-673fbbe5d063
# ╠═15aebe76-9e6d-4085-bd2c-a8d5d25ea111
# ╟─4cddabf9-412a-4290-8210-3b7504c7d134
# ╠═67bd5e06-164f-4ca7-8a69-21bb1b2ca3d6
# ╟─029312c5-192a-420a-ac85-83a1dbd82711
# ╠═efc7d7d7-dd9f-447b-832f-04ee5d1a1bc1
# ╠═c5a83862-062a-4472-80ff-8a291a122fa1
# ╟─e8bd4fac-df94-4025-922b-27f7567e75d1
# ╠═cd4eb269-73a7-437a-bff7-2d453558b925
# ╠═1986fd1b-30cd-4db2-9376-ca6afab1ba1e
# ╠═534e7124-9889-443d-aef4-0d7d82d35019
