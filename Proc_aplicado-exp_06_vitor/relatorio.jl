### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 680e5c90-fc8e-11ec-1304-0deea9b53834
begin
	using MAT
	using DSP
	using Statistics
	using Plots
	plotly()
end

# ╔═╡ fb2fa972-a156-4980-9d15-0be562958957
md"# Arranjos de antenas e ganho espacial"

# ╔═╡ 12c8cb4e-bc06-44c3-975b-74d1b9323e29
md" ## 01 - Projeto de antena DAS

Projeto de arranjo de antena DAS miradas para 20º usando os coeficientes de cada antena definida por:

$a_m = \frac{1}{M}e^{j\Omega\tau(\theta)}$

E com os coeficientes das antenas definimos o ganho espacial como 

$B(\theta) = \sum a_m e^{-j\Omega\tau(\theta)}$

"

# ╔═╡ e6d1fc78-8723-418d-be3b-8eddc3430be0
md" ## 02 - Diferentes ganhos

Vamos calcular o ganho das antenas para diferentes valores de distâncias.

Distâncias : λ/4,  λ/2,  3λ/4,  λ 
"

# ╔═╡ 539fca48-8236-4734-a7bc-4049bc8d347b
md" ## 03

Dado um sinal x1 vindo de θ1 conhecido e um sinal x2 vindo de θ2, queremos saber a amplitude máxima do sinal a depender da direção θ2 da interferência"

# ╔═╡ a9338f7b-e329-4bd5-a576-5da4896b6ba4
md" ## 1 - Processamento DAS do sinal "

# ╔═╡ 138691ba-b2ec-4bf8-ac09-0d0cc8ac3bcb
md" # Sinal Real"

# ╔═╡ f9f27eba-5e67-4f7e-bb87-932a4273db0c
md"#### Filtro Passa-Baixas

Implementação do filtro Passa-Biaxas com frequÊncia de corte igual a portadora.

Usaremos um filtro Butterworth de 6 coeficientes.
"

# ╔═╡ 74784532-bf4a-4d84-8bf1-3f1b279a8567
md" #### Divisão dos sinais em componente imaginária e real"

# ╔═╡ 41ce1800-51b8-4424-a01f-2d7ebf5d2117
md" ## 2 - Projeto dos coeficientes das antenas"

# ╔═╡ 7e4fb1df-3d37-4290-b0d3-3b2f985ae234
md" ## 3 - Decodificação"

# ╔═╡ 6e695e78-b2db-4128-bdca-fe244d4667f5
matread("mensagem.mat")["mensagem"]

# ╔═╡ 58b09e7c-9bd0-459e-8077-97c020758932
md" ## Functions"

# ╔═╡ f41cf0d5-af4b-49ee-b9ce-917cd6e57866
function B(θin, θ0in, M, d, Ω0)
	θ = deg2rad(θin)
	θ0 = deg2rad(θ0in)
	c = 3e8
	exponencial = exp( im*(M-1)/2 * Ω0*d/c *( sin(θ)-sin(θ0) ) )
	if (sin(θ)-sin(θ0))==0.0
		senos = M
	else
		senos = sin(M*Ω0*d/(2c) * (sin(θ)-sin(θ0)) ) / ( sin( Ω0*d/(2c) * (sin(θ)-sin(θ0) )) )
	end
	valor =exponencial *1/M* senos
	return valor
end

# ╔═╡ 72f16a75-fc9d-442c-be94-f443db657efa
function τ(θ, m, d; c = 3e8)
	return m*d*sind(θ)/c
end

# ╔═╡ 199b6831-1cf9-44fa-9afc-f1bc838ce090
function B(a, θ::Vector{Float64}, d ;c = 3e8)

	N = length(θ) #numero de graus do ganho
	B_ = zeros(Complex, N)
	
	for i in 1:N
		B_[i] = B(a, θ[i], d; c = c)	
	end

	return B_
end

# ╔═╡ edc6e700-01ee-4fa9-8eec-b98d323e1fa2
function B(a, θ::StepRangeLen{Float64}, d; c = 3e8)
	return B(a, collect(θ),d; c=c)
end

# ╔═╡ 36ef6932-be01-4760-a104-6764d4746787
noprint = md""

# ╔═╡ 4d2fb09b-36a4-474b-a908-05b174b68c4f
begin
	M = 8
	F = 60e9 #Hz
	Ω = 2*π*F
	c = 3e8 #m/s
	λ = c/F #m
	j = im 
	θ0 =  20 #graus
	noprint
end

# ╔═╡ 9562b052-e711-4aaa-8aec-072f83a70de3
F

# ╔═╡ 9ea61348-08e5-4070-812b-a013478afeb8
function B(a, θ, d ;c = 3e8)
	M = length(a)
	B_ = 0 + 0*j
	
	for m in 1:M
		B_ += a[m]*exp(-j*Ω*τ(θ,m, d; c = c))
	end
	return B_
end

# ╔═╡ fecad3c1-f01d-49a0-9141-77577f29c6e3
function DAS(θ, M, d; Ω = Ω, c=3e8)
	a = zeros(Complex , M )
	for m in 1:M
		a[m] = 1/M * exp(j * Ω*τ(θ, m-1, d; c=c))
	end
	return a
end

# ╔═╡ b62af42c-3682-4787-90e2-718e41ac084e
begin
	a = DAS(θ0, M, λ/4)
	θ = -90:0.5:90 
	B20 = B(a, θ, λ/4)
	noprint
end

# ╔═╡ 31f41ff9-0ab3-4b66-bd01-a4f5daac6375
begin
	plot(θ, abs.(B20), legend = false)
	plot!(title = "B($(θ0)º, θ)")
	plot!(xlabel = "graus")
	vline!([20], linestyle = :dash)
end

# ╔═╡ 4c047c22-e13d-4298-a5ea-78c0fb39cf45
begin
	plot(θ, abs.(B(DAS(θ0, M, λ/4), θ, λ/4)), label = "λ/4")
	plot!(θ, abs.(B(DAS(θ0, M, λ/2), θ, λ/2)), label = "λ/2")
	plot!(θ, abs.(B(DAS(θ0, M, 3λ/4), θ, 3λ/4)), label = "3λ/4")
	plot!(θ, abs.(B(DAS(θ0, M, λ), θ, λ)), label = "λ")
	plot!(title = "Ganhos para diferentes valores de distâncias de antenas")
	plot!(legend = (0.9,0.3))
end

# ╔═╡ 089f1d31-0dff-4bae-8526-190b83f24189
begin
	A1=1
	A2=0.5
	a3 = DAS(20, M, λ/2)
	A_ =  A1*abs(B(a3, 20, λ/2)) .+ A2.*abs.(B(a3, θ, λ/2))
	noprint
end

# ╔═╡ 6d723d92-a4cd-44b3-8c1e-79bdd219d8bc
begin
	plot(θ, A_)
	plot!(title = "Amplitude máxima", legend = false, xlabel = "graus")
end

# ╔═╡ db7b6207-6b79-45c3-9047-072a83d0e2c8
begin
	fa = 1e12
	A0 = 3
	θ1 = 45 #direção da primeira interferência
	θ2 = -15 #direção da segunda interferência
	x = matread("sinais.mat")["sinal_recebido"] # sinais de cada antenas
	noprint
end

# ╔═╡ d5815751-c3f5-4ae0-91cb-b15554823a92
fa/2

# ╔═╡ 5894017c-47dc-4122-8e3d-0357089b55db
ω2 = Ω/fa

# ╔═╡ 074c5f54-13f0-45cb-81b2-a7111e476eb6
begin
	N = round(Int, fa*1e-9)
end

# ╔═╡ f43b4fca-03ca-4d45-abe5-5fd899e7b059
begin
	# pb = digitalfilter(Lowpass(Ω; fs=fa) ,Butterworth(6))
	pb = digitalfilter(Lowpass(Ω/fa) ,Butterworth(6))
	PB, ω = freqresp(pb)
	noprint
end

# ╔═╡ 031f83f8-8687-43c5-b54d-695467a2fb41
begin
	plot(ω, abs.(PB))
	plot!(title= "Resposta do filtro", xlabel = "rad/amostra")
	vline!([Ω/fa*π], label = "Freq corte")
end
	

# ╔═╡ b1c5dba2-8e82-4127-a917-b71c9e48ac08
begin
	coss = 2*cos.(Ω* range(0, length = size(x)[1], step = 1/fa))
	senn = -2*sin.(Ω* range(0, length = size(x)[1], step = 1/fa))

	
	x̃i = zeros(size(x))
	x̃q = zeros(size(x))
	
	for m in 1:M
		x̃i[:, m] = x[:,m].*coss
		x̃q[:, m] = x[:,m].*senn
	end
	
	xi = zeros(size(x))
	xq = zeros(size(x))
	
	for m in 1:M
		xi[:, m] = filt(pb, x̃i[:,m])
		xq[:, m] = filt(pb, x̃q[:,m])
	end
end

# ╔═╡ 864a1b91-b216-49e3-b26d-827a53bec374
begin
	plot(xi[:,1][1:500], label = "xi")
	plot!(xq[:,1][1:500], label = "xq")
	plot!(title = "Sinal imaginário na antena 1")
	
end

# ╔═╡ 8530f05a-9bd5-49a9-8e86-affd823221ad
begin
	w = DAS(θ0, M, λ/4)
	xim = xi + j*xq
	y = xim *w
	noprint
end

# ╔═╡ 7646d2ca-dd7e-40d8-b675-7a5c3f62b9a0
begin
	plot(real.(y), label = "Re{y}")
	plot!(imag.(y), label = "Im{y}", alpha = 0.5)
	plot!(title= "sinal resultado das antenas DAS")
end

# ╔═╡ 9eb5f1da-cda4-4ae5-af60-5197215d412b
Ω/fa*pi

# ╔═╡ b3fbd957-7764-4663-bce9-54bd6f5518c2
function QAM(real, imaginario)
	
	local a =  real
	local b = imaginario
	
	if 		a == -3 && b == 3
		return 0
	elseif  a == -3 && b == 1
		return 1
	elseif  a == -3 && b == -3
		return 2
	elseif  a == -3 && b == -1
		return 3
	elseif  a == -1 && b == 3
		return 4
	elseif  a == -1 && b == 1
		return 5
	elseif  a == -1 && b == -3
		return 6
	elseif  a == -1 && b == -1
		return 7
	elseif  a == 3 && b == 3
		return 8
	elseif  a == 3 && b == 1
		return 9
	elseif  a == 3 && b == -3
		return 10
	elseif  a == 3 && b == -1
		return 11
	elseif  a == 1 && b == 3
		return 12
	elseif  a == 1 && b == 1
		return 13
	elseif  a == 1 && b == -3
		return 14
	elseif  a == 1 && b == -1
		return 15
	else
		return -1
	end
	
end
	
		

# ╔═╡ f3161d89-dc95-47cf-8ab2-c46496923313
begin
	i = 3
	
	rea = round(Int, mean(real.(y[i*N + 1:(i+1)*N])))
	ima = round(Int, mean(imag.(y[i*N + 1:(i+1)*N])))
	QAM(rea, ima)
end

# ╔═╡ ee596dc8-574c-4ba3-b720-5d7aa416372c
rea

# ╔═╡ 314c138a-be11-482e-91dc-273909c1f3cf
ima

# ╔═╡ Cell order:
# ╠═680e5c90-fc8e-11ec-1304-0deea9b53834
# ╟─fb2fa972-a156-4980-9d15-0be562958957
# ╟─12c8cb4e-bc06-44c3-975b-74d1b9323e29
# ╠═4d2fb09b-36a4-474b-a908-05b174b68c4f
# ╠═b62af42c-3682-4787-90e2-718e41ac084e
# ╟─31f41ff9-0ab3-4b66-bd01-a4f5daac6375
# ╟─e6d1fc78-8723-418d-be3b-8eddc3430be0
# ╟─4c047c22-e13d-4298-a5ea-78c0fb39cf45
# ╟─539fca48-8236-4734-a7bc-4049bc8d347b
# ╠═089f1d31-0dff-4bae-8526-190b83f24189
# ╟─6d723d92-a4cd-44b3-8c1e-79bdd219d8bc
# ╟─a9338f7b-e329-4bd5-a576-5da4896b6ba4
# ╟─138691ba-b2ec-4bf8-ac09-0d0cc8ac3bcb
# ╠═db7b6207-6b79-45c3-9047-072a83d0e2c8
# ╠═9562b052-e711-4aaa-8aec-072f83a70de3
# ╠═d5815751-c3f5-4ae0-91cb-b15554823a92
# ╠═5894017c-47dc-4122-8e3d-0357089b55db
# ╟─f9f27eba-5e67-4f7e-bb87-932a4273db0c
# ╠═f43b4fca-03ca-4d45-abe5-5fd899e7b059
# ╠═031f83f8-8687-43c5-b54d-695467a2fb41
# ╠═74784532-bf4a-4d84-8bf1-3f1b279a8567
# ╠═b1c5dba2-8e82-4127-a917-b71c9e48ac08
# ╠═864a1b91-b216-49e3-b26d-827a53bec374
# ╠═41ce1800-51b8-4424-a01f-2d7ebf5d2117
# ╠═8530f05a-9bd5-49a9-8e86-affd823221ad
# ╠═7646d2ca-dd7e-40d8-b675-7a5c3f62b9a0
# ╠═7e4fb1df-3d37-4290-b0d3-3b2f985ae234
# ╠═074c5f54-13f0-45cb-81b2-a7111e476eb6
# ╠═f3161d89-dc95-47cf-8ab2-c46496923313
# ╠═ee596dc8-574c-4ba3-b720-5d7aa416372c
# ╠═314c138a-be11-482e-91dc-273909c1f3cf
# ╠═6e695e78-b2db-4128-bdca-fe244d4667f5
# ╠═58b09e7c-9bd0-459e-8077-97c020758932
# ╠═f41cf0d5-af4b-49ee-b9ce-917cd6e57866
# ╠═72f16a75-fc9d-442c-be94-f443db657efa
# ╠═9ea61348-08e5-4070-812b-a013478afeb8
# ╠═199b6831-1cf9-44fa-9afc-f1bc838ce090
# ╠═edc6e700-01ee-4fa9-8eec-b98d323e1fa2
# ╠═fecad3c1-f01d-49a0-9141-77577f29c6e3
# ╠═36ef6932-be01-4760-a104-6764d4746787
# ╠═9eb5f1da-cda4-4ae5-af60-5197215d412b
# ╠═b3fbd957-7764-4663-bce9-54bd6f5518c2
