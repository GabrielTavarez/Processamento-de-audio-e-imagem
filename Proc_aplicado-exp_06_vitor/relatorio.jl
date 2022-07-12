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

# ╔═╡ 0c2dc34f-d88e-4a2c-a767-d9e06d02071a
md" # ExperiÊncia 06

Gabriel Tavares - 10773801

Guilherme Reis - 10773700"

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
md" ## 03 - Ganho da interferência

Dado um sinal x1 vindo de θ1 conhecido e um sinal x2 vindo de θ2, queremos saber a amplitude máxima do sinal a depender da direção θ2 da interferência"

# ╔═╡ 138691ba-b2ec-4bf8-ac09-0d0cc8ac3bcb
md" # Sinal Real"

# ╔═╡ a9338f7b-e329-4bd5-a576-5da4896b6ba4
md" ## 1 - Processamento DAS do sinal "

# ╔═╡ f9f27eba-5e67-4f7e-bb87-932a4273db0c
md"#### Filtro Passa-Baixas

Implementação do filtro Passa-Biaxas com frequÊncia de corte igual a portadora.

Usaremos um filtro Butterworth de 6 coeficientes.
"

# ╔═╡ 74784532-bf4a-4d84-8bf1-3f1b279a8567
md" #### Divisão dos sinais em componente imaginária e real

Nessa etapa multiplicamos o sinal por $$2cos(Ω)$$ e $$-2sen(Ω)$$ para reinterpretar o sinal real como uma parte imaginária e uma parte real do sinal modulado."

# ╔═╡ 41ce1800-51b8-4424-a01f-2d7ebf5d2117
md" ## 2 - Projeto dos coeficientes das antenas

Tendo o sinal reinterpretado em cada antena, vamos multiplicar cada sinal por um coeficiente definido da mesma maneira que o item anterior"

# ╔═╡ 7e4fb1df-3d37-4290-b0d3-3b2f985ae234
md" ## 3 - Decodificação

A decodificação irá determinar o tamanho do intervalo de transmissão de cada símbolo calcula o valor médio do número nesse intervalo e aproxima pra algum valor da tabela de 16QAM."

# ╔═╡ 0e0adc5c-14fe-4e76-96c2-8ab4b1372664
md" ## 4 - Decodificação do sinal completo

Agora repetimos o processo a cima em todos os 10 trechos de N amostras"

# ╔═╡ 300585fb-1ec8-4f86-9669-ea23b3db04a0
md" ## 05 - Adição de ruído ao sinal


Depois iremos repetir o processo a cima de leitura das antenas mas adicionando um ruído gaussiano nos sinais de entrada para verificar quão robusto é a modulação a ruidos."

# ╔═╡ 58b09e7c-9bd0-459e-8077-97c020758932
md" # Functions

Definições das funções usadas nesse código"

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

# ╔═╡ db7b6207-6b79-45c3-9047-072a83d0e2c8
begin
	fa = 1e12
	A0 = 3
	θ1 = 45 #direção da primeira interferência
	θ2 = -15 #direção da segunda interferência
	x = matread("sinais.mat")["sinal_recebido"] # sinais de cada antenas
	noprint
end

# ╔═╡ f43b4fca-03ca-4d45-abe5-5fd899e7b059
begin
	fc = 2*F/fa
	pb = digitalfilter(Lowpass(F, fs = fa) ,Butterworth(6))
	PB, ω = freqresp(pb)
	noprint
end

# ╔═╡ 031f83f8-8687-43c5-b54d-695467a2fb41
begin
	plot(ω, abs.(PB), label =false)
	plot!(title= "Resposta do filtro", xlabel = "rad/amostra")
	vline!([fc*π], label = "Freq corte")
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

# ╔═╡ 72f16a75-fc9d-442c-be94-f443db657efa
function τ(θ, m, d; c = 3e8)
	return m*d*sind(θ)/c
end

# ╔═╡ 9ea61348-08e5-4070-812b-a013478afeb8
function B(a, θ, d ;c = 3e8)
	M = length(a)
	B_ = 0 + 0*j
	
	for m in 1:M
		B_ += a[m]*exp(-j*Ω*τ(θ,m, d; c = c))
	end
	return B_
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

# ╔═╡ 8530f05a-9bd5-49a9-8e86-affd823221ad
begin
	w = DAS(θ0, M, λ/2)
	xim = xi + j*xq
	y = xim * conj.(w)
	noprint
end

# ╔═╡ 7646d2ca-dd7e-40d8-b675-7a5c3f62b9a0
begin
	index = 1
	plot(real.(y), label = "Re{y}")
	plot!(imag.(y), label = "Im{y}", alpha = 0.5)
	plot!(title= "sinal resultado das antenas DAS")
end

# ╔═╡ e73d9f93-bde8-4925-92da-2011759dc0f8
begin
	A_ruido = 30 #Amplitude do ruído
	x_ruido = x + A_ruido*randn(size(x))
	
	x̃i_ruido = zeros(size(x))
	x̃q_ruido = zeros(size(x))
	
	for m in 1:M
		x̃i_ruido[:, m] = x_ruido[:,m].*coss
		x̃q_ruido[:, m] = x_ruido[:,m].*senn
	end
	
	xi_ruido = zeros(size(x))
	xq_ruido = zeros(size(x))
	
	for m in 1:M
		xi_ruido[:, m] = filt(pb, x̃i_ruido[:,m])
		xq_ruido[:, m] = filt(pb, x̃q_ruido[:,m])
	end
	
	# w = DAS(θ0, M, λ/2)
	xim_ruido = xi_ruido + j*xq_ruido
	y_ruido = xim_ruido * conj.(w)
	
	noprint
end

# ╔═╡ e99c41cd-aa41-49d3-a976-e06e7ce0bbd9
md" Adicionando um ruído de amplitude $A_ruido vemos que ainda temos uma taxa de acerto bastante alta. Isso ocorre porque a média de um ruído tende a zero quando temos pontos suficientes"

# ╔═╡ ac88b1f6-7ad4-4798-855d-82725ba7119c
function QAM(numero::Complex{})
	
	QAM_table = [
		-3+3j, # 0
		-3+1j, # 1
		-3-3j, # 2
		-3-1j, # 3
		-1+3j, # 4
		-1+1j, # 5
		-1-3j, # 6
		-1-1j, # 7
		3+3j,  # 8
		3+1j,  # 9
		3-3j,  # 10
		3-1j,  # 11
		1+3j,  # 12
		1+1j,  # 13
		1-3j,  # 14
		1-1j   # 15
		]
	
	min_dist = Inf
	dist = Inf
	melhor_indice = 0
	
	for indice in 1:16
		dist = abs(numero - QAM_table[indice])
		if dist < min_dist
			min_dist = dist
			melhor_indice = indice-1
		end	
	end
	
	return melhor_indice
end		

# ╔═╡ 6e695e78-b2db-4128-bdca-fe244d4667f5
begin
	N = round(Int, fa*1e-9)
	
	codigo = QAM(mean(real.(y[1:N]))+ j* mean(imag.(y[1:N])))
	
	noprint
end

# ╔═╡ eeb5e002-a2ba-4ea9-8cd0-6dfdd16f4260
md" O simbolo do primeiro intervalo lido é **$codigo**"

# ╔═╡ 1230836e-2560-4514-a9b3-cd9ce599d0f7
begin
	mensagem = matread("mensagem.mat")["mensagem"]
	
	num_simbolos = round(Int,size(x)[1]/N)
	mensagem_lida = zeros(num_simbolos)
	QAM_lido = zeros(num_simbolos,2)
	
	for i in 1:num_simbolos
		QAM_lido[i,1] = mean(real.(y[(i-1)*N + 1 : (i)*N]))
		QAM_lido[i,2] = mean(imag.(y[(i-1)*N + 1 : (i)*N]))
		
		mensagem_lida[i] = QAM(QAM_lido[i,1]+ j* QAM_lido[i,2]) 
	end
	
	mensagem_lida
end

# ╔═╡ 4244cb05-047e-43bb-a4d2-91d6d6c0e7c7
begin
	scatter(QAM_lido[:,1], QAM_lido[:,2])
	plot!(title = "Mapa complexo da mensagem lida em 16QAM")
	plot!(legend = false)
end

# ╔═╡ bb3a275c-f813-45f1-8f65-3a61cb7a38b2
begin
	mensagem_lida_ruido = zeros(num_simbolos)
	QAM_lido_ruido = zeros(num_simbolos,2)
	
	for i in 1:num_simbolos
		QAM_lido_ruido[i,1] = round(mean(real.(y_ruido[(i-1)*N + 1 : (i)*N])))
		QAM_lido_ruido[i,2] = round(mean(imag.(y_ruido[(i-1)*N + 1 : (i)*N])))
		
		mensagem_lida_ruido[i] = QAM(QAM_lido_ruido[i,1] + j*QAM_lido_ruido[i,2]) 
	end
	
	mensagem_lida_ruido
end

# ╔═╡ d2cf837e-f131-447f-996b-f850132b3b10
mensagem_lida_ruido - mensagem

# ╔═╡ Cell order:
# ╟─0c2dc34f-d88e-4a2c-a767-d9e06d02071a
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
# ╟─138691ba-b2ec-4bf8-ac09-0d0cc8ac3bcb
# ╠═db7b6207-6b79-45c3-9047-072a83d0e2c8
# ╟─a9338f7b-e329-4bd5-a576-5da4896b6ba4
# ╟─f9f27eba-5e67-4f7e-bb87-932a4273db0c
# ╠═f43b4fca-03ca-4d45-abe5-5fd899e7b059
# ╠═031f83f8-8687-43c5-b54d-695467a2fb41
# ╟─74784532-bf4a-4d84-8bf1-3f1b279a8567
# ╠═b1c5dba2-8e82-4127-a917-b71c9e48ac08
# ╠═864a1b91-b216-49e3-b26d-827a53bec374
# ╟─41ce1800-51b8-4424-a01f-2d7ebf5d2117
# ╠═8530f05a-9bd5-49a9-8e86-affd823221ad
# ╟─7646d2ca-dd7e-40d8-b675-7a5c3f62b9a0
# ╟─7e4fb1df-3d37-4290-b0d3-3b2f985ae234
# ╠═6e695e78-b2db-4128-bdca-fe244d4667f5
# ╟─eeb5e002-a2ba-4ea9-8cd0-6dfdd16f4260
# ╟─0e0adc5c-14fe-4e76-96c2-8ab4b1372664
# ╠═1230836e-2560-4514-a9b3-cd9ce599d0f7
# ╟─4244cb05-047e-43bb-a4d2-91d6d6c0e7c7
# ╟─300585fb-1ec8-4f86-9669-ea23b3db04a0
# ╠═e73d9f93-bde8-4925-92da-2011759dc0f8
# ╠═bb3a275c-f813-45f1-8f65-3a61cb7a38b2
# ╠═d2cf837e-f131-447f-996b-f850132b3b10
# ╟─e99c41cd-aa41-49d3-a976-e06e7ce0bbd9
# ╟─58b09e7c-9bd0-459e-8077-97c020758932
# ╠═36ef6932-be01-4760-a104-6764d4746787
# ╠═72f16a75-fc9d-442c-be94-f443db657efa
# ╠═9ea61348-08e5-4070-812b-a013478afeb8
# ╠═199b6831-1cf9-44fa-9afc-f1bc838ce090
# ╠═edc6e700-01ee-4fa9-8eec-b98d323e1fa2
# ╠═fecad3c1-f01d-49a0-9141-77577f29c6e3
# ╠═ac88b1f6-7ad4-4798-855d-82725ba7119c
