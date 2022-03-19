### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 57ddaf90-1a1e-11ec-30d9-47139b4361f4
begin
	using Plots
	using DSP
	using FFTW
	using WAV
	using SampledSignals
	using LaTeXStrings
	using FixedPointNumbers
	using Statistics
	plotly()
	md"Bibliotécas usadas: Plots, DSP, FFTW, WAV, SampledSignals, FixedPointNumbers,Statistcs"
end

# ╔═╡ a0be60d9-5d35-4342-8cfe-b9b28e717af6
md" # Experiência 1

Gabriel Tavares 10773801

Guilhemre Reis 10773700

Diego Hidek 10336622

"

# ╔═╡ 20b24c92-b9a1-49b8-8aa8-8f7f6f244ced
   md" # Redução de taxa de amostragem"

# ╔═╡ dc3fb707-66b4-459e-b62a-7826b27a5dbf
md" ## Sinal enunciado

Usaremos o sinal $x(t) = sin (1000.2\pi t) - 1/3 sin(21000.2\pi t)$ para realizar o experimento. Inicialmente a frequência de amostragem de 48kHz

"

# ╔═╡ bbcffcd2-0aef-4c4b-a586-ba17581b31d9
begin
	fa1 = 48e3 #48kHz	
	fa2 = 24e3 #24kHz
end

# ╔═╡ 248ad4a8-991d-4b2a-b5ec-9044c87c662b
begin 
	range_plot = 100:300
	plot_range = range_plot
	t = 0:1/fa1:0.05
	t_ms = t*1000
	x_48 = sin.(1e3*2π * t) - 1/3*sin.(21e3*2π * t)
	plot(t[range_plot]*1000,x_48[range_plot], label = "x(t)")
	plot!(title = "Sinal original", xlabel = "t(ms)")
end

# ╔═╡ bfa602aa-eedf-4fa2-9099-c2fa22914939
md" ## 1) Projeto de filtro

Projeto de filtro com frequência de passagem de 11kHz, frequência de rejeição de 13kHz e tolerância de 0.01

"

# ╔═╡ 3e85535b-9e1a-4457-8421-952bea01cd01
md" ## 2) Filtragem do sinal "

# ╔═╡ bb5bd0dc-901c-439f-86b4-08a2b9bcfbe6
md" ## 3) Calculando espectros e respostas do filtro

Nesta secção iremos comparar o ganho do filtro do sinal em coparação com a resposta em frequência, e a razão entre as TDFs dos sinais de entrada e saída

"

# ╔═╡ e0f03303-5ced-4651-bd3f-a8ddb1479385
md" ## 4) Redução da taxa de amostragem

Nesta seção iremos reduzir a taxa de amostragem e comparar a redução do sinal filtrado e não filtrado

"

# ╔═╡ 0e663744-3401-428a-9e95-d24938584c42
md" # Aumento da Taxa de Amostragem

Nesta etapa iremos fazer uma reamostragem de um sinal que foi amostrado a 44kHz para uma taxa L vezes maior."

# ╔═╡ 0f8dad6c-a251-406c-9cd5-1097f93ce63f
md"
Inicialmente usaremos o mesmo sinal  $x(t) = sin (1000.2\pi t) - 1/3 sin(21000.2\pi t)$ amostrado a 48kHz e reamostraremos ele a 144kHz.

$144kHz = 48kHz*L= 48kHz*3$
"

# ╔═╡ cad2b83d-4f12-4382-baa3-363d3ebcd928
begin
	L = 3
	fa_48 = 48_000
	fa_144 = fa_48*L
end

# ╔═╡ 9a9f35ec-aa90-45b5-848f-d2d236e209c0
md" ## 1) Filtro Passa Baixas

Primeiramente é necessário criar um filtro que servirá como interpolador ideal. O interpolador terá a frequência de corte de $\omega_c = \pi/3$ e ganho de passagem de 3.

Os parâmetros da másara do filtro estão descritos na célula a baixo
"

# ╔═╡ acf76b2a-8ae8-4b45-bcb8-7fec77f2b3a5
md" ## 2) Criação do sinal super amostrado

O sinal superamostrado será criado colocando 0's entre as amostras e depois interpolando este sinal com um filtro (sinc).
"

# ╔═╡ 7d8d9ca6-c878-479f-a10a-663195db9174
begin
	t_144 = range(0,t[end],step = 1/fa_144)
	t_144_ms = t_144.*1000
	y = zeros(length(x_48)*3)
	for i  in 1:length(x_48)
		y[1+(i-1)*3] = x_48[i]
	end
end

# ╔═╡ 11c794dd-45f3-480d-aaf1-e5b2df8e875d
begin
	plot(t_144_ms[range_plot],y[range_plot])#, line=:stem, marker=:circle)
	title!("Sinal superamostrado com zeros", legend = false, xlabel ="t(ms)" )
end

# ╔═╡ b0efe5ba-7f2e-40d1-9e24-223d113734be
md" Interpolação/Filtragem do sinal amostrado com zeros"

# ╔═╡ 17e7ad63-05e0-4ca3-98f0-3e704502b3aa
md" ## 3) 	Comparação do sinal ideal com o reamostrado

Para verificar a qualidade dessa super amostragem, iremos amostrar o sinal original a 144kHz e comparar os resultados.

Obs: iremos atrasar o sinal  ideal para compensar o efeito da filtragem do sinal"

# ╔═╡ 40185644-29c2-423f-9cfc-b546d0351c57
begin 
	x_144_ideal = sin.(1e3*2π * t_144) - 1/3*sin.(21e3*2π * t_144)
end

# ╔═╡ 11004fa9-d2f8-4a66-8ec2-416f2dd2653e
md"É possível ver que a reamostragem e interpolação com o filtro Passa-Baixas gerá um sinal muito parecido com o sinal originl amostrado a 144kHz"

# ╔═╡ 729d9438-3c69-4596-b4e0-4d1b1f6227b2
md" ## 4)Espectro do sinal reamostrado
"

# ╔═╡ 1b1e4af6-379e-4d1b-9097-5fb29ea5c505
md" # Aumento da taxa de amostragem usando iterpolação linear

Nesta etapa iremos comparar a interpolação com o filtro Passa-Baixas e uma interpolação linear simples
"

# ╔═╡ 8a6ffc9a-1845-41d7-b3ae-f3fa69d6b6ec
md"
Para fazer a interpolação linear usaremos um filtro triangular como interpolador para poder analisar sua resposta em frequência e netender a origem da distorção do sinal 
"

# ╔═╡ 28bc738b-7e46-4d84-a0d8-969150975f4d
md" Observa-se que o sinal com interpolação linear tem diferenças consideráveis em relação ao sinal original. Portanto conclui-se que a interpolação com o filtro Passa-Baixas tem um resultado muito melhor."

# ╔═╡ 0d610ccf-c261-4db1-b20a-d9aabb1e7d96
md" # Conversão A/D com sobreamostragem

Nesta etapa iremos fazer a comparação entre a uma conversão analógico digital simples (com taxa de conversã de 40kHz) e uma conersão com sobre amostragem (com fator de M=19)."

# ╔═╡ 075e5e24-ad15-44da-86fe-1eef1a4c5ed0
md" ## Convesão simples

Aqui teremos uma conversão direta de um sinal analógico para um sinal digital de 5 bits a 40kHz
"

# ╔═╡ 2b95cdbf-c860-477a-a178-1eb39b0e093e
md" ### Sinal a 40kHz"

# ╔═╡ 0bc85a2d-bc2e-48f5-b5e5-13044f6c8673
md" ### Relação Sinal Ruído

Aqui iremos modelar o ruído observado como uma variável aleatória uniforme para ter uma noção teórico do ruído esperado.

"

# ╔═╡ d8ff7a10-3e34-4ccb-a318-71c6032265b5
md" ## Conversão com sobre amostragem

Aqui usaremos um arranjo de convesão com sobreamostragem, onde o sinal é amostrado a uma taxa alta, depois filtrado e reamostrado à taxa desejada.
"

# ╔═╡ 0c06bc9a-349a-439d-b682-ed7720fdb98b
md" ### Criação do sinal"

# ╔═╡ 11583038-72cd-435c-9f8b-bbbe949ff8d5
md" ### Filtro Passa-Baixas"

# ╔═╡ 99357a1b-e2fc-4894-8a52-5cf9cd0fbf2c
md" ### Filtragem e reamostragem"

# ╔═╡ 7bc5802b-7519-4394-a630-f34f377186e2
md" ### Relação Sinal-Ruído
"

# ╔═╡ 6b2ddf23-4b75-40b5-a4ff-bfb696a8bb68
md" ### Número de bits relativo"

# ╔═╡ e03d504d-cc90-4aef-856b-28ad65f1d500
md" # Anexos

Funções usadas nessa experiência
"

# ╔═╡ 149d7c6e-5c53-4847-b989-e2e30bb5b9bb
DSP.freqz(filt_fir::Vector, ω) = DSP.freqz(PolynomialRatio(filt_fir,[1]),ω)

# ╔═╡ 22c75efc-0d44-475a-846a-f48045862e2c
function kaiserbeta(δp, δr, Δω)
	 δ = min(δp, δr)
	 A = -20log10(δ)
	 if A < 21
		return 0.0
	 elseif A ≤ 50
		return 0.5842(A-21)^0.4 + 0.07886(A-21)
	 else
		return 0.1102(A-8.7)
	 end
end	


# ╔═╡ f1235f75-71ad-457b-8937-ac6ea8c3acba
"""
    kaiser_filter_lowpass{δp, δr, ωp, ωr}
Retorna os coeficientes de um filtro FIR passa baixas de ganho 1 definido com os parâmetros:

`δp` : Atenuação linear da banda passante em escala linear 

`δr` : Atenuação linear da banda de rejeição em escala linear 

`ωp` : Frequência passante normalizada (0-π) 

`ωr` : Frequência de rejeição normalizada (0-π)
"""
function kaiser_filter_lowpass(δp, δr, ωp, ωr)
	#retorna um filtro passa baixas apenas
	Δω = ωr - ωp
	A= -20log10(min(δp, δr))
	Nk_aux = ceil(Int,(A - 8) / (2.285 * Δω) + 1)
	Nk = ifelse(iseven(Nk_aux), Nk_aux+1, Nk_aux)
	β = kaiserbeta(δp, δr, Δω)
	nk = 0:Nk-1
	kaiser_window = kaiser(Nk,β/π)

	Lk = (Nk -1)÷2
	ωc = (ωr+ωp)/2
	h = ωc/π * sinc.(ωc/π .* (nk.-Lk))
	hk = h.*kaiser_window
	return hk
end

# ╔═╡ 86695cac-aa9c-494a-a43b-cdbdd5ea3f22
begin
	δp = 0.02
	δr = 0.01
	
	fp = 11e3
	fr = 13e3
	
	ωp = fp/fa1 *2π
	ωr = fr/fa1 * 2π
	
	pb_12k_kaiser = kaiser_filter_lowpass(δp, δr, ωp, ωr)
end

# ╔═╡ aacf5e9c-6cf3-44c9-9656-1e3952e81129
begin
	ω = range(0,π, length = 500)
	PB_12k_kaiser = freqz(pb_12k_kaiser,ω)
	f = ω/(2π)*fa1
	plot(f/1000, abs.(PB_12k_kaiser))
	plot!([0,fp]./1000, [1-δp,1-δp], color =:red)
	plot!([0,fp]./1000, [1+δp,1+δp], color =:red)
	plot!([fr,fa1/2]./1000, [δr,δr], color =:red)
	plot!(title = "Resposta em frequência", xlabel = "f(kHz)", legend = false)
end

# ╔═╡ eadfca3f-4c38-4fc3-aef9-a6eef900bb52
begin
	filt_interp_linear = [1/3,2/3,1,2/3,1/3]
	y_linear = filt(filt_interp_linear,y)
	
	FILT_interp_linear = freqz(filt_interp_linear,ω)
end

# ╔═╡ d6b9635b-43cf-476a-9a91-b070b67332f5
begin
	plot(filt_interp_linear)
	plot!(title="Resposta ao impulso da interpolação linear", legend = false)
end

# ╔═╡ 46f5473c-d6f9-48e6-bffa-c88ddf63a604
begin
	plot(ω/π,abs.(FILT_interp_linear))
	plot!(title="Resposta em frquência do interpolador linear", legend = false)
end

# ╔═╡ 687b32c0-4f09-4ef7-9a72-59916ae8af60
x_filtrado = filt(pb_12k_kaiser,x_48)

# ╔═╡ 93c5cb04-d476-46f7-b9ad-763f207cc534
begin	
	plot(t_ms[range_plot], x_48[range_plot.+21], label = "Sinal entrada e deslocado", alpha = 0.5)
	plot!(t_ms[range_plot], x_filtrado[range_plot], label = "Sinal filtrado")
	plot!(title= "Sinal filtrado", xlabel = "t(ms)")
end

# ╔═╡ b15444df-5c5f-40ac-9c80-ddfc5bde5aba
begin
	X =  fft(x_48[1:48])
	X_saida = fft(x_filtrado[2*length(pb_12k_kaiser):2*length(pb_12k_kaiser)+48-1])
	
	fk = range(0,48,length = 48)
end

# ╔═╡ 10a0943d-c215-4b6d-85f8-ca5c8fd04745
begin 
	plot(fk, abs.(X)/48, line =:stem , marker=:circle , label = "TDF entrada")
	plot!(fk, abs.(X_saida)./48, line =:stem, marker=:circle, label = "TDF saída")
	plot!(title="TDF do sinal", xlabel ="f(kHz)", legend = (0.75,1))
end

# ╔═╡ 40993b17-b9ce-44a2-a15f-2162178a7e38
begin 
	ω_1kHz = 1e3*2π/fa1
	ω_21kHz = 21e3*2π/fa1
	indice_1k_ω = findall(x->isapprox(x,ω_1kHz,atol = 0.003),ω)
	indice_21k_ω = findall(x->isapprox(x,ω_21kHz,atol = 0.003),ω)
	
	PB_1k_filtro = abs.(PB_12k_kaiser[indice_1k_ω])[]
	PB_21k_filtro = abs.(PB_12k_kaiser[indice_21k_ω])[]
	
	indice_1k_fk = findall(x-> isapprox(x,1,atol = 0.1),fk)
	indice_21k_fk = findall(x-> isapprox(x,21,atol = 0.5),fk)
	
	PB_1k_fft = abs(X_saida[indice_1k_fk][])/abs(X[indice_1k_fk][])
	PB_21k_fft =abs(X_saida[indice_21k_fk][])/abs(X[indice_21k_fk][])
end

# ╔═╡ 4f3e7be3-807e-4856-a2a8-53a80b9abd7c
md"

Tabela de ganhos calculados pela razão entre a fft na entrada e saída do sinal; e ganho da resposta em frequência do sinal.

Frequência (kHz) | Ganho Teórico ($H(e^{j\omega})$) | Ganho Calculado (FFT)
:--------- | :------------------------: | :---:
$1$   | $(round(abs(PB_1k_filtro),digits=4)) | $(round(abs(PB_1k_fft),digits=4))
$21$  | $(round(abs(PB_21k_filtro),digits=4)) | $(round(abs(PB_21k_fft),digits=4))

"

# ╔═╡ 49a1a852-8f8f-4414-a9b9-0257680a08a6
begin
	ωp3 =2π*21_000/fa_144
	ωr3 = 2π/L - 2π*21_000/fa_144
	δp3 = 0.005
	δr3 = 0.01
	pb_interpolador = kaiser_filter_lowpass(δp3,δr3,ωp3,ωr3)*3
end

# ╔═╡ 19ab55de-04ba-487b-8260-2c1b3165054f
begin
	PB_interpolador = freqz(pb_interpolador,ω)
	plot(ω/π , abs.(PB_interpolador))
	plot!([0,ωp3/π], [3*(1-δp3),3*(1-δp3)], color =:red)
	plot!([0,ωp3/π], [3*(1+δp3),3*(1+δp3)], color =:red)
	plot!([ωr3/π,1], [3*δr,3*δr], color =:red)
	plot!(title = "Resposta em frequência", xlabel = "f(kHz)", legend = false)
end

# ╔═╡ c1a60f10-2f7e-4915-bdc2-ae7641f07312
begin
	x_144 = filt(pb_interpolador,y)
end

# ╔═╡ a1ce07ff-4f39-4603-a091-bc62d44cfd78
begin
	plot(t_144_ms[range_plot],x_144[range_plot])
	plot!(title="Sinal super amostrado e filtrado", legend = false, xlabel = "t(ms)")
end

# ╔═╡ ecbd8b83-bddb-41d5-8e50-55fca4356e92
begin
	plot(x_144_ideal[range_plot.-32], width = 6, alpha = 0.5, label = "Sinal amostrado ideal")
	plot!(x_144[range_plot], width = 2, label = "Sinal reamostrado com filtro")
	plot!(title = "Comparaçõ entre sinal ideal e reamostrado" ,legend = (0.68,0.2), xlabel = "t(ms)")
end

# ╔═╡ b5052fd0-d8d9-465e-90ed-94f9cec385c4
begin
	X_144 = fft(x_144[200:1200])/length(x_144[200:1200])
end

# ╔═╡ 7b164f17-30a6-4820-8419-7e16b2ff3d79
begin
	f_144 = range(-fa_144/2, fa_144/2, length = length(X_144))
	plot(f_144/1000, abs.(fftshift(X_144)))
	plot!(title="Espectro do sinal reamostrado", xlabel = "f(kHz)", legend = false)
	#colocar as frequências
end

# ╔═╡ 919d504d-b67e-46a1-9883-5675e0d5b9b6
begin
	plot(x_144_ideal[range_plot], label = "Sinal ideal", width = 4, alpha = 0.5, color = :blue)
	
	plot!(x_144[range_plot.+32], color =:yellow, label = "Interpolação com filtro Passa-Baixas")
	
	plot!(y_linear[range_plot.+2], label = "Sinal reamostrado com interpolação linear", width = 2, color =:red)
	
	plot!(title = "Comparação entre sinal ideal e com intepolação linear", legend = (0.55,0.2), xlabel = "t(ms)")
	
end

# ╔═╡ 0d1eab07-8971-4ea2-ae14-4a2dcf6dd3be
begin
	δp5 = 0.0001
	δr5 = 0.0001
	ωp5 = 6/40
	Δω5 = π/100
	ωr5 = ωp5 + Δω5
	
	M5 = Int16(ceil(π/(ωp5 + Δω5/2) ))
	
	pb5 = kaiser_filter_lowpass(δp5, δr5, ωp5, ωr5)
	PB5 = freqz(pb5, ω)
	length(pb5)
end

# ╔═╡ 8c5058c3-a03f-4fb6-a8a9-e131b4059d90
begin
	#parâmetros do exercício
	Ω0 = 3_000
	ΔΩ = 3_000
	Ω1 = 2π*750
	
	
	B = 5 #bits
	
	
	f_40 = 40_000
	fa_40 = 40_000
	fa_5 = f_40*M5 #760kHz
	t_40 = 0:1/fa_40:2
	t_40_ms = t_40*1000
	t_5 = 0:1/fa_5:2
	t_5_ms = t_5*1000	
end

# ╔═╡ 44a2dce8-b991-450d-85a9-90f49406f709
begin
	s_40_simples = 0.7sin.((Ω0 .+ 0.5*ΔΩ*t_40).*t_40) + 0.3cos.(Ω1*t_40)
	sq_40_simples = Fixed{Int16,B-1}.(s_40_simples)
end

# ╔═╡ ae4d01f1-1648-4b83-8288-75b1c18bf934
begin
	p_s_40 = plot(t_40_ms[range_plot], s_40_simples[range_plot])
	plot!(title="Sinal amostrado a 40kHz")
	
	p_sq_40 = plot(t_40_ms[range_plot], sq_40_simples[range_plot])
	plot!(title = "Sinal quatizado amostrado a 40 kHz")
	plot(p_s_40,p_sq_40,layout=(2,1))
end

# ╔═╡ a3f7b9a6-5eca-46af-b754-da2edaadd9a7
md" Som original

$(SampleBuf(s_40_simples, fa_40))

Som quantizado

$(SampleBuf(sq_40_simples, fa_40))


Vemos que há um ruído na quantização do sinal, que será modelado como um ruído uniforma para efeito de cálculos teóricos.
"

# ╔═╡ c72d4c62-39a9-4e60-9973-87013e8ca386
begin
	S_40_teo = 0.7^2.0/2 + 0.3^2.0/2
	N_40_teo = 2^(-2.0*B)/3
	σ5 = N_40_teo

	ϵ_simples = s_40_simples-sq_40_simples
	
	S_40_real = mean(s_40_simples.^2)
	N_40_real = mean(ϵ_simples.^2)
	
	SN_teo_simples = pow2db(S_40_teo/N_40_teo)
	SN_real_simples = pow2db(S_40_real/N_40_real)	
end

# ╔═╡ 491a73a4-27bd-4a3d-9254-d12d1a96a885
md" Temos uma relação **SNR = $(round(SN_real_simples, digits = 3))dB** na conversão simples a 40kHz de 5 bits. "

# ╔═╡ 924f8479-bb95-47f7-83ab-478bd8f93923
begin
	s_5 = 0.7sin.((Ω0 .+ 0.5*ΔΩ*t_5).*t_5) + 0.3cos.(Ω1*t_5)
	sq_5 = Fixed{Int16,B-1}.(s_5)
end	

# ╔═╡ 29df134f-9aa2-4583-bb42-c0d1c5c58ed0
begin
	p_s_5 = plot(t_5_ms[range_plot], s_5[range_plot])
	plot!(title="Sinal amostrado a 40kHz")
	
	p_sq_5 = plot(t_5_ms[range_plot], sq_5[range_plot])
	plot!(title = "Sinal quatizado amostrado a 40 kHz")
	plot(p_s_40,p_sq_40,layout=(2,1))
end

# ╔═╡ e8b9651e-8b73-4fd6-b0a2-080f29873133
begin
	plot(ω/π , pow2db.(abs.(PB5)))
	plot!([0,ωp5/π], [pow2db((1-δp5)),pow2db((1-δp5))], color =:red)
	plot!([0,ωp5/π], [(1+δp5),(1+δp5)], color =:red)
	plot!([ωr5/π,1], [pow2db(δr5),pow2db(δr5)], color =:red)
	plot!(title = "Resposta em frequência", xlabel = "f(kHz)", legend = false)
end

# ╔═╡ 1ce6cdab-6454-477c-98d8-9d0aa15835bc
function reduz_amostragem(sinal, M)
	#reduz a taxa de amostragem do sinal para a taxa M
	saida = sinal[1:M:end]
end

# ╔═╡ 17f85898-748a-490f-a277-b5a445177d85
begin
	M = 2
	x_reamostrado = reduz_amostragem(x_filtrado,M)
end

# ╔═╡ 2318e5b5-b666-44a0-97c7-be9860bf7865
begin
	plot(t_ms[1:M:end], x_reamostrado, label="Sinal filtrado reamostrado")
	plot!(t_ms[1:M:end], reduz_amostragem(x_48,M)[1:end], label="Sinal original reamostrado")
	
	plot!(xlims = (10,15), xlabel = "t(ms)", title = "Comparação dos sinais reamostrados", legend = (0,0.9) )
end
	


# ╔═╡ de21a988-515b-4e77-898f-f4363aeeaec1
begin
	sq_5_filtrado = filt(pb5, Float64.(sq_5))
	s_5_filtrado = filt(pb5, Float64.(s_5))
	sq_5_40 = reduz_amostragem(sq_5_filtrado,M5)
	s_5_40 = reduz_amostragem(s_5_filtrado,M5)
end

# ╔═╡ d0127175-bb82-44ff-9a58-2aa09848ab91
begin
	plot(t_40_ms[range_plot], sq_40_simples[range_plot], label = "Conversão AD simples a 40kHz", width = 5, alpha = 0.5)
	
	plot!(t_40_ms[range_plot], sq_5_40[range_plot.+26], label = "Conversão AD com sobreamostragem", width = 2)
	
end

# ╔═╡ ca5d533a-4da5-4186-9397-a23388305060
begin
	ϵ = s_5_40-sq_5_40
	
	S_real = mean(s_5_40.^2)
	N_real = mean(ϵ.^2)
	
	SN_real = pow2db(S_real/N_real)
end

# ╔═╡ 89d1f69c-8daa-4312-802e-68bc7d35a026
md" Temos agora uma relação sinal ruído de  **SNR = $(round(SN_real, digits = 3))dB** usando a conversão com sobreamostragem"

# ╔═╡ cbab2437-39e5-4612-b362-ebc8d96c7db4
begin
	B_rel =-log2(3*N_real)/2
end

# ╔═╡ 8e14edcf-ff1c-4218-bc63-e0dba3e167e9
md"O número de bits relativo a essa conversão é de **B' = $(round(B_rel, digits = 3))**"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FixedPointNumbers = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
SampledSignals = "bd7594eb-a658-542f-9e75-4c4d8908c167"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
WAV = "8149f6b0-98f6-5db9-b78f-408fbbb8ef88"

[compat]
DSP = "~0.6.10"
FFTW = "~1.4.5"
FixedPointNumbers = "~0.8.4"
LaTeXStrings = "~1.2.1"
Plots = "~1.22.1"
SampledSignals = "~2.1.2"
WAV = "~1.1.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "b8d49c34c3da35f220e7295659cd0bab8e739fed"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.33"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4ce9393e871aca86cc457d9f66976c3da6902ea7"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.4.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "9995eb3977fbf67b86d0a0a0508e83017ded03f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.14.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "4866e381721b30fac8dda4c8cb1d9db45c8d2994"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.37.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DSP]]
deps = ["FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "2a63cb5fc0e8c1f0f139475ef94228c7441dc7d0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.6.10"

[[DataAPI]]
git-tree-sha1 = "bec2532f8adb82005476c141ec23e921fc20971b"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.8.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "8041575f021cba5a099a456b4163c9a08b566a02"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.1.0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "3c041d2ac0a52a12a27af2782b34900d9c3ee68c"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c2178cfbc0a5a552e16d097fae508f2024de61a3"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.59.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "ef49a187604f865f4708c90e3f431890724e9012"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.59.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "60ed5f1643927479f845b0135bb369b031b541fa"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.14"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "8a954fed8ac097d5be04921d595f741115c1b2ad"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+0"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[IrrationalConstants]]
git-tree-sha1 = "f76424439413893a832026ca355fe273e93bce94"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "34dc30f868e368f8a17b728a1238f3fcda43931a"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.3"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "748f6e1e4de814b101911e64cc12d83a6af66782"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.2"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "c0e9e582987d36d5a61e650e6e543b9e44d9914b"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.7"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "2537ed3c0ed5e03896927187f5f2ee6a4ab342db"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.14"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "4c2637482176b1c2fb99af4d83cb2ff0328fc33c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.22.1"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "OffsetArrays", "RecipesBase"]
git-tree-sha1 = "0b15f3597b01eb76764dd03c3c23d6679a3c32c8"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "1.2.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SampledSignals]]
deps = ["Base64", "Compat", "DSP", "FFTW", "FixedPointNumbers", "IntervalSets", "LinearAlgebra", "Random", "TreeViews", "Unitful"]
git-tree-sha1 = "4b7e413f20fa56fa47b8433c96f96a1acfe372a6"
uuid = "bd7594eb-a658-542f-9e75-4c4d8908c167"
version = "2.1.2"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "a322a9493e49c5f3a10b50df3aedaf1cdb3244b7"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.1"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "a8f30abc7c64a39d389680b74e749cf33f872a70"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.3"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8cbbc098554648c84f79a463c9ff0fd277144b6c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.10"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "1162ce4a6c4b7e31e0e6b14486a6986951c73be9"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.2"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "6c9040665b2da00d30143261aea22c7427aada1c"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.7"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a981a8ef8714cba2fd9780b22fd7a469e7aaf56d"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.9.0"

[[WAV]]
deps = ["Base64", "FileIO", "Libdl", "Logging"]
git-tree-sha1 = "1d5dc6568ab6b2846efd10cc4d070bb6be73a6b8"
uuid = "8149f6b0-98f6-5db9-b78f-408fbbb8ef88"
version = "1.1.1"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─57ddaf90-1a1e-11ec-30d9-47139b4361f4
# ╟─a0be60d9-5d35-4342-8cfe-b9b28e717af6
# ╟─20b24c92-b9a1-49b8-8aa8-8f7f6f244ced
# ╟─dc3fb707-66b4-459e-b62a-7826b27a5dbf
# ╠═bbcffcd2-0aef-4c4b-a586-ba17581b31d9
# ╠═248ad4a8-991d-4b2a-b5ec-9044c87c662b
# ╟─bfa602aa-eedf-4fa2-9099-c2fa22914939
# ╠═86695cac-aa9c-494a-a43b-cdbdd5ea3f22
# ╟─aacf5e9c-6cf3-44c9-9656-1e3952e81129
# ╟─3e85535b-9e1a-4457-8421-952bea01cd01
# ╠═687b32c0-4f09-4ef7-9a72-59916ae8af60
# ╟─93c5cb04-d476-46f7-b9ad-763f207cc534
# ╟─bb5bd0dc-901c-439f-86b4-08a2b9bcfbe6
# ╠═b15444df-5c5f-40ac-9c80-ddfc5bde5aba
# ╟─10a0943d-c215-4b6d-85f8-ca5c8fd04745
# ╠═40993b17-b9ce-44a2-a15f-2162178a7e38
# ╟─4f3e7be3-807e-4856-a2a8-53a80b9abd7c
# ╟─e0f03303-5ced-4651-bd3f-a8ddb1479385
# ╠═17f85898-748a-490f-a277-b5a445177d85
# ╟─2318e5b5-b666-44a0-97c7-be9860bf7865
# ╟─0e663744-3401-428a-9e95-d24938584c42
# ╟─0f8dad6c-a251-406c-9cd5-1097f93ce63f
# ╠═cad2b83d-4f12-4382-baa3-363d3ebcd928
# ╟─9a9f35ec-aa90-45b5-848f-d2d236e209c0
# ╠═49a1a852-8f8f-4414-a9b9-0257680a08a6
# ╟─19ab55de-04ba-487b-8260-2c1b3165054f
# ╟─acf76b2a-8ae8-4b45-bcb8-7fec77f2b3a5
# ╠═7d8d9ca6-c878-479f-a10a-663195db9174
# ╟─11c794dd-45f3-480d-aaf1-e5b2df8e875d
# ╟─b0efe5ba-7f2e-40d1-9e24-223d113734be
# ╠═c1a60f10-2f7e-4915-bdc2-ae7641f07312
# ╟─a1ce07ff-4f39-4603-a091-bc62d44cfd78
# ╟─17e7ad63-05e0-4ca3-98f0-3e704502b3aa
# ╠═40185644-29c2-423f-9cfc-b546d0351c57
# ╟─ecbd8b83-bddb-41d5-8e50-55fca4356e92
# ╟─11004fa9-d2f8-4a66-8ec2-416f2dd2653e
# ╟─729d9438-3c69-4596-b4e0-4d1b1f6227b2
# ╠═b5052fd0-d8d9-465e-90ed-94f9cec385c4
# ╟─7b164f17-30a6-4820-8419-7e16b2ff3d79
# ╟─1b1e4af6-379e-4d1b-9097-5fb29ea5c505
# ╟─8a6ffc9a-1845-41d7-b3ae-f3fa69d6b6ec
# ╠═eadfca3f-4c38-4fc3-aef9-a6eef900bb52
# ╟─d6b9635b-43cf-476a-9a91-b070b67332f5
# ╟─46f5473c-d6f9-48e6-bffa-c88ddf63a604
# ╟─919d504d-b67e-46a1-9883-5675e0d5b9b6
# ╟─28bc738b-7e46-4d84-a0d8-969150975f4d
# ╟─0d610ccf-c261-4db1-b20a-d9aabb1e7d96
# ╠═8c5058c3-a03f-4fb6-a8a9-e131b4059d90
# ╟─075e5e24-ad15-44da-86fe-1eef1a4c5ed0
# ╟─2b95cdbf-c860-477a-a178-1eb39b0e093e
# ╠═44a2dce8-b991-450d-85a9-90f49406f709
# ╠═ae4d01f1-1648-4b83-8288-75b1c18bf934
# ╟─a3f7b9a6-5eca-46af-b754-da2edaadd9a7
# ╟─0bc85a2d-bc2e-48f5-b5e5-13044f6c8673
# ╠═c72d4c62-39a9-4e60-9973-87013e8ca386
# ╟─491a73a4-27bd-4a3d-9254-d12d1a96a885
# ╟─d8ff7a10-3e34-4ccb-a318-71c6032265b5
# ╟─0c06bc9a-349a-439d-b682-ed7720fdb98b
# ╠═924f8479-bb95-47f7-83ab-478bd8f93923
# ╠═29df134f-9aa2-4583-bb42-c0d1c5c58ed0
# ╟─11583038-72cd-435c-9f8b-bbbe949ff8d5
# ╠═0d1eab07-8971-4ea2-ae14-4a2dcf6dd3be
# ╟─e8b9651e-8b73-4fd6-b0a2-080f29873133
# ╟─99357a1b-e2fc-4894-8a52-5cf9cd0fbf2c
# ╠═de21a988-515b-4e77-898f-f4363aeeaec1
# ╟─d0127175-bb82-44ff-9a58-2aa09848ab91
# ╟─7bc5802b-7519-4394-a630-f34f377186e2
# ╠═ca5d533a-4da5-4186-9397-a23388305060
# ╟─89d1f69c-8daa-4312-802e-68bc7d35a026
# ╟─6b2ddf23-4b75-40b5-a4ff-bfb696a8bb68
# ╠═cbab2437-39e5-4612-b362-ebc8d96c7db4
# ╟─8e14edcf-ff1c-4218-bc63-e0dba3e167e9
# ╟─e03d504d-cc90-4aef-856b-28ad65f1d500
# ╠═149d7c6e-5c53-4847-b989-e2e30bb5b9bb
# ╠═f1235f75-71ad-457b-8937-ac6ea8c3acba
# ╠═22c75efc-0d44-475a-846a-f48045862e2c
# ╠═1ce6cdab-6454-477c-98d8-9d0aa15835bc
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
