### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 662faab4-c00c-11ec-0d8b-b597dda91e9a
																					begin
	using Pkg
	using DSP
	using Plots
	using WAV
	using FFTW
	using SampledSignals
	using PlutoUI
	plotly()
end

# ╔═╡ 979a7ac8-8940-4ff3-9a6f-3e300c347529
md" ## Leitura do sinal


"

# ╔═╡ 9d31b63a-ecd6-4e83-9268-24ece795133d
begin
	sinal, fs = wavread("antarctica.wav")
	plot(sinal)
end

# ╔═╡ fc7430c9-83b7-43ad-a63e-b644329f099a
SampleBuf(sinal,fs)

# ╔═╡ f3d9110c-eb91-4514-ab65-f849de960aec
md" ## Estimação do filtro de trato vocal"

# ╔═╡ 89bea7bb-22d4-473a-9050-23c293e9dd68
begin
	#Estimação dos coeficientes do filtro
	trecho = sinal[200:439]
	ak10, pot_erro = lpc(trecho.*hamming(240), 10)
	
	filtro10 = PolynomialRatio([1],[1;ak10]) 
	ω = range(0,π, length= 512)
	H = freqz(filtro10, ω)
	
	#peridiograma do trecho
	per = periodogram(trecho; fs = fs, nfft = 512)
end

# ╔═╡ 86779651-51e3-4f48-85a3-f12665891430
md"

$(@bind G_ Slider(0:0.01:0.2, default = 0.1))
"

# ╔═╡ 9d653a46-9e2b-4c6f-a5da-e5352d71f124
md"Ganho = $(G_)"

# ╔═╡ b643e40e-ab39-4f1d-b3a8-a3aee4774ec0
G_

# ╔═╡ b3a12d13-1def-4b46-ae9b-ce8133a9064c
md" #### Filtro com mais coeficientes"

# ╔═╡ d11ad970-ba00-46f1-9744-a40709f478a5
begin
	ak20, p = lpc(trecho.*hamming(240), 10)
	ak40, p = lpc(trecho.*hamming(240), 40)
	ak80, p = lpc(trecho.*hamming(240), 80)
	
	filtro20 = PolynomialRatio([1],[1;ak20]) 
	filtro40 = PolynomialRatio([1],[1;ak40]) 
	filtro80 = PolynomialRatio([1],[1;ak80]) 
	
	H20 = freqz(filtro20, ω)
	H40 = freqz(filtro40, ω)
	H80 = freqz(filtro80, ω)
end

# ╔═╡ 351d6990-3589-4705-978c-b3dccdf702ef
md" ### Testando estimador de Pitch"

# ╔═╡ 776174d2-0974-420b-950e-151d775a0399
SampleBuf(trecho, fs)

# ╔═╡ 1c567043-99bf-43f5-a9ae-29b0e81e539d
md" ## Reconstrução do sinal completo"

# ╔═╡ 90506c2e-5435-497a-8846-5c454b745f00
begin
	plot(sinal)
	plot!(title = "Sinal original")
end

# ╔═╡ f1f8eead-3d53-4114-b596-99d77bec2872
md" # Functions"

# ╔═╡ 8cd21f1d-b1c7-44ab-a726-0a0a41908a56
"""
    pitch(quadro)
Adaptado de T. Dutoit (2007)
Calcula o período fundamental de um trecho de voz de 30ms amostrado a 8kHz.  Retorna 0 se o trecho é surdo (não-vozeado).

"""
function pitch(quadro)
    ai,ξ = lpc(quadro,10)
    h = PolynomialRatio([1],[1;ai])
    lpc_residual=filt(h,quadro)
    M = length(lpc_residual)
    C=xcorr(lpc_residual, lpc_residual; padmode = :longest)
    Cxx=C[M:end] / C[M] # normaliza a autocorrelação pela variância
    Cxx[1:26] .= 0
    Amax,Imax = findmax(Cxx[1:min(133,length(quadro))]) # Limita T para 60-300Hz
    if (Amax > 0.20) # teste bem simples para sinal sonoro/surdo
       T0 = Imax - 1
    else 
       T0 = 0
    end

    return T0
end

# ╔═╡ 9d136466-a9b8-40e4-add4-34c28c7ded54
begin
	Tp = pitch(trecho)
	pulsos = zeros(length(trecho))
	pulsos[1:Tp:end] .= 1
	plot(pulsos, marker=:circle, line= :stem)
	plot!(title = "Gerador de Pulsos")
end

# ╔═╡ 5bd05fc1-3bc4-43a0-95d7-61a6a067d9b2
SampleBuf(filt(filtro10,pulsos), fs)

# ╔═╡ 9ded7d0f-d68b-4b4e-8795-a33fd8940137
begin
	
	tamanho_analise = 240 
	tamanho_sintese = 80 
	off_set_analise = (Int)((tamanho_analise - tamanho_sintese)/2) #80 amostras pra tras e pra frente
	
	
	
	sinal_analise = vcat(zeros(off_set_analise),sinal, zeros(tamanho_analise - mod(length(sinal),tamanho_analise)), zeros(off_set_analise)) 
	
	sinal_sintese =  zeros(length(sinal_analise))
	
	for i in range(1,length = (Int)(length(sinal_analise)/tamanho_sintese)-2)
		
		trecho_analise = sinal_analise[i*tamanho_sintese+1-off_set_analise:(i+1)*tamanho_sintese+off_set_analise]
		
		if sum(trecho_analise) != 0
			ak, pot = lpc(trecho_analise.*hamming(tamanho_analise), 10)
			tp = pitch(trecho_analise)

			filtro = PolynomialRatio([1],[1;ak]) 

			if tp == 0 
				G = sqrt(pot)
				excitacao = randn(tamanho_sintese)
			else
				G = sqrt(pot*tp)
				excitacao = zeros(tamanho_sintese)
				excitacao[1:Tp:end] .= 1
			end
			trecho_sintese = filt(G*filtro,excitacao)

			sinal_sintese[i*tamanho_sintese+1:(i+1)*tamanho_sintese] = trecho_sintese
		else
			sinal_sintese[i*tamanho_sintese+1:(i+1)*tamanho_sintese] .= 0
		end

	end
	
end

# ╔═╡ b3c1f863-c2dd-44d9-a09f-37af7f2720d8
begin
	plot(sinal_sintese)
	plot!(title = "Sinal reconstruido")
end

# ╔═╡ 5b15f3bd-7de5-48dc-a7e7-a2ca49c9ed67
begin 
	plot(sinal, label = "Sinal original")
	plot!(sinal_sintese, label = "Sinal reconstruido")
	plot!(title = "Sinais sobrepostos")
end

# ╔═╡ 4b0da572-7724-41ef-b862-c8956a7f10e6
md"
##### Sinal Original
$(SampleBuf(sinal, fs))

##### Sinal reconstruído

$(SampleBuf(sinal_sintese, fs))

"

# ╔═╡ efe4b6ab-1c2f-4222-a033-e56a4d17686c
wavwrite(sinal_sintese, "sintese_LPC.wav", Fs = fs)

# ╔═╡ d45b68d0-963f-4c09-a611-2d7475d912e9
function pow2db(x)
	
	return 20log10(x)
end

# ╔═╡ 545c4793-a255-457f-9a5f-80549b992de3
begin
	plot(ω/π*fs/2, pow2db.(G_*abs.(H).^2), label = "Filtro vocal")
	plot!(per.freq, pow2db.(per.power*fs/π), label = "Periodograma do sinal", color = :green)
	plot!(title="Filtro do trato vocal")
	
end

# ╔═╡ c0f03a26-d2ef-459d-8143-7e454e984e5d
begin
	plot(ω/π*fs/2, pow2db.(G_*abs.(H20).^2), label = "Filtro vocal 20 coeficiente", color="#de00de")
	plot!(ω/π*fs/2, pow2db.(G_*abs.(H40).^2), label = "Filtro vocal 40 coeficiente", color="#00dede")
	plot!(ω/π*fs/2, pow2db.(G_*abs.(H80).^2), label = "Filtro vocal 80 coeficiente", color="#dede00")
	plot!(per.freq, pow2db.(per.power*fs/π), label = "Periodograma do sinal", color = :green)
	plot!(title="Filtro do trato vocal")
	
end

# ╔═╡ Cell order:
# ╠═662faab4-c00c-11ec-0d8b-b597dda91e9a
# ╟─979a7ac8-8940-4ff3-9a6f-3e300c347529
# ╠═9d31b63a-ecd6-4e83-9268-24ece795133d
# ╠═fc7430c9-83b7-43ad-a63e-b644329f099a
# ╟─f3d9110c-eb91-4514-ab65-f849de960aec
# ╠═89bea7bb-22d4-473a-9050-23c293e9dd68
# ╟─9d653a46-9e2b-4c6f-a5da-e5352d71f124
# ╟─86779651-51e3-4f48-85a3-f12665891430
# ╠═b643e40e-ab39-4f1d-b3a8-a3aee4774ec0
# ╠═545c4793-a255-457f-9a5f-80549b992de3
# ╟─b3a12d13-1def-4b46-ae9b-ce8133a9064c
# ╠═d11ad970-ba00-46f1-9744-a40709f478a5
# ╟─c0f03a26-d2ef-459d-8143-7e454e984e5d
# ╟─351d6990-3589-4705-978c-b3dccdf702ef
# ╠═9d136466-a9b8-40e4-add4-34c28c7ded54
# ╠═5bd05fc1-3bc4-43a0-95d7-61a6a067d9b2
# ╠═776174d2-0974-420b-950e-151d775a0399
# ╟─1c567043-99bf-43f5-a9ae-29b0e81e539d
# ╠═9ded7d0f-d68b-4b4e-8795-a33fd8940137
# ╟─90506c2e-5435-497a-8846-5c454b745f00
# ╟─b3c1f863-c2dd-44d9-a09f-37af7f2720d8
# ╟─5b15f3bd-7de5-48dc-a7e7-a2ca49c9ed67
# ╟─4b0da572-7724-41ef-b862-c8956a7f10e6
# ╠═efe4b6ab-1c2f-4222-a033-e56a4d17686c
# ╟─f1f8eead-3d53-4114-b596-99d77bec2872
# ╠═8cd21f1d-b1c7-44ab-a726-0a0a41908a56
# ╠═d45b68d0-963f-4c09-a611-2d7475d912e9
