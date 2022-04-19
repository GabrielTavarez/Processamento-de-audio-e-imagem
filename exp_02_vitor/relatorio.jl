### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 662faab4-c00c-11ec-0d8b-b597dda91e9a
begin
	using Pkg
	using DSP
	using Plots
	using WAV
	using FFTW
	plotly()
end

# ╔═╡ 9d31b63a-ecd6-4e83-9268-24ece795133d
begin
	sinal, fs = wavread("antarctica.wav")
	
end

# ╔═╡ 89bea7bb-22d4-473a-9050-23c293e9dd68
begin
	#Estimação dos coeficientes do filtro
	trecho = sinal[200:439]
	ak, pot_erro = lpc(trecho.*hamming(240), 10)
end

# ╔═╡ 30751353-0395-4590-a811-4af7da0da3e8
begin
	filtro = PolynomialRatio([1],[1;ak]) 
	ω = range(0,π, length= 512)
	H = freqz(filtro, ω)
end

# ╔═╡ 2a22f0a7-f831-46ac-8123-9e3d102bbaca
begin
	#peridiograma do trecho
	per = periodogram(trecho; fs = fs, nfft = 512)
	
end

# ╔═╡ 86779651-51e3-4f48-85a3-f12665891430
G=sqrt(pot_erro)

# ╔═╡ afd543ea-1a6c-43e0-83a1-d2461908b885
sum((G*abs.(H).^2))/π

# ╔═╡ d301cb9e-e5b1-4801-8dd3-41f39d4266f9
sum(per.power/π)

# ╔═╡ 7d79fb6e-2988-48f4-a7f1-d80cb07aeede
pot_erro

# ╔═╡ ecd464d8-481b-462d-bcb8-eac240a719ac
fft(trecho)

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

# ╔═╡ d45b68d0-963f-4c09-a611-2d7475d912e9
function pow2db(x)
	
	return 20log10(x)
end

# ╔═╡ 545c4793-a255-457f-9a5f-80549b992de3
begin
	plot(ω/π*fs/2, pow2db.(G*abs.(H).^2))
	plot!(per.freq, pow2db.(per.power*fs/π))
end

# ╔═╡ Cell order:
# ╠═662faab4-c00c-11ec-0d8b-b597dda91e9a
# ╠═9d31b63a-ecd6-4e83-9268-24ece795133d
# ╠═89bea7bb-22d4-473a-9050-23c293e9dd68
# ╠═30751353-0395-4590-a811-4af7da0da3e8
# ╠═2a22f0a7-f831-46ac-8123-9e3d102bbaca
# ╠═86779651-51e3-4f48-85a3-f12665891430
# ╠═545c4793-a255-457f-9a5f-80549b992de3
# ╠═afd543ea-1a6c-43e0-83a1-d2461908b885
# ╠═d301cb9e-e5b1-4801-8dd3-41f39d4266f9
# ╠═7d79fb6e-2988-48f4-a7f1-d80cb07aeede
# ╠═ecd464d8-481b-462d-bcb8-eac240a719ac
# ╟─8cd21f1d-b1c7-44ab-a726-0a0a41908a56
# ╠═d45b68d0-963f-4c09-a611-2d7475d912e9
