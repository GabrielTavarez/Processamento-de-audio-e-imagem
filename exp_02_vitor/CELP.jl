### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 03a4d57c-d709-11ec-223b-1fed26b50381
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

# ╔═╡ 20ff7aa3-4256-4673-9142-cbc8df19fc1d
begin
	sinal, fs = wavread("antarctica.wav")
	plot(sinal)
end

# ╔═╡ 91579668-2c98-4269-ba3c-a157eec9c6b9
SampleBuf(sinal,fs)

# ╔═╡ 5ba134f1-75dc-4677-a6e3-9e0f0b1207fc
begin
	tamanho_analise = 240 
	tamanho_sintese = 80 
	off_set_analise = (Int)((tamanho_analise - tamanho_sintese)/2) #80 amostras pra tras e pra frente
end

# ╔═╡ aba080c7-314a-498d-b007-338ed8af36af


# ╔═╡ 6b9fac24-54cd-4864-a523-9ff5b544b0da
begin
	Q = randn(tamanho_sintese, 512) #sinais são as colunas/5
	zs = zeros(10)
end

# ╔═╡ 74a87a9f-6637-410e-a931-242f09887bc0
begin	
	sinal_analise = vcat(zeros(off_set_analise),sinal, zeros(tamanho_analise - mod(length(sinal),tamanho_analise)), zeros(off_set_analise)) 
	
	sinal_sintese =  zeros(length(sinal_analise))
	
	for i in range(1,length = (Int)(length(sinal_analise)/tamanho_sintese)-2)
		
		trecho_analise = sinal_analise[i*tamanho_sintese+1-off_set_analise:(i+1)*tamanho_sintese+off_set_analise]
		
		if sum(trecho_analise) != 0
			
			ak, G = lpc(trecho_analise .* hamming(tamanho_analise),10)

			subquadro = trecho_analise[off_set_analise+1:off_set_analise+1+tamanho_sintese]
			
			
			filtro = PolynomialRatio([1],[1;ak])
			Qfilt = filt(filtro,Q)
			
			
			
		else

		end

	end
	
end

# ╔═╡ Cell order:
# ╠═03a4d57c-d709-11ec-223b-1fed26b50381
# ╠═20ff7aa3-4256-4673-9142-cbc8df19fc1d
# ╠═91579668-2c98-4269-ba3c-a157eec9c6b9
# ╠═5ba134f1-75dc-4677-a6e3-9e0f0b1207fc
# ╠═74a87a9f-6637-410e-a931-242f09887bc0
# ╠═aba080c7-314a-498d-b007-338ed8af36af
# ╠═6b9fac24-54cd-4864-a523-9ff5b544b0da
