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
	using LinearAlgebra
	using FixedPoint
	plotly()
	
	include("fxfilt.jl")
end

# ╔═╡ 747c21d1-51e2-4d15-bd13-48665637b1c4
Pkg.add("FixedPoint")

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

# ╔═╡ a8ecb7fc-4626-41d8-8d14-aaffe14155ad
begin
	K = 2
	Q = 512
	N = tamanho_sintese
	func_base = randn(N,Q)
	func_filtradas = zeros(N,Q)
	
	zs = zeros(10)
end

# ╔═╡ ae069d6a-6148-4528-99ae-7d8058e8c073
tamanho_analise - mod(length(sinal),tamanho_analise) + 80

# ╔═╡ 16771c79-3b50-40ea-b029-6f80310a0aac
minimum(func_base)

# ╔═╡ aba080c7-314a-498d-b007-338ed8af36af


"""
    function find_Nbest_components(s, codebook_vectors, N)
Adaptado de T. Dutoit (2009)
Acha os N melhores componentes de s a partir dos vetores no livro-código codebook_vectors, minimizando o quadrado do erro erro = s - codebook_vectors[indices]*ganhos.
Retorna (ganhos, indices)
"""
function find_Nbest_components(signal, codebook_vectors, N)

    M, L = size(codebook_vectors)
    codebook_norms = zeros(L)
    
    for j = 1:L
            codebook_norms[j] = norm(codebook_vectors[:,j])
    end

    gains = zeros(N)
    indices = ones(Int, N)

    for k = 1:N
        max_norm = 0.0
        for j = 1:L
            beta = codebook_vectors[:,j] ⋅ signal
            if codebook_norms[j] != 0
                component_norm = abs(beta)/codebook_norms[j]
            else
                component_norm = 0.0
            end
            if component_norm  > max_norm
                gains[k] = beta/(codebook_norms[j]^2)
                indices[k] = j
                max_norm = component_norm 
            end
        end
        signal = signal - gains[k]*codebook_vectors[:,indices[k]]   
    end
    return gains, indices
end

# ╔═╡ 74a87a9f-6637-410e-a931-242f09887bc0
begin	
	sinal_analise = vcat(zeros(off_set_analise), sinal, zeros(tamanho_analise - mod(length(sinal),tamanho_analise)), zeros(off_set_analise)) 
	
	sinal_sintese =  zeros(length(sinal_analise))
	
	for i in range(1,length = (Int)(length(sinal_analise)/tamanho_sintese)-2)
		
		trecho_analise = sinal_analise[i*tamanho_sintese+1-off_set_analise:(i+1)*tamanho_sintese+off_set_analise]
		
		if sum(trecho_analise) != 0
			
			ak, G = lpc(trecho_analise .* hamming(tamanho_analise),10)

			# aq = Fixed{Int16, 7-1}.(ak)
			
			subquadro = trecho_analise[off_set_analise+1:off_set_analise+tamanho_sintese]
			
			
			filtro = PolynomialRatio([1],[1;ak])
			
			for coluna in 1:Q
				func_filtradas[:, coluna] = filt(filtro, func_base[: ,coluna])
			end
			
			y0 = filt(filtro,zeros(tamanho_sintese))

			e0 = subquadro .- y0
			
			ganhos, indices = find_Nbest_components(e0, func_filtradas, K);
			
			trecho_sintese = func_filtradas[:, indices[1]] * ganhos[1] + func_filtradas[:, indices[2]] * ganhos[2]
			
			sinal_sintese[i*tamanho_sintese+1:(i+1)*tamanho_sintese] = trecho_sintese
		else
			
		end

	end
	
end

# ╔═╡ 8324941b-a582-4e42-992f-2ab7150d6e0c
SampleBuf(sinal_sintese, fs)

# ╔═╡ 7b108c8b-9b05-4454-9340-0f4c0ccaae5d
wavwrite(sinal_sintese, "sintese_CELP.wav", Fs = fs)

# ╔═╡ a26d407d-5f31-4a44-8a1b-8fdd1b7d7219
begin
	ak = 0
	i = 20
	trecho_analise = sinal_analise[i*tamanho_sintese+1-off_set_analise:(i+1)*tamanho_sintese+off_set_analise]

		if sum(trecho_analise) != 0

			ak, G = lpc(trecho_analise .* hamming(tamanho_analise),10)

			# aq = Fixed{Int16, 7-1}.(ak)

			subquadro = trecho_analise[off_set_analise+1:off_set_analise+tamanho_sintese]


			filtro = PolynomialRatio([1],[1;ak])

			for coluna in 1:Q
				func_filtradas[:, coluna] = filt(filtro, func_base[: ,coluna])
			end

			y0 = filt(filtro,zeros(tamanho_sintese))

			e0 = subquadro .- y0

			ganhos, indices = find_Nbest_components(e0, func_filtradas, K);

			trecho_sintese = func_filtradas[:, indices[1]] * ganhos[1] + func_filtradas[:, indices[2]] * ganhos[2]

			sinal_sintese[i*tamanho_sintese+1:(i+1)*tamanho_sintese] = trecho_sintese
		else

		end
end

# ╔═╡ Cell order:
# ╠═03a4d57c-d709-11ec-223b-1fed26b50381
# ╠═747c21d1-51e2-4d15-bd13-48665637b1c4
# ╠═20ff7aa3-4256-4673-9142-cbc8df19fc1d
# ╠═91579668-2c98-4269-ba3c-a157eec9c6b9
# ╠═5ba134f1-75dc-4677-a6e3-9e0f0b1207fc
# ╠═a8ecb7fc-4626-41d8-8d14-aaffe14155ad
# ╠═74a87a9f-6637-410e-a931-242f09887bc0
# ╠═ae069d6a-6148-4528-99ae-7d8058e8c073
# ╠═8324941b-a582-4e42-992f-2ab7150d6e0c
# ╠═7b108c8b-9b05-4454-9340-0f4c0ccaae5d
# ╠═a26d407d-5f31-4a44-8a1b-8fdd1b7d7219
# ╠═16771c79-3b50-40ea-b029-6f80310a0aac
# ╠═aba080c7-314a-498d-b007-338ed8af36af
