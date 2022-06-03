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

# ╔═╡ b61539e5-72e7-430a-91d2-20599314ec01
md" # Codificado de voz CELP

Gabriel Tavares 10773801

Guilherme Reis 10773700

"

# ╔═╡ bf5b4d85-52f1-4bd3-9b60-807d7a21fcff
md" ## Leitura do Sinal
"

# ╔═╡ 20ff7aa3-4256-4673-9142-cbc8df19fc1d
begin
	sinal, fs = wavread("antarctica.wav")
	plot(sinal)
	plot!(title="Sinal de som")
end

# ╔═╡ 91579668-2c98-4269-ba3c-a157eec9c6b9
SampleBuf(sinal,fs)

# ╔═╡ 25e2bde5-4092-4e5a-8c6c-db62953f9e5e
md" ## Análise e Sintese do sinal

Para fazer a reconstrução do sinal iremos dividí-lo em trechos de 240 amostras andando a passos de 80 amostras. A cada trecho de 240 amoostras iremos reconstruir as 80 amostras centrais desse trecho.

Em todos os trechos faremos o seguinte processo:

* Fazer a análise LPC do filtro de trato vocal daquele trecho
* Filtrar a família de funções bases por esse filtro 
* Passar as funções base pela função FindBest que acha a melhor correspondência entre as funções e o centro do trecho

Dessa forma para cada trecho teremos os coeficientes de trato vocal, os indices das melhores K funções da família de funções bases e os ganhos dessas funções.

"

# ╔═╡ 5ba134f1-75dc-4677-a6e3-9e0f0b1207fc
begin
	tamanho_analise = 240 
	tamanho_sintese = 80 
	off_set_analise = (Int)((tamanho_analise - tamanho_sintese)/2) #80 amostras pra tras e pra frente
end

# ╔═╡ e804c053-2059-4c01-a7b3-ef1038943eee
md" ## Resultados

Após a análise podemos fazer a reconstrução desse sinal com os parâmetros achados em cada trecho.
"

# ╔═╡ f019b8d7-12b4-46dd-9b63-8cdf976ce04b
md" Podemos ver que o sinal sintetizado se assemelha bastante ao sinal original como o desejado, e o tamanho do sinal é bastante menor, pois só possui os parâmtros da sintese"

# ╔═╡ ddf2261a-50e1-450f-9864-6ade9998ba41
md" # Funções"

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

# ╔═╡ 8bb89636-dea1-4140-a3a5-598d7e3dd2a1
noprint = md""

# ╔═╡ a8ecb7fc-4626-41d8-8d14-aaffe14155ad
begin
	K = 2
	Q = 512
	N = tamanho_sintese
	func_base = randn(N,Q)
	func_filtradas = zeros(N,Q)
	noprint
end

# ╔═╡ 74a87a9f-6637-410e-a931-242f09887bc0
begin	
	sinal_analise = vcat(zeros(off_set_analise), sinal, zeros(tamanho_analise - mod(length(sinal),tamanho_analise)), zeros(off_set_analise)) 
	
	sinal_sintese =  zeros(length(sinal_analise))
	
	for i in range(1,length = (Int)(length(sinal_analise)/tamanho_sintese)-2)
		
		trecho_analise = sinal_analise[i*tamanho_sintese+1-off_set_analise:(i+1)*tamanho_sintese+off_set_analise]
		
		if sum(trecho_analise) != 0
			
			ak, G = lpc(trecho_analise .* hamming(tamanho_analise),10)

			
			subquadro = trecho_analise[off_set_analise+1:off_set_analise+tamanho_sintese]
			
			
			filtro = PolynomialRatio([1],[1;ak])
			
			for coluna in 1:Q
				func_filtradas[:, coluna] = filt(filtro, func_base[: ,coluna])
			end

			ganhos, indices = find_Nbest_components(subquadro, func_filtradas, K);
			
			trecho_sintese = func_filtradas[:, indices[1]] * ganhos[1] + func_filtradas[:, indices[2]] * ganhos[2]
			
			sinal_sintese[i*tamanho_sintese+1:(i+1)*tamanho_sintese] = trecho_sintese
		else
			
		end

	end
	
end

# ╔═╡ 8324941b-a582-4e42-992f-2ab7150d6e0c
SampleBuf(sinal_sintese, fs)

# ╔═╡ 597e71b3-de00-4582-8760-3d8e385cb3b2
begin 
	plot(sinal, label= "Sinal Original")
	plot!(sinal_sintese, label ="Sinal Sintetizado", alpha = 0.7)
	plot!(title="Comparação dos sinais")
end

# ╔═╡ 7b108c8b-9b05-4454-9340-0f4c0ccaae5d
wavwrite(sinal_sintese, "sintese_CELP.wav", Fs = fs)

# ╔═╡ Cell order:
# ╟─b61539e5-72e7-430a-91d2-20599314ec01
# ╟─03a4d57c-d709-11ec-223b-1fed26b50381
# ╟─bf5b4d85-52f1-4bd3-9b60-807d7a21fcff
# ╟─20ff7aa3-4256-4673-9142-cbc8df19fc1d
# ╟─91579668-2c98-4269-ba3c-a157eec9c6b9
# ╟─25e2bde5-4092-4e5a-8c6c-db62953f9e5e
# ╠═5ba134f1-75dc-4677-a6e3-9e0f0b1207fc
# ╠═a8ecb7fc-4626-41d8-8d14-aaffe14155ad
# ╠═74a87a9f-6637-410e-a931-242f09887bc0
# ╟─e804c053-2059-4c01-a7b3-ef1038943eee
# ╠═8324941b-a582-4e42-992f-2ab7150d6e0c
# ╟─597e71b3-de00-4582-8760-3d8e385cb3b2
# ╠═7b108c8b-9b05-4454-9340-0f4c0ccaae5d
# ╟─f019b8d7-12b4-46dd-9b63-8cdf976ce04b
# ╟─ddf2261a-50e1-450f-9864-6ade9998ba41
# ╠═aba080c7-314a-498d-b007-338ed8af36af
# ╠═8bb89636-dea1-4140-a3a5-598d7e3dd2a1
