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

# ╔═╡ 10d943de-ddc6-11ec-1c31-adfb7f99248c
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
end

# ╔═╡ 2bcb6054-2b66-43e1-818c-0b94d9e218cc
md"
Gabriel Tavares

nusp : 10773801
"

# ╔═╡ 2e5988b3-5b2e-46e7-91df-154c48c744e3
md" ## Leitura do sinal"

# ╔═╡ d472d46a-66ae-4b51-978a-12652cdb191a
begin
	sinal, fs = wavread("antarctica.wav")
	SampleBuf(sinal, fs)
end

# ╔═╡ 69d9b4fc-a951-49d1-a210-f06eb6e3ef0f
@bind vel_reproducao Slider(0.5:0.1:3.0, default = 2.0)

# ╔═╡ 47e07e28-b664-46bd-bbaf-31351034899b
md"Velocidade = **$vel_reproducao**"

# ╔═╡ af80d64e-eca7-4847-8dfc-3c1d17cb62fa
md" # Mudança de taxa de amostragem


Aqui iremos alterar a velocidade de reprodução a partir da taxa de amostragem do sinal. Dessa forma o sinal irá tocar mais devagar ou mais rápido conforme a velocidade desejada"

# ╔═╡ 730c96a2-deba-411b-a1a3-99668cd98afe
begin
	wavwrite(sinal, "sinal_x05_taxaamostrage.wav", Fs = 0.5fs)
	SampleBuf(sinal, vel_reproducao*fs)
end

# ╔═╡ 29016173-efb9-4bcc-8232-109308377996
md" Vemos que a velocidade do sinal se altera, mas a tonalidade do sinal também é modificada. Quando aumentamos a velocidade, a voz fica mais aguda e quando diminuimos a velocidade, ela fica mais grave"

# ╔═╡ 4c6895da-95ed-46cb-af64-2721db2affe9
md" # Mudança do codificação de voz


Para corrigir o problema a cima, a estratégia é usar um algoritmo de reconstrução de voz.

O algoritimo irá reconstruir trechos maiores ou menos de voz, dependendo da velocidade desejada. 

Por exemplo, se queremos uma velocidade de x2.0, o algoritimo irá analisar trechos de 240 amostras a passos de 80 amostras, mas irá reconstruir apenas 40 amostras por passo (e não 80). Dessa forma o tamanho do sinal irá diminuir (ele ira tocar mais rápido) mas a tonalidade da voz não irá se alterar, porque os fonemas ainda estão sendo analisados e reconstruídos conforme a velocidade original.
"

# ╔═╡ 54724dc7-8711-467b-b852-3e1de5825c87
begin
	tamanho_analise = 240 
	tamanho_sintese = floor(Int,80/vel_reproducao)
	passo_analise = 80 #quantas amostras andam a cada trecho
	sobreposicao = 80  #numero de amostras de sobreposicao entre analises
						 # significa 80 amostras do trecho anterior e do futuro
end

# ╔═╡ 6fbbfb97-1286-4c62-886c-4fee20d9faf2
begin
	K = 2
	Q = 512
	N = tamanho_sintese
	func_base = randn(N,Q)
	func_filtradas = zeros(N,Q)
	md"" #noprint
end

# ╔═╡ 6683e551-1634-4e3d-82a3-5eb91282165b
begin
	#adiciona zeros no sinal para deixar o número certo de analises no sinal
	sinal_analise = vcat(zeros(passo_analise), sinal, zeros(passo_analise - mod(length(sinal),passo_analise)), zeros(passo_analise)) 

	num_trechos = floor(Int,length(sinal_analise)/(tamanho_analise - 2*sobreposicao))

	sinal_sintese =  zeros(tamanho_sintese*(num_trechos-2))
	md""
end

# ╔═╡ 3a05cdad-690e-423b-a83b-f53360ffcee2
begin
	wavwrite(sinal_sintese, "sinal_x2_sintese.wav", Fs= fs)
	SampleBuf(sinal_sintese, fs)
end	

# ╔═╡ 54557982-16e4-4c33-a2b0-476b8a2dd3c3
md"Observamos que o sinal foi acelerado, mas manteve a tonalidade da voz do locutor, como desejado."

# ╔═╡ 83d95f65-cee9-499d-90c4-6cc4a56560be
md"# Functions"

# ╔═╡ a814c82a-5cc2-4f5a-a04f-52a4dfe7d2bc
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

# ╔═╡ 4866125f-b023-4b4f-8b5d-612f1577e1aa
begin	
	for i in 1:num_trechos-3
		
		trecho_analise = sinal_analise[i*passo_analise+1-sobreposicao:(i+1)*passo_analise+sobreposicao]
		
		if sum(trecho_analise) != 0
			
			ak, G = lpc(trecho_analise .* hamming(tamanho_analise),10)

			
			subquadro = 
trecho_analise[floor(Int,(tamanho_analise -tamanho_sintese)/2) + 1 : floor(Int,(tamanho_analise -tamanho_sintese)/2)+tamanho_sintese]
			
			
			filtro = PolynomialRatio([1],[1;ak])
			
			for coluna in 1:Q
				func_filtradas[:, coluna] = filt(filtro, func_base[: ,coluna])
			end
						
			ganhos, indices = find_Nbest_components(subquadro, func_filtradas, K);
			
			trecho_sintese = func_filtradas[:, indices[1]] * ganhos[1] + func_filtradas[:, indices[2]] * ganhos[2]
			
			sinal_sintese[i*tamanho_sintese+1:(i+1)*tamanho_sintese] = trecho_sintese
			
		end

	end
	
end

# ╔═╡ Cell order:
# ╟─2bcb6054-2b66-43e1-818c-0b94d9e218cc
# ╠═10d943de-ddc6-11ec-1c31-adfb7f99248c
# ╠═2e5988b3-5b2e-46e7-91df-154c48c744e3
# ╠═d472d46a-66ae-4b51-978a-12652cdb191a
# ╟─47e07e28-b664-46bd-bbaf-31351034899b
# ╟─69d9b4fc-a951-49d1-a210-f06eb6e3ef0f
# ╟─af80d64e-eca7-4847-8dfc-3c1d17cb62fa
# ╠═730c96a2-deba-411b-a1a3-99668cd98afe
# ╟─29016173-efb9-4bcc-8232-109308377996
# ╟─4c6895da-95ed-46cb-af64-2721db2affe9
# ╠═54724dc7-8711-467b-b852-3e1de5825c87
# ╠═6fbbfb97-1286-4c62-886c-4fee20d9faf2
# ╠═6683e551-1634-4e3d-82a3-5eb91282165b
# ╠═4866125f-b023-4b4f-8b5d-612f1577e1aa
# ╠═3a05cdad-690e-423b-a83b-f53360ffcee2
# ╟─54557982-16e4-4c33-a2b0-476b8a2dd3c3
# ╟─83d95f65-cee9-499d-90c4-6cc4a56560be
# ╠═a814c82a-5cc2-4f5a-a04f-52a4dfe7d2bc
