### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 5394d850-dddc-11ec-36df-09b9e40c37f7
begin
	using Plots
	using DSP
	using MAT
	using SampledSignals
	using WAV
	using Pkg
	using Statistics
	plotly()
end

# ╔═╡ 10f8b443-519a-4699-8712-ed4ebdf763ad
md" 
Gabriel Tavares Ferrarez

Nusp: 10773801
"

# ╔═╡ 371033e4-459f-4f73-b237-ee75cec05576
md" # Exercicio 1"

# ╔═╡ 1868856d-de56-4f0f-86c5-884ddf436b87
begin
	desejado1_file = matopen("desejado1.mat");
	desejado1 = read(desejado1_file, "d1")[1,:];
	close(desejado1_file);
	
	entrada1_file = matopen("entrada1.mat");
	entrada1 = read(entrada1_file, "u1")[1,:];
	close(entrada1_file);
end

# ╔═╡ 3615b7d8-30bb-484b-a7c2-b251bac43ce7
md" Primeiro vamos testar o algoritimo com um valor grande de coeficientes e com um valor pequeno de passo de adaptação.

M = 200 coeficientes

μ = 0.0001 

Vamos analisar se os coeficientes estão convergindo e quantos dos coeficiente são muito pequenos e vão influenciar menos na estimação do filtro.
"

# ╔═╡ 8d114f66-37ff-48f2-998e-f7cf3cf5adab
md" #### Estimação do passo de convergência

Nessa etapa vamos plotar 10 dos coeficientes e ver como eles se comportam com a mudança de μ.
"

# ╔═╡ 193383c7-79ca-4190-a5ff-8ed8908aef01
md"Nos gráficos vemos que os 3 comportamentos possíveis a partir do passo de convergência.

* No primeiro gráfico o valor de μ era muito pequeno e o algoritimo não chegou a convergir para o valor ótimo

* No último gráfico, o valor de μ era muito grande e o algoritimo não consegue convergir para nenhum valor

* No gráfico do meio o  algoritimo conseguiu convergir para valores ótimos do filtro, portanto o valor idela deve estar próximo dessa ordem de grandeza

Com isso em mente, vamos alterar o valor de μ manualmente até chegar em um valor ótimo.
"

# ╔═╡ dfcdea6c-9122-4cc6-8baa-3d0cac68b8a4
md" Após testes manuais, o valor escolhido para μ foi **μ = 0.0004**, dessa forma o algoritimo tem tempo de convergir e e mantém uma certa estabilidade"

# ╔═╡ 500f65db-f1a0-42f5-bfbb-1b245845154f
md" #### Estimação do número de coeficientes

Para estimar o número ideal de coeficientes, vamos plotar o valor final de todos os coeficientes e analisar quantos deles são pequenos e podem ser desconsiderados

"

# ╔═╡ c0b15d35-e917-404a-9214-86bbd0b20bcb
md"No gráfico vemos que parte dos coeficientes tem um valor bastante baixo e pode ser desconsiderado a fim de diminuir o tamanho do filtro. Os primeiros coeficientes estão em terno de 0.08, portanto iremos desconsiderar os coeficientes uma ordem de grandeza menor (<0.01). Isso nos deixa com algo em torno de **100 coeficientes**.

Para verificar a influência desses coeficientes, vamos plotar a potência do erro em trechos ao longo do tempo dos dois filtros e verificar que a potência dos dois diminui"

# ╔═╡ 95c405a3-52a1-490c-8ade-cac8e5196a87
md" No gráfico podemos ver a influência do número de coeficientes na saída do filtro. 

* Com 200 coeficientes, o gráfico convege para uma potência de erro baixa  e constante.

* Com 50 coeficientes, o gráfico chega a convergia para uma potência, mas com mais instabilidade.

Usando **100 coeficientes** como apontado antes, é possível reduzir o tamanho do filtro e ainda ter uma qualidade suficientemente boa para o filtro.
"

# ╔═╡ 10ac890e-6e45-4795-b32d-b89072735a52
md"
#### Filtro final"

# ╔═╡ 54ad474d-8c20-4385-816e-0d589b443e42
md" 
No final dessa etapa o filtro que simula o ambiente real é:"

# ╔═╡ dfd3a8d9-c986-4fcb-b7cf-b762acf81c88
md" #### Potência do ruido v[n]"

# ╔═╡ 9c01e7c1-5cd2-4b03-a7cd-1143f1c090a3
md"Para calcular a potência do ruído usaremos o fato que nosso filtro se aproxima do sistema real H.

Com isso temos que a saída Y do filtro adaptativo é muito proxima da saída X do sistema real. Portanto:

$E = Y - D =  Y - (V+X) = V$
portanto:

$V = Y - D$

Portanto para calcular a potência do ruído **v[n]** basta passarmos o sinal desejado pelos  coeficientes finais do filtro adptativo e calcular sua potência.
"

# ╔═╡ 4d952870-d87a-45fc-a480-8fd6c613a151
md" # Exercicio 2"

# ╔═╡ 3b8150c3-6260-4528-aa7b-332b2105c7b1
begin
	fa = 8_000
	desejado2_file = matopen("desejado2.mat");
	desejado2 = read(desejado2_file, "d2")[:,1]
	close(desejado2_file);
	
	entrada2_file = matopen("entrada2.mat");
	entrada2 = read(entrada2_file, "u2")[:,1]
	close(entrada2_file);
end

# ╔═╡ 8d2fac19-e6a8-430d-af5e-1d23cb026b48
md"
**Entrada**

$(SampleBuf(entrada2, fa))

**Eco**

$(SampleBuf(desejado2, fa))

"

# ╔═╡ a017765c-c31a-4e52-b2a7-7a8af0282b9c
md" #### a) Eliminação do eco. 

Para eliminar o eco iremos usar o LMS para tentar convergir o sistema real de que gerou o eco da voz. Com esse sistema estimado, iremos passar o sinal original da voz por ele, e subtrair do sinal com eco. Com isso, esperamos eliminar o eco e ficar apenas com o ruído do som desejado.


Usaremos o mesmo número de coeficientes do ultimo exercício, mas será necessário alterar o passo de adaptação. 

Teremos que modificar o passo de adaptação, porque os sinais de entrada mudaram. Isso faz com que o valores de correlação entre sinais que geram o filtro ótimo e consequentemente o valor máximo e ideal de μ também se alteram.
"

# ╔═╡ eef91f61-6ed9-4635-be7c-9bffcd37b994
md" Da mesma forma que o exercício anterior, se nosso filtro se aproxima do sistema real, temos:

$X = Y$

Potanto 

$E = Y-D = Y-(X+V) = V$

Sendo **X** o sinal de eco e **V** o sinal que queremos obter. Portanto

$V =  Y - D$
"

# ╔═╡ d70e75d0-f811-4bce-afff-d2b5a26caa0a
md" #### b) Potência do ruído de fundo "

# ╔═╡ 8acce428-b2e9-4bff-a298-53128a5d2c1f
md"Após eliminar boa parte do eco, o que sobre é aproximadamente apenas o ruído de fundo, então apenar precisamos calcular a potência do sinal n[n]"

# ╔═╡ c50a0aaf-9329-46fa-873e-e20aa9cd0562
md" #### c) Comparação dos sinais

Aqui iremos comparar o sinal de saída. Como o esperado, o sinal ainda possui um pouco da fala do narrador, mas essa fala está bastante atenuada devido a eliminação do filtro adaptativo. 

Somado a isso, ouvimo o ruído que calculamos a potência.
"

# ╔═╡ 3309d953-e26f-4a7e-a709-c359f2927950
md" #### d) Curva ERLE do sinal filtrado

Aqui podemos ver a eficiência do nosso eliminador de eco.  No gráfico abaixo vemos a quantidade de eco eliminado em cada trecho de audio. Podemos ver que nos trechos de som, a curva ELRE tem um valor alto indicando que eliminou ruído.
"

# ╔═╡ 04e0025d-5fe4-42cc-92e8-114c071c1722
md" #### e) Convergência dos coeficientes"

# ╔═╡ 7a0066ca-62fe-484e-be73-b14c4731ec7a
md" Vemos que os coeficientes realmente convergem para algo próximo. 

O objetivo desse filtro é simular um sistema real, portanto a razão mais provável do porque os dois filtros se aproximam é que o sistema real que os dois coeficientes tentam simular é o mesmo. Com isso o algoritimo irá fazer com que os coeficientes dos dois sistemas convirjam."

# ╔═╡ ed9b0a72-ec16-442a-b0b6-bb4206322c32
md" # Functions"

# ╔═╡ df70a5ce-c529-4a05-94d1-8305c64155e5
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

# ╔═╡ 04e8e0d7-2ae2-4f57-8c5a-8973b93f13e2
function pow2db(x)	
	20*log10(x)
end

# ╔═╡ 365f993f-7652-422f-864c-0087c9c3e0f3
noprint = md""

# ╔═╡ 5dbb787a-c7db-45bf-8be0-936d0350c761
begin
	M200 = 200 #num filtros
	N = length(entrada1) #num interações
	noprint
end

# ╔═╡ 8419582b-9bdf-4f71-8175-707c4bccd961
begin
	μ_peq = 1e-4
	μ_med = 1e-3
	μ_grande = 1e-2
	W_peq, erro_peq =  LMS(entrada1, desejado1, M200,N,μ_peq)
	W_med, erro_med =  LMS(entrada1, desejado1, M200,N,μ_med)
	W_grande, erro_grande =  LMS(entrada1, desejado1, M200,N,μ_grande)
	noprint
end

# ╔═╡ f8b4d2fc-b8d4-4c75-b436-f2a6430421cb
begin
	p_peq = plot(W_peq[1:100:end,1:10], legend = false, title = "μ = 0.0001")
	p_med = plot(W_med[1:100:end,1:10], legend = false, title = "μ = 0.001")
	p_grande = plot(W_grande[1:100:end,1:10], legend = false,title = "μ = 0.01")
	plot(p_peq, p_med, p_grande, layout=(3,1))
end

# ╔═╡ e95fd133-af21-4bfb-ae41-2630d6b52db6
begin
	μ = 4e-4

	W_μ_otimo, erro_μ_otimo =  LMS(entrada1, desejado1, M200,N,μ)
	noprint
end

# ╔═╡ 4dfcaeaf-4671-4382-82f4-50d2739df6ab
begin
	#plots feitos a cada 10 interações, para não sobrecarregar o plot desnecessariamente
	plot(W_μ_otimo[1:10:end,1:10], legend = false, title = "μ = $μ")
end

# ╔═╡ 5da49bc3-a442-4df5-bd99-7a79a71baa27
begin
	plot(sort(abs.(W_μ_otimo[end,:]), rev = true), marker=:circle, line=:stem)
	plot!(title = "Coeficientes ordenados", legend = false)
end

# ╔═╡ e7364c37-a550-4829-abf0-b4e5952b1314
begin
	M = 100
	M50 = 50
	W1, erro1 =  LMS(entrada1, desejado1, M,N,μ)
	W_50, erro_50 =  LMS(entrada1, desejado1, M50,N,μ)
	noprint	
end

# ╔═╡ 4c8334ae-24ce-42b7-99ee-63088f1c62cd
begin
	
	filtro = PolynomialRatio(W1[end, :],[1])

	ω = range(1, π, length = 1000)
	resp_freq = freqz(filtro, ω)
	
	p_filtro = plot(W1[end,:], line =:stem, marker = :circle, title =  "Resposta ao impulso")
	p_resp = plot(ω/π, pow2db.(abs.(resp_freq)), title = "Resp em frequência")
	
	
	plot(p_filtro, p_resp, layout=(2,1))
end

# ╔═╡ 4e118182-b843-47ff-ac33-0ab19e241e6c
begin
	filter = PolynomialRatio(W1[end, :],[1])
	y = filt(filter, entrada1)
	v = y - desejado1
	pot_v = sum(v.^2)/length(v)
end

# ╔═╡ e32611e0-52d8-4caa-ae99-73f84b53fc5e
md" Pot de v[n] = **$(round(pot_v, digits=4))**"

# ╔═╡ e057970c-89cb-4d0d-8fe2-8ad00e600b45
begin
	SampleBuf(v,8_000)
end

# ╔═╡ 7ce7bc6e-60bc-49da-8b20-9838bdf4aad7
begin
	
	pot_erro = []
	for i in 0:49
		tamanho_trecho = 200
		pot_trecho = sum(erro1[i*tamanho_trecho+1:(i+1)*tamanho_trecho].^2)/tamanho_trecho
		append!(pot_erro,pot_trecho)
	end
	
	pot_erro_50 = []
	for i in 0:49
		tamanho_trecho = 200
		pot_trecho = sum(erro_50[i*tamanho_trecho+1:(i+1)*tamanho_trecho].^2)/tamanho_trecho
		append!(pot_erro_50,pot_trecho)
	end
	
	pot_erro_200 = []
	for i in 0:49
		tamanho_trecho = 200
		pot_trecho = sum(erro_μ_otimo[i*tamanho_trecho+1:(i+1)*tamanho_trecho].^2)/tamanho_trecho
		append!(pot_erro_200,pot_trecho)
	end
	noprint
end

# ╔═╡ f842eda6-d38d-4b18-90fc-0b28281cfee9
begin
	plot(pot_erro_200, label= "M = 200 coeficientes")
	plot!(pot_erro, label= "M = 100 coeficientes")
	plot!(pot_erro_50, label= "M = 500 coeficientes")
	plot!(title = "Potência do erro ao longo do tempo para diferentes coeficientes")
	plot!(legend = (0,0))
end

# ╔═╡ 96d95e66-9c9d-421b-8ded-1c70b3f09009
begin
	N2 = length(entrada2) #190_000 pontors
	μ2 = 1e-2
	W2, erro2 = LMS(entrada2, desejado2, M, N2, μ2)
	noprint
end

# ╔═╡ 662f022d-32e3-4b64-b390-74f20c23a190
begin
	p_μ_peq = plot(LMS(entrada2, desejado2, M, N2, 1e-3)[1][1:1000:end,1], title= "μ = 0.001")
	p_μ_medio = plot(LMS(entrada2, desejado2, M, N2, 1e-2)[1][1:1000:end,1], title= "μ = 0.01")
	p_μ_grande = plot(LMS(entrada2, desejado2, M, N2, 1e-1)[1][1:1000:end,1], title= "μ = 0.1")
	plot(p_μ_peq, p_μ_medio, p_μ_grande, layout = (3,1))
end

# ╔═╡ 2905d396-c0e1-4a73-a680-bc2ba1faba75
md" Empiricamente observamos que com 

* μ = 0.001, o sistema não tem tempo de convergir
* μ = 0.1, o sistema não converge
* μ = 0.01, o sistema converge bem, portanto o valor ótimo de μ deve estar nessa ordem de grandeza

O μ escolhido foi **μ = $μ2**
"

# ╔═╡ 09d54b52-66b5-491c-a540-537ffe5cbecc
begin
	plot(W2[1:1000:end,:], title= "Convergência dos coeficientes com μ = $μ2", legend = false)
end

# ╔═╡ 2c56f1d2-8e07-447b-a9d8-953ced9bd12d
begin
	plot(W1[end,:], line=:stem, marker=:circle, label="Coeficientes ex1")
	plot!(W2[end,:], line=:stem, marker=:circle, label="Coeficientes ex2")
end

# ╔═╡ 61f97b3a-c503-43aa-91f5-ac21cd505625
begin
	filtro2 = PolynomialRatio(W2[170_000,:], [1])
	y2 = filt(filtro2, entrada2)
	v2 = desejado2 .- y2
	noprint
end

# ╔═╡ 1e7f4850-059d-4cbe-82da-3f8c107f41e6
begin
	plot(desejado2, label = "sinal com eco", alpha = 0.8)
	plot!(v2, label = "sinal sem eco")
	plot!(title = "Comparação do sinal filtrado")
end

# ╔═╡ ca7be7c0-6798-4af4-a879-6e4edf073093
begin
	pot_ruido2 = sum(v2.^2)/length(v2)
end

# ╔═╡ 7a46af1e-51e6-4830-aa77-5a93893de833
md"pot ruído = **$(round(pot_ruido2, digits = 6))**"

# ╔═╡ 5ba0b9d7-caaa-4c2b-b43e-ae1fcbf1bd98
md"
**Sinal com eco**

$(SampleBuf(desejado2, fa))

**Sinal sem eco**

$(SampleBuf(v2, fa))"

# ╔═╡ 46dc0fae-26ec-4a57-9077-1f2421b39cfc
wavwrite(v2, "sinal_sem_eco.wav", Fs = fa)

# ╔═╡ 7baedb85-1b9a-449e-bb73-0c5f1068da3a
function ERLE(d,e,Nw,fa)
# %
# % Função para cálcudo do ERLE 
# % (echo return loss enhancement)
# % após o cancelador de eco
# % 
# % Saída 
# % Estimativa do ERLE em dB
# %
# % Entradas
# % d: sinal de eco (desejado)
# % e: sinal de eco residual (erro do filtro adaptativo)
# % Nw: número de amostras na janela para estimativa do eco ao longo do tempo
# % fa: frequência de amostragem 
# %
# %

# % MTMS, 08/2018

N= length(d);
Nb= (Int)(floor(N/Nw));
ERLEx = zeros(1,Nb);
Ta= 1/fa ;

for i in 1:Nb
   l=Nw*(i-1)+1:Nw*i;
   ERLEx[i] = mean(d[l].^2)/mean((e[l].+eps()).^2);
end
ERLEdB=10*log10.(ERLEx);

return ERLEdB[1,:]

# figure
# t=0:(N-1)*Ta/(N/Nw-1):(N-1)*Ta;
# plot(t,ERLEdB)
# grid on
# hold on
# ylabel('     Eco                                      ERLE (dB)')
# xlabel('tempo (s)')
# td=t(1):(t(end)-t(1))/(length(d)-1):t(end);
# plot(td,10*d/max(abs(d))-11)
# plot([td(1) td(end)],[-0.5 -0.5],'k','LineWidth',1)
# set(gca,'YTick',[0:10:max(ERLEdB), ceil(max(ERLEdB))]');
# axis([td(1) td(end) -21 ceil(max(ERLEdB))])
end

# ╔═╡ bd27bd51-1d78-4980-8de3-cfa3a6edd133
begin
	p_erle = plot(ERLE(desejado2, v2, 1024, fa), title = "ERLE(s)")
	p_som = plot(desejado2, title= "Som com eco")
	plot(p_erle, p_som, layout = (2,1), legend = false)
end

# ╔═╡ 1f27650d-02b9-42ea-9fb4-e4d2e7f49cb7


# ╔═╡ Cell order:
# ╠═5394d850-dddc-11ec-36df-09b9e40c37f7
# ╟─10f8b443-519a-4699-8712-ed4ebdf763ad
# ╟─371033e4-459f-4f73-b237-ee75cec05576
# ╠═1868856d-de56-4f0f-86c5-884ddf436b87
# ╟─3615b7d8-30bb-484b-a7c2-b251bac43ce7
# ╠═5dbb787a-c7db-45bf-8be0-936d0350c761
# ╟─8d114f66-37ff-48f2-998e-f7cf3cf5adab
# ╠═8419582b-9bdf-4f71-8175-707c4bccd961
# ╟─f8b4d2fc-b8d4-4c75-b436-f2a6430421cb
# ╟─193383c7-79ca-4190-a5ff-8ed8908aef01
# ╠═e95fd133-af21-4bfb-ae41-2630d6b52db6
# ╟─4dfcaeaf-4671-4382-82f4-50d2739df6ab
# ╟─dfcdea6c-9122-4cc6-8baa-3d0cac68b8a4
# ╟─500f65db-f1a0-42f5-bfbb-1b245845154f
# ╟─5da49bc3-a442-4df5-bd99-7a79a71baa27
# ╟─c0b15d35-e917-404a-9214-86bbd0b20bcb
# ╠═e7364c37-a550-4829-abf0-b4e5952b1314
# ╠═7ce7bc6e-60bc-49da-8b20-9838bdf4aad7
# ╟─f842eda6-d38d-4b18-90fc-0b28281cfee9
# ╟─95c405a3-52a1-490c-8ade-cac8e5196a87
# ╟─10ac890e-6e45-4795-b32d-b89072735a52
# ╟─54ad474d-8c20-4385-816e-0d589b443e42
# ╟─4c8334ae-24ce-42b7-99ee-63088f1c62cd
# ╟─dfd3a8d9-c986-4fcb-b7cf-b762acf81c88
# ╟─9c01e7c1-5cd2-4b03-a7cd-1143f1c090a3
# ╠═4e118182-b843-47ff-ac33-0ab19e241e6c
# ╟─e32611e0-52d8-4caa-ae99-73f84b53fc5e
# ╟─e057970c-89cb-4d0d-8fe2-8ad00e600b45
# ╟─4d952870-d87a-45fc-a480-8fd6c613a151
# ╠═3b8150c3-6260-4528-aa7b-332b2105c7b1
# ╟─8d2fac19-e6a8-430d-af5e-1d23cb026b48
# ╟─a017765c-c31a-4e52-b2a7-7a8af0282b9c
# ╟─662f022d-32e3-4b64-b390-74f20c23a190
# ╟─2905d396-c0e1-4a73-a680-bc2ba1faba75
# ╠═96d95e66-9c9d-421b-8ded-1c70b3f09009
# ╟─09d54b52-66b5-491c-a540-537ffe5cbecc
# ╟─eef91f61-6ed9-4635-be7c-9bffcd37b994
# ╠═61f97b3a-c503-43aa-91f5-ac21cd505625
# ╟─1e7f4850-059d-4cbe-82da-3f8c107f41e6
# ╟─d70e75d0-f811-4bce-afff-d2b5a26caa0a
# ╟─8acce428-b2e9-4bff-a298-53128a5d2c1f
# ╠═ca7be7c0-6798-4af4-a879-6e4edf073093
# ╟─7a46af1e-51e6-4830-aa77-5a93893de833
# ╟─c50a0aaf-9329-46fa-873e-e20aa9cd0562
# ╟─5ba0b9d7-caaa-4c2b-b43e-ae1fcbf1bd98
# ╠═46dc0fae-26ec-4a57-9077-1f2421b39cfc
# ╟─3309d953-e26f-4a7e-a709-c359f2927950
# ╟─bd27bd51-1d78-4980-8de3-cfa3a6edd133
# ╟─04e0025d-5fe4-42cc-92e8-114c071c1722
# ╟─2c56f1d2-8e07-447b-a9d8-953ced9bd12d
# ╟─7a0066ca-62fe-484e-be73-b14c4731ec7a
# ╟─ed9b0a72-ec16-442a-b0b6-bb4206322c32
# ╠═df70a5ce-c529-4a05-94d1-8305c64155e5
# ╠═04e8e0d7-2ae2-4f57-8c5a-8973b93f13e2
# ╠═365f993f-7652-422f-864c-0087c9c3e0f3
# ╠═7baedb85-1b9a-449e-bb73-0c5f1068da3a
# ╠═1f27650d-02b9-42ea-9fb4-e4d2e7f49cb7
