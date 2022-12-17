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

# ╔═╡ 33968390-052e-11ed-3e8c-5174d9b8b5e9
begin
	using Pkg
	using PlutoUI
	using Images
	using DSP
	using FileIO
	using FFTW
	using ImageShow
	using Statistics
	using ImageCore
	using Plots
	plotly()
end

# ╔═╡ ffbb31d9-b981-49e9-9855-5eb9beb8de55
md" # Questão 1

Gabriel Tavares - 10773801
"

# ╔═╡ b929e304-0c84-4e3e-a5f9-4c9effe59f54
begin
	imagem_original = load("cat.png")[1:end-1,:]
end

# ╔═╡ 7d4cc985-0005-4c3f-a00f-5db518371098
md" ## Conversão RGB -> YCbCr


Aqui os canais de cor da imagem são convertidos de RGB para YCbCr seguindo a relação


$Y= \alpha_r R + \alpha_b B + \alpha_g G$
$Cb=  \frac{1}{2(1-\alpha_b)}(B-Y)$
$Cr=  \frac{1}{2(1-\alpha_r)}(R-Y)$
"

# ╔═╡ d7a8f6f3-1d06-4da2-85eb-84f4f31ff30c
md" ## Subamostragem Cb e Cr

Nos novos canais de cor, iremos subamostrar pela metade os canais de crominância. Para isso, iremos pular linhas e colunas da imagem nesses canais"

# ╔═╡ b3a24803-fefc-455a-83ee-bcc5f51d8ffa
md" ## DCT blocos

Agora em cada canal iremos iremos divir a imagem em bloco de 4x4, 8x8 e 16x16. Em cada bloco iremos aplicar a transformada discreta de cossenos. 

Para isso, iremos priemeiro subtrair os blocos por 0.5 para a imagem variar entre **-0.5** até **0.5** e multiplicar os blocos por **255** para podermos quantizar esses valores no item seguinte. A quantização considera valores de pixels que variam de **0** até **255**, e não **0.0** até **1.0** como é nessa implementação de imagens.

Para poder dividir a imagem em blocos sem perda, vamos adicionar linhas e colunas na borda da imagem para termos um número inteiro de blocos.


"

# ╔═╡ de11aa58-a703-4b51-963b-b7e028bf5790
md" #### Blocos DCT - 4x4"

# ╔═╡ d551e0e2-1661-4e66-9fb6-84738c54bf87
md" #### Blocos DCT - 8x8"

# ╔═╡ e74fd7c4-7700-4a5e-a125-9a81bcf327ee
md" #### Blocos DCT - 16x16"

# ╔═╡ e504814e-c731-4584-ab4b-36ab8daad80b
@bind k Slider(1:7, default = 3, show_value = true)

# ╔═╡ 365553ab-4910-4685-9c17-43172c28fdce
md" ## Quantização

Com os canais divididos em blocos de DCTs, iremos quantizar os coeficientes de cada bloco para eliminar as informações menos importantes. Usaremos a matriz de quantização:

$$Q = k\left[\begin{array}{cc} 
8 & 11 & 10 & 16 & 24 & 40 & 51 & 61 \\
12 & 12 & 14 & 19 & 26 & 58 & 60 & 55 \\
14 & 13 & 16 & 24 & 40 & 57 & 69 & 56 \\
14 & 17 & 22 & 29 & 51 & 87 & 80 & 62 \\
18 & 22 & 37 & 56 & 68 & 109 & 103 & 77 \\
24 & 35 & 55 & 64 & 81 & 104 & 113 & 92 \\
49 & 64 & 78 & 87 & 103 & 121 & 120 & 101 \\
72 & 92 & 95 & 98 & 112 & 100 & 103 & 99
\end{array}\right]$$ 

Em cada bloco iremos realizar uma divisão inteira entre o bloco e a matriz para ter a quantização.

Essa é uma matriz 8x8, portanto terá que ser adaptada pros blocos 4x4 e 16x16.

No bloco **4x4** iremos apenas subamostrar a matriz, pulando linhas e colunas.

No bloco **16x16** iremos superamostrar a matriz colocando zeros entre as linhas e colunas e depois passar por um filtro de interpolação linear. A última linha e coluna serão apenas replicadas das anteriores por causa do efeito de borda da interpolação.

Além disso a matriz **Q_n** é multiplicada por um fator **k = $k**, para termos diferentes níveis de qualidade na quantização.
"

# ╔═╡ 735181d5-1d0e-4f16-94a5-aea5892fab20
md" #### Quantização - 4x4"

# ╔═╡ 5bb28256-ea66-4e3d-8637-1d5046fa6f2f
md" #### Quantização - 8x8 "

# ╔═╡ c2587a9e-ee35-4fea-a887-f3cab92a8669
md" #### Quantização - 16x16"

# ╔═╡ de58696c-0a39-4b5e-a59b-8d6e8390f479
md" ## Reconstrução


Para reconstruir a imagems, iremos fazer o inverso dos itens anteriores.

* Em cada bloco iremos multiplicar o bloco pela matriz de quantização para recuperar os valores quantizados.
* Depois faremos a Transformada Discreta de Cossenos Inversa para ter a imagem no domínio do espaço
* Dividimos essa imagem por **255** para recuperar a forma que julia interpreta as imagens (0.0 até 1.0)
* Por fim somamos **0.5** a todos os pixels da imagem

Com isso temos os 3 canais da imagem recuprados no domínio espacial e na escala correta. 

Por fim, é necessário superamostrar os canais de cromância que estão subamostrados. Para isso, preenchemos as linhas e colunas com zeros e fazemos uma interpolação linear usando um filtro interpolador.
"

# ╔═╡ 128b32a0-d957-4f22-89f6-3d7693a390cc
md" #### Reconstrução - 4x4"

# ╔═╡ ae3dc871-899b-4192-b7d2-3d77ea98498e
md" #### Reconstrução - 8x8 "

# ╔═╡ 1ee91a69-c81d-4e6d-bbbf-bc23f836bdc1
md" #### Reconstrução - 16x16"

# ╔═╡ 8c56ee11-1cdc-461d-a3f0-f3e91a4f31d5
md" ## Análise 

Aqui iremos analisar a qualidade da imagem de forma qualitativa e de forma analítica"

# ╔═╡ 6de48ba3-da92-4fa9-89e5-a59726817fe6
md" #### Qualitativamente" 

# ╔═╡ 72cecc71-379c-4f9f-b185-f747cabec38f
md" Qualitativamente, usando **k = $k** temos

* **16x16** : a coloração está correta, mas vemos grandes blocos se formando, o que deixa a imagem bastante pixelada. Isso ocorre por causa da divisão das DCTs em grandes blocos
* **8x8** : a coloração está correta, e vemos a formações de blocos, mas esses blocos são menores e mais discretos.
* **4x4** : a coloração começa a ter erro, e puxar um pouco para o vermelho. Além disso, é possível ver os blocos se formando, mas por serem pequenos, a sensação de pixelamento da imagem é menor"

# ╔═╡ d294309e-0dfa-48de-80c4-0409f0c417b9
md" #### PSNR

Aqui vamos cacular o **PSNR** de cada imagem usandos as diferentes janelas. Esse calculo ve a PSNR de cada canal e faz a média entre os 3 canais para ter o valor final.
"

# ╔═╡ 959e2df7-7c27-40c9-8f8b-1ea9057caeec
md" #### Taxa de números diferentes de zero

Aqui calculamos a quantidade de valores diferentes de zero para analisar a quantidade de coeficientes que carregam informações que ainda temos.


"

# ╔═╡ 6461b6fc-8136-4a83-80f8-d2fcce279dac
md" ## Conclusão

Com essa simulação, vemos a importância da escolha de uma janela ideal no padrão JPEG. 

A janela de **4x4** é uma janela que tem menor compressão da imagem (mais valores diferentes de zero) e traz algumas distorçoes na imagem final (canais de cor principalmente), mas não tem muito efeito do pixelamento visual. 

A janela **16x16** apresenta maior compressão e maior PSNR, mas ao vermos a imagem, fica bastante visível a formação de grandes blocos que deixa a imagem visualmente pior.

Já a janela **8x8** acaba equilibrando essas relações. Ela tem uma PSNR alta, uma taxa de compressão de dados alta e o efeito de pixelamente é menos notável na imagem final.


"

# ╔═╡ 8f1e1677-dedc-4888-8236-862412a5d8df
md" # Functions"

# ╔═╡ 0d332e67-5a18-40b2-991f-36dbcbdb4c65
noprint = md""

# ╔═╡ f4a254ef-5f61-404f-9b7d-12630219874d
begin
	Q = 
	[8 11 10 16 24 40 51 61;
	12 12 14 19 26 58 60 55;
	14 13 16 24 40 57 69 56;
	14 17 22 29 51 87 80 62;
	18 22 37 56 68 109 103 77;
	24 35 55 64 81 104 113 92;
	49 64 78 87 103 121 120 101;
	72 92 95 98 112 100 103 99]
	
	Q = Q*k

	noprint
end

# ╔═╡ 4c20f6c0-eb5b-4ebd-a7bd-c3afd567418e
function n_element(matriz)
	return size(matriz)[1]*size(matriz)[2]
end

# ╔═╡ 44a5ebcf-0e49-4309-9d17-47a2143eec20
md" ### Gray matrix images functions"

# ╔═╡ 0bb36cc1-09eb-4199-9c58-801ac78d8ed4
function mat_to_image(mat)
	return RGB.(mat,mat,mat)	
end

# ╔═╡ e15b85cf-1c09-46ab-9374-1eef22531076
function image(mat)
	return RGB.(mat,mat,mat)
end

# ╔═╡ fffe60c9-e44e-4930-82b1-373c42631002
md" ### RGB and YCbCr functions

"

# ╔═╡ d1a0d53f-59dd-4d9c-9a25-9f5bf615fbb3
struct YCbCr_pixel
	Y
	Cb
	Cr
end

# ╔═╡ 2cc4635e-1ef3-4569-9a5e-1d1adad47798
function luminancia(pixel::YCbCr_pixel)
	return pixel.Y
end

# ╔═╡ 5613743d-d9b3-4605-8c30-bd5d5bed4031
function croma_blue(pixel::YCbCr_pixel)
	return pixel.Cb
end

# ╔═╡ 3bdf35aa-c9a1-4a66-b9b6-214845b6f400
function croma_red(pixel::YCbCr_pixel)
	return pixel.Cr
end

# ╔═╡ 4df910c9-021a-4fb1-947c-783ac606a247
function get_Y(pixel::RGB)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	Y = αr*red(pixel) +αg*green(pixel) + αb*blue(pixel)
	
	return Y
end

# ╔═╡ 744cc57d-f1d7-492a-9c63-0e30b58f89b1
function get_Cb(pixel::RGB)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	Y = get_Y(pixel)
	Cb = 1/(2*(1-αb))*(blue(pixel)-Y)
end

# ╔═╡ c175bc56-0eb7-491f-a1ee-ba4a6d2d39e8
function get_Cr(pixel::RGB)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	Y = get_Y(pixel)
	Cr = 1/(2*(1-αr))*(red(pixel)-Y)
end

# ╔═╡ d6a9d811-330d-4fd6-abcc-6e019a39e0ad
begin
	Y_original = get_Y.(imagem_original)
	Cb_original = get_Cb.(imagem_original)
	Cr_original = get_Cr.(imagem_original)
	noprint
end

# ╔═╡ 03b031db-c076-4cec-bacf-8d1fa7d46a55
begin
	Cr_sub = Cr_original[1:2:end, 1:2:end]
	Cb_sub = Cb_original[1:2:end, 1:2:end]
	noprint
end

# ╔═╡ ecc162e1-c58b-4239-a7e6-57c4bf0967ac
begin
	#Tamanhos finais da imagem
	linhas_original = size(imagem_original)[1]
	colunas_original = size(imagem_original)[2]
	
	#Tamanhos finais da imagem
	linhas_original_Y = size(Y_original)[1]
	colunas_original_Y = size(Y_original)[2]
	
	linhas_original_C = size(Cb_sub)[1]
	colunas_original_C = size(Cb_sub)[2]
	noprint
end

# ╔═╡ 89b990d9-9d6c-4e3d-94b8-1b49a12f716a
function get_YCbCr(pixel::RGB)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	Y = get_Y(pixel)
	Cb = get_Cb(pixel)
	Cr = get_Cr(pixel)
	return Y, Cb, Cr
end

# ╔═╡ c237e2bf-c35a-468e-aab1-bdc1bf04559c
function get_YCbCr(R,G,B)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	pixel = RGB(R,G,B)
	Y = get_Y(pixel)
	Cb = get_Cb(pixel)
	Cr = get_Cr(pixel)
	return Y, Cb, Cr
end

# ╔═╡ d92ecfb8-5357-4f3a-a022-b8d7eb2464d8
function get_R(Y, Cb, Cr)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	return Y + (2-2*αr)*Cr
end

# ╔═╡ 9b4f4fd7-fdbe-482b-a066-1c25cac2977f
function get_G(Y, Cb, Cr)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	return Y - αb/αg*(2-2*αb)*Cb -αr/αg*(2-2*αr)*Cr
end

# ╔═╡ 9488c2f2-296b-4944-8e2c-93a123080657
function get_B(Y, Cb, Cr)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	R = get_R(Y,Cb,Cr)
	G = get_G(Y, Cb, Cr)
	
	return Y + (2-2*αb)*Cb
end

# ╔═╡ 3a3cf5a1-c3e8-4688-bd4a-2f726442d1f6
function get_RGB(Y, Cb, Cr)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	R = get_R(Y,Cb,Cr)
	G = get_G(Y, Cb, Cr)
	B = get_B(Y, Cb, Cr)
	
	return RGB(R,G,B)
end

# ╔═╡ 60410518-0037-4ac3-88ce-87b6ba10718c
md" ## PSNR"

# ╔═╡ 625577ec-9f9e-467d-abae-bfb44a78c475
function MSE(imagem_exata, imagem)
	MSE = 0
	linhas, colunas = size(imagem)
	
	for lin in 1:linhas
		for col in 1:colunas
			MSE += (imagem_exata[lin, col] - imagem[lin, col])^2
		end
	end
	
	MSE = MSE/(linhas*colunas)
		
	return MSE
end

# ╔═╡ d14449d6-128e-4259-87b4-db19fa0f71f1
function PSNR(imagem_exata, imagem)
	PSNR = 0
	maxi = 1.0
	MSE_ = MSE(imagem_exata, imagem)
	PSNR = 10*log10(maxi/ MSE_)
	return PSNR
end

# ╔═╡ 759323db-7d17-482a-a8f6-ef32524f8126
function PSNR(imagem_exata::Matrix{RGB{N0f8}}, imagem::Matrix{RGB{Float64}})
	PSNR_r = PSNR(red.(imagem_exata), red.(imagem))
	PSNR_g = PSNR(green.(imagem_exata), green.(imagem))
	PSNR_b = PSNR(blue.(imagem_exata), blue.(imagem))	
	PSNR_ = mean([PSNR_r, PSNR_g, PSNR_b])
	
	
	return PSNR_

end

# ╔═╡ 999a96e4-58bb-4b15-96cb-1d189ea9c1dd
md" ## Interpolador"

# ╔═╡ 898459dd-3a04-46cc-b647-3dd7ef2664c7
int_filter = [0.25  0.5 	0.25;
			  0.5  	1.0 	0.5 ;
			  0.25 	0.5 	0.25]

# ╔═╡ 0a694982-3836-4034-9659-961043e57e4b
md" ## Padarray"

# ╔═╡ aa6cefe0-d6a1-4745-b77a-9ad7327c0c22
#Adiciona linhas colunas no canto direito e linhas na parte inferior da matrix
# Só faz o espelhamento da imagem

#exemplo
#A = [ 1 2 3
#	   4 5 6
#	   7 8 9]
# 
#padarray(A, 2, 2)
#[ 1 2 3 3 2
#  4 5 6 6 5
#  7 8 9 9 8
#  7 8 9 9 6
#  4 5 6 8 5

function my_padarray(imagem, linhas_add, colunas_add)
	linhas_original = size(imagem)[1]
	colunas_original = size(imagem)[2]
	
	imagem_out = zeros(typeof(imagem[1]), linhas_original + linhas_add, colunas_original+colunas_add)
	
	imagem_out[1:linhas_original, 1:colunas_original] = imagem
	
	
	linhas_preenche = imagem[ end:-1:(end-linhas_add+1) , :]
	colunas_preenche = imagem[ :, end:-1:(end-colunas_add+1)]
	quadrado_preenche = imagem[ end:-1:(end-linhas_add+1), end:-1:(end-colunas_add+1)]
	
	
	imagem_out[linhas_original + 1: end, 1:colunas_original] = linhas_preenche
	imagem_out[1:linhas_original, colunas_original + 1 : end] = colunas_preenche
	imagem_out[linhas_original+1:end, colunas_original+1:end] = quadrado_preenche
	
	
	return imagem_out
end
	

# ╔═╡ 417ef625-6735-4bf9-85cd-4db85639e4c0
begin
	#ADIÇÃO DE LINHAS PRA TER BLOCOS COMPLETOS DE 8X8
	linhas_add_Y_4 = linhas_original_Y % 4 == 0 ? 0 : 4 - linhas_original_Y % 4
	colunas_add_Y_4 = colunas_original_Y % 4 == 0 ? 0 : 4 - colunas_original_Y % 4
	
	linhas_expand_Y_4 = linhas_original_Y + linhas_add_Y_4
	colunas_expand_Y_4 = colunas_original_Y + colunas_add_Y_4
	

	linhas_add_C_4 = linhas_original_C % 4
	colunas_add_C_4 = colunas_original_C % 4
	
	linhas_expand_C_4 = linhas_original_C + linhas_add_C_4
	colunas_expand_C_4 = colunas_original_C + colunas_add_C_4

	
	Y_expand_4= my_padarray(Y_original, linhas_add_Y_4, colunas_add_Y_4) .- 0.5
	Cb_expand_4 = my_padarray(Cb_sub, linhas_add_C_4, colunas_add_C_4) .- 0.5
	Cr_expand_4 = my_padarray(Cr_sub, linhas_add_C_4, colunas_add_C_4) .- 0.5
	
	noprint
end

# ╔═╡ 9ae5b74d-12ea-4dcb-9ae9-e36cdc8effce
begin
	#DCT por bloco
	
	#Luminancia=======================
	numx_blocos_Y_4 =size(Y_expand_4)[1]÷4
	numy_blocos_Y_4 =size(Y_expand_4)[2]÷4
	Y_dct_4 = zeros(size(Y_expand_4))
	
	for i in 1:numx_blocos_Y_4
		for j in 1:numy_blocos_Y_4
			bloco = Y_expand_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4] 
			Y_dct_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4]  = dct(bloco*255)
		end
	end
	
	#Cromância=======================
	numx_blocos_C_4 =size(Cb_expand_4)[1]÷4
	numy_blocos_C_4 =size(Cb_expand_4)[2]÷4
	Cb_dct_4 = zeros(size(Cb_expand_4))
	
	for i in 1:numx_blocos_C_4
		for j in 1:numy_blocos_C_4
			bloco = Cb_expand_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4] 
			Cb_dct_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4]  = dct(bloco*255)
		end
	end
	
	Cr_dct_4 = zeros(size(Cr_expand_4))
	
	for i in 1:numx_blocos_C_4
		for j in 1:numy_blocos_C_4
			bloco = Cr_expand_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4] 
			Cr_dct_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4]  = dct(bloco*255)
		end
	end
end

# ╔═╡ 74ee7277-514c-4a40-8c17-380f64849bd9
begin
	Q_4 = Q[1:2:end, 1:2:end]
	#Luminancia======================	
	Y_dct_q_4 = zeros(size(Y_expand_4))
	
	for i in 1:numx_blocos_Y_4
		for j in 1:numy_blocos_Y_4
			bloco = Y_dct_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4] 
			Y_dct_q_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4]  = (bloco).÷Q_4
		end
	end
	
	Cb_dct_q_4 = zeros(size(Cb_expand_4))
	
	for i in 1:numx_blocos_C_4
		for j in 1:numy_blocos_C_4
			bloco = Cb_dct_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4] 
			Cb_dct_q_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4]  = (bloco).÷Q_4
		end
	end
	
	Cr_dct_q_4 = zeros(size(Cr_expand_4))
	
	for i in 1:numx_blocos_C_4
		for j in 1:numy_blocos_C_4
			bloco = Cr_dct_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4] 
			Cr_dct_q_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4]  = (bloco).÷Q_4
		end
	end
end

# ╔═╡ 639fa9c2-a6e9-4b22-bb48-b0ea93ad1010
begin
	#Luminancia======================
	Y_rec_4 = zeros(size(Y_expand_4))
	
	for i in 1:numx_blocos_Y_4
		for j in 1:numy_blocos_Y_4
			bloco = Y_dct_q_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4] 
			bloco = bloco.*Q_4
			bloco = idct(bloco)./255
			Y_rec_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4]  = bloco .+ 0.5
		end
	end
	Y_rec_4 = Y_rec_4[1:linhas_original, 1:colunas_original]
	
	Cb_rec_sub_4 = zeros(size(Cb_expand_4))
	
	for i in 1:numx_blocos_C_4
		for j in 1:numy_blocos_C_4
			bloco = Cb_dct_q_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4] 
			bloco = bloco.*Q_4
			bloco = idct(bloco)./255
			Cb_rec_sub_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4]  = bloco .+ 0.5
		end
	end
	
	Cb_rec_4 = zeros(size(Cb_rec_sub_4).*2)
	Cb_rec_4[1:2:end, 1:2:end] = Cb_rec_sub_4
	Cb_rec_4 = imfilter(Cb_rec_4, int_filter)
	Cb_rec_4 = Cb_rec_4[1:linhas_original, 1:colunas_original]
	
	
	
	Cr_rec_sub_4 = zeros(size(Cr_expand_4))
	
	for i in 1:numx_blocos_C_4
		for j in 1:numy_blocos_C_4
			bloco = Cr_dct_q_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4] 
			bloco = bloco.*Q_4
			bloco = idct(bloco)./255
			Cr_rec_sub_4[(i-1)*4 + 1:i*4, (j-1)*4+1:j*4]  = bloco .+ 0.5
		end
	end
	
	Cr_rec_4 = zeros(size(Cr_rec_sub_4).*2)
	Cr_rec_4[1:2:end, 1:2:end] = Cr_rec_sub_4
	Cr_rec_4 = imfilter(Cr_rec_4, int_filter)
	Cr_rec_4 = Cr_rec_4[1:linhas_original, 1:colunas_original]
	noprint
	
end

# ╔═╡ 2744e922-0c65-4bfd-bbfa-7810ec85180e
begin
	#ADIÇÃO DE LINHAS PRA TER BLOCOS COMPLETOS DE 8X8
	linhas_add_Y_8 = linhas_original_Y % 8 == 0 ? 0 : 8 - linhas_original_Y % 8
	colunas_add_Y_8 = colunas_original_Y % 8 == 0 ? 0 : 8 - colunas_original_Y % 8
	
	linhas_expand_Y_8 = linhas_original_Y + linhas_add_Y_8
	colunas_expand_Y_8 = colunas_original_Y + colunas_add_Y_8
	

	linhas_add_C_8 = linhas_original_C % 8
	colunas_add_C_8 = colunas_original_C % 8
	
	linhas_expand_C_8 = linhas_original_C + linhas_add_C_8
	colunas_expand_C_8 = colunas_original_C + colunas_add_C_8

	
	Y_expand_8= my_padarray(Y_original, linhas_add_Y_8, colunas_add_Y_8) .- 0.5
	Cb_expand_8 = my_padarray(Cb_sub, linhas_add_C_8, colunas_add_C_8) .- 0.5
	Cr_expand_8 = my_padarray(Cr_sub, linhas_add_C_8, colunas_add_C_8).- 0.5
	
	noprint
end

# ╔═╡ 34faf8c1-43d8-4721-86da-ef5d09d00461
begin
	#DCT por bloco
	#Luminancia=======================
	numx_blocos_Y_8 =size(Y_expand_8)[1]÷8
	numy_blocos_Y_8 =size(Y_expand_8)[2]÷8
	Y_dct_8 = zeros(size(Y_expand_8))
	
	for i in 1:numx_blocos_Y_8
		for j in 1:numy_blocos_Y_8
			bloco = Y_expand_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8] 
			Y_dct_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8]  = dct(bloco*255)
		end
	end
	
	numx_blocos_C_8 =size(Cb_expand_8)[1]÷8
	numy_blocos_C_8 =size(Cb_expand_8)[2]÷8
	Cb_dct_8 = zeros(size(Cb_expand_8))
	
	for i in 1:numx_blocos_C_8
		for j in 1:numy_blocos_C_8
			bloco = Cb_expand_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8] 
			Cb_dct_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8]  = dct(bloco*255)
		end
	end
	
	Cr_dct_8 = zeros(size(Cr_expand_8))
	
	for i in 1:numx_blocos_C_8
		for j in 1:numy_blocos_C_8
			bloco = Cr_expand_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8] 
			Cr_dct_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8]  = dct(bloco*255)
		end
	end
end

# ╔═╡ 71718be8-c46d-4619-9c84-8315f6e42019
begin
	Q_8 = Q
	
	#Luminancia======================	
	Y_dct_q_8 = zeros(size(Y_expand_8))
	
	for i in 1:numx_blocos_Y_8
		for j in 1:numy_blocos_Y_8
			bloco = Y_dct_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8] 
			Y_dct_q_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8]  = (bloco).÷Q_8
		end
	end
	
	Cb_dct_q_8 = zeros(size(Cb_expand_8))
	
	for i in 1:numx_blocos_C_8
		for j in 1:numy_blocos_C_8
			bloco = Cb_dct_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8] 
			Cb_dct_q_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8]  = (bloco).÷Q_8
		end
	end
	
	Cr_dct_q_8 = zeros(size(Cr_expand_8))
	
	for i in 1:numx_blocos_C_8
		for j in 1:numy_blocos_C_8
			bloco = Cr_dct_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8] 
			Cr_dct_q_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8]  = (bloco).÷Q_8
		end
	end
end

# ╔═╡ 4ad88d56-6f9c-4572-b6e5-8e8d77f39bb5
begin
	#Luminancia======================
	Y_rec_8 = zeros(size(Y_expand_8))
	
	for i in 1:numx_blocos_Y_8
		for j in 1:numy_blocos_Y_8
			bloco = Y_dct_q_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8] 
			bloco = bloco.*Q_8
			bloco = idct(bloco)./255
			Y_rec_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8]  = bloco .+ 0.5
		end
	end
	Y_rec_8 = Y_rec_8[1:linhas_original, 1:colunas_original]
	
	Cb_rec_sub_8 = zeros(size(Cb_expand_8))
	
	for i in 1:numx_blocos_C_8
		for j in 1:numy_blocos_C_8
			bloco = Cb_dct_q_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8] 
			bloco = bloco.*Q_8
			bloco = idct(bloco)./255
			Cb_rec_sub_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8]  = bloco .+ 0.5
		end
	end
	
	Cb_rec_8 = zeros(size(Cb_rec_sub_8).*2)
	Cb_rec_8[1:2:end, 1:2:end] = Cb_rec_sub_8
	Cb_rec_8 = imfilter(Cb_rec_8, int_filter)
	Cb_rec_8 = Cb_rec_8[1:linhas_original, 1:colunas_original]
	
	
	#---
	
	Cr_rec_sub_8 = zeros(size(Cr_expand_8))
	
	for i in 1:numx_blocos_C_8
		for j in 1:numy_blocos_C_8
			bloco = Cr_dct_q_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8] 
			bloco = bloco.*Q_8
			bloco = idct(bloco)./255
			Cr_rec_sub_8[(i-1)*8 + 1:i*8, (j-1)*8+1:j*8]  = bloco .+ 0.5
		end
	end
	
	Cr_rec_8 = zeros(size(Cr_rec_sub_8).*2)
	Cr_rec_8[1:2:end, 1:2:end] = Cr_rec_sub_8
	Cr_rec_8 = imfilter(Cr_rec_8, int_filter)
	Cr_rec_8 = Cr_rec_8[1:linhas_original, 1:colunas_original]
	noprint
	
end

# ╔═╡ 20750b35-0f91-402a-80b1-1c31a8a5eb0d
begin
	#ADIÇÃO DE LINHAS PRA TER BLOCOS COMPLETOS DE 8X8
	linhas_add_Y_16 = linhas_original_Y % 16 == 0 ? 0 : 16 - linhas_original_Y % 16
	colunas_add_Y_16 = colunas_original_Y % 16 == 0 ? 0 : 16 - colunas_original_Y % 16
	
	linhas_expand_Y_16 = linhas_original_Y + linhas_add_Y_16
	colunas_expand_Y_16 = colunas_original_Y + colunas_add_Y_16
	

	linhas_add_C_16 = linhas_original_C % 16
	colunas_add_C_16 = colunas_original_C % 16
	
	linhas_expand_C_16 = linhas_original_C + linhas_add_C_16
	colunas_expand_C_16 = colunas_original_C + colunas_add_C_16

	
	Y_expand_16= my_padarray(Y_original, linhas_add_Y_16, colunas_add_Y_16) .- 0.5
	Cb_expand_16 = my_padarray(Cb_sub, linhas_add_C_16, colunas_add_C_16) .- 0.5
	Cr_expand_16 = my_padarray(Cr_sub, linhas_add_C_16, colunas_add_C_16).- 0.5
	
	noprint
end

# ╔═╡ c577ca22-dfe4-4f5e-ae05-bda144de173e
begin
	#DCT por bloco
	#Luminancia=======================
	numx_blocos_Y_16 =size(Y_expand_16)[1]÷16
	numy_blocos_Y_16 =size(Y_expand_16)[2]÷16
	Y_dct_16 = zeros(size(Y_expand_16))
	
	for i in 1:numx_blocos_Y_16
		for j in 1:numy_blocos_Y_16
			bloco = Y_expand_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16] 
			Y_dct_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16]  = dct(bloco*255)
		end
	end
	
	numx_blocos_C_16 =size(Cb_expand_16)[1]÷16
	numy_blocos_C_16 =size(Cb_expand_16)[2]÷16
	Cb_dct_16 = zeros(size(Cb_expand_16))
	
	for i in 1:numx_blocos_C_16
		for j in 1:numy_blocos_C_16
			bloco = Cb_expand_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16] 
			Cb_dct_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16]  = dct(bloco*255)
		end
	end
	
	Cr_dct_16 = zeros(size(Cr_expand_16))
	
	for i in 1:numx_blocos_C_16
		for j in 1:numy_blocos_C_16
			bloco = Cr_expand_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16] 
			Cr_dct_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16]  = dct(bloco*255)
		end
	end
end

# ╔═╡ 98740d4f-db45-4b16-a7ae-a0daa9f4e1c6
begin
	Q_16 = zeros(16,16)
	Q_16[1:2:end, 1:2:end] = Q
	Q_16 = imfilter(Q_16, int_filter, Fill(0)) #interpolação das frequencias altas
	Q_16[end-1, end] = Q_16[end-1, end - 1] 
	Q_16[end, end-1] = Q_16[end-1, end-1] 
	Q_16[:, end] = Q_16[:, end-1] 
	Q_16[end, :] = Q_16[end-1, :] 
	
	
	#Luminancia======================	
	Y_dct_q_16 = zeros(size(Y_expand_16))
	
	for i in 1:numx_blocos_Y_16
		for j in 1:numy_blocos_Y_16
			bloco = Y_dct_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16] 
			Y_dct_q_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16]  = (bloco).÷Q_16
		end
	end
	
	Cb_dct_q_16 = zeros(size(Cb_expand_16))
	
	for i in 1:numx_blocos_C_16
		for j in 1:numy_blocos_C_16
			bloco = Cb_dct_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16] 
			Cb_dct_q_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16]  = (bloco).÷Q_16
		end
	end
	
	Cr_dct_q_16 = zeros(size(Cr_expand_16))
	
	for i in 1:numx_blocos_C_16
		for j in 1:numy_blocos_C_16
			bloco = Cr_dct_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16] 
			Cr_dct_q_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16]  = (bloco).÷Q_16
		end
	end
end

# ╔═╡ 0cc3b73b-db04-4906-b55b-50d5b1200138
begin
	num_zeros_4 = sum(Y_dct_q_4 .!= 0) + sum(Cb_dct_q_4 .!= 0) + sum(Cr_dct_q_4 .!= 0)
	not_zeros_ratio_4 = num_zeros_4/(n_element(Y_dct_q_4)+n_element(Cb_dct_q_4)+n_element(Cr_dct_q_4))
	
	num_zeros_8 = sum(Y_dct_q_8 .!= 0) + sum(Cb_dct_q_8 .!= 0) + sum(Cr_dct_q_8 .!= 0)
	not_zeros_ratio_8 = num_zeros_8/(n_element(Y_dct_q_8)+n_element(Cb_dct_q_8)+n_element(Cr_dct_q_8))
	
	num_zeros_16 = sum(Y_dct_q_8 .!= 0) + sum(Cb_dct_q_8 .!= 0) + sum(Cr_dct_q_8 .!= 0)
	not_zeros_ratio_16 = num_zeros_16/(n_element(Y_dct_q_16)+n_element(Cb_dct_q_16)+n_element(Cr_dct_q_16))
	noprint
end

# ╔═╡ 1fcb68e3-dbb4-4251-9344-5a6ac95e7bd7
md" 

Com **k =  $k** os diferentes padrões de JPEG tem as seguintes taxas de números diferentes de zero

**Taxa de numeros:**
* **4x4** : $(round(not_zeros_ratio_4*100, digits = 2))%
* **8x8** : $(round(not_zeros_ratio_8*100, digits = 2))%
* **16x16** : $(round(not_zeros_ratio_16*100, digits = 2))%
"

# ╔═╡ 5a4b6bf3-67b7-4850-8f8b-36f3b4c74909
begin
	#Luminancia======================
	Y_rec_16 = zeros(size(Y_expand_16))
	
	for i in 1:numx_blocos_Y_16
		for j in 1:numy_blocos_Y_16
			bloco = Y_dct_q_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16] 
			bloco = bloco.*Q_16
			bloco = idct(bloco)./255
			Y_rec_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16]  = bloco .+ 0.5
		end
	end
	Y_rec_16 = Y_rec_16[1:linhas_original, 1:colunas_original]
	
	Cb_rec_sub_16 = zeros(size(Cb_expand_16))
	
	for i in 1:numx_blocos_C_16
		for j in 1:numy_blocos_C_16
			bloco = Cb_dct_q_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16] 
			bloco = bloco.*Q_16
			bloco = idct(bloco)./255
			Cb_rec_sub_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16]  = bloco .+ 0.5
		end
	end
	
	Cb_rec_16 = zeros(size(Cb_rec_sub_16).*2)
	Cb_rec_16[1:2:end, 1:2:end] = Cb_rec_sub_16
	Cb_rec_16 = imfilter(Cb_rec_16, int_filter)
	Cb_rec_16 = Cb_rec_16[1:linhas_original, 1:colunas_original]
	
	
	#---
	
	Cr_rec_sub_16 = zeros(size(Cr_expand_16))
	
	for i in 1:numx_blocos_C_16
		for j in 1:numy_blocos_C_16
			bloco = Cr_dct_q_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16] 
			bloco = bloco.*Q_16
			bloco = idct(bloco)./255
			Cr_rec_sub_16[(i-1)*16 + 1:i*16, (j-1)*16+1:j*16]  = bloco .+ 0.5
		end
	end
	
	Cr_rec_16 = zeros(size(Cr_rec_sub_16).*2)
	Cr_rec_16[1:2:end, 1:2:end] = Cr_rec_sub_16
	Cr_rec_16 = imfilter(Cr_rec_16, int_filter)
	Cr_rec_16 = Cr_rec_16[1:linhas_original, 1:colunas_original]
	noprint
	
end

# ╔═╡ 3f0685b7-556e-4de2-873f-4483d3b0732a
begin
	#visualização das imagens
	imagem_recuparada_4 = get_RGB.(Y_rec_4, Cb_rec_4, Cr_rec_4)
	imagem_recuparada_8 = get_RGB.(Y_rec_8, Cb_rec_8, Cr_rec_8)
	imagem_recuparada_16 = get_RGB.(Y_rec_16, Cb_rec_16, Cr_rec_16)
	
	mosaico = hcat(imagem_recuparada_4, imagem_recuparada_8, imagem_recuparada_16)
	noprint
end

# ╔═╡ b076a993-e31a-4e78-a897-7673ec644a36
mosaico

# ╔═╡ 743d04ba-eafd-473e-9fc1-a14d17ae4b50
begin
	PSNR_4 = PSNR(imagem_original, imagem_recuparada_4)
	PSNR_8 = PSNR(imagem_original, imagem_recuparada_8)
	PSNR_16 = PSNR(imagem_original, imagem_recuparada_16)
	noprint
end

# ╔═╡ f64b222a-359c-401b-aa48-95e7c0bcfff3
md" 

Com **k =  $k** os diferentes padrões de JPEG tem as seguintes PSNRs

**PSNR:**
* **4x4** : $(round(PSNR_4, digits = 2))dB
* **8x8** : $(round(PSNR_8, digits = 2))dB
* **16x16** : $(round(PSNR_16, digits = 2))dB
"

# ╔═╡ Cell order:
# ╟─ffbb31d9-b981-49e9-9855-5eb9beb8de55
# ╠═33968390-052e-11ed-3e8c-5174d9b8b5e9
# ╠═b929e304-0c84-4e3e-a5f9-4c9effe59f54
# ╟─7d4cc985-0005-4c3f-a00f-5db518371098
# ╠═d6a9d811-330d-4fd6-abcc-6e019a39e0ad
# ╟─d7a8f6f3-1d06-4da2-85eb-84f4f31ff30c
# ╠═03b031db-c076-4cec-bacf-8d1fa7d46a55
# ╟─b3a24803-fefc-455a-83ee-bcc5f51d8ffa
# ╠═ecc162e1-c58b-4239-a7e6-57c4bf0967ac
# ╟─de11aa58-a703-4b51-963b-b7e028bf5790
# ╠═417ef625-6735-4bf9-85cd-4db85639e4c0
# ╠═9ae5b74d-12ea-4dcb-9ae9-e36cdc8effce
# ╟─d551e0e2-1661-4e66-9fb6-84738c54bf87
# ╠═2744e922-0c65-4bfd-bbfa-7810ec85180e
# ╠═34faf8c1-43d8-4721-86da-ef5d09d00461
# ╟─e74fd7c4-7700-4a5e-a125-9a81bcf327ee
# ╠═20750b35-0f91-402a-80b1-1c31a8a5eb0d
# ╠═c577ca22-dfe4-4f5e-ae05-bda144de173e
# ╟─365553ab-4910-4685-9c17-43172c28fdce
# ╠═e504814e-c731-4584-ab4b-36ab8daad80b
# ╠═f4a254ef-5f61-404f-9b7d-12630219874d
# ╟─735181d5-1d0e-4f16-94a5-aea5892fab20
# ╠═74ee7277-514c-4a40-8c17-380f64849bd9
# ╟─5bb28256-ea66-4e3d-8637-1d5046fa6f2f
# ╠═71718be8-c46d-4619-9c84-8315f6e42019
# ╟─c2587a9e-ee35-4fea-a887-f3cab92a8669
# ╠═98740d4f-db45-4b16-a7ae-a0daa9f4e1c6
# ╟─de58696c-0a39-4b5e-a59b-8d6e8390f479
# ╟─128b32a0-d957-4f22-89f6-3d7693a390cc
# ╠═639fa9c2-a6e9-4b22-bb48-b0ea93ad1010
# ╟─ae3dc871-899b-4192-b7d2-3d77ea98498e
# ╠═4ad88d56-6f9c-4572-b6e5-8e8d77f39bb5
# ╟─1ee91a69-c81d-4e6d-bbbf-bc23f836bdc1
# ╠═5a4b6bf3-67b7-4850-8f8b-36f3b4c74909
# ╟─8c56ee11-1cdc-461d-a3f0-f3e91a4f31d5
# ╠═3f0685b7-556e-4de2-873f-4483d3b0732a
# ╟─6de48ba3-da92-4fa9-89e5-a59726817fe6
# ╟─b076a993-e31a-4e78-a897-7673ec644a36
# ╟─72cecc71-379c-4f9f-b185-f747cabec38f
# ╟─d294309e-0dfa-48de-80c4-0409f0c417b9
# ╠═743d04ba-eafd-473e-9fc1-a14d17ae4b50
# ╟─f64b222a-359c-401b-aa48-95e7c0bcfff3
# ╟─959e2df7-7c27-40c9-8f8b-1ea9057caeec
# ╠═0cc3b73b-db04-4906-b55b-50d5b1200138
# ╟─1fcb68e3-dbb4-4251-9344-5a6ac95e7bd7
# ╟─6461b6fc-8136-4a83-80f8-d2fcce279dac
# ╟─8f1e1677-dedc-4888-8236-862412a5d8df
# ╠═0d332e67-5a18-40b2-991f-36dbcbdb4c65
# ╠═4c20f6c0-eb5b-4ebd-a7bd-c3afd567418e
# ╠═44a5ebcf-0e49-4309-9d17-47a2143eec20
# ╠═0bb36cc1-09eb-4199-9c58-801ac78d8ed4
# ╠═e15b85cf-1c09-46ab-9374-1eef22531076
# ╠═fffe60c9-e44e-4930-82b1-373c42631002
# ╠═d1a0d53f-59dd-4d9c-9a25-9f5bf615fbb3
# ╠═2cc4635e-1ef3-4569-9a5e-1d1adad47798
# ╠═5613743d-d9b3-4605-8c30-bd5d5bed4031
# ╠═3bdf35aa-c9a1-4a66-b9b6-214845b6f400
# ╠═4df910c9-021a-4fb1-947c-783ac606a247
# ╠═744cc57d-f1d7-492a-9c63-0e30b58f89b1
# ╠═c175bc56-0eb7-491f-a1ee-ba4a6d2d39e8
# ╠═89b990d9-9d6c-4e3d-94b8-1b49a12f716a
# ╠═c237e2bf-c35a-468e-aab1-bdc1bf04559c
# ╠═d92ecfb8-5357-4f3a-a022-b8d7eb2464d8
# ╠═9b4f4fd7-fdbe-482b-a066-1c25cac2977f
# ╠═9488c2f2-296b-4944-8e2c-93a123080657
# ╠═3a3cf5a1-c3e8-4688-bd4a-2f726442d1f6
# ╠═60410518-0037-4ac3-88ce-87b6ba10718c
# ╠═625577ec-9f9e-467d-abae-bfb44a78c475
# ╠═d14449d6-128e-4259-87b4-db19fa0f71f1
# ╠═759323db-7d17-482a-a8f6-ef32524f8126
# ╠═999a96e4-58bb-4b15-96cb-1d189ea9c1dd
# ╠═898459dd-3a04-46cc-b647-3dd7ef2664c7
# ╠═0a694982-3836-4034-9659-961043e57e4b
# ╠═aa6cefe0-d6a1-4745-b77a-9ad7327c0c22
