### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 5cd0e390-e111-11ec-31c9-d93bec27fd1c
begin
	using Pkg
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

# ╔═╡ aa423a97-fa9b-4378-b2ae-88102430a63f
md" # JPEG

Gabriel Tavares 10773801

Guilherme Reis 10773700
"

# ╔═╡ 5855e0be-5416-4ed9-824a-d9ae85e4e977
md" ## Leitura da imagem"

# ╔═╡ 195e65e3-8c48-422a-8f37-8435d65b03b6
begin
	imagem_original = load("cat.png")
end

# ╔═╡ 6b4a91cb-1a3b-4af5-9915-4e966a77dc0f
md" # Subamostragem de Cb e Cr

A imagem será convertida do espaço RGB para o espaço YCbCr e depois os canis Cb e Cr serão subamostrados."

# ╔═╡ 712bc6c3-8928-43cf-8ffb-16895f078d97
md" #### Taxa de compresão com a subamostragem


Até agora diminuios os canais de cores 4 vezes e mantivemos o canal de luminancia, portanco reduzimos pela **metade o tamanho da imagem**"

# ╔═╡ 0ccd9de5-8554-4195-af60-ea52a682087f
md" ## Reconstrução da imagem

Com os canais Cb e Cr subamostrados, iremos reconstruir a imagem usando um interpolador de média de 2 pontos na imagem para recriar os canais de cor no tamanho original.
"

# ╔═╡ 12d12111-2ce6-41b0-b337-9587c9b70307
int_filter = [0.25  0.5 	0.25;
			  0.5  	1.0 	0.5 ;
			  0.25 	0.5 	0.25]

# ╔═╡ 483e3ae6-9f83-4017-afe2-e03bd19691c7
md" # DCT"

# ╔═╡ 1f687d20-37cb-4c0c-b3e7-54b08a814f1f
md" #### Expansão da imagem 

Aqui a imagens devem ter um tamanho multiplo de 8 para definirmos os blocos de DCT na imagem"

# ╔═╡ 7fc069f9-0592-4ad2-8260-244163953a52
md" #### DCT do blocos

A imagem é divida em blocos de 8x8 e em cada bloco é feita a DCT desse bloco.

"

# ╔═╡ 036e142f-9562-4784-ab44-33837b1d0d74
md" ## Quantização das DCTS


Depois de ter o blocos transformados, iremos quatizar os coeficientes da DCT seguindo a tabela de quantização :

$$Q = \left[\begin{array}{cc} 
8 & 11 & 10 & 16 & 24 & 40 & 51 & 61 \\
12 & 12 & 14 & 19 & 26 & 58 & 60 & 55 \\
14 & 13 & 16 & 24 & 40 & 57 & 69 & 56 \\
14 & 17 & 22 & 29 & 51 & 87 & 80 & 62 \\
18 & 22 & 37 & 56 & 68 & 109 & 103 & 77 \\
24 & 35 & 55 & 64 & 81 & 104 & 113 & 92 \\
49 & 64 & 78 & 87 & 103 & 121 & 120 & 101 \\
72 & 92 & 95 & 98 & 112 & 100 & 103 & 99
\end{array}\right]$$ 
"

# ╔═╡ bc9097c6-0533-40e6-98b0-895d2bba2c80
md" ## Reconstrução

Com os coeficientes quantizados, iremos reconstruir a imagem multiplicado os coneficientes quantizados pela matriz de quantização novamente e analizar a diferença entre as imagens."

# ╔═╡ 6eb66248-558f-4394-b507-a792e42aaa5e
md" Abaixo temos a imagem original à esquerda e a imagem reconstruída à direita"

# ╔═╡ d9b6eb4f-82ae-495d-a2d3-56c12ef9f9e5
md" #### Número de coeficientes nulos"

# ╔═╡ 97462fbe-3c1e-49d5-ae6b-c264f320b6c9
md" ## Entropia"

# ╔═╡ be853fe4-851d-4b12-ac90-e53b0a130ae0
md" #### Coeficientes DC"

# ╔═╡ 353b5a34-3f66-42e7-a49a-2bdfdfaf0c7a
md" #### Análise de correlação dos coeficientes DC

Na imagem abaixo estamos plotando os coeficientes DC da transformada de luminância para observarmos que eles praticamente formam a imagem de volta subamostrada. Isso é um indicativo da correlação existente entre trechos da imagem.

"

# ╔═╡ 5d8b495c-1120-4d5f-835b-5a60a5bccfc0
md" #### AC"

# ╔═╡ 29f033be-9a57-40f2-9429-7a658a585cd8
md" # Functions"

# ╔═╡ 2fbdc13c-3ca1-4faa-8353-ba3afccd251b
noprint = md""

# ╔═╡ 944839eb-9d6d-439b-8a69-4ca18eb8698a
begin
	αr = 0.299
	αg = 0.587
	αb = 0.114
	noprint
end

# ╔═╡ 3f1f0ad4-bbc9-4fc7-9ecd-d6ec324461e1
begin
	R = red.(imagem_original)
	G = green.(imagem_original)
	B = blue.(imagem_original)
	noprint
end

# ╔═╡ 954f37df-a215-4c6a-897d-0e37fa88f9f4
begin
	linhas_original = size(imagem_original)[1]
	colunas_original = size(imagem_original)[2]
	
	linhas_add = 8 - (linhas_original - (linhas_original ÷  8)*8)
	colunas_add = 8 - (colunas_original - (colunas_original ÷  8)*8)
	
	linhas_expand = linhas_original + linhas_add
	colunas_expand = colunas_original + colunas_add
	noprint
end

# ╔═╡ 4c8da3f3-157e-4d05-99ea-6f134261a827
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
	
	k = 3
	
	noprint
end

# ╔═╡ 94a7f969-f06e-4b53-b99d-c4cdfeea6a35
begin
	filt_descorrelacao = 
	   [0.25  	0.5  	0.25
		0.5   	1    	0.5
		0.25 	0.5  	0.25]
	noprint
end

# ╔═╡ 51dce6e1-712f-4061-8f54-b6be746ee431
function n_element(matriz)
	return size(matriz)[1]*size(matriz)[2]
end

# ╔═╡ 5cfd422e-7411-485a-a513-7b3312f25dba
md" ### Gray matrix images functions"

# ╔═╡ 871c8016-4807-4f8d-8194-11f54b59f47e
function mat_to_image(mat)
	return RGB.(mat,mat,mat)	
end

# ╔═╡ 2a7afdf6-81fb-4a01-b072-42f041bc712b
function image(mat)
	return RGB.(mat,mat,mat)
end

# ╔═╡ dc11f115-1cef-4a1a-96d6-7844e73eeeef
md" ### RGB and YCbCr functions

"

# ╔═╡ 075a4d40-e26a-41ff-b879-f0cbdee3ac90
struct YCbCr_pixel
	Y
	Cb
	Cr
end

# ╔═╡ eb3e8ef2-eea1-429d-88c1-a6da088b841b
function luminancia(pixel::YCbCr_pixel)
	return pixel.Y
end

# ╔═╡ f539b11c-d588-436d-9ff9-8ef74c29e355
function croma_blue(pixel::YCbCr_pixel)
	return pixel.Cb
end

# ╔═╡ 35244517-13fc-4861-bf06-51c0f626951c
function croma_red(pixel::YCbCr_pixel)
	return pixel.Cr
end

# ╔═╡ 818376c8-844f-497d-bea3-238898e6b4aa
function get_Y(pixel::RGB)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	Y = αr*red(pixel) +αg*green(pixel) + αb*blue(pixel)
	
	return Y
end

# ╔═╡ 1da88d32-e54a-487d-ac35-28d3f6fde9d8
function get_Cb(pixel::RGB)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	Y = get_Y(pixel)
	Cb = 1/(2*(1-αb))*(blue(pixel)-Y)
end

# ╔═╡ d61d27fc-2b19-4113-adca-10d06fee748e
function get_Cr(pixel::RGB)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	Y = get_Y(pixel)
	Cr = 1/(2*(1-αr))*(red(pixel)-Y)
end

# ╔═╡ d41027c2-769d-4597-b203-fb6dca67037e
begin
	Y = get_Y.(imagem_original)
	Cb = get_Cb.(imagem_original)
	Cr = get_Cr.(imagem_original)
	
	Cbsub = Cb[1:2:end, 1:2:end] # Subamostragem de Cb
	Crsub = Cr[1:2:end, 1:2:end] # Subamostragem de Cr
	noprint
end

# ╔═╡ 254611e5-fbac-4dfc-9433-e68b808bfa1e
begin
	tamanho_original_imagem = n_element(imagem_original) * 3
	tamanho_int = n_element(Y) + n_element(Cbsub) + n_element(Crsub)
	compresao_int = round(tamanho_int/tamanho_original_imagem, digits = 2)
end

# ╔═╡ 39f85367-dbe5-4422-8731-d48023465c66
begin
	Y_expand=padarray(Y,Pad(:symmetric, linhas_add, colunas_add))[1:end, 1:end] .- 0.5
	Cb_expand = padarray(Cb, Pad(:symmetric, linhas_add, colunas_add))[1:end, 1:end]
	Cr_expand = padarray(Cr, Pad(:symmetric, linhas_add, colunas_add))[1:end, 1:end]
	
	Cb_sub = Cb_expand[1:2:end, 1:2:end] .- 0.5 #expandido e subamostrado
	Cr_sub = Cr_expand[1:2:end, 1:2:end] .- 0.5 #expandido e subamostrado
	
	noprint
end

# ╔═╡ b4fb610a-faaa-499d-8cb1-a0dc9b555d11
begin
	numx_blocos_Y =size(Y_expand)[1]÷8
	numy_blocos_Y =size(Y_expand)[2]÷8
	Y_dct = zeros(size(Y_expand))
	
	for i in 1:numx_blocos_Y
		for j in 1:numy_blocos_Y
			bloco = Y_expand[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8] 
			Y_dct[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]  = dct(bloco)
		end
	end
	
	numx_blocos_Cb =size(Cb_sub)[1]÷8
	numy_blocos_Cb =size(Cb_sub)[2]÷8
	Cb_dct = zeros(size(Cb_sub))
	
	for i in 1:numx_blocos_Cb
		for j in 1:numy_blocos_Cb
			bloco = Cb_sub[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8] 
			Cb_dct[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]  = dct(bloco)
		end
	end
	
	numx_blocos_Cr =size(Cr_sub)[1]÷8
	numy_blocos_Cr =size(Cr_sub)[2]÷8
	Cr_dct = zeros(size(Cr_sub))
	
	for i in 1:numx_blocos_Cr
		for j in 1:numy_blocos_Cr
			bloco = Cr_sub[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8] 
			Cr_dct[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]  = dct(bloco)
		end
	end
end

# ╔═╡ 85ba0d5e-e999-4b6b-a394-9159bced3d80
begin
	probabilidades_Y_AC = zeros(2047*2 +1)
	for i in 1:length(probabilidades_Y_AC)
		probabilidades_Y_AC[i] = sum(Y_dct .== i-2048)
	end
	probabilidades_Y_AC = probabilidades_Y_AC/n_element(Y_dct)
	
	
	probabilidades_Cb_AC = zeros(2047*2 +1)
	for i in 1:length(probabilidades_Cb_AC)
		probabilidades_Cb_AC[i] = sum(Cb_dct .== i-2048)
	end
	probabilidades_Cb_AC = probabilidades_Cb_AC./n_element(Cb_dct)
	
	probabilidades_Cr_AC = zeros(2047*2 +1)
	for i in 1:length(probabilidades_Cr_AC)
		probabilidades_Cr_AC[i] = sum(Cr_dct .== i-2048)
	end
	probabilidades_Cr_AC = probabilidades_Cr_AC./n_element(Cr_dct)
	noprint
end

# ╔═╡ aa633ca2-dbad-4081-8553-5bc51f123eac
begin
	#Quantização de Y
	
	Y_dct_q = zeros(size(Y_expand))
	
	for i in 1:numx_blocos_Y
		for j in 1:numy_blocos_Y
			bloco = Y_dct[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]*255
			Y_dct_q[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]  = round.(bloco./(k*Q))/255
		end
	end
	
	Cb_dct_q = zeros(size(Cb_sub))
	
	for i in 1:numx_blocos_Cb
		for j in 1:numy_blocos_Cb
			bloco = Cb_dct[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]*255
			Cb_dct_q[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]  = round.(bloco./(k*Q))/255
		end
	end
	
	Cr_dct_q = zeros(size(Cb_sub))
	
	for i in 1:numx_blocos_Cr
		for j in 1:numy_blocos_Cr
			bloco = Cr_dct[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]*255
			Cr_dct_q[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]  = round.(bloco./(k*Q))/255
		end
	end
	
end

# ╔═╡ c97fac72-2daa-489e-8273-22a4bd21fe70
begin
	nulos_Cb = sum(Cb_dct_q .== 0)
	nulos_Cr = sum(Cr_dct_q .== 0)	
	nulos_Y  = sum(Y_dct_q .== 0)
	
	razao_nulos = (sum([nulos_Cb, nulos_Cr, nulos_Y]))/
		sum([n_element(Cb_dct_q),n_element(Cr_dct_q),n_element(Y_dct_q)])
	noprint
end

# ╔═╡ 54a7f7cb-448b-4beb-8d33-cc3a71795de0
md" Com isso vemos que **$(round(razao_nulos*100; digits = 2))%** dos coeficientes são nulos, o que permite uma codificação dos valores que diminua muito o tamanho da imagem."

# ╔═╡ 01e2e642-4406-4335-8986-6ba3ee92923b
begin
	Y_dct_DC = round.(Y_dct_q[1:8:end, 1:8:end] * 2047)
	Cb_dct_DC = round.(Cb_dct_q[1:8:end, 1:8:end] * 2047)
	Cr_dct_DC = round.(Cr_dct_q[1:8:end, 1:8:end] * 2047)
	
	noprint
end

# ╔═╡ 77af673d-d4ae-4328-aa5c-557d534e84dc
begin
	probabilidades_Y_DC = zeros(2047*2 +1)
	for i in 1:length(probabilidades_Y_DC)
		probabilidades_Y_DC[i] = sum(Y_dct_DC .== i-2048)
	end
	probabilidades_Y_DC = probabilidades_Y_DC/n_element(Y_dct_DC)
	
	
	probabilidades_Cb_DC = zeros(2047*2 +1)
	for i in 1:length(probabilidades_Cb_DC)
		probabilidades_Cb_DC[i] = sum(Cb_dct_DC .== i-2048)
	end
	probabilidades_Cb_DC = probabilidades_Cb_DC./n_element(Cb_dct_DC)
	
	probabilidades_Cr_DC = zeros(2047*2 +1)
	for i in 1:length(probabilidades_Cr_DC)
		probabilidades_Cr_DC[i] = sum(Cr_dct_DC .== i-2048)
	end
	probabilidades_Cr_DC = probabilidades_Cr_DC./n_element(Cr_dct_DC)
	noprint
end

# ╔═╡ 4223c82b-5170-4676-a4b4-4df0554581df
begin
	H_DC_Y = 0
	for k in 1: length(probabilidades_Y_DC)
		pk = probabilidades_Y_DC[k]
		H_DC_Y += pk !=0 ? - log2(pk) * pk : 0
	end
	
	H_DC_Cb = 0
	for k in 1: length(probabilidades_Cb_DC)
		pk = probabilidades_Cb_DC[k]
		H_DC_Cb += pk!=0 ? - log2(pk) * pk : 0
	end
	
	H_DC_Cr = 0
	for k in 1: length(probabilidades_Cr_DC)
		pk = probabilidades_Cr_DC[k]
		H_DC_Cr += pk!=0 ? - log2(pk) * pk : 0
	end
	
end

# ╔═╡ 7611b105-0287-4da5-8f8e-68ee9d9a8eb2
md" Entropia de Y é $(round(H_DC_Y, digits = 2)), portanto o número mínimo de bits para a transmissão do canal de luminância sem erros é $(ceil(H_DC_Y)) bits

Entropia de Cb é $(round(H_DC_Cb, digits = 2)), portanto o número mínimo de bits para a transmissão do canal de luminância sem erros é $(ceil(H_DC_Cb)) bits

Entropia de Cr é $(round(H_DC_Cr, digits = 2)), portanto o número mínimo de bits para a transmissão do canal de luminância sem erros é $(ceil(H_DC_Cr)) bits
"

# ╔═╡ 8465d820-557d-469d-af3b-1a4270f2513b
begin
	mat_to_image(Y_dct_DC./50)
end

# ╔═╡ de976d83-b9e6-44a7-91f3-a857c1c5a9cd
begin
	Y_dct_DC_desc = round.(imfilter(Y_dct_DC, filt_descorrelacao))
	Cb_dct_DC_desc = round.(imfilter(Cb_dct_DC, filt_descorrelacao))
	Cr_dct_DC_desc = round.(imfilter(Cr_dct_DC, filt_descorrelacao))
	noprint
end

# ╔═╡ c5dfc3ae-c8b6-4806-b2c7-cf79edec8db0
begin
	num_bit_Y = ceil(H_DC_Y)
	num_bit_Cb = ceil(H_DC_Cb)	
	num_bit_Cr = ceil(H_DC_Cr)
	
	tamanho_Y = num_bit_Y*n_element(Y_dct_q)
	tamanho_Cb = num_bit_Cb*n_element(Cr_dct_q)
	tamanho_Cr = num_bit_Cr*n_element(Cb_dct_q)
	
	tamanho_total = tamanho_Y +tamanho_Cb + tamanho_Cr
	
end

# ╔═╡ 2b6fa440-767b-46c4-9c60-fe47be186c26
md" Por fim, o tamanho total da imagem seria de $(round(tamanho_total/1000_000, digits = 2))Mbytes sem considerar uma compactação de agrupamentos"

# ╔═╡ b80899fd-9fcf-4441-bd2b-f450b4742a73
begin
	#DCT inversa
	
	Y_rec_expand = zeros(size(Y_expand))
	
	for i in 1:numx_blocos_Y
		for j in 1:numy_blocos_Y
			bloco = Y_dct_q[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]*255
			bloco = bloco .* (k*Q)
			bloco = bloco./255
			bloco = idct(bloco)
			bloco = bloco .+ 0.5
			Y_rec_expand[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]  .= bloco
		end
	end
	
	Cb_rec_expand = zeros(size(Cb_sub))
	
	for i in 1:numx_blocos_Cb
		for j in 1:numy_blocos_Cb
			bloco = Cb_dct_q[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]*255
			bloco = bloco .* (k*Q)
			bloco = bloco./255
			bloco = idct(bloco)
			bloco = bloco .+ 0.5
			Cb_rec_expand[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]  .= bloco
		end
	end
	
	Cr_rec_expand = zeros(size(Cr_sub))
	
	for i in 1:numx_blocos_Cr
		for j in 1:numy_blocos_Cr
			bloco = Cr_dct_q[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]*255
			bloco = bloco .* (k*Q)
			bloco = bloco./255
			bloco = idct(bloco)
			bloco = bloco .+ 0.5
			Cr_rec_expand[(i-1)*8 + 1:(i)*8, (j-1)*8+1:(j)*8]  .= bloco
		end
	end
	
end

# ╔═╡ 420dfd0d-c5e1-4f8e-b7b0-2b20de182cba
begin
	#Interpolação de Cb e Cr
	Cb_rec_int_expand = zeros(size(Cb_expand))
	Cb_rec_int_expand[1:2:end, 1:2:end] = Cb_rec_expand
	Cb_rec_int_expand = imfilter(Cb_rec_int_expand, int_filter)
	
	Cr_rec_int_expand = zeros(size(Cr_expand))
	Cr_rec_int_expand[1:2:end, 1:2:end] = Cr_rec_expand
	Cr_rec_int_expand = imfilter(Cr_rec_int_expand, int_filter)
	noprint
end

# ╔═╡ 5e857a99-9f2f-470e-97e4-fc4a9f4854a4
begin 
	#Recuperação do tamanho original da imagem
	Y_rec = Y_rec_expand[1:linhas_original, 1:colunas_original]
	Cb_rec =  Cb_rec_int_expand[1:linhas_original, 1:colunas_original]
	Cr_rec =  Cr_rec_int_expand[1:linhas_original, 1:colunas_original]
	noprint
	
end

# ╔═╡ 6f001133-37b1-47f2-862e-00d3cdb51615
function get_YCbCr(pixel::RGB)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	Y = get_Y(pixel)
	Cb = get_Cb(pixel)
	Cr = get_Cr(pixel)
	return Y, Cb, Cr
end

# ╔═╡ aa81d420-7bf5-4e59-8a31-ad4bf46535c3
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

# ╔═╡ 521314f0-c659-4e6f-8f29-747b25cadb69
function get_R(Y, Cb, Cr)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	return Y + (2-2*αr)*Cr
end

# ╔═╡ baf0e37e-976b-4d89-835c-f852aba1e42e
function get_G(Y, Cb, Cr)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	return Y - αb/αg*(2-2*αb)*Cb -αr/αg*(2-2*αr)*Cr
end

# ╔═╡ 2ace941c-b958-4ba2-a80f-6225890b2ad2
function get_B(Y, Cb, Cr)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	R = get_R(Y,Cb,Cr)
	G = get_G(Y, Cb, Cr)
	
	return Y + (2-2*αb)*Cb
end

# ╔═╡ 384b7a5b-ebe2-4113-8caa-8b35fe47f9f8
function get_RGB(Y, Cb, Cr)
	αr = 0.299
	αg = 0.587
	αb = 0.114
	
	R = get_R(Y,Cb,Cr)
	G = get_G(Y, Cb, Cr)
	B = get_B(Y, Cb, Cr)
	
	return RGB(R,G,B)
end

# ╔═╡ 64db704f-ad48-4246-9d18-5ac92733a13f
begin
	Cb_int = zeros(size(Y))
	Cb_int[1:2:end, 1:2:end] = Cbsub
	Cb_int = imfilter(Cb_int, int_filter)
	
	Cr_int = zeros(size(Y))
	Cr_int[1:2:end, 1:2:end] = Crsub
	Cr_int = imfilter(Cr_int, int_filter)
	
	RGB_int = get_RGB.(Y, Cb_int, Cr_int)
	
	imagem_reconstruida = RGB_int
end

# ╔═╡ 87a128dc-33b6-4edd-8f6f-6b5191dd8432
hcat(imagem_original, get_RGB.(Y_rec,Cb_rec, Cr_rec))

# ╔═╡ 7ea72842-901f-4093-ab3c-00ccd894e23b
md" ## PSNR"

# ╔═╡ fe61f917-e874-40c7-8761-32e3f7c1ab73
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

# ╔═╡ 04eb8d3c-6da1-4d24-a99f-af85e8256c6a
begin
	PSNR_red = MSE(red.(imagem_original), red.(imagem_reconstruida))
	PSNR_green = MSE(green.(imagem_original), green.(imagem_reconstruida))
	PSNR_blue = MSE(blue.(imagem_original), blue.(imagem_reconstruida))
	PSNR_medio = mean([PSNR_red, PSNR_green, PSNR_blue])
	noprint
end

# ╔═╡ e4a1a8c8-652b-4d85-9c7a-805457da6909
md" #### Análise da reconstrução

Qualitativamente vemos que a qualidade da imagem permanece muito próxima da original, já que o olho humano é menos sensível a cor do que a luz.

Analiticamente o PSNR médio dos 3 canais da imagem reconstruída : **$(round(PSNR_medio, digits = 5))dB**, o que explica a alta qualidade da imagem em compração com a a original"

# ╔═╡ 0a051a2d-83bf-4595-9c75-5cc6a016661e
begin
	imagem_rec = get_RGB.(Y_rec,Cb_rec, Cr_rec)
	
	PSNR_red_ = MSE(red.(imagem_original), red.(imagem_rec))
	PSNR_green_ = MSE(green.(imagem_original), green.(imagem_rec))
	PSNR_blue_ = MSE(blue.(imagem_original), blue.(imagem_rec))
	PSNR_medio_ = mean([PSNR_red_, PSNR_green_, PSNR_blue_])
	noprint
end

# ╔═╡ 15359b2f-a203-491d-8aeb-b53515c907b5
md" #### Análise

Qualitativamente observamos que as imagens são relativamente parecidas. É possível ver a formação de blocos quadrados mais no fundo da imagem, mas o conteúdo principal ainda é bom. Essa qualidade é bastante dependente do valor de **k**.

* k = 1 : a imagem permanece com uma qualidade muito boa e os erros são pouco perceptíveis.
* k = 3 : começamos a reparar os blocos de na imagem, principalemnte ao fundo
* k = 5 : os blocos ficam muito visíveis no fundo e começam a aparecer no gato


Analiticamente o PSNR médio dos 3 canais da imagem reconstruída com **k = $k** é: **$(round(PSNR_medio_, digits = 5))dB**"

# ╔═╡ 8440deed-d2ed-433e-bde9-73aa1acfca91
function PSNR(imagem_exata, imagem)
	PSNR = 0
	maxi = 1
	MSE_ = MSE(imagem_exata, imagem)
	PSNR = 10*log10(maxi/ MSE_)
	return PSNR
end

# ╔═╡ Cell order:
# ╟─aa423a97-fa9b-4378-b2ae-88102430a63f
# ╠═5cd0e390-e111-11ec-31c9-d93bec27fd1c
# ╟─5855e0be-5416-4ed9-824a-d9ae85e4e977
# ╠═195e65e3-8c48-422a-8f37-8435d65b03b6
# ╟─6b4a91cb-1a3b-4af5-9915-4e966a77dc0f
# ╠═944839eb-9d6d-439b-8a69-4ca18eb8698a
# ╠═3f1f0ad4-bbc9-4fc7-9ecd-d6ec324461e1
# ╠═d41027c2-769d-4597-b203-fb6dca67037e
# ╟─712bc6c3-8928-43cf-8ffb-16895f078d97
# ╠═254611e5-fbac-4dfc-9433-e68b808bfa1e
# ╟─0ccd9de5-8554-4195-af60-ea52a682087f
# ╠═12d12111-2ce6-41b0-b337-9587c9b70307
# ╠═64db704f-ad48-4246-9d18-5ac92733a13f
# ╠═04eb8d3c-6da1-4d24-a99f-af85e8256c6a
# ╟─e4a1a8c8-652b-4d85-9c7a-805457da6909
# ╟─483e3ae6-9f83-4017-afe2-e03bd19691c7
# ╟─1f687d20-37cb-4c0c-b3e7-54b08a814f1f
# ╠═954f37df-a215-4c6a-897d-0e37fa88f9f4
# ╠═39f85367-dbe5-4422-8731-d48023465c66
# ╟─7fc069f9-0592-4ad2-8260-244163953a52
# ╠═b4fb610a-faaa-499d-8cb1-a0dc9b555d11
# ╠═036e142f-9562-4784-ab44-33837b1d0d74
# ╠═4c8da3f3-157e-4d05-99ea-6f134261a827
# ╠═aa633ca2-dbad-4081-8553-5bc51f123eac
# ╟─bc9097c6-0533-40e6-98b0-895d2bba2c80
# ╠═b80899fd-9fcf-4441-bd2b-f450b4742a73
# ╠═420dfd0d-c5e1-4f8e-b7b0-2b20de182cba
# ╠═5e857a99-9f2f-470e-97e4-fc4a9f4854a4
# ╟─6eb66248-558f-4394-b507-a792e42aaa5e
# ╟─87a128dc-33b6-4edd-8f6f-6b5191dd8432
# ╠═0a051a2d-83bf-4595-9c75-5cc6a016661e
# ╟─15359b2f-a203-491d-8aeb-b53515c907b5
# ╟─d9b6eb4f-82ae-495d-a2d3-56c12ef9f9e5
# ╠═c97fac72-2daa-489e-8273-22a4bd21fe70
# ╟─54a7f7cb-448b-4beb-8d33-cc3a71795de0
# ╟─97462fbe-3c1e-49d5-ae6b-c264f320b6c9
# ╟─be853fe4-851d-4b12-ac90-e53b0a130ae0
# ╠═01e2e642-4406-4335-8986-6ba3ee92923b
# ╠═77af673d-d4ae-4328-aa5c-557d534e84dc
# ╠═85ba0d5e-e999-4b6b-a394-9159bced3d80
# ╠═4223c82b-5170-4676-a4b4-4df0554581df
# ╟─7611b105-0287-4da5-8f8e-68ee9d9a8eb2
# ╠═c5dfc3ae-c8b6-4806-b2c7-cf79edec8db0
# ╟─2b6fa440-767b-46c4-9c60-fe47be186c26
# ╟─353b5a34-3f66-42e7-a49a-2bdfdfaf0c7a
# ╟─8465d820-557d-469d-af3b-1a4270f2513b
# ╟─5d8b495c-1120-4d5f-835b-5a60a5bccfc0
# ╠═94a7f969-f06e-4b53-b99d-c4cdfeea6a35
# ╠═de976d83-b9e6-44a7-91f3-a857c1c5a9cd
# ╟─29f033be-9a57-40f2-9429-7a658a585cd8
# ╠═2fbdc13c-3ca1-4faa-8353-ba3afccd251b
# ╠═51dce6e1-712f-4061-8f54-b6be746ee431
# ╟─5cfd422e-7411-485a-a513-7b3312f25dba
# ╠═871c8016-4807-4f8d-8194-11f54b59f47e
# ╠═2a7afdf6-81fb-4a01-b072-42f041bc712b
# ╠═dc11f115-1cef-4a1a-96d6-7844e73eeeef
# ╠═075a4d40-e26a-41ff-b879-f0cbdee3ac90
# ╠═eb3e8ef2-eea1-429d-88c1-a6da088b841b
# ╠═f539b11c-d588-436d-9ff9-8ef74c29e355
# ╠═35244517-13fc-4861-bf06-51c0f626951c
# ╠═818376c8-844f-497d-bea3-238898e6b4aa
# ╠═1da88d32-e54a-487d-ac35-28d3f6fde9d8
# ╠═d61d27fc-2b19-4113-adca-10d06fee748e
# ╠═6f001133-37b1-47f2-862e-00d3cdb51615
# ╠═aa81d420-7bf5-4e59-8a31-ad4bf46535c3
# ╠═521314f0-c659-4e6f-8f29-747b25cadb69
# ╠═baf0e37e-976b-4d89-835c-f852aba1e42e
# ╠═2ace941c-b958-4ba2-a80f-6225890b2ad2
# ╠═384b7a5b-ebe2-4113-8caa-8b35fe47f9f8
# ╠═7ea72842-901f-4093-ab3c-00ccd894e23b
# ╠═fe61f917-e874-40c7-8761-32e3f7c1ab73
# ╠═8440deed-d2ed-433e-bde9-73aa1acfca91
