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
end

# ╔═╡ c28ca588-b092-4dd5-b96b-7692ebd9a540
Pkg.add("ImageShow")

# ╔═╡ 5855e0be-5416-4ed9-824a-d9ae85e4e977
md" ## Leitura da imagem"

# ╔═╡ 195e65e3-8c48-422a-8f37-8435d65b03b6
begin
	imagem_original = load("cat.png")
end

# ╔═╡ 6b4a91cb-1a3b-4af5-9915-4e966a77dc0f
md" ## Conversão RGB to YCC"

# ╔═╡ 712bc6c3-8928-43cf-8ffb-16895f078d97
md" ## Taxa de compresão com a subamostragem


Até agora diminuios os canais de cores 4 vezes e mantivemos o canal de luminancia, portanco reduzimos pela metade o tamanho da imagem"

# ╔═╡ 0ccd9de5-8554-4195-af60-ea52a682087f
md" ## Reconstrução da imagem
"

# ╔═╡ 12d12111-2ce6-41b0-b337-9587c9b70307
int_filter = [0.25 0.5 0.25 ;
				0.5 1.0 0.5;
				0.25 0.5 0.25]

# ╔═╡ 483e3ae6-9f83-4017-afe2-e03bd19691c7
md" # DCT"

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
	
	Cbsub = Cb[1:2:end, 1:2:end]
	Crsub = Cr[1:2:end, 1:2:end]
	noprint
end

# ╔═╡ 254611e5-fbac-4dfc-9433-e68b808bfa1e
begin
	tamanho_original_imagem = n_element(imagem_original)
	tamanho_int = n_element(Y) + n_element(Cbsub) + n_element(Crsub)
	compresao_int = tamanho_int/tamanho_original_imagem
end

# ╔═╡ 39f85367-dbe5-4422-8731-d48023465c66
begin
	Y_expand=padarray(Y,Pad(:symmetric, linhas_add, colunas_add))[1:end, 1:end] .- 0.5
	Cb_expand = padarray(Cb, Pad(:symmetric, linhas_add, colunas_add))[1:end, 1:end]
	Cr_expand = padarray(Cr, Pad(:symmetric, linhas_add, colunas_add))[1:end, 1:end]
	
	Cb_sub = Cb_expand[1:2:end, 1:2:end] .- 0.5
	Cr_sub = Cr_expand[1:2:end, 1:2:end] .- 0.5
end

# ╔═╡ b4fb610a-faaa-499d-8cb1-a0dc9b555d11
begin
	numx_blocos_Y =size(Y_expand)[1]÷8
	numy_blocos_Y =size(Y_expand)[2]÷8
	Y_dct = zeros(size(Y_expand)[1], size(Y_expand)[2])
	
	for i in 1:numx_blocos_Y
		for j in 1:numy_blocos_Y
			bloco = Y_expand[(i-1)*8 + 1:(i)*8 + 1, (j-1)*8+1:(j)*8+1] 
		end
	end
	
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
	maxi = 2^-8 -1
	PSNR_red = MSE(red.(imagem_original), red.(imagem_reconstruida))
	PSNR_green = MSE(green.(imagem_original), green.(imagem_reconstruida))
	PSNR_blue = MSE(blue.(imagem_original), blue.(imagem_reconstruida))
	PSNR_medio = mean([PSNR_red, PSNR_green, PSNR_blue])
end

# ╔═╡ e4a1a8c8-652b-4d85-9c7a-805457da6909
md" PSNR médio dos 3 canais da imagem : **$(round(PSNR_medio, digits = 5))**"

# ╔═╡ 8440deed-d2ed-433e-bde9-73aa1acfca91
function PSNR(imagem_exata, imagem)
	PSNR = 0
	maxi = 2^8 -1
	MSE_ = MSE(imagem_exata, imagem)
	PSNR = 10*log10(maxi/ MSE_)
	return PSNR
end

# ╔═╡ Cell order:
# ╠═5cd0e390-e111-11ec-31c9-d93bec27fd1c
# ╠═c28ca588-b092-4dd5-b96b-7692ebd9a540
# ╟─5855e0be-5416-4ed9-824a-d9ae85e4e977
# ╠═195e65e3-8c48-422a-8f37-8435d65b03b6
# ╟─6b4a91cb-1a3b-4af5-9915-4e966a77dc0f
# ╠═944839eb-9d6d-439b-8a69-4ca18eb8698a
# ╠═3f1f0ad4-bbc9-4fc7-9ecd-d6ec324461e1
# ╠═d41027c2-769d-4597-b203-fb6dca67037e
# ╠═712bc6c3-8928-43cf-8ffb-16895f078d97
# ╠═254611e5-fbac-4dfc-9433-e68b808bfa1e
# ╠═0ccd9de5-8554-4195-af60-ea52a682087f
# ╠═12d12111-2ce6-41b0-b337-9587c9b70307
# ╠═64db704f-ad48-4246-9d18-5ac92733a13f
# ╠═04eb8d3c-6da1-4d24-a99f-af85e8256c6a
# ╠═e4a1a8c8-652b-4d85-9c7a-805457da6909
# ╠═483e3ae6-9f83-4017-afe2-e03bd19691c7
# ╠═954f37df-a215-4c6a-897d-0e37fa88f9f4
# ╠═39f85367-dbe5-4422-8731-d48023465c66
# ╠═b4fb610a-faaa-499d-8cb1-a0dc9b555d11
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
