### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 6662db60-23de-11ed-3ee1-83ae6c1eb90c
begin
	using DSP
	using Plots
	plotly()
end

# ╔═╡ 85b8d338-200f-4e50-b35f-f755ccc3470c
md" # 1 Criação dos sinais"

# ╔═╡ 1bee1342-76d4-46b8-989f-fc355e500916
md" # 2 Algoritmo LMS"

# ╔═╡ 04ca8de8-1a0d-4b6c-9a21-543f7dd402c5
begin 
	"""
		LMS(x, d, H0, α)
	
	* x : sinal original
	* d : sinal com interferência
	* H0 : coeficientes iniciais do filtro adaptativo
	* α : passo de adptação do filtro
	
	- Retorna uma matriz com a evolução do filtro ao longo do tempo. Cada linha dessa matriz é o estado do filtro no instante **n**
	
	"""
	function LMS(x, d, H0, α)
		if length(x) != length(d)
			throw(DimensionMismatch("sinais x e d devem ser do mesmo tamanho"))
		end
		
		M = length(H0) #num coeficientes do filtro
		
		x_filtro = zeros(M) #vetor que guarda o trecho de sinal que será filtrado 
		H = zeros(length(x), M)
		erros = zeros(length(x))
		
		for i in 1:length(x)-1
			x_filtro = [x[i];x_filtro[1:M-1]]
			y = H[i,:]'*x_filtro
			erro = d[i] - y
			erros[i] = erro
			H[i+1, :] = H[i,:] + α*erro*x_filtro
		end
		
		return H, erros
	end
end

# ╔═╡ c14b1233-2d3c-4d67-855a-2d92fea14c4f


# ╔═╡ d0b34c2d-f53c-4587-a144-3fc394610782
md" # Functions"

# ╔═╡ b9798d57-8ce0-4119-bbf5-c69a07b9bf48
noprint = md""

# ╔═╡ cb98290f-9897-4ab7-a489-c2a1cb1e8a8b
begin
	ω0 = π/6
	ω1 = π/20
	
	H = [1, -0.5]
	n = 1:10_000
	
	x = sin.(ω0*n)
	v = sin.(ω1*n)
	
	noprint
end

# ╔═╡ f3183376-8af7-4b24-badc-95f6d101f2bd
begin
	p1 = plot(x)
	p2 = plot(v)
	plot(p1,p2, layout=(2,1))
	plot!(title = ["x[n]" "v[n]"])
end

# ╔═╡ 15047666-5e20-44d7-94c9-1070d0cd542f
begin
	filt_adap, erro = LMS(x, x+v, H, 0.001)
	plot(erro)
	
end

# ╔═╡ 637ccfd8-ffbb-4a79-b763-99dcf20f7b18
erro

# ╔═╡ fd1f16ab-3dcc-4320-bd83-2316445d0214


# ╔═╡ f0410113-f959-44c9-bafd-fdd500e5895b
"""

# ╔═╡ e072d0b4-e184-4d31-8b6b-d5d4acffb674
bar([1, 2], [1, 2])

# ╔═╡ Cell order:
# ╠═6662db60-23de-11ed-3ee1-83ae6c1eb90c
# ╠═85b8d338-200f-4e50-b35f-f755ccc3470c
# ╠═cb98290f-9897-4ab7-a489-c2a1cb1e8a8b
# ╟─f3183376-8af7-4b24-badc-95f6d101f2bd
# ╟─1bee1342-76d4-46b8-989f-fc355e500916
# ╠═04ca8de8-1a0d-4b6c-9a21-543f7dd402c5
# ╠═15047666-5e20-44d7-94c9-1070d0cd542f
# ╠═637ccfd8-ffbb-4a79-b763-99dcf20f7b18
# ╠═c14b1233-2d3c-4d67-855a-2d92fea14c4f
# ╠═d0b34c2d-f53c-4587-a144-3fc394610782
# ╠═b9798d57-8ce0-4119-bbf5-c69a07b9bf48
# ╠═fd1f16ab-3dcc-4320-bd83-2316445d0214
# ╠═f0410113-f959-44c9-bafd-fdd500e5895b
# ╠═e072d0b4-e184-4d31-8b6b-d5d4acffb674
