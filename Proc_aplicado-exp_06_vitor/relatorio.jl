### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 680e5c90-fc8e-11ec-1304-0deea9b53834
begin
	using MAT
	using DSP
	using Plots
	plotly()
end

# ╔═╡ 12c8cb4e-bc06-44c3-975b-74d1b9323e29
md" ## 01"

# ╔═╡ 4d2fb09b-36a4-474b-a908-05b174b68c4f
begin
	M = 8
	F = 60e9 #Hz
	Ω = 2*π*F
	c = 3e8 #m/s
	λ = c/F #m
	# d = λ/4 #m
	j = im
	θ0 =  20 #graus
end

# ╔═╡ 5579a80a-978b-4050-ae57-9123080fc3a1
md" Calculo para ter ganho 1 em θ = 20°"

# ╔═╡ e6d1fc78-8723-418d-be3b-8eddc3430be0
md" ## 02"

# ╔═╡ 539fca48-8236-4734-a7bc-4049bc8d347b
md" ## 03"

# ╔═╡ 138691ba-b2ec-4bf8-ac09-0d0cc8ac3bcb
md" ## 04"

# ╔═╡ db7b6207-6b79-45c3-9047-072a83d0e2c8


# ╔═╡ ec484da9-740f-4792-8525-ea5d10da527b


# ╔═╡ 58b09e7c-9bd0-459e-8077-97c020758932
md" ## Functions"

# ╔═╡ f41cf0d5-af4b-49ee-b9ce-917cd6e57866
function B(θin, θ0in, M, d, Ω0)
	θ = deg2rad(θin)
	θ0 = deg2rad(θ0in)
	c = 3e8
	exponencial = exp( im*(M-1)/2 * Ω0*d/c *( sin(θ)-sin(θ0) ) )
	if (sin(θ)-sin(θ0))==0.0
		senos = M
	else
		senos = sin(M*Ω0*d/(2c) * (sin(θ)-sin(θ0)) ) / ( sin( Ω0*d/(2c) * (sin(θ)-sin(θ0) )) )
	end
	valor =exponencial *1/M* senos
	return valor
end

# ╔═╡ 72f16a75-fc9d-442c-be94-f443db657efa
function τ(θ, m, d; c = 3e8)
	return m*d*sind(θ)/c
end

# ╔═╡ 9ea61348-08e5-4070-812b-a013478afeb8
function B(a, θ, d ;c = 3e8)
	M = length(a)
	B_ = 0 + 0*j
	
	for m in 1:M
		B_ += a[m]*exp(-j*Ω*τ(θ,m, d; c = c))
	end
	return B_
end

# ╔═╡ 199b6831-1cf9-44fa-9afc-f1bc838ce090
function B(a, θ::Vector{Float64}, d ;c = 3e8)

	N = length(θ) #numero de graus do ganho
	B_ = zeros(Complex, N)
	
	for i in range(1, N)
		B_[i] = B(a, θ[i], d; c = c)	
	end

	return B_
end

# ╔═╡ edc6e700-01ee-4fa9-8eec-b98d323e1fa2
function B(a, θ::StepRangeLen{Float64}, d; c = 3e8)
	return B(a, collect(θ),d; c=c)
end

# ╔═╡ fecad3c1-f01d-49a0-9141-77577f29c6e3
function DAS(θ, M, d; Ω = Ω, c=3e8)
	a = zeros(Complex , M )
	for m in 1:M
		a[m] = 1/M * exp(j * Ω*τ(θ, m-1, d; c=c))
	end
	return a
end

# ╔═╡ b62af42c-3682-4787-90e2-718e41ac084e
begin
	a = DAS(θ0, M, λ/4)
	θ = -90:0.1:90 
	B20 = B(a, θ, λ/4)
end

# ╔═╡ 31f41ff9-0ab3-4b66-bd01-a4f5daac6375
begin
	plot(θ, abs.(B20))
end

# ╔═╡ 9d39c9db-db58-4c9a-af6e-1f829ad6ab0a
begin
	plot(θ, abs.(B.(θ0, θ, M, λ/4, Ω)) )
	plot!(θ, abs.(B.(θ0, θ, M, λ/2, Ω)) )	
	plot!(θ, abs.(B.(θ0, θ, M, 3λ/4, Ω)) )	
	plot!(θ, abs.(B.(θ0, θ, M, λ, Ω)) )	
end

# ╔═╡ 089f1d31-0dff-4bae-8526-190b83f24189
begin
	A1=1
	A2=0.5
	a3 = DAS(20, M, λ/2)
	A =  A1*abs(B(a3, 20, λ/2)) .+ A2.*abs.(B(a3, θ, λ/2))

end

# ╔═╡ 6d723d92-a4cd-44b3-8c1e-79bdd219d8bc
begin
	plot(θ, A)
end

# ╔═╡ Cell order:
# ╠═680e5c90-fc8e-11ec-1304-0deea9b53834
# ╠═12c8cb4e-bc06-44c3-975b-74d1b9323e29
# ╠═4d2fb09b-36a4-474b-a908-05b174b68c4f
# ╠═5579a80a-978b-4050-ae57-9123080fc3a1
# ╠═b62af42c-3682-4787-90e2-718e41ac084e
# ╠═31f41ff9-0ab3-4b66-bd01-a4f5daac6375
# ╠═e6d1fc78-8723-418d-be3b-8eddc3430be0
# ╠═9d39c9db-db58-4c9a-af6e-1f829ad6ab0a
# ╠═539fca48-8236-4734-a7bc-4049bc8d347b
# ╠═089f1d31-0dff-4bae-8526-190b83f24189
# ╠═6d723d92-a4cd-44b3-8c1e-79bdd219d8bc
# ╠═138691ba-b2ec-4bf8-ac09-0d0cc8ac3bcb
# ╠═db7b6207-6b79-45c3-9047-072a83d0e2c8
# ╠═ec484da9-740f-4792-8525-ea5d10da527b
# ╠═58b09e7c-9bd0-459e-8077-97c020758932
# ╠═f41cf0d5-af4b-49ee-b9ce-917cd6e57866
# ╠═72f16a75-fc9d-442c-be94-f443db657efa
# ╠═9ea61348-08e5-4070-812b-a013478afeb8
# ╠═199b6831-1cf9-44fa-9afc-f1bc838ce090
# ╠═edc6e700-01ee-4fa9-8eec-b98d323e1fa2
# ╠═fecad3c1-f01d-49a0-9141-77577f29c6e3
