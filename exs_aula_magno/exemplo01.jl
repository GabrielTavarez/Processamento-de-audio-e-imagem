### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 621aae7e-aae5-11ec-22ba-9f5fdecb4088
begin
	using Plots
	using DSP
end

# ╔═╡ 9ddfd5da-4dc3-4002-96b7-91732746be30
begin
	num = [1, -1, 1]
	den = [1, -0.9, 0.81]
	h1 = PolynomialRatio(num,[1])
	h2 = PolynomialRatio(0.91*num, den)
end

# ╔═╡ 2abf7644-38bd-4962-bdc9-e267d9d687ca
begin
	ω = range(0,π, length = 500)
	H1 = freqz(h1, ω)
	H2 = freqz(h2, ω)
end

# ╔═╡ 51d4eda7-6b9c-494a-8f79-a6a6edf7f185
begin
	plot(ω/π, real.(H1))
	plot!(ω/π, real.(H2))
	plot!([1/3,1/3],[-1,3])
end

# ╔═╡ Cell order:
# ╠═621aae7e-aae5-11ec-22ba-9f5fdecb4088
# ╠═9ddfd5da-4dc3-4002-96b7-91732746be30
# ╠═2abf7644-38bd-4962-bdc9-e267d9d687ca
# ╠═51d4eda7-6b9c-494a-8f79-a6a6edf7f185
