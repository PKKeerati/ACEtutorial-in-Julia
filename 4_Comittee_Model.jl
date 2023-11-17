### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 15830db7-7ef6-4222-838d-1a86eff4ca85
begin
	using Pkg
	Pkg.activate(".")
	Pkg.add("LaTeXStrings")
	Pkg.add("MultivariateStats")
	Pkg.add("Plots")
	Pkg.add("Suppressor")
	using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor
end

# ╔═╡ 4a5d9b66-6f39-4c6b-856a-1f4cf873ed29
begin
	Pkg.activate(".")
	Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add("ACEpotentials")
	using ACEpotentials
end

# ╔═╡ dc0eccf0-8541-11ee-2421-0713757269ef
begin
	model = acemodel(elements = [:Si,],
	                 order = 3,
	                 totaldegree = 12,
	                 Eref = [:Si => -158.54496821]);
	
	acefit!(model, Si_tiny_dataset;
	        solver = ACEfit.BLR(committee_size=50, factorization=:svd),
	        energy_key = "dft_energy", force_key = "dft_force",
	        verbose = false);
end

# ╔═╡ Cell order:
# ╠═15830db7-7ef6-4222-838d-1a86eff4ca85
# ╠═4a5d9b66-6f39-4c6b-856a-1f4cf873ed29
# ╠═dc0eccf0-8541-11ee-2421-0713757269ef
