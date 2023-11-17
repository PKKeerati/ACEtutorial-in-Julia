### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 81f97bc0-8546-11ee-3ef6-414aa5615b41
begin
	using Pkg
	Pkg.activate(".")
	Pkg.add("LaTeXStrings")
	Pkg.add("MultivariateStats")
	Pkg.add("Plots")
	Pkg.add("Suppressor")
	using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor
end

# ╔═╡ 449fdd2e-170c-4463-b3ef-b34a2251df56
begin
	Pkg.activate(".")
	Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add("ACEpotentials")
	using ACEpotentials
end

# ╔═╡ 5e3d4522-4959-4e6e-bb77-134b0ce23a47
tial_data, _, _ = ACEpotentials.example_dataset("TiAl_tutorial");

# ╔═╡ 0c7fb646-5750-4647-89b0-bee88bf34c97
begin
	r_cut = 6.0
	rdf = ACEpotentials.get_rdf(tial_data, r_cut)
	plt_TiTi = histogram(rdf[(:Ti, :Ti)], bins=100, xlabel = "", c = 1,
	         ylabel = "RDF - TiTi", label = "", yticks = [], xlims = (0, r_cut) )
	plt_TiAl = histogram(rdf[(:Ti, :Al)], bins=100, xlabel = "", c = 2,
	         ylabel = "RDF - TiAl", label = "", yticks = [], xlims = (0, r_cut) )
	plt_AlAl = histogram(rdf[(:Al, :Al)], bins=100, xlabel = L"r [\AA]", c = 3,
	         ylabel = "RDF - AlAl", label = "", yticks = [], xlims = (0, r_cut), )
	plot(plt_TiTi, plt_TiAl, plt_AlAl, layout = (3,1), size = (500, 500), left_margin = 6Plots.mm)
end

# ╔═╡ 0b5cee76-1772-49a4-ab4e-35d7e4698455
begin
	model = acemodel(elements = [:Ti, :Al],
	                 order = 3,
	                 totaldegree = 6,
	                 rcut = 5.5,
	                 Eref = [:Ti => -1586.0195, :Al => -105.5954])
	@show length(model.basis);
end

# ╔═╡ c56e1ee1-bee6-4fcc-9a15-15b6ff7e32af
begin
	acefit!(model, tial_data[1:5:end]);
	ACEpotentials.linear_errors(tial_data[1:5:end], model);
end

# ╔═╡ Cell order:
# ╠═81f97bc0-8546-11ee-3ef6-414aa5615b41
# ╠═449fdd2e-170c-4463-b3ef-b34a2251df56
# ╠═5e3d4522-4959-4e6e-bb77-134b0ce23a47
# ╠═0c7fb646-5750-4647-89b0-bee88bf34c97
# ╠═0b5cee76-1772-49a4-ab4e-35d7e4698455
# ╠═c56e1ee1-bee6-4fcc-9a15-15b6ff7e32af
