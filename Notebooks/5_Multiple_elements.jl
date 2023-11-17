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
	         ylabel = "RDF - TiTi", label = "rdf", yticks = [], xlims = (0, r_cut) )
	vline!(rnn(:Ti)*[1.0, 1.333, 1.633, 1.915, 2.3, 2.5], label = "r1, r2, ...", lw=5, color = "black")
	
	plt_TiAl = histogram(rdf[(:Ti, :Al)], bins=100, xlabel = "", c = 2,
	         ylabel = "RDF - TiAl", label = "rdf", yticks = [], xlims = (0, r_cut) )
	vline!((rnn(:Al)+rnn(:Ti))/2*[1.0, 1.633, 2.3, 2.5], label = "r1, r2, ...", lw=5, color = "black")
	
	plt_AlAl = histogram(rdf[(:Al, :Al)], bins=100, xlabel = L"r [\AA]", c = 3,
	         ylabel = "RDF - AlAl", label = "rdf", yticks = [], xlims = (0, r_cut), )
	vline!(rnn(:Al)*[1.0, 1.333, 1.633, 1.915, 2.3, 2.5], label = "r1, r2, ...", lw=5, color = "black" )

	
	plot(plt_TiTi, plt_TiAl, plt_AlAl, layout = (3,1), size = (600, 500), left_margin = 6Plots.mm)
end

# ╔═╡ 9276daa7-113a-4c06-aaa9-9787b7d71bfc
begin
	r_cut_adf_Ti = 1.25 * rnn(:Ti) #
	r_cut_adf_Al = 1.25 * rnn(:Al) #
	r_cut_adf_TiAl = 1.25 * (rnn(:Al)+rnn(:Ti))/2 #
	
	eq_angle = [60, 90, 120]* pi/180 # radians
	
	adf_tiny_Ti = ACEpotentials.get_adf(tial_data, r_cut_adf_Ti);
	plt_adf_1 = histogram(adf_tiny_Ti, bins=50, label = "adf", yticks = [], c = 3,
	                    title = "Ti", titlefontsize = 10,
	                    xlabel = L"\theta", ylabel = "ADF",
	                    xlims = (0, π), size=(400,200), left_margin = 2Plots.mm)
	vline!([ eq_angle,], label = "60˚, 90˚, 120˚", lw=5, color = "black")


	adf_tiny_Al = ACEpotentials.get_adf(tial_data, r_cut_adf_Al);
	plt_adf_2= histogram(adf_tiny_Al, bins=50, label = "adf", yticks = [], c = 3,
	                    title = "Al", titlefontsize = 10,
	                    xlabel = L"\theta", ylabel = "ADF",
	                    xlims = (0, π), size=(400,200), left_margin = 2Plots.mm)
	vline!([ eq_angle,], label = "60˚, 90˚, 120˚", lw=5, color = "black")
	
	adf_tiny_TiAl = ACEpotentials.get_adf(tial_data, r_cut_adf_TiAl);
	plt_adf_3= histogram(adf_tiny_TiAl, bins=50, label = "adf", yticks = [], c = 3,
	                    title = "TiAl", titlefontsize = 10,
	                    xlabel = L"\theta", ylabel = "ADF",
	                    xlims = (0, π), size=(400,200), left_margin = 2Plots.mm)
	vline!([ eq_angle,], label = "60˚, 90˚, 120˚", lw=5, color = "black")

	plot(plt_adf_1, plt_adf_2, plt_adf_3, layout=(3,1), size=(600,600))
end

# ╔═╡ 88a9dbf0-ee53-4b04-b36f-598e9c124f02
60* pi/180

# ╔═╡ 6e9ee5ba-95ac-4588-875d-7b26799424f3
109.5* pi/180

# ╔═╡ 0b5cee76-1772-49a4-ab4e-35d7e4698455
begin
	model = acemodel(elements = [:Ti, :Al],
	                 order = 3,
	                 totaldegree = 6,
	                 rcut = 5.5,
	                 Eref = [:Ti => -1586.0195, :Al => -105.5954])
	@show length(model.basis);
end

# ╔═╡ 29e644d3-0d44-45b6-b728-87c45bc50c88
config_types = [at.data["config_type"].data for at in tial_data]

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
# ╠═9276daa7-113a-4c06-aaa9-9787b7d71bfc
# ╠═88a9dbf0-ee53-4b04-b36f-598e9c124f02
# ╠═6e9ee5ba-95ac-4588-875d-7b26799424f3
# ╠═0b5cee76-1772-49a4-ab4e-35d7e4698455
# ╠═29e644d3-0d44-45b6-b728-87c45bc50c88
# ╠═c56e1ee1-bee6-4fcc-9a15-15b6ff7e32af
