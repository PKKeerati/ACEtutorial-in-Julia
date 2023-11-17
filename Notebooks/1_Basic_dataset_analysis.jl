### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ da3158ce-84dd-11ee-0e81-3516ae8782d9
begin
	using Pkg
	Pkg.activate(".")
	Pkg.add("LaTeXStrings")
	Pkg.add("MultivariateStats")
	Pkg.add("Plots")
	Pkg.add("Suppressor")
	using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor
end

# ╔═╡ 1b5cb9ff-8487-4e75-b17b-3a373bb8a1b0
begin
	Pkg.activate(".")
	Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add("ACEpotentials")
	using ACEpotentials
end

# ╔═╡ 95379cd9-37fa-4450-ab28-55cf4ac92e7c
md"""
##### Introduction
The ACEpotentials.jl documentation contains a number of short, focused tutorials on key topics. This tutorial is longer and has a single narrative. While it is not specifically intended as a Julia-language primer, many Julia commands are introduced by example.
"""

# ╔═╡ a2f772d5-7735-4b3d-a8ff-9f31b1eb3337
ACEpotentials.list_example_datasets()

# ╔═╡ 9a13a2f0-94dd-44d7-89c7-5080460ffa3c
Si_tiny_dataset, _, _ = ACEpotentials.example_dataset("Si_tiny");

# ╔═╡ bf1b81ac-97eb-46ee-abfb-f0b19cad2a3c
begin
	download("https://www.dropbox.com/scl/fi/mzd7zcb1x1l4rw5eswxcd/gp_iter6_sparse9k.xml.xyz?rlkey=o4avtpkka6jnqn7qg375vg7z0&dl=0",
	         "Si_dataset.xyz");
	
	Si_dataset = read_extxyz("Si_dataset.xyz");
end

# ╔═╡ d61552f2-bf22-4fca-8622-8a001698aa41
begin
	println("The tiny dataset has ", length(Si_tiny_dataset), " structures.")
	println("The large dataset has ", length(Si_dataset), " structures.")
end

# ╔═╡ 8ee22d81-db9e-4985-a0fe-bf2f4184e10d
begin
	config_types_tiny = [at.data["config_type"].data for at in Si_tiny_dataset]
	config_types = [at.data["config_type"].data for at in Si_dataset]
	
	function count_configs(config_types)
	    config_counts = [sum(config_types.==ct) for ct in unique(config_types)]
	    config_dict = Dict([ct=>cc for (ct,cc) in zip(unique(config_types), config_counts)])
	end;
end

# ╔═╡ 8a33d459-6fd1-44e4-986f-e51b4c852fba
begin
	println("There are ", length(unique(config_types_tiny)), " unique config_types "*
	        "in the tiny dataset:")
	display(count_configs(config_types_tiny))
end

# ╔═╡ cb3f1c8d-1042-43bf-a357-9562b140f743
begin
	println("There are ", length(unique(config_types)), " unique config_types "*
	        "in the full dataset:")
	display(count_configs(config_types))
end

# ╔═╡ ab88d602-778a-44a4-9273-9a9f08865626
rnn(:Si)

# ╔═╡ 0b79e0b5-3a83-4431-ae9a-81950c457c60
begin
	r_cut = 6.0
	
	rdf_tiny = ACEpotentials.get_rdf(Si_tiny_dataset, r_cut; rescale = true)
	plt_rdf_1 = histogram(rdf_tiny[(:Si, :Si)], bins=150, label = "rdf",
	                      title="Si_tiny_dataset", titlefontsize=10,
	                      xlabel = L"r[\AA]", ylabel = "RDF", yticks = [],
	                      xlims=(1.5,6), size=(400,200), left_margin = 2Plots.mm)
	vline!(rnn(:Si)*[1.0, 1.633, 1.915, 2.3, 2.5], label = "r1, r2, ...", lw=3)
	
	rdf = ACEpotentials.get_rdf(Si_dataset, r_cut; rescale = true);
	plt_rdf_2 = histogram(rdf[(:Si, :Si)], bins=150, label = "rdf",
	                      title="Si_dataset", titlefontsize=10,
	                      xlabel = L"r[\AA]", ylabel = "RDF", yticks = [],
	                      xlims=(1.5,6), size=(400,200), left_margin = 2Plots.mm)
	vline!(rnn(:Si)*[1.0, 1.633, 1.915, 2.3, 2.5], label = "r1, r2, ...", lw=3)
	
	plot(plt_rdf_1, plt_rdf_2, layout=(2,1), size=(400,400))
end

# ╔═╡ f3cba55e-19ee-445f-8caf-167483f3bdc0
begin
	r_cut_adf = 1.25 * rnn(:Si)
	eq_angle = 1.91 # radians
	adf_tiny = ACEpotentials.get_adf(Si_tiny_dataset, r_cut_adf);
	plt_adf_1 = histogram(adf_tiny, bins=50, label = "adf", yticks = [], c = 3,
	                    title = "Si_tiny_dataset", titlefontsize = 10,
	                    xlabel = L"\theta", ylabel = "ADF",
	                    xlims = (0, π), size=(400,200), left_margin = 2Plots.mm)
	vline!([ eq_angle,], label = "109.5˚", lw=3)
	
	adf = ACEpotentials.get_adf(Si_dataset, r_cut_adf);
	plt_adf_2 = histogram(adf, bins=50, label = "adf", yticks = [], c = 3,
	                    title = "Si_dataset", titlefontsize = 10,
	                    xlabel = L"\theta", ylabel = "ADF",
	                    xlims = (0, π), size=(400,200), left_margin = 2Plots.mm)
	vline!([ eq_angle,], label = "109.5˚", lw=3)
	
	plot(plt_adf_1, plt_adf_2, layout=(2,1), size=(400,400))
end

# ╔═╡ 5b4944fd-411a-45ec-a97b-8abc38febf50
begin
	function extract_energies(dataset)
	    energies = []
	    for atoms in dataset
	        for key in keys(atoms.data)
	            if lowercase(key) == "dft_energy"
	                push!(energies, atoms.data[key].data/length(atoms))
	            end
	        end
	    end
	    return energies
	end;
	
	Si_dataset_energies = extract_energies(Si_dataset)
	
	GC.gc()
end

# ╔═╡ Cell order:
# ╠═95379cd9-37fa-4450-ab28-55cf4ac92e7c
# ╠═da3158ce-84dd-11ee-0e81-3516ae8782d9
# ╠═1b5cb9ff-8487-4e75-b17b-3a373bb8a1b0
# ╠═a2f772d5-7735-4b3d-a8ff-9f31b1eb3337
# ╠═9a13a2f0-94dd-44d7-89c7-5080460ffa3c
# ╠═bf1b81ac-97eb-46ee-abfb-f0b19cad2a3c
# ╠═d61552f2-bf22-4fca-8622-8a001698aa41
# ╠═8ee22d81-db9e-4985-a0fe-bf2f4184e10d
# ╠═8a33d459-6fd1-44e4-986f-e51b4c852fba
# ╠═cb3f1c8d-1042-43bf-a357-9562b140f743
# ╠═ab88d602-778a-44a4-9273-9a9f08865626
# ╠═0b79e0b5-3a83-4431-ae9a-81950c457c60
# ╠═f3cba55e-19ee-445f-8caf-167483f3bdc0
# ╠═5b4944fd-411a-45ec-a97b-8abc38febf50
