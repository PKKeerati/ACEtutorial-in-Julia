### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 9f0c04f0-8535-11ee-2a38-b553837a87a0
begin
	using Pkg
	Pkg.activate(".")
	Pkg.add("LaTeXStrings")
	Pkg.add("MultivariateStats")
	Pkg.add("Plots")
	Pkg.add("Suppressor")
	using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor
end

# ╔═╡ 82a9fcc1-735e-48ec-9b69-d2d7c4ab8d92
begin
	Pkg.activate(".")
	Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add("ACEpotentials")
	using ACEpotentials
end

# ╔═╡ 0d64c44f-8acb-412e-a97c-443887971ab9
begin
	model = acemodel(elements = [:Si,],
	                 order = 3,
	                 totaldegree = 8,
	                 rcut = 5.0,
	                 Eref = [:Si => -158.54496821])
	@show length(model.basis);

	#model = acemodel(elements = [:Ti, :Al],
	#                 order = 3,
	#                 totaldegree = 12,
	#                 rcut = 5.5,
	#                 Eref = [:Ti => -1586.0195, :Al => -105.5954])
	#@show length(model.basis);
end

# ╔═╡ 8782a242-7b79-47a0-a08c-81159e3e8310
begin
	Si_tiny_dataset, _, _ = ACEpotentials.example_dataset("Si_tiny");
	config_types_tiny = [at.data["config_type"].data for at in Si_tiny_dataset]
	
	
	download("https://www.dropbox.com/scl/fi/mzd7zcb1x1l4rw5eswxcd/gp_iter6_sparse9k.xml.xyz?rlkey=o4avtpkka6jnqn7qg375vg7z0&dl=0",
	         "Si_dataset.xyz");
	
	Si_dataset = read_extxyz("Si_dataset.xyz");
	config_types = [at.data["config_type"].data for at in Si_dataset]
end

# ╔═╡ 82b8e75e-795e-43e1-a951-2785be68f212
begin
	solver = ACEfit.QR(lambda=1e-1)
	data_keys = (energy_key = "dft_energy", force_key = "dft_force", virial_key = "dft_virial")
	acefit!(model, Si_tiny_dataset;
	        solver=solver, data_keys...);
end

# ╔═╡ 0d6b1720-e89f-4112-812d-f7c4a246ecfa
acefit!(model, Si_dataset;
	        solver=solver, data_keys...);

# ╔═╡ 86011b92-5e5f-46a6-a425-e2df1d554728
begin
	@info("Training Errors")
	ACEpotentials.linear_errors(Si_tiny_dataset, model; data_keys...);
	
	@info("Test Error")
	ACEpotentials.linear_errors(Si_dataset, model; data_keys...);
end

# ╔═╡ c97410c0-6f63-4929-b74c-0737a48ec3fe
begin
	export2json("model.json", model)
	export2lammps("model.yace", model)
end

# ╔═╡ Cell order:
# ╠═9f0c04f0-8535-11ee-2a38-b553837a87a0
# ╠═82a9fcc1-735e-48ec-9b69-d2d7c4ab8d92
# ╠═0d64c44f-8acb-412e-a97c-443887971ab9
# ╠═8782a242-7b79-47a0-a08c-81159e3e8310
# ╠═82b8e75e-795e-43e1-a951-2785be68f212
# ╠═0d6b1720-e89f-4112-812d-f7c4a246ecfa
# ╠═86011b92-5e5f-46a6-a425-e2df1d554728
# ╠═c97410c0-6f63-4929-b74c-0737a48ec3fe
