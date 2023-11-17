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

	Pkg.activate(".")
	Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add("ACEpotentials")
	using ACEpotentials
end

# ╔═╡ 3ff5eddb-6436-4a89-8894-7f0a7a841fba
md"""
### Downloading dataset
"""

# ╔═╡ 8782a242-7b79-47a0-a08c-81159e3e8310
begin
	Si_tiny_dataset, _, _ = ACEpotentials.example_dataset("Si_tiny");
	config_types_tiny = [at.data["config_type"].data for at in Si_tiny_dataset]
	
	
	download("https://www.dropbox.com/scl/fi/mzd7zcb1x1l4rw5eswxcd/gp_iter6_sparse9k.xml.xyz?rlkey=o4avtpkka6jnqn7qg375vg7z0&dl=0",
	         "Si_dataset.xyz");
	
	Si_dataset = read_extxyz("Si_dataset.xyz");
	config_types = [at.data["config_type"].data for at in Si_dataset]
end;

# ╔═╡ d2ccaf17-c762-4a5f-a3f7-594a7b6a2230
md"""
## Part 3: Basic model fitting
We begin by defining an (extremely simple) ```acemodel```. The parameters have the same meaning as for ```ace_basis``` in the previous notebook, with an additional Eref providing a reference energy.
"""

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

# ╔═╡ c09391c7-3972-4d3a-b575-27444442de33
md"""
Next, we fit determine the model parameters using the tiny dataset and ridge regression via the ```QR``` solver which is one of the four different solvers as following
* ```QR``` (QR decomposition)
* ```LSQR``` (Krylov Method)
* ```RRQR``` (Rank-revealing QR decomposition)
* ```BLR``` (Bayesian Linear Regression).
More detail illustrated in Table1 in ```ACEpotential.jl``` paper
"""

# ╔═╡ 82b8e75e-795e-43e1-a951-2785be68f212
begin
	solver = ACEfit.QR(lambda=1e-1)
	data_keys = (energy_key = "dft_energy", force_key = "dft_force", virial_key = "dft_virial")
	acefit!(model, Si_tiny_dataset;
	        solver=solver, data_keys...);
end

# ╔═╡ 86011b92-5e5f-46a6-a425-e2df1d554728
begin
	@info("Training Errors")
	ACEpotentials.linear_errors(Si_tiny_dataset, model; data_keys...);
	
	@info("Test Error")
	ACEpotentials.linear_errors(Si_dataset, model; data_keys...);
end

# ╔═╡ 4995a9a9-eb40-4514-8621-973f024ab3c1
md"""
###### A model may be exported to JSON or LAMMPS formats with the following.
"""

# ╔═╡ c97410c0-6f63-4929-b74c-0737a48ec3fe
begin
	export2json("Si_model.json", model)
	export2lammps("Si_model.yace", model)
end

# ╔═╡ 7e7b7a3e-90cf-4bdc-871a-598391ff8adf
potential = load_potential("Si_model.yace")

# ╔═╡ Cell order:
# ╠═9f0c04f0-8535-11ee-2a38-b553837a87a0
# ╠═3ff5eddb-6436-4a89-8894-7f0a7a841fba
# ╠═8782a242-7b79-47a0-a08c-81159e3e8310
# ╠═d2ccaf17-c762-4a5f-a3f7-594a7b6a2230
# ╠═0d64c44f-8acb-412e-a97c-443887971ab9
# ╠═c09391c7-3972-4d3a-b575-27444442de33
# ╠═82b8e75e-795e-43e1-a951-2785be68f212
# ╠═86011b92-5e5f-46a6-a425-e2df1d554728
# ╠═4995a9a9-eb40-4514-8621-973f024ab3c1
# ╠═c97410c0-6f63-4929-b74c-0737a48ec3fe
# ╠═7e7b7a3e-90cf-4bdc-871a-598391ff8adf
