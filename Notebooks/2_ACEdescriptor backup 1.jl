### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ e5af15c0-84e8-11ee-11a4-7be33f619d04
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

# ╔═╡ 252323d6-ca55-44aa-94e3-dd232c780b83
md"""
### Installing ACEpotentials
"""

# ╔═╡ f2ee641e-76af-4267-adc6-f554974bf0bd
basis = ACE1x.ace_basis(elements = [:Si],
                        rcut = 5.5,
                        order = 3,        # body-order - 1
                        totaldegree = 8);

# ╔═╡ 13a20f55-3a1a-4a5a-80e9-acdb8b7fae8d
md"""
### Downloading dataset
"""

# ╔═╡ eaa670b3-2058-4a72-8fd9-15ccc3ff8825
begin
	Si_tiny_dataset, _, _ = ACEpotentials.example_dataset("Si_tiny");
	
	download("https://www.dropbox.com/scl/fi/mzd7zcb1x1l4rw5eswxcd/gp_iter6_sparse9k.xml.xyz?rlkey=o4avtpkka6jnqn7qg375vg7z0&dl=0",
	         "Si_dataset.xyz");
	
	Si_dataset = read_extxyz("Si_dataset.xyz");
	
	config_types_tiny = [at.data["config_type"].data for at in Si_tiny_dataset]
	config_types = [at.data["config_type"].data for at in Si_dataset]
end;

# ╔═╡ a313e117-a7fa-4e2a-afce-be192484931a
begin
	config_types_tiny = [at.data["config_type"].data for at in Si_tiny_dataset]
	config_types = [at.data["config_type"].data for at in Si_dataset]
end

# ╔═╡ 3225d9a5-db98-41fb-96fc-44be39b74f7b
begin
	descriptorsTiny = []
	for atoms in Si_tiny_dataset
	    struct_descriptor = sum(site_descriptors(basis, atoms)) / length(atoms)
		#print(length(atoms))
	    push!(descriptorsTiny, struct_descriptor)
	end
end

# ╔═╡ 73cfabfb-1c9d-45b7-a03b-03231dd5e938
begin
	descriptorsM = hcat(descriptors...)  # convert to matrix
	M = fit(PCA, descriptorsM; maxoutdim=3, pratio=1)
	descriptors_trans = transform(M, descriptorsM)
	p = scatter(
	     descriptors_trans[1,:], descriptors_trans[2,:], descriptors_trans[3,:],
	     marker=:circle, linewidth=0, group=config_types_tiny, legend=:right)
	plot!(p, xlabel="PC1", ylabel="PC2", zlabel="PC3", camera=(20,10))
end

# ╔═╡ 24e3e640-6dc3-493e-b144-d034a75617f7
begin
	descriptors2 = []
	for atoms in Si_dataset
	    struct_descriptor = sum(site_descriptors(basis, atoms)) / length(atoms)
	    push!(descriptors2, struct_descriptor)
	end
	
	descriptors2M = hcat(descriptors2...)  # convert to matrix
	M2 = fit(PCA, descriptors2M; maxoutdim=3, pratio=1)
	descriptors_trans2 = transform(M2, descriptors2M)
	p2 = scatter(
	     descriptors_trans2[1,:], descriptors_trans2[2,:], descriptors_trans2[3,:],
	     marker=:circle, linewidth=0, group=config_types, legend=:right)
	plot!(p2, xlabel="PC1", ylabel="PC2", zlabel="PC3", camera=(10,10))
end

# ╔═╡ 6798ab9c-e280-402b-b303-8ad6caf97a0e
GC.gc()

# ╔═╡ Cell order:
# ╠═252323d6-ca55-44aa-94e3-dd232c780b83
# ╠═e5af15c0-84e8-11ee-11a4-7be33f619d04
# ╠═f2ee641e-76af-4267-adc6-f554974bf0bd
# ╠═13a20f55-3a1a-4a5a-80e9-acdb8b7fae8d
# ╠═eaa670b3-2058-4a72-8fd9-15ccc3ff8825
# ╠═8781b613-0041-44ba-964f-9c207c1bd82c
# ╠═3225d9a5-db98-41fb-96fc-44be39b74f7b
# ╠═73cfabfb-1c9d-45b7-a03b-03231dd5e938
# ╠═24e3e640-6dc3-493e-b144-d034a75617f7
# ╠═6798ab9c-e280-402b-b303-8ad6caf97a0e
