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

# ╔═╡ b4423e07-1b6e-4a77-ae50-d1f2ecd70d1a
begin
	Si_tiny_dataset, _, _ = ACEpotentials.example_dataset("Si_tiny");
	config_types_tiny = [at.data["config_type"].data for at in Si_tiny_dataset]
	
	
	download("https://www.dropbox.com/scl/fi/mzd7zcb1x1l4rw5eswxcd/gp_iter6_sparse9k.xml.xyz?rlkey=o4avtpkka6jnqn7qg375vg7z0&dl=0",
	         "Si_dataset.xyz");
	
	Si_dataset = read_extxyz("Si_dataset.xyz");
	config_types = [at.data["config_type"].data for at in Si_dataset]
end

# ╔═╡ d6c59890-d1e2-4bb2-b281-db1ca4b667e0
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

# ╔═╡ c2481c27-3bf4-474e-b0fa-4932ddcea13b
function assess_model(model, train_dataset)

    plot([-164,-158], [-164,-158]; lc=:black, label="")

    model_energies = []
    model_std = []
    for atoms in Si_dataset
        ene, co_ene = ACE1.co_energy(model.potential, atoms)
        push!(model_energies, ene/length(atoms))
        push!(model_std, std(co_ene/length(atoms)))
    end
    rmse = sqrt(sum((model_energies-Si_dataset_energies).^2)/length(Si_dataset))
    mae = sum(abs.(model_energies-Si_dataset_energies))/length(Si_dataset)
    scatter!(Si_dataset_energies, model_energies;
             label="full dataset",
             title = @sprintf("Structures Used In Training:  %i out of %i\n", length(train_dataset), length(Si_dataset)) *
                     @sprintf("RMSE (MAE) For Entire Dataset:  %.0f (%.0f) meV/atom", 1000*rmse, 1000*mae),
             titlefontsize = 8,
             yerror = model_std,
             xlabel="Energy [eV/atom]", xlims=(-164,-158),
             ylabel="Model Energy [eV/atom]", ylims=(-164,-158),
             aspect_ratio = :equal, color=1)

    model_energies = [energy(model.potential,atoms)/length(atoms) for atoms in train_dataset]
    scatter!(extract_energies(train_dataset), model_energies;
             label="training set", color=2)

end;

# ╔═╡ 4c8e1eb4-aedc-4344-808f-5cfb6e37d607
assess_model(model, Si_tiny_dataset)

# ╔═╡ adf9c84d-21df-44c1-bb57-f80e8d09fb72
function augment(old_dataset, old_model; num=5)

    new_dataset = deepcopy(old_dataset)
    new_model = deepcopy(old_model)

    model_std = []
    for atoms in Si_dataset
        ene, co_ene = ACE1.co_energy(new_model.potential, atoms)
        push!(model_std, std(co_ene/length(atoms)))
    end
    for atoms in Si_dataset[sortperm(model_std; rev=true)[1:num]]
        push!(new_dataset, atoms)
    end
    @suppress acefit!(new_model, new_dataset;
            solver = ACEfit.BLR(committee_size=50, factorization=:svd),
            energy_key = "dft_energy", force_key = "dft_force",
            verbose = false);

    return new_dataset, new_model
end;

# ╔═╡ 7c387fde-81df-415e-863d-ec585d328900
begin
	new_dataset, new_model = augment(Si_tiny_dataset, model; num=5);
	assess_model(new_model, new_dataset)
end

# ╔═╡ 6a99ca1d-6f23-4ba5-8dc3-f23601b4b67b
begin
	for i in 1:5
	    @show i
	    new_dataset, new_model = augment(new_dataset, new_model; num=5);
	end
	assess_model(new_model, new_dataset)
end

# ╔═╡ 81992eb6-8fa0-473d-aa34-d8b54045b17d
GC.gc()

# ╔═╡ Cell order:
# ╠═15830db7-7ef6-4222-838d-1a86eff4ca85
# ╠═4a5d9b66-6f39-4c6b-856a-1f4cf873ed29
# ╠═b4423e07-1b6e-4a77-ae50-d1f2ecd70d1a
# ╠═d6c59890-d1e2-4bb2-b281-db1ca4b667e0
# ╠═dc0eccf0-8541-11ee-2421-0713757269ef
# ╠═c2481c27-3bf4-474e-b0fa-4932ddcea13b
# ╠═4c8e1eb4-aedc-4344-808f-5cfb6e37d607
# ╠═adf9c84d-21df-44c1-bb57-f80e8d09fb72
# ╠═7c387fde-81df-415e-863d-ec585d328900
# ╠═6a99ca1d-6f23-4ba5-8dc3-f23601b4b67b
# ╠═81992eb6-8fa0-473d-aa34-d8b54045b17d
