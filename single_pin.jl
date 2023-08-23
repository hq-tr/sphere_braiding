include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("/home/trung/_qhe-julia/Potentials.jl")

using .FQH_states
using .Potentials
using LinearAlgebra
using SparseArrays

function main()
    println("Input root file name:")
    fname = readline()
    roots = readwf(fname)
    jack_list = Vector{FQH_state}()

    for root in roots.basis
        rootstring = prod(string.((Int.(root))))
        print("\r$(rootstring)\t")
        jack = sphere_normalize(readwf("jacks/J_$(rootstring)"))
        push!(jack_list, jack)
    end

    Ne = count(jack_list[1].basis[1])
    No = length(jack_list[1].basis[1])
    println("\n\n$(Ne) electrons and $(No) orbitals\n")


    all_basis, all_coef = collate_many_vectors(jack_list; separate_out=true, collumn_vector=true)

    println("Orthonormalizing basis using QR decomposition")
    @time all_coef_ortho = collect(transpose(Matrix(qr(all_coef).Q)))
    println("------")

    north_pin = sphere_point_matrix(all_basis, 0.0, 0.0)
    mat_MR = conj.(all_coef_ortho) * north_pin * transpose(all_coef_ortho)

    println("Executing ED...   ")
    @time ED = eigen(mat_MR)

    println()

    

    println("Eigenvalues = ")
    if length(ED.values) > 3
        display(ED.values[1:5])
        println("\t.\n\t.\n\t.")
    else
        display(ED.values)
    end

    gap = abs(ED.values[2])-abs(ED.values[1])
    println("Energy gap = $(gap)")


    vecs = ED.vectors[:,sortperm(abs.(ED.values))] # Sort by increasing abs(E)


    gs_coef = transpose(all_coef_ortho) * vecs[:,1]
    gs = FQH_state(all_basis, gs_coef)
    printwf(gs;fname="ground_states/$(Ne)e$(No)_gs_1pin_0")

    
    ex1_coef = transpose(all_coef_ortho) * vecs[:,2]
    ex1 = FQH_state(all_basis, ex1_coef) # 1st excited state
    printwf(ex1;fname="ground_states/$(Ne)e$(No)_gs_1pin_1")

end

@time main()