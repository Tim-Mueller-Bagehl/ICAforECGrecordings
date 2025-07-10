using Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()
import ICAforECGrecordings
using ICAforECGrecordings: whiten, plot_dataset, read_dataset_from_dat, shibbs, jade
using Plots: savefig


# load data
data = read_dataset_from_dat("./data/foetal_ecg.dat")
println("Data loaded successfully.")

# whiten data
whitened_data = whiten(data)
println("Data whitened successfully.")

# do shibbs on it
shibbs_result = shibbs(whitened_data)
println("Shibbs decomposition completed successfully.")

# do jade on it
# jade_result = jade(whitened_data)
# println("Jade decomposition completed successfully.")

# create plots dir
try 
    mkdir("plots")
catch e
    println(e)
end


# save plots 
savefig(plot_dataset(whitened_data), "plots/whiten.png")
savefig(plot_dataset(shibbs_result), "plots/shibbs.png")
# savefig(plot_dataset(jade_result), "plots/jade.png")