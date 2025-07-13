using Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()
Pkg.add("Plots")
import ICAforECGrecordings
using ICAforECGrecordings: whiten, plot_dataset, read_dataset_from_dat, solve, JadeSeperator, ShibbsSeperator, PicardoSeperator
using Plots: savefig


# load data
data = read_dataset_from_dat("./data/foetal_ecg.dat")
println("Data loaded successfully.")

# create plots dir
try 
    mkdir("plots")
catch e
    println(e)
end


# save plots 
savefig(plot_dataset(data), "plots/original.png")
savefig(plot_dataset(solve(ShibbsSeperator(), data)), "plots/shibbs.png")
savefig(plot_dataset(solve(JadeSeperator(), data)), "plots/jade.png")
savefig(plot_dataset(solve(PicardoSeperator(), data)), "plots/picard.png")