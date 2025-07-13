using Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()
Pkg.add("Plots")
import ICAforECGrecordings
using ICAforECGrecordings: whiten, plot_dataset, read_dataset_from_dat, solve, JadeSeperator, ShibbsSeperator, PicardSeperator
using Plots: savefig


data = read_dataset_from_dat("./data/foetal_ecg.dat")

try 
    mkdir("plots")
catch e
    println(e)
end


# save plots 
savefig(plot_dataset(data), "plots/original.png")
savefig(plot_dataset(solve(ShibbsSeperator(), data)), "plots/shibbs.png")
savefig(plot_dataset(solve(JadeSeperator(), data)), "plots/jade.png")
savefig(plot_dataset(solve(PicardSeperator(), data)), "plots/picard.png")