# Getting Started Guide

### 1. Read the data

```julia
data = read_dataset_from_dat("./data/foetal_ecg.dat")
```

### 2. Call select the Seperator and call solve()
```julia
shibbs_signals = solve(ShibbsSeperator(), data)
jade_signals = solve(JadeSeperator(), data)
```

### 3. Plot 
```julia
plot_dataset(shibbs_signals)
plot_dataset(jade_signals)
```

----
### Full Example

```julia
using Pkg
Pkg.activate(temp=true)
Pkg.add(url="https://github.com/Tim-Mueller-Bagehl/ICAforECGrecordings")
Pkg.resolve()
Pkg.instantiate()
Pkg.add("Plots")
import ICAforECGrecordings
using ICAforECGrecordings: whiten, plot_dataset, read_dataset_from_dat, solve, JadeSeperator, ShibbsSeperator
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
savefig(plot_dataset(solve(ShibbsSeperator(), data)), "plots/shibbs.png")
savefig(plot_dataset(solve(JadeSeperator(), data)), "plots/jade.png")
savefig(plot_dataset(solve(PicardoSeperator(), data)), "plots/picardo.png")
```