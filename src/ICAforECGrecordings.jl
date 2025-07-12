module ICAforECGrecordings

using LinearAlgebra: transpose, sqrt, I, norm, svd, Diagonal, pinv, eigen, Symmetric, dot, inv
using DelimitedFiles: readdlm
using Statistics: norm, mean
using Plots: plot, plot!, xlabel!


# Write your package code here.
include("Preprocessing.jl")
include("Visualization.jl")
include("Parser.jl")
include("Utils.jl")

# Algos
include("Shibbs.jl")
include("Jade.jl")


abstract type AbstractSeperator end

struct JadeSeperator <: AbstractSeperator end
struct ShibbsSeperator <: AbstractSeperator end
struct PicardoSeperator <: AbstractSeperator end


"""
    solve(seperator::AbstractSeperator, data::AbstractMatrix) -> AbstractMatrix
Performs Independent Component Analysis (ICA) using the specified algorithm on the provided data.
# Arguments
- `seperator`: An instance of `AbstractSeperator`, which can be `JadeSeperator`, `ShibbsSeperator`, or `PicardoSeperator`.
- `data`: A matrix of size `(n_{samples}, n_{signals} + 1)` where the first column is time in seconds and `data[:, 2:end]` contains the signal measurements.
# Returns
- A matrix of size `(n_{samples}, n_{signals} + 1)` where the first column is time and `[:, 2:end]` are the separated signals.  
"""
solve(seperator::JadeSeperator, data::AbstractMatrix) = 
begin 
    data_w, W_white = whiten(data)
    signals = jade(data_w, W_white)
    return signals = signals[:, 1:3] 
end 


"""
    solve(seperator::ShibbsSeperator, data::AbstractMatrix) -> AbstractMatrix

Performs Independent Component Analysis (ICA) using the Shibbs algorithm on the provided data.

# Arguments
- `seperator`: An instance of `ShibbsSeperator`.
- `data`: A matrix of size `(n_{samples}, n_{signals} + 1)` where the first column is time in seconds and `data[:, 2:end]` contains the signal measurements.

# Returns
- A matrix of size `(n_{samples}, n_{signals} + 1)` where the first column is time and `[:, 2:end]` are the separated signals.
"""
solve(seperator::ShibbsSeperator, data::AbstractMatrix) = 
begin 
    t = data[:,1]
    X, _ = whiten(data) 
    B = shibbs(X,2) 
    signals = X[:, 2:end] 
    X = Matrix(signals') 
    S = B * X
    S = Matrix(S') 
    S = hcat(t,S)  
    return S
end 

solve(seperator::PicardoSeperator, data::AbstractMatrix) = 
begin 
    data = white(data)
    signal1 = picardo(data)
    signal2 = hcat(data[:,1], data[:,2:end]-signal1[:,2:end])
    return signal1, signal2
end 





@doc """
    load_example_data()
Load example data from a .dat file for testing purposes.
Returns a matrix containing the data.
"""
function load_example_data()
    datpath = joinpath(@__DIR__, "..", "data", "foetal_ecg.dat")
    data = read_dataset_from_dat(datpath)
    return data
end

export whiten, plot_dataset, read_dataset_from_dat, shibbs, jade, load_example_data, solve, PicardoSeperator, 
       JadeSeperator, ShibbsSeperator, AbstractSeperator

end