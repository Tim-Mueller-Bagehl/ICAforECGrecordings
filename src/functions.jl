using Plots

"""
    plot_dataset(data::Matrix{Float64})
Plot a dataset with time on the x-axis and multiple signals on the y-axis.
# Arguments
- `data::Matrix{Float64}`: Matrix where the first column represents time and the remaining columns represent different signals.
"""
function plot_dataset(data::Matrix{Float64})
    time = data[:, 1]
    signals = data[:, 2:end]
    n_signals = size(signals, 2)

    plt = plot(layout=(n_signals, 1), link=:x, size=(800, 200 * n_signals))

    for i in 1:n_signals
        plot!(plt[i], time, signals[:, i], label="Signal $i", ylabel="Value")
        if i == n_signals
            xlabel!(plt[i], "Time")
        end
    end

    savefig(plt, "dataset_plot.png")
    return plt
end