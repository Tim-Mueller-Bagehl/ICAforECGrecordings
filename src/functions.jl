using Plots


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