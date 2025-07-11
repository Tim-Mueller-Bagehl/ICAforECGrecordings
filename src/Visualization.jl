"""
    plot_dataset(data::AbstractMatrix{Float64})
Plot a dataset with time on the x-axis and multiple signals on the y-axis.
# Arguments
- `data::AbstractMatrix{Float64}`: Matrix where the first column represents time and the remaining columns represent different signals.
"""
function plot_dataset(data::AbstractMatrix{Float64})
    time = data[:, 1]
    signals = data[:, 2:end]
    n_signals = size(signals, 2)

    plt = plot(size=(1600, 400), legend=:topright, title="ECG Signals Over Time",
               xlabel="Time", ylabel="Voltage", grid=true, background_color=:white)

    for i in 1:n_signals
        plot!(plt, time, signals[:, i], label="Signal $i", lw=3)
    end


    return plt
end