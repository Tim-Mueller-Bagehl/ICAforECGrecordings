"""
    plot_dataset(data::AbstractMatrix{Float64})
Plot a dataset with time on the x-axis and multiple signals on the y-axis.
# Arguments
- `data::AbstractMatrix{Float64}`: Matrix where the first column represents time and the remaining columns represent different signals.
"""
function plot_dataset(data::AbstractMatrix{Float64})
    if isempty(data)
        @warn "Input data is empty. Returning nothing."
        return nothing
    end
    time = data[:, 1]
    signals = data[:, 2:end]
    n_signals = size(signals, 2)

    if n_signals > 2
        plt = plot(size=(1600, 400* n_signals), layout = (n_signals, 1), legend=:topright, 
                xlabel="Time", grid=true, ylabel="Signal Value", 
                background_color=:white, link=:x,)
    else
        plt = plot(size=(1600, 400), legend=:topright, 
                xlabel="Time", grid=true, background_color=:white)
    end

    for i in 1:n_signals
        if n_signals > 2
            plot!(plt[i], time, signals[:, i], label="Signal $i", lw=2)
        else
            plot!(plt, time, signals[:, i], label="Signal $i", lw=3)
        end
    end

    return plt
end



