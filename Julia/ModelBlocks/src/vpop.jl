# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import BlackBoxOptim
import NearestNeighbors
import SpecialFunctions
using LinearAlgebra

function generatePPop(blockWithOutputs::AbstractBlock,
    timeRange::AbstractRange,
    parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}},
    outputsAndStds::AbstractDict{String, Tuple{Float64, Float64}}, # Name to value and std
    count::Int;
    kwargs...)
    generatePPop([blockWithOutputs], timeRange, parametersAndBounds, outputsAndStds, count; kwargs...)
end

function generatePPop(blocksWithOutputs::Vector{T},
                      timeRange::AbstractRange,
                      parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}},
                      outputsAndStds::AbstractDict{String, Tuple{Float64, Float64}}, # Indexed to match blocks.  Name to value and std.
                      count::Int;
                      MaxTime = 1000) where T <: AbstractBlock
    boundBlocks = [BlockWithBindings(block, [name for (name, lower, upper) in parametersAndBounds]) for block in blocksWithOutputs];
    lower = [lower for (name, lower, upper) in parametersAndBounds];
    upper = [upper for (name, lower, upper) in parametersAndBounds];

    outputs = [getOutputs(block, timeRange) for block in blocksWithOutputs];

    expectedValuesArray = [];
    stdsArray = [];
    for output in outputs
        expectedValues = [];
        stds = [];
        for variable in output.variables
            push!(expectedValues, outputsAndStds[variable.name][1]);
            push!(stds, outputsAndStds[variable.name][2]);
        end
        push!(expectedValuesArray, expectedValues);
        push!(stdsArray, stds);
    end

    objective = x -> begin
        obj = 0.;
        for i in 1:length(boundBlocks)
            setParameters!(boundBlocks[i], x);
            outputs = getOutputs(boundBlocks[i], timeRange);
            outputs = (outputs.values - expectedValuesArray[i]) ./ stdsArray[i];
            obj += sum(max.(outputs .^ 2, 1) .- 1)
        end
        return obj
    end

    vpop = Matrix{Float64}(undef, length(lower), count);

    for i in 1:count
        result = BlackBoxOptim.bboptimize(objective; SearchRange = collect(zip(lower, upper)), MaxTime = MaxTime / count);
        vpop[:, i] = BlackBoxOptim.best_candidate(result);
    end
    vpop;
end

const NEIGHBORS = 5;

# outputs are columns of outputs.  Returns a vector of indices.
function subsamplePPop(outputs::Matrix, mean::Vector, covariance::Matrix, count::Integer)
    tree = NearestNeighbors.KDTree(outputs; leafsize = 10);
    idxs, dists = NearestNeighbors.knn(tree, outputs, NEIGHBORS + 1, true);
    radii = [norm(outputs[idx[end]] - outputs[idx[1]]) for idx in idxs];
    densities = NEIGHBORS ./ nBallVolume(radii, size(outputs, 1));
    scaledProbabilities = normalPDF(outputs, mean, covariance) ./ densities;
    scaledProbabilities = scaledProbabilities / sum(scaledProbabilities);
    return weightedSample(scaledProbabilities, count);
end

function weightedSample(weights::Vector, count::Integer)
    fAndIdx = [(weights[i] < 1e-14 ? -1e100 : log(rand()) / weights[i], i) for i in 1:length(weights)];
    fAndIdx = sort(fAndIdx, rev = true);
    count = min(count, length(weights));
    return sort([s[2] for s in fAndIdx[1:count]]);
end

function normalPDF(points::Matrix, mean::Vector, covariance::Matrix)
    points = points .- mean;
    result = exp.(-0.5 * sum(points .* (inv(covariance) * points), dims=1)[1,:]) / sqrt((2 * pi) ^ length(mean) * abs(det(covariance)));
    return result;
end

function nBallVolume(radius, d::Number)
    (pi ^ (d / 2) / SpecialFunctions.gamma(d / 2 + 1)) * radius .^ d 
end
