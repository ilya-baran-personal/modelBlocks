# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import BlackBoxOptim
import NearestNeighbors
import SpecialFunctions
using LinearAlgebra

# Generate plausible population
function generatePPop(blockWithOutputs::AbstractBlock,
    parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}},
    outputsAndStds::AbstractDict{String, Tuple{S, T}}, # Name to value and std
    count::Int;
    kwargs...) where {S, T}
    generatePPop([blockWithOutputs], parametersAndBounds, outputsAndStds, count; kwargs...)
end

# Generate plausible population with multiple block runs
function generatePPop(blocksWithOutputs::Vector{<:AbstractBlock},
                      parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}},
                      outputsAndStds::AbstractDict{String, Tuple{S, T}}, # Name to value and std.
                      count::Int;
                      MaxTime = 1000) where {S, T}
    boundBlocks = [BlockWithBindings(block, [name for (name, lower, upper) in parametersAndBounds]) for block in blocksWithOutputs];
    lower = [lower for (_, lower, _) in parametersAndBounds];
    upper = [upper for (_, _, upper) in parametersAndBounds];

    outputs = [computeOutputs(block) for block in blocksWithOutputs];

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
            blockCopy = deepcopy(boundBlocks[i]);
            setParameters!(blockCopy, x);
            outputs = computeOutputs(blockCopy);
            outputs = unfold(outputs.values - expectedValuesArray[i]) ./ unfold(stdsArray[i]);
            obj += sum(max.(outputs .^ 2, 1) .- 1)
        end
        return obj
    end

    ppop = Matrix{Float64}(undef, length(lower), count);

    Threads.@threads for i in 1:count
        result = BlackBoxOptim.bboptimize(objective; SearchRange = collect(zip(lower, upper)),
                                          MaxTime = Threads.nthreads() * MaxTime / count, TargetFitness = 0.);
        ppop[:, i] = BlackBoxOptim.best_candidate(result);
    end

    # Print results
    outputResults = Vector(undef, count);
    Threads.@threads for i in 1:count
        allOutputs = "";
        for j in 1:length(boundBlocks)
            blockCopy = deepcopy(boundBlocks[j]);
            setParameters!(blockCopy, ppop[:, i]);
            local outputs = computeOutputs(blockCopy);
            blockOutputs = "";
            for k = 1:length(expectedValuesArray[j])
                v = outputs.values[k];
                min = expectedValuesArray[j][k] - stdsArray[j][k];
                max = expectedValuesArray[j][k] + stdsArray[j][k];
                if min < v < max
                    blockOutputs *= "   $(outputs.variables[k].name): $v in ($min, $max) \n";
                else
                    blockOutputs *= " * $(outputs.variables[k].name): $v OUT ($min, $max) \n";
                end
            end
            allOutputs = allOutputs * blockOutputs;
        end
        outputResults[i] = allOutputs;
    end

    for i in 1:count
        println("Plausible patient $i");
        println(outputResults[i]);
    end

    ppop;
end

function subsamplePPop(ppop::Matrix, block::AbstractBlock, parameterNames::AbstractArray{String}, mean::Vector, covariance::Matrix, count::Integer)
    if length(parameterNames) != 0
        block = BlockWithBindings(block, parameterNames);
    end
    outputs = Array{Variables}(undef, size(ppop, 2));
    for i in 1:length(outputs)
        setParameters!(block, ppop[:, i]);
        outputs[i] = computeOutputs(block);
    end
    outputsMatrix = [outputs[j].values[i] for i=1:length(outputs[1].values), j=1:length(outputs)];
    return subsamplePPop(outputsMatrix, mean, covariance, count);
end

const NEIGHBORS = 5;

# outputs are columns of outputs.  Returns a vector of indices.
function subsamplePPop(outputs::Matrix, mean::Vector, covariance::Matrix, count::Integer)
    tree = NearestNeighbors.KDTree(outputs; leafsize = 10);
    idxs, _ = NearestNeighbors.knn(tree, outputs, NEIGHBORS + 1, true);
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
