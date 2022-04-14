# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import BlackBoxOptim
import NearestNeighbors
import SpecialFunctions
import Statistics
import Distributions
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

    (expectedValuesArray, stdsArray) = makeExpectedValuesAndStdsArrays(outputs, outputsAndStds);

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

# Generate plausible population with multiple block runs
function generatePPopFarthest(blocksWithOutputs::Vector{<:AbstractBlock},
                              parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}},
                              outputsAndStds::AbstractDict{String, Tuple{S, T}}, # Name to value and std.
                              count::Int;
                              MaxTime = 1000, DistanceFactor = 1., threads = 0) where {S, T}
    boundBlocks = [BlockWithBindings(block, [name for (name, lower, upper) in parametersAndBounds]) for block in blocksWithOutputs];
    lower = [lower for (_, lower, _) in parametersAndBounds];
    upper = [upper for (_, _, upper) in parametersAndBounds];

    outputs = [computeOutputs(block) for block in blocksWithOutputs];

    (expectedValuesArray, stdsArray) = makeExpectedValuesAndStdsArrays(outputs, outputsAndStds);

    metricPPopSoFar::Vector{Vector{Float64}} = [];
    ppop = Matrix{Float64}(undef, length(lower), count);

    toMetric = p -> begin
        result = deepcopy(p);
        for i in 1:length(p)
            if lower[i] <= 0
                result[i] = (p[i] - lower[i]) / (upper[i] - lower[i]);
            else
                result[i] = log(p[i]);
            end
        end
        return result;
    end

    if threads == 0
        threads = Threads.nthreads();
    end
    rounds::Integer = ceil(count / threads);

    objective = (x, idx) -> begin
        obj = 0.;
        for i in 1:length(boundBlocks)
            blockCopy = deepcopy(boundBlocks[i]);
            setParameters!(blockCopy, x);
            outputs = computeOutputs(blockCopy);
            outputs = unfold(outputs.values - expectedValuesArray[i]) ./ unfold(stdsArray[i]);
            obj += sum(max.(outputs .^ 2, 1) .- 1);
        end
        if length(metricPPopSoFar) > 0
            m = toMetric(x);
            closest = Inf64;
            for i in 1:length(metricPPopSoFar)
                dist = sum((m - metricPPopSoFar[i]) .^ 2);
                if mod(i + idx, threads) == 1
                    dist = dist * 2;
                end
                closest = min(closest, dist);
            end
            obj -= closest * DistanceFactor;
        end
        return obj;
    end

    for i in 1:rounds
        println("Starting round $i of $rounds");
        batch::Integer = min(threads, count - (i - 1) * threads);
        batchPpop = Matrix{Float64}(undef, length(lower), batch);
        Threads.@threads for j in 1:batch
            result = BlackBoxOptim.bboptimize((x) -> objective(x, j); SearchRange = collect(zip(lower, upper)),
                                              MaxTime = MaxTime / rounds);
            batchPpop[:, j] = BlackBoxOptim.best_candidate(result);
        end
        for j in 1:batch
            push!(metricPPopSoFar, toMetric(batchPpop[:, j]));
            ppop[:, (i - 1) * threads + j] = batchPpop[:, j];
        end
    end

    return ppop;
end

function computeDistanceCurve(blocksWithOutputs::Vector{<:AbstractBlock},
                              parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}},
                              outputsAndStds::AbstractDict{String, Tuple{S, T}}, # Name to value and std.
                              ppop::Matrix{Float64},
                              maxDistance::Float64;
                              samples = 1000) where {S, T}
    boundBlocks = [BlockWithBindings(block, [name for (name, lower, upper) in parametersAndBounds]) for block in blocksWithOutputs];
    lower = [lower for (_, lower, _) in parametersAndBounds];
    upper = [upper for (_, _, upper) in parametersAndBounds];

    outputs = [computeOutputs(block) for block in blocksWithOutputs];

    (expectedValuesArray, stdsArray) = makeExpectedValuesAndStdsArrays(outputs, outputsAndStds);

    radii = 0:(maxDistance / 500):maxDistance;
    result = zeros(length(radii));
    
    n = size(ppop)[2];
    logMin = log.(lower);
    logMax = log.(upper);

    successful = 0;
    lk = ReentrantLock();
    
    Threads.@threads for i = 1:samples
        sample = rand(size(ppop)[1]) .* (logMax - logMin) + logMin;
        # Check outputs
        ok = true;
        for i in 1:length(boundBlocks)
            blockCopy = deepcopy(boundBlocks[i]);
            setParameters!(blockCopy, exp.(sample));
            outputs = computeOutputs(blockCopy);
            outputs = unfold(outputs.values - expectedValuesArray[i]) ./ unfold(stdsArray[i]);
            if any(abs.(outputs) .> 1)
                ok = false;
            end
        end
        if !ok
            continue
        end

        distance = sum((log.(ppop) - repeat(sample, 1, n)) .^ 2, dims = 1);
        closest = sqrt(minimum(distance));
        lock(lk) do 
            result[radii .> closest] .+= 1;
            successful += 1;
        end
    end

    return (radii, result / successful);
end

# Take a population and a set of blocks and output ranges.  Optionally filter the population to just those
# whose outputs are in the ranges.  Generate random new patients by sampling in a normal distribution around
# random patients and keeping those that satisfy the output constraints.  The covariance matrix is the
# covariance matrix of the input ppop, scaled by ratio.  We generate and try count new patients (and a subset of
# those are kept).  The return value includes the original (possibly filtered) ppop.
function expandPPop(blocksWithOutputs::Vector{<:AbstractBlock},
                    parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}},
                    outputsAndStds::AbstractDict{String, Tuple{S, T}}, # Name to value and std.
                    ppop::Matrix{Float64},
                    ratio::Number, # What fraction of the covariance to use
                    count::Int;
                    filter = false) where {S, T}
                    
    boundBlocks = [BlockWithBindings(block, [name for (name, lower, upper) in parametersAndBounds]) for block in blocksWithOutputs];
    lower = [lower for (_, lower, _) in parametersAndBounds];
    upper = [upper for (_, _, upper) in parametersAndBounds];

    outputs = [computeOutputs(block) for block in blocksWithOutputs];

    (expectedValuesArray, stdsArray) = makeExpectedValuesAndStdsArrays(outputs, outputsAndStds);

    check = function(parameters)
        for j in 1:length(boundBlocks)
            blockCopy = deepcopy(boundBlocks[j]);
            setParameters!(blockCopy, parameters);
            outputs = computeOutputs(blockCopy);
            outputs = unfold(outputs.values - expectedValuesArray[j]) ./ unfold(stdsArray[j]);
            if any(abs.(outputs) .> 1)
                return false;
            end
        end
        return true;
    end

    isOk = [true for i = 1:size(ppop, 2)];
    if filter
        Threads.@threads for i in 1:size(ppop, 2)
            isOk[i] = check(ppop[:, i]);
        end
        ppop = ppop[:, isOk];
        println("Kept $(size(ppop, 2)) patients");
    end

    covariance = Statistics.cov(ppop, dims = 2);
    distribution = Distributions.MvNormal(covariance * ratio);

    candidates = Matrix{Float64}(undef, length(lower), count);
    for i in 1:count
        candidates[:, i] = ppop[:, rand(1:size(ppop, 2))] + rand(distribution);
        candidates[:, i] = max.(lower, min.(upper, candidates[:, i]));
    end

    isOk = [true for i = 1:count];
    
    Threads.@threads for i in 1:count
        isOk[i] = check(candidates[:, i]);
    end

    for i in 1:count
        if isOk[i]
            ppop = hcat(ppop, candidates[:, i]);
        end
    end

    return ppop;
end

function makeExpectedValuesAndStdsArrays(outputs::Vector{Variables}, outputsAndStds::AbstractDict{String, Tuple{S, T}}) where {S, T}
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
    return (expectedValuesArray, stdsArray);
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
