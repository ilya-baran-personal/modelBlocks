# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using ModelBlocks
using Test
using LinearAlgebra
using Interpolations
using Plots


# Part 1: define model - variables, parameters, reactions, block, setTimeRange, outputs

println("Land model force transient Test starting");

#include("testLandForceTransient.jl");
# OR
variables = Variables([
    Variable("TRPN", 0.0752, (0, 100), "kg", ""),
    Variable("XB", 0.00046, (0, 100), "kg", ""),
]);

parameters = Variables([
    Variable("Tref", 156.067, (0, 1000), "kPa", "Maximum force that can be generated"),
    Variable("beta1", -1.5, (-100, 100), "", ""),
    Variable("lambda_const", 1.0, (0, 100), "", ""),
    Variable("n_xb", 5.0, (0, 100), "", ""),
    Variable("n_TRPN", 2.0, (0, 100), "", ""),
    Variable("TRPN_50", 0.35, (0, 100), "", ""),
    Variable("k_on", 0.1, (0, 100), "", ""),
    Variable("k_off", 0.0515, (0, 100), "", ""),
    Variable("k_xb", .0172, (0, 100), "", ""),
    Variable("CaT50", 1.0, (0, 100), "", ""),
]);

Ca50 = 2.1708;
parameters.CaT50 = Ca50 *(1 + parameters.beta1 * (parameters.lambda_const - 1));

Cai_transient = [0.14615806, 0.16809044, 0.23975627, 0.24182846, 0.26624959, 0.32315128, 0.40846157, 0.50964583, 0.6141502, 0.71285975, 0.80032514, 0.87399915, 0.93330914, 0.9789494, 1.01229084, 1.03501069, 1.04879625, 1.05545079, 1.0565415, 1.05305274, 1.0466294, 1.03790691, 1.02759521, 1.01665548, 1.00477103, 0.99288059, 0.9808191, 0.96862228, 0.95646617, 0.94456834, 0.93265957, 0.92062926, 0.90892036, 0.89721261, 0.88539131, 0.87372928, 0.86214212, 0.85059505, 0.83911839, 0.82774447, 0.81647891, 0.8053252, 0.79424877, 0.78327627, 0.7724101, 0.76160429, 0.7508795, 0.74023296, 0.72966264, 0.71916982, 0.70875044, 0.69840241, 0.68812979, 0.67794566, 0.66784322, 0.65782355, 0.64789118, 0.63805268, 0.62830954, 0.61866595, 0.60912213, 0.59967838, 0.59033832, 0.58110354, 0.57197444, 0.56295247, 0.5540384, 0.54523488, 0.5365426, 0.52796248, 0.51949613, 0.51114537, 0.50291124, 0.49479501, 0.48679746, 0.47892009, 0.47116413, 0.46352924, 0.45601652, 0.44862701, 0.44136046, 0.43421679, 0.42719693, 0.42030099, 0.41352811, 0.40687818, 0.40035057, 0.39394477, 0.38766014, 0.38149583, 0.37545106, 0.3695245, 0.36371534, 0.35802232, 0.352444, 0.34697929, 0.34162643, 0.33638408, 0.33125046, 0.32622404, 0.32130302, 0.31648575, 0.31177038, 0.30715516, 0.30263818, 0.29821769, 0.29389171, 0.28965844, 0.28551596, 0.28146237, 0.2774958, 0.27361434, 0.26981612, 0.26609934, 0.26246209, 0.2589025, 0.25541869, 0.25200905, 0.24867187, 0.2454052, 0.24220723, 0.23907637, 0.23601143, 0.23301034, 0.23007141, 0.22719294, 0.22437432, 0.22161345, 0.21890874, 0.21625863, 0.21366243, 0.21111866, 0.20862562, 0.20618187, 0.20378597, 0.20143794, 0.19913583, 0.19687828, 0.19466402, 0.19249196, 0.19036213, 0.18827266, 0.18622241, 0.18421026, 0.18223555, 0.18029796, 0.17839595, 0.17652854, 0.17469475, 0.17289523, 0.17112907, 0.16939451, 0.16769064, 0.16601657, 0.16437141, 0.16275424, 0.16116418, 0.15960033, 0.1580647, 0.15655702, 0.15507466, 0.15361691, 0.15218304, 0.15077233, 0.14938403, 0.14801744, 0.14667182];
Cai_interp = interpolate(Cai_transient, BSpline(Constant()));
Cai = extrapolate(Cai_interp, Flat());

reactions = [
    GeneralRateReaction("TRPN IN", [], ["TRPN"], (t, v, p) -> p.k_on * (1 - v.TRPN) * (Cai(t + 0.5) / p.CaT50) ^ p.n_TRPN),
    SimpleReaction("TRPN OUT", ["TRPN"], [], "k_off"),
    GeneralRateReaction("XB", [], ["XB"], (t, v, p) -> begin
        permtot = max(0, v.TRPN / p.TRPN_50) ^ (p.n_xb / 2);
        p.k_xb * (permtot * (1 - v.XB) - v.XB / permtot);
    end),
];

block = Block(variables, parameters, reactions);
setTimeRange!(block, 0:5:334);

outputs = Variables([
    Variable("f_act_vel", 0.0, (0, 100), "?", ""),
    Variable("f_relax_vel", 0.0, (0, 100), "?", ""),
    Variable("time_to_70perForce", 0.0, (0, 100), "?", ""),
    Variable("f_act_vel_norm", 0.0, (0, 100), "?", ""),
    Variable("f_relax_vel_norm", 0.0, (0, 100), "?", ""),
]);

setOutputDefinition!(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
    force = block.parameters.Tref * solution.XB;
    idx_maxF = min(length(force) - 1, max(2, argmax(force)));
    maxF = force[idx_maxF];
    t_at_maxF = timeRange[idx_maxF];
    
    dy = diff(force[1:idx_maxF]) / step(timeRange);
    outputs.f_act_vel = maximum(dy);
    outputs.f_act_vel_norm = maximum(dy) / maxF;

    dy = diff(force[idx_maxF:end]) / step(timeRange);
    outputs.f_relax_vel = maximum(-dy);
    outputs.f_relax_vel_norm = maximum(-dy) / maxF;

    outputs.time_to_70perForce = 3; #TODO: interpolate

    return outputs;
end);

println("Land model force transient Test defined");


# compute outputs
@time outputs = computeOutputs(block)

# Part 2: generate PPop
rangeMin = [120., 0.5];
rangeMax = [160., 1.5];

parameterBounds = [
    ("Tref", rangeMin[1], rangeMax[1]),
    ("CaT50", rangeMin[2], rangeMax[2])
];

outputBounds = Dict(
    # All outputs ranges for PPs presented as (MEAN, STD)
     "f_act_vel" => (50, 50),
     "f_relax_vel" => (50, 50),
     "time_to_70perForce" => (50, 50),
     "f_act_vel_norm" =>  (0.1, 0.04),
     "f_relax_vel_norm" => (0.02, 0.01),
 );

println("Generate PPops using Farthest Point optimization");
@time ppopFarthest = generatePPopFarthest([block], parameterBounds, outputBounds, 100;
                                          MaxTime = 120, DistanceFactor = 0.1, threads = 1); # num of pts, total max time for all PPs generation, distance to nearest existing pt is weighted relative to staying wihtin the output bounds
# ideally, would want at least 30s per pt, may need to add multithreading

plt = plot(ppopFarthest[1,:], ppopFarthest[2,:], seriestype = :scatter, legend = false, size = (600, 600));

# log plot
#plt = plot(log.(ppopFarthest[1,:]), log.(ppopFarthest[2,:]), seriestype = :scatter, legend = false, aspect_ratio = :equal, size = (600, 600));
display(plt);


println("Generate PPops");

@time ppop = generatePPop(block, parameterBounds, outputBounds, 100; MaxTime = 120);

plt = plot(ppop[1,:], ppop[2,:], seriestype = :scatter, legend = false, size = (600, 600));

# log plot
#plt = plot(log.(ppop[1,:]), log.(ppop[2,:]), seriestype = :scatter, legend = false, aspect_ratio = :equal, size = (600, 600));
display(plt);

(radii, areas) = computeDistanceCurve([block], parameterBounds, outputBounds, ppop, .1; samples = 10000);
(radii, areasFarthest) = computeDistanceCurve([block], parameterBounds, outputBounds, ppopFarthest, .1; samples = 10000);

plt = plot(radii, hcat(areas, areasFarthest), legend = true);
display(plt);
