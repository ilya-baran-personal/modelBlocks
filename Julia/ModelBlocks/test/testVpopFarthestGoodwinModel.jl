# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

# Reproducing Goodwin model of multi-component negative feedback oscillator
# Nature, 2008. https://www.nature.com/articles/nrm2530
# Figure 2d–f: Negative feedback loop with nuclear transport
# mRNA is synthesized in the nucleus (xn) and transported into the cytoplasm (xc)
# where it gets translated into protein (yc) which is tranlocated into the nucleus (yn) #
# eps = Vnuc/Vcyt
# half-life of mRNA in nucleus = 0.693/kdxn
# half-life of prot in cytoplasm = 0.693/kdyc
# X = mRNA concentration 
# Y = protein concentration
# n - in nucleus, c - in cytoplasm
# kdx and ksy are degradation and synthesis constants respectively
#dxn/dt = kdxn*(sig/(1 + yn^p) - xn) - kexport*xn dxc/dt = eps*kexport*xn - kdxc*xc
#dyc/dt = kdyc*(xc - yc) - eps*kimport*yc dyn/dt = kimport*yc - kdyn*yn/(Km + yn)
#p Sig=1000, p=2, kdxn=10, kexport=0.2,kdxc=0.2, eps=1 p kdyn=8, kdyc=0.1, Km=0.1, kimport=0.1
#@ XP=t, YP=xn, TOTAL=100, METH=stiff, XLO=0, XHI=100, YLO=0, YHI=1000, bounds=10000

using ModelBlocks
using Test
using Plots
using Interpolations
using LinearAlgebra

println("Goodwin model Test starting");

variables = Variables([
    Variable("xn", 100, (0, 500), "", ""),
    Variable("xc", 70, (0, 500), "", ""),
    Variable("yn", 50, (0, 500), "", ""),
    Variable("yc",100, (0, 500), "", ""),
]);

parameters = Variables([
    Variable("sig", 1000.0, (0, 2000), "", ""),
    Variable("pparam", 2.0, (0, 100), "", ""),
    Variable("kdxn", 10.0, (0, 100), "", ""), # half-life of mRNA in nucleus = 0.693/kdxn
    Variable("kexport", 0.2, (0, 100), "", ""),
    Variable("kdxc", 0.2, (0, 100), "", ""), 
    Variable("epsparam", 1.0, (0, 100), "", ""),  # eps = Vnuc/Vcyt
    Variable("kdyn", 8.0, (0, 100), "", ""),
    Variable("kdyc", 0.1, (0, 100), "", ""), # half-life of prot in cytoplasm = 0.693/kdyc
    Variable("Km", 0.1, (0, 100), "", ""),
    Variable("kimport", 0.1, (0, 100), "", ""),
]);

reactions = [
     GeneralReaction("Goodwinmodel", (dvdt, t, v, p, e) -> begin
        dvdt.xn = p.kdxn*(p.sig/(1.0 + max(0, v.yn)^p.pparam) - v.xn) - p.kexport*v.xn;
        dvdt.xc = p.epsparam*p.kexport*v.xn - p.kdxc*v.xc;
        dvdt.yn = p.kimport*v.yc - p.kdyn*v.yn/(p.Km + v.yn);
        dvdt.yc = p.kdyc*(v.xc - v.yc) - p.epsparam*p.kimport*v.yc;
    end),
];

block = Block(variables, parameters, reactions);
setTimeRange!(block, 0:1:100);

@time solution = runBlock(block)


println("Goodwin model ODEs solved. Solution (hopefully) plotted");
display(plot(solution))

function isOscillatory(solution)::Bool
    local m = hcat(solution.u...);
    local numOscillatory = 0;
    for v in 1:size(m, 1)
        local f = m[v,:];
        local thresh = maximum(f) * 0.9 + minimum(f) * 0.1;
        # count how many times f crosses the threshold
        local count = sum(f[1:end-1] .< thresh .&& f[2:end] .>= thresh);
        numOscillatory += (count >= 3 ? 1 : 0);
    end
    return numOscillatory > 1;
end

outputs = Variables([
   Variable("mRNAnucleus", 0.0, (0, 1000), "", ""),
   Variable("mRNAcytoplasm", 0.0, (0, 1000), "", ""),
   Variable("Protnucleus", 0.0, (0, 1000), "", ""),
   Variable("Protcytoplasm", 0.0, (0, 1000), "", ""),
   Variable("oscillationmRNA", 0.0, (0, 1000), "", ""),
   Variable("oscilationsProt", 0.0, (0, 1000), "", ""),
   ]);

setOutputDefinition!(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
    outputs.mRNAnucleus = 4;    # min(solution.xn);
    outputs.mRNAcytoplasm = 2;  # min(solution.xc);
    outputs.Protnucleus = 2;    # maximum(solution.yn);
    outputs.Protcytoplasm = 3;  # maximum(solution.yc);
    outputs.oscillationmRNA = 1;
    outputs.oscilationsProt = 1;
    return outputs;
end);

println("Goodwin model outputs specified");

# compute outputs
@time outputs = computeOutputs(block)


## Part 2: generate PPop

# Vary 2 parameters to generate plausible patients
parameterBounds = [
    ("pparam", 0.0, 15.0),
    ("kdyc", 0.001, 100.0),
];

outputBounds = Dict(
    # All outputs ranges for PPs presented as (MEAN, STD)
     "mRNAnucleus" => (500, 500),
     "mRNAcytoplasm" => (500, 500),
     "Protnucleus" => (500, 500),
     "Protcytoplasm" =>  (500, 500),
     "oscillationmRNA" => (1, 0.5),
     "oscilationsProt" => (1,0.5),
 );

display(isOscillatory(solution))
 # Generate PPop using FPO
println("Generating PPs for Goodwin model using FPO");
println("Generate PPops using Farthest Point optimization");
@time ppopFarthest = generatePPopFarthest([block], parameterBounds, outputBounds, 20;
                                          MaxTime = 60, DistanceFactor = 0.1, threads = 1); # num of pts, total max time for all PPs generation, distance to nearest existing pt is weighted relative to staying wihtin the output bounds
# ideally, would want at least 30s per pt, may need to add multithreading

plt = plot(ppopFarthest[1,:], ppopFarthest[2,:], seriestype = :scatter, legend = false, size = (600, 600));

display(plt);

 # Generate PPop using traditional approach
println("Generate PPops");
@time ppop = generatePPop(block, parameterBounds, outputBounds, 20; MaxTime = 60);

plt = plot(ppop[1,:], ppop[2,:], seriestype = :scatter, legend = false, size = (600, 600));

display(plt);

# # Compute CDF 
# (radii, areas) = computeDistanceCurve([block], parameterBounds, outputBounds, ppop, .1; samples = 10000);
# (radii, areasFarthest) = computeDistanceCurve([block], parameterBounds, outputBounds, ppopFarthest, .1; samples = 10000);

# plt = plot(radii, hcat(areas, areasFarthest), legend = false, label = ["SA" "FP"]);
# plt = plot(plt, thickness_scaling=2, tickfontsize=10/2, labelfontsize=14/2, colorbar_tickfontsize=8/2, reuse=false);
# display(plt);

# # Print the integral of the difference
# display((sum(areasFarthest) - sum(areas)) * (radii[2] - radii[1]));

