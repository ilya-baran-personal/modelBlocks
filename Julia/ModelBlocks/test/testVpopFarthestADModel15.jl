# module admodel
using Test
using LinearAlgebra
using Interpolations
using Plots

using ModelBlocks

# Part 1: define model - variables, parameters, reactions, block, setTimeRange, outputs

println("AD model starting");
# include("module_admodel.jl")
# Read text file with nominal parameter values
using CSV
using DataFrames

x1 = CSV.read("Julia/ModelBlocks/test/data/mu.csv", DataFrame, header = false)
x1 = Array{Float64}(x1)
x = exp.(x1)

# Define variables
variables = Variables([
    Variable("s", 0.5931, (0, 100 * 0.5931), "a.u.", "1"),
    Variable("path", 0.4069, (0, 100 * 0.4069), "a.u.", "2"),
    Variable("Th1", 3.1, (0, 100 * 3.1), "pg/mL", "3"),
    Variable("Th2", 8.7, (0, 100 * 8.7), "pg/mL", "4"),
    Variable("Th17", 2.0, (0, 100 * 2.0), "pg/mL", "5"),
    Variable("Th22", 21.0, (0, 100 * 21.0), "pg/mL", "6"),
    Variable("IL4", 38.0, (0, 100 * 38.0), "pg/mL", "7"),
    Variable("IL13", 40.5, (0, 100 * 40.5), "pg/mL", "8"),
    Variable("IL17", 5.4, (0, 100 * 5.4), "pg/mL", "9"),
    Variable("IL22", 3.0, (0, 100 * 3.0), "pg/mL", "10"),
    Variable("IL31", 2.0, (0, 100 * 2.0), "pg/mL", "11"),
    Variable("IFNg", 1.5, (0, 100 * 1.5), "pg/mL", "12"),
    Variable("TSLP", 4.4, (0, 100 * 4.4), "pg/mL", "13"),
    Variable("OX40", 9.7, (0, 100 * 9.7), "pg/mL", "14")])

# Define parameters 
parameters = Variables([
    Variable("k1", x[1, 1], (0, 100 * x[1, 1]), "a.u.", "p1"),
    Variable("k2", x[2, 1], (0, 100 * x[2, 1]), "a.u.", "p2"),
    Variable("b1", x[4, 1], (0, 100 * x[4, 1]), "a.u.", "p4"),
    Variable("b2", x[5, 1], (0, 100 * x[5, 1]), "a.u.", "p5"),
    Variable("b3", x[6, 1], (0, 100 * x[6, 1]), "a.u.", "p6"),
    Variable("b4", x[7, 1], (0, 100 * x[7, 1]), "a.u.", "p7"),
    Variable("b5", x[8, 1], (0, 100 * x[8, 1]), "a.u.", "p8"),
    Variable("d1", x[9, 1], (0, 100 * x[9, 1]), "a.u.", "p9"),
    Variable("d2", x[10, 1], (0, 100 * x[10, 1]), "a.u.", "p10"),
    Variable("d3", x[11, 1], (0, 100 * x[11, 1]), "a.u.", "p11"),
    Variable("b6", x[12, 1], (0, 100 * x[12, 1]), "a.u.", "p12"),
    Variable("d4", x[13, 1], (0, 100 * x[13, 1]), "a.u.", "p13"),
    Variable("d5", x[14, 1], (0, 100 * x[14, 1]), "a.u.", "p14"),
    Variable("d6", x[15, 1], (0, 100 * x[15, 1]), "a.u.", "p15"),
    Variable("d7", x[16, 1], (0, 100 * x[16, 1]), "a.u.", "p16"),
    Variable("b7", x[17, 1], (0, 100 * x[17, 1]), "a.u.", "p17"),
    Variable("b8", x[18, 1], (0, 100 * x[18, 1]), "a.u.", "p18"),
    Variable("d8", x[19, 1], (0, 100 * x[19, 1]), "a.u.", "p19"),
    Variable("k5", x[20, 1], (0, 100 * x[20, 1]), "a.u.", "p20"),
    Variable("k9", x[21, 1], (0, 100 * x[21, 1]), "a.u.", "p21"),
    Variable("d9", x[22, 1], (0, 100 * x[22, 1]), "a.u.", "p22"),
    Variable("b9", x[23, 1], (0, 100 * x[23, 1]), "a.u.", "p23"),
    Variable("k6", x[24, 1], (0, 100 * x[24, 1]), "a.u.", "p24"),
    Variable("k10", x[25, 1], (0, 100 * x[25, 1]), "a.u.", "p25"),
    Variable("k7", x[26, 1], (0, 100 * x[26, 1]), "a.u.", "p26"),
    Variable("k8", x[27, 1], (0, 100 * x[27, 1]), "a.u.", "p27"),
    Variable("k11", x[28, 1], (0, 100 * x[28, 1]), "a.u.", "p28"),
    Variable("k12", x[29, 1], (0, 100 * x[29, 1]), "a.u.", "p29"),
    Variable("d10", x[30, 1], (0, 100 * x[30, 1]), "a.u.", "p30"),
    Variable("k13", x[31, 1], (0, 100 * x[31, 1]), "a.u.", "p31"),
    Variable("k14", x[32, 1], (0, 100 * x[32, 1]), "a.u.", "p32"),
    Variable("d11", x[33, 1], (0, 100 * x[33, 1]), "a.u.", "p33"),
    Variable("k15", x[34, 1], (0, 100 * x[34, 1]), "a.u.", "p34"),
    Variable("k16", x[35, 1], (0, 100 * x[35, 1]), "a.u.", "p35"),
    Variable("d12", x[36, 1], (0, 100 * x[36, 1]), "a.u.", "p36"),
    Variable("k17", x[37, 1], (0, 100 * x[37, 1]), "a.u.", "p37"),
    Variable("k18", x[38, 1], (0, 100 * x[38, 1]), "a.u.", "p38"),
    Variable("d13", x[39, 1], (0, 100 * x[39, 1]), "a.u.", "p39"),
    Variable("k19", x[40, 1], (0, 100 * x[40, 1]), "a.u.", "p40"),
    Variable("k20", x[41, 1], (0, 100 * x[41, 1]), "a.u.", "p41"),
    Variable("d14", x[42, 1], (0, 100 * x[42, 1]), "a.u.", "p42"),
    Variable("k21", x[43, 1], (0, 100 * x[43, 1]), "a.u.", "p43"),
    Variable("k22", x[44, 1], (0, 100 * x[44, 1]), "a.u.", "p44"),
    Variable("d15", x[45, 1], (0, 100 * x[45, 1]), "a.u.", "p45"),
    Variable("k23", x[46, 1], (0, 100 * x[46, 1]), "a.u.", "p46"),
    Variable("k24", x[47, 1], (0, 100 * x[47, 1]), "a.u.", "p47"),
    Variable("d16", x[48, 1], (0, 100 * x[48, 1]), "a.u.", "p48"),
    Variable("k25", x[49, 1], (0, 100 * x[49, 1]), "a.u.", "p49"),
    Variable("k26", x[50, 1], (0, 100 * x[50, 1]), "a.u.", "p50"),
    Variable("d17", x[51, 1], (0, 100 * x[51, 1]), "a.u.", "p51"),
    Variable("de1", 1e20, (0, 1e20), "a.u.", "p51"),
    Variable("de2", 0, (0, 1), "a.u.", "p51"),
    Variable("de3", 0, (0, 1), "a.u.", "p51"),
    Variable("de4", 0, (0, 1), "a.u.", "p51"),
    Variable("de5", 0, (0, 1), "a.u.", "p51"),
    Variable("de6", 0, (0, 1), "a.u.", "p51"),
    Variable("de7", 0, (0, 1), "a.u.", "p51"),
    Variable("de8", 0, (0, 1), "a.u.", "p51"),
    Variable("de9", 0, (0, 1), "a.u.", "p51"),
    Variable("de10", 0, (0, 1), "a.u.", "p51"),
    Variable("k3a", x[3, 1], (0, 100 * x[3, 1]), "a.u.", "p3"),
])

reactions = [
    GeneralReaction("ADmodel", (dvdt, t, v, p, e) -> begin
        de_il4 = (1 - p.de2)
        ea2 = max(0.4396, p.de10)
        de_il13 = (1 - p.de3 * ea2)
        de_il17 = (1 - p.de4)
        de_il22 = (1 - p.de5)
        de_il31 = (1 - p.de6)
        de_tslp = (1 - p.de7)
        de_ox40 = (1 - p.de8)
        de_ifng = p.de9

        k3 = min(p.k3a, p.de1)
        k4 = p.d8
        #   ODEs
        # skin barrier integrity        
        dvdt.s = (1 - v.s) * (p.k1 + p.k2 * de_il22 * v.IL22 + k3) / ((1 + p.b1 * de_il4 * v.IL4) * (1 + p.b2 * de_il13 * v.IL13) * (1 + p.b3 * de_il17 * v.IL17) * (1 + p.b4 * de_il22 * v.IL22) * (1 + p.b5 * de_il31 * v.IL31)) - v.s * (p.d1 * (1 + p.d3 * v.path) + p.d2 * de_il31 * v.IL31)

        # infiltrated pathogens        
        dvdt.path = k4 / (1 + p.b6 * v.s) - v.path * (((1 + p.d4 * v.path) * (1 + p.d5 * de_il17 * v.IL17) * (1 + p.d6 * de_il22 * v.IL22) * (1 + p.d7 * (v.IFNg + de_ifng))) / ((1 + p.b7 * de_il4 * v.IL4) * (1 + p.b8 * de_il13 * v.IL13)) + p.d8)

        # Th cells
        dvdt.Th1 = p.k5 * v.path * (1 + p.k9 * (v.IFNg + de_ifng)) / (4 + p.k9 * (v.IFNg + de_ifng) + p.k10 * de_il4 * v.IL4) - p.d9 * v.Th1 / (1 + p.b9 * de_ox40 * v.OX40) # Th1
        dvdt.Th2 = p.k6 * v.path * (1 + p.k10 * de_il4 * v.IL4) / (4 + p.k9 * (v.IFNg + de_ifng) + p.k10 * de_il4 * v.IL4) - p.d9 * v.Th2 / (1 + p.b9 * de_ox40 * v.OX40)  # Th2
        dvdt.Th17 = p.k7 * v.path / (4 + p.k9 * (v.IFNg + de_ifng) + p.k10 * de_il4 * v.IL4) - p.d9 * v.Th17 / (1 + p.b9 * de_ox40 * v.OX40) # Th17
        dvdt.Th22 = p.k8 * v.path / (4 + p.k9 * (v.IFNg + de_ifng) + p.k10 * de_il4 * v.IL4) - p.d9 * v.Th22 / (1 + p.b9 * de_ox40 * v.OX40) # Th22

        # cytokines        
        dvdt.IL4 = p.k11 * v.Th2 + p.k12 - p.d10 * v.IL4  # v.IL4
        dvdt.IL13 = p.k13 * v.Th2 + p.k14 - p.d11 * v.IL13  # v.IL13
        dvdt.IL17 = p.k15 * v.Th17 + p.k16 - p.d12 * v.IL17  # v.IL17
        dvdt.IL22 = p.k17 * v.Th22 + p.k18 - p.d13 * v.IL22  # v.IL22
        dvdt.IL31 = p.k19 * v.Th2 + p.k20 - p.d14 * v.IL31 # v.IL31
        dvdt.IFNg = p.k21 * v.Th1 + p.k22 - p.d15 * v.IFNg # v.IFNg
        dvdt.TSLP = p.k23 * v.path + p.k24 - p.d16 * v.TSLP # TSLP
        dvdt.OX40 = p.k25 * de_tslp * v.TSLP + p.k26 - p.d17 * v.OX40 # OX40L
    end),
]


block = Block(variables, parameters, reactions)
setTimeRange!(block, 0:100) # can specify output points (pts at which differential equation is computed) as 0:0.5:100

@time solution = runBlock(block)

plot(solution)



# For case of varying just 2 parameters
# outputs = Variables([
#     Variable("s_out", 0.5931, (0, 100 * 0.5931), "a.u.", "1"),
#     Variable("path_out", 0.4069, (0, 100 * 0.4069), "a.u.", "2"),
#     Variable("Th1_out", 3.1, (0, 100 * 3.1), "pg/mL", "3"),
# ])

# setOutputDefinition!(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
#     # outputs.vector = [solution.s, solution.path, solution.Th1, solution.Th2, solution.Th17, solution.Th22, solution.IL4, solution.IL13, solution.IL17, solution.IL22, solution.IL31, solution.IFNg, solution.TSLP, solution.OX40]
#     outputs.Th1_out = solution.Th1[end]
#     outputs.s_out = solution.s[end]
#     outputs.path_out = solution.path[end]

#     return outputs
# end)


# For the case of varying 5 parameters
outputs = Variables([
    Variable("Th2_out", 0.5931, (0, 100 * 0.5931), "a.u.", "1"),
    Variable("path_out", 0.4069, (0, 100 * 0.4069), "a.u.", "2"),
    Variable("Th1_out", 3.1, (0, 100 * 3.1), "cells/uL", "3"),
    Variable("IL13_out", 40, (0, 100 * 40), "pg/mL", "4"),
    Variable("IL4_out", 35, (0, 100 * 40), "pg/mL", "5"),
])

setOutputDefinition!(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
    # outputs.vector = [solution.s, solution.path, solution.Th1, solution.Th2, solution.Th17, solution.Th22, solution.IL4, solution.IL13, solution.IL17, solution.IL22, solution.IL31, solution.IFNg, solution.TSLP, solution.OX40]
    outputs.Th1_out = solution.Th1[end]
    outputs.Th2_out = solution.Th2[end]
    outputs.path_out = solution.path[end]
    outputs.IL13_out = solution.IL13[end]
    outputs.IL4_out = solution.IL4[end]

    return outputs
end)

computedOutputs = computeOutputs(block)
@time outputs = computeOutputs(block)


# For the case of varying just 2 parameters
# # Part 2: generate PPop
# rangeMin = [x[1,1]/10, x[19,1]/10, x[20,1]/10];
# rangeMax = [x[1, 1]*10., x[19, 1]*10, x[20,1]*10];
# parameterBounds = [
#     ("k1", rangeMin[1], rangeMax[1]),
#     ("d8", rangeMin[2], rangeMax[2]),
#     # ("k5", rangeMin[3], rangeMax[3])
# ];

# outputBounds = Dict(
#     # All outputs ranges for PPs presented as (MEAN, STD)
#      "s_out" => (0.5, 1.0),
#      "path_out" => (0.1, 1.0),
#      "Th1_out" => (0.5, 10.0)
#  );

# Move to logarithm space
extraParameters = Variables([
    Variable("logK1", 1.0, (0, 100), "kg/s", ""),
    Variable("logD8", 1.0, (0, 100), "kg/s", ""),
    Variable("logK5", 1.0, (0, 100), "kg/s", ""),
    Variable("logK6", 1.0, (0, 100), "kg/s", ""),
    Variable("logK11", 1.0, (0, 100), "kg/s", ""),
    Variable("logK12", 1.0, (0, 100), "kg/s", ""),
    Variable("logK13", 1.0, (0, 100), "kg/s", ""),
    Variable("logK14", 1.0, (0, 100), "kg/s", ""),
    Variable("logK15", 1.0, (0, 100), "kg/s", ""),
    Variable("logK17", 1.0, (0, 100), "kg/s", ""),
    Variable("logK18", 1.0, (0, 100), "kg/s", ""),
    Variable("logK19", 1.0, (0, 100), "kg/s", ""),
    Variable("logK20", 1.0, (0, 100), "kg/s", ""),
    Variable("logK21", 1.0, (0, 100), "kg/s", ""),
    Variable("logK22", 1.0, (0, 100), "kg/s", ""),
]);

block = BlockCombo(
    [("first", block)],
    [
        ("first", "k1", (t, v, p) -> exp(p.logK1)),
        ("first", "d8", (t, v, p) -> exp(p.logD8)),
        ("first", "k5", (t, v, p) -> exp(p.logK5)),
        ("first", "k6", (t, v, p) -> exp(p.logK6)),
        ("first", "k11", (t, v, p) -> exp(p.logK11)),
        ("first", "k12", (t, v, p) -> exp(p.logK12)),
        ("first", "k13", (t, v, p) -> exp(p.logK13)),
        ("first", "k14", (t, v, p) -> exp(p.logK14)),
        ("first", "k15", (t, v, p) -> exp(p.logK15)),
        ("first", "k17", (t, v, p) -> exp(p.logK17)),
        ("first", "k18", (t, v, p) -> exp(p.logK18)),
        ("first", "k19", (t, v, p) -> exp(p.logK19)),
        ("first", "k20", (t, v, p) -> exp(p.logK20)),
        ("first", "k21", (t, v, p) -> exp(p.logK21)),
        ("first", "k22", (t, v, p) -> exp(p.logK22)),
    ],
    extraParameters
);

# For the case of varying 5 parameters
# Part 2: generate PPop
indices = [1, 19, 20, 24, 28, 29, 31, 32, 34, 37, 38, 40, 41, 43, 44];
 
rangeMin = [log(x[idx, 1] / 10) for idx in indices];
rangeMax = [log(x[idx, 1] * 10) for idx in indices];

parameterBounds = [
    ("logK1", rangeMin[1], rangeMax[1]),
    ("logD8", rangeMin[2], rangeMax[2]),
    ("logK5", rangeMin[3], rangeMax[3]),
    ("logK6", rangeMin[4], rangeMax[4]),
    ("logK11", rangeMin[5], rangeMax[5]),
    ("logK12", rangeMin[6], rangeMax[6]),
    ("logK13", rangeMin[7], rangeMax[7]),
    ("logK14", rangeMin[8], rangeMax[8]),
    ("logK15", rangeMin[9], rangeMax[9]),
    ("logK17", rangeMin[10], rangeMax[10]),
    ("logK18", rangeMin[11], rangeMax[11]),
    ("logK19", rangeMin[12], rangeMax[12]),
    ("logK20", rangeMin[13], rangeMax[13]),
    ("logK21", rangeMin[14], rangeMax[14]),
    ("logK22", rangeMin[15], rangeMax[15]),
    ];

outputBounds = Dict(
    # All outputs ranges for PPs presented as (MEAN, STD)
    "Th2_out" => (0.5, 1.0),
    "path_out" => (0.1, 1.0),
    "Th1_out" => (0.5, 5.0),
    "IL4_out" => (25, 50),
    "IL13_out" => (35, 60)
);


println("Generate PPops");

@time ppop = generatePPop(block, parameterBounds, outputBounds, 500; MaxTime = 500);

plt = plot(ppop[1,:], ppop[2,:], seriestype = :scatter, legend = false, size = (600, 600),title = "non-farthest point");

# log plot
#plt = plot(log.(ppop[1,:]), log.(ppop[2,:]), seriestype = :scatter, legend = false, aspect_ratio = :equal, size = (600, 600),title = "log non-farthest point (BB)");
display(plt);

println("Generate PPops using Farthest Point optimization");
@time ppopFarthest = generatePPopFarthest([block], parameterBounds, outputBounds, 500;
                                          MaxTime = 500, DistanceFactor = 0.1, threads = 1); # num of pts, total max time for all PPs generation, distance to nearest existing pt is weighted relative to staying wihtin the output bounds
# ideally, would want at least 30s per pt, may need to add multithreading

plt = plot(ppopFarthest[1,:], ppopFarthest[2,:], seriestype = :scatter, legend = false, size = (600, 600),title = "Farthest point");

# log plot
#plt = plot(log.(ppopFarthest[1,:]), log.(ppopFarthest[2,:]), seriestype = :scatter, legend = false, aspect_ratio = :equal, size = (600, 600),title = "log farthest point");
display(plt);


# Growth volume curve plotting
(radii, areas) = computeDistanceCurve([block], parameterBounds, outputBounds, ppop, 1.2; samples = 10000);
(radii, areasFarthest) = computeDistanceCurve([block], parameterBounds, outputBounds, ppopFarthest, 1.2; samples = 10000);

#plt = plot(radii, hcat(areas, areasFarthest), legend = false,title = "increased to 6.0");
#display(plt);


plt = plot(radii, hcat(areas, areasFarthest), label =["SA" "FP"]);#, title = "CDF or Volume-Growth curves: AD model- 5 params");
plt = plot(plt, thickness_scaling=2, tickfontsize=10/2, labelfontsize=14/2, colorbar_tickfontsize=8/2, reuse=false);
display(plt);

# Print the integral of the difference
display((sum(areasFarthest) - sum(areas)) * (radii[2] - radii[1]));
