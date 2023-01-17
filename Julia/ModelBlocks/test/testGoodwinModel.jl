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


