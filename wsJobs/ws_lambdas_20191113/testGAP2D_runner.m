clear all;
close all;
clc;

try
	pc = parcluster('local');
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
end

testGAP2D();