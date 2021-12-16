close all 
clear all
%331.379e+3;
w = 2*pi*200e+3;
l = 4.7e-3;
r = 2.2e+3;

vsal = abs ((r*10) / (r+j*w*l));
vsaldb = abs ((r) / (r+j*w*l));
potenciaV = 20*log10(vsaldb);