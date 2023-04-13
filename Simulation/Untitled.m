clc;
clear;
close all;
global dt;
dt = 0.001;
Ts = 10;
times = Ts / dt;
length = times + 1;
t = linspace(0, Ts, length);

f = sin(t);
fhat = zeros(length, 1);
for i = 1 : times
[fhat(i+1), fhat1d] = NonlinearBoundedErrorApproximator(f(i), fhat(i));
end
plot(t, f, t, fhat);

function [uhatNext, uhat1d] = NonlinearBoundedErrorApproximator(u, uhat)
    global dt;
    persistent epsilonu;
    if isempty(epsilonu)
        epsilonu = 0.01;
    end
    uhat1d = atanh(epsilonu^-1 * (u - uhat));
    uhatNext = dt * uhat1d + uhat;
end