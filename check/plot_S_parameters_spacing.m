close all; clear; clc;

data  =   load("../data/S_matrix_spacing.txt");

figure()
hold on
plot(data(:, 1)/1E-2/20, 20*log10(abs(data(:, 2))))
plot(data(:, 1)/1E-2/20, 20*log10(abs(data(:, 3))))
hold off
##ylim([-20 0])