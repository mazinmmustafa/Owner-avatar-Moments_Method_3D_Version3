close all; clear; clc;

data  =   load("../data/S_matrix_freq.txt");

figure()
hold on
plot(data(:, 1)/1E6, 20*log10(abs(data(:, 2))))
hold off
##ylim([-20 0])