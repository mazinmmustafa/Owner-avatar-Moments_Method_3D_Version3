close all; clear; clc;

data  =   load("../data/S_matrix_cap.txt");

figure()
hold on
plot(data(:, 1)/1E-12, 20*log10(abs(data(:, 2))))
hold off
##ylim([-20 0])
[V, I] = min(data(:, 2))
data(I, 1)/1E-12