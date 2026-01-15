clc; clear; close all;

% Define M-PSK order (M = 2 for BPSK, M = 4 for QPSK, M = 8 for 8-PSK, etc.)
M = 8;  % Change this value to see different modulations

% Number of symbols to generate
numSymbols = 1000; 

% Generate random bits
numBitsPerSymbol = log2(M);
bitStream = randi([0 1], numSymbols * numBitsPerSymbol, 1);

% Convert bits to decimal symbol indices
symbolIndices = bi2de(reshape(bitStream, numBitsPerSymbol, []).', 'left-msb');

% Perform M-PSK Modulation
pskSymbols = exp(1j * (2 * pi * symbolIndices / M));

% Plot constellation diagram
figure;
scatter(real(pskSymbols), imag(pskSymbols), 'bo'); 
grid on;
hold on;
title([num2str(M) '-PSK Constellation']);
xlabel('In-Phase');
ylabel('Quadrature');
axis([-1.2 1.2 -1.2 1.2]);

% Plot ideal constellation points
idealSymbols = exp(1j * (2 * pi * (0:M-1) / M));
scatter(real(idealSymbols), imag(idealSymbols), 'rx', 'LineWidth', 2);

legend('Transmitted Symbols', 'Ideal Constellation Points');
hold off;
