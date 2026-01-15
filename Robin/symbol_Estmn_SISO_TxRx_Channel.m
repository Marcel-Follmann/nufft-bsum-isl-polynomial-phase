clc; clear; close all;

% Define M-PSK order (e.g., M=4 for QPSK, M=8 for 8-PSK)
M = 4;  
numSymbols = 10;  % Number of symbols

% Generate random bit stream
numBitsPerSymbol = log2(M);
bitStream = randi([0 1], numSymbols * numBitsPerSymbol, 1);

% Convert bits to decimal symbol indices
symbolIndices = bi2de(reshape(bitStream, numBitsPerSymbol, []).', 'left-msb');

% Perform M-PSK Modulation
txSymbols = exp(1j * (2 * pi * symbolIndices / M));  % Transmitted symbols

% Simulate Rayleigh fading channel (random complex gain per symbol)
h = (randn(numSymbols,1) + 1j*randn(numSymbols,1)) / sqrt(2);  % Rayleigh channel

% Additive White Gaussian Noise (AWGN)
SNR_dB = 10; % Signal-to-noise ratio in dB
SNR_linear = 10^(SNR_dB/10);
noisePower = 1 / (2 * SNR_linear);
noise = sqrt(noisePower) * (randn(numSymbols,1) + 1j*randn(numSymbols,1));

% Transmit symbols through channel
rxSymbols_no_eq = h .* txSymbols + noise;  % Received symbols before equalization

% Channel Estimation & Equalization
h_est = h;  % Assuming perfect channel estimation
rxSymbols_eq = rxSymbols_no_eq ./ h_est;  % Equalized symbols

% Plot results
figure;

subplot(1,3,1);
scatter(real(txSymbols), imag(txSymbols), 'bo');
title('Transmitted Symbols');
xlabel('In-Phase');
ylabel('Quadrature');
axis equal; grid on;

subplot(1,3,2);
scatter(real(rxSymbols_no_eq), imag(rxSymbols_no_eq), 'ro');
title('Received Symbols (Before Equalization)');
xlabel('In-Phase');
ylabel('Quadrature');
axis equal; grid on;

subplot(1,3,3);
scatter(real(rxSymbols_eq), imag(rxSymbols_eq), 'go');
title('Received Symbols (After Equalization)');
xlabel('In-Phase');
ylabel('Quadrature');
axis equal; grid on;

legend('Symbols');
