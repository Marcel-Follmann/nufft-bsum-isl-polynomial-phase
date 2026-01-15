clc; clear; close all;

% MIMO dimensions
Nt = 2;  % Number of transmit antennas
Nr = 4;  % Number of receive antennas
Ns = 1000;  % Number of symbols per antenna

% Generate random QPSK symbols for 2 Tx antennas
txSymbols = (randi([0, 1], Nt, Ns) + 1j * randi([0, 1], Nt, Ns)) * 2 - (1+1j);

% Rayleigh fading channel (Nr x Nt)
H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2);

% Add AWGN noise
SNR_dB = 10;
SNR_linear = 10^(SNR_dB/10);
noisePower = 1 / (2 * SNR_linear);
noise = sqrt(noisePower) * (randn(Nr, Ns) + 1j * randn(Nr, Ns));

% Received signal at 4 Rx antennas
rxSignal = H * txSymbols + noise;

% Zero Forcing Equalization
H_inv_ZF = pinv(H);  % Moore-Penrose Pseudo-Inverse
zfEqualized = H_inv_ZF * rxSignal;

% MMSE Equalization
H_MMSE = (H' * H + (1/SNR_linear) * eye(Nt)) \ (H');  
mmseEqualized = H_MMSE * rxSignal;

% Plot Original and Equalized Symbols
figure;
subplot(1,3,1);
scatter(real(txSymbols(:)), imag(txSymbols(:)), 'b');
title('Transmitted Symbols'); grid on; axis equal;

subplot(1,3,2);
scatter(real(zfEqualized(:)), imag(zfEqualized(:)), 'r');
title('Zero Forcing Equalized'); grid on; axis equal;

subplot(1,3,3);
scatter(real(mmseEqualized(:)), imag(mmseEqualized(:)), 'g');
title('MMSE Equalized'); grid on; axis equal;
