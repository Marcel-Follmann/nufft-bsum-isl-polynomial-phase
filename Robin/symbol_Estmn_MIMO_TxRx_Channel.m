clc; clear; close all;

% MIMO Configuration
Nt = 2;  % Number of transmit antennas
Nr = 3;  % Number of receive antennas
M = 4;   % QPSK (M = 4 for QPSK, 8 for 8-PSK, etc.)
numSymbols = 1000; % Number of symbols per antenna

% Generate Random Bits
numBitsPerSymbol = log2(M);
bitStream = randi([0 1], numSymbols * numBitsPerSymbol * Nt, 1);

% Convert bits to decimal symbol indices
symbolIndices = bi2de(reshape(bitStream, numBitsPerSymbol, []).', 'left-msb');

% Perform M-PSK Modulation (QPSK)
txSymbols = exp(1j * (2 * pi * symbolIndices / M));  % QPSK symbols

% Reshape to (Nt x numSymbols)
txMatrix = reshape(txSymbols, Nt, []);

% MIMO Rayleigh Fading Channel (Nr x Nt matrix per symbol)
H = (randn(Nr, Nt, numSymbols) + 1j * randn(Nr, Nt, numSymbols)) / sqrt(2);

% Add AWGN
SNR_dB = 10;  % SNR in dB
SNR_linear = 10^(SNR_dB/10);
noisePower = 1 / (2 * SNR_linear);
noise = sqrt(noisePower) * (randn(Nr, numSymbols) + 1j * randn(Nr, numSymbols));

% Transmit through MIMO Channel
rxMatrix = zeros(Nr, numSymbols);
for k = 1:numSymbols
    rxMatrix(:,k) = H(:,:,k) * txMatrix(:,k) + noise(:,k); % Received signal
end

% MIMO Channel Estimation (Assume perfect CSI)
H_est = H;

% Zero-Forcing (ZF) Equalization
rxSymbols_eq = zeros(Nt, numSymbols);
for k = 1:numSymbols
    rxSymbols_eq(:,k) = pinv(H_est(:,:,k)) * rxMatrix(:,k);  % Pseudo-inverse (ZF detection)
end

% Plot Constellations
figure; 
subplot(1,3,1);box on;
scatter(real(txMatrix(:)), imag(txMatrix(:)), 'bo');
title('Transmitted Symbols');
xlabel('In-Phase'); ylabel('Quadrature');
axis equal; 
xlim([-5,5]); ylim([-5,5]);
grid on;

subplot(1,3,2);box on;
scatter(real(rxMatrix(:)), imag(rxMatrix(:)), 'ro');
title('Received Symbols (Before Equalization)');
xlabel('In-Phase'); ylabel('Quadrature');
% axis equal;
xlim([-5,5]); ylim([-5,5]);
grid on;
legend('Received Signal');

subplot(1,3,3);box on;
scatter(real(rxSymbols_eq(:)), imag(rxSymbols_eq(:)), 'go');
title('Received Symbols (After Equalization)');
xlabel('In-Phase'); ylabel('Quadrature');
axis equal; 
xlim([-5,5]); ylim([-5,5]);
grid on;
legend('Symbols');
