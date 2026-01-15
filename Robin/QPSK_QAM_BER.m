clc; clear; close all;

% Simulation Parameters
M_qpsk = 4;        % QPSK
M_qam = 16;        % 16-QAM
numSymbols = 10000; % Number of symbols
SNR_dB = 0:2:20;   % SNR range (0 dB to 20 dB)
numBitsPerSymbol_qpsk = log2(M_qpsk);
numBitsPerSymbol_qam = log2(M_qam);

% Generate Random Bits
bits_qpsk = randi([0 1], numSymbols * numBitsPerSymbol_qpsk, 1);
bits_qam = randi([0 1], numSymbols * numBitsPerSymbol_qam, 1);

% Modulation
symbols_qpsk = pskmod(bi2de(reshape(bits_qpsk, numBitsPerSymbol_qpsk, []).'), M_qpsk, pi/M_qpsk);
symbols_qam = qammod(bi2de(reshape(bits_qam, numBitsPerSymbol_qam, []).'), M_qam, 'UnitAveragePower', true);

% Plot Constellation Diagrams
figure;
subplot(1,2,1);
scatter(real(symbols_qpsk), imag(symbols_qpsk), 'bo');
title('QPSK Constellation');
xlabel('In-Phase'); ylabel('Quadrature');
% axis equal; 
xlim([-5,5]); ylim([-5,5]);
grid on;

subplot(1,2,2);
scatter(real(symbols_qam), imag(symbols_qam), 'ro');
title('16-QAM Constellation');
xlabel('In-Phase'); ylabel('Quadrature');
% axis equal; 
xlim([-5,5]); ylim([-5,5]);
grid on;

% BER Analysis
ber_qpsk = zeros(size(SNR_dB));
ber_qam = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    SNR_linear = 10^(SNR_dB(i)/10);
    
    % Add AWGN Noise
    noisy_qpsk = awgn(symbols_qpsk, SNR_dB(i), 'measured');
    noisy_qam = awgn(symbols_qam, SNR_dB(i), 'measured');

    % Demodulation
    receivedBits_qpsk = de2bi(pskdemod(noisy_qpsk, M_qpsk, pi/M_qpsk), numBitsPerSymbol_qpsk, 'left-msb');
    receivedBits_qam = de2bi(qamdemod(noisy_qam, M_qam, 'UnitAveragePower', true), numBitsPerSymbol_qam, 'left-msb');

    % BER Calculation
    ber_qpsk(i) = sum(bits_qpsk ~= receivedBits_qpsk(:)) / length(bits_qpsk);
    ber_qam(i) = sum(bits_qam ~= receivedBits_qam(:)) / length(bits_qam);
end

% Plot BER vs. SNR
figure;
semilogy(SNR_dB, ber_qpsk, 'b-o', 'LineWidth', 2); hold on;
semilogy(SNR_dB, ber_qam, 'r-s', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR for QPSK and 16-QAM');
legend('QPSK', '16-QAM');
