clc; clear; close all;

% MIMO Configuration
Nt = 2;                 % Number of transmit antennas
Nr = 3;                 % Number of receive antennas (try 2 vs. 8)
M_values = [2, 4, 8];      % QPSK (M=4) and 8-PSK (M=8)
numSymbols = 10000;     % Increase for smoother BER curves
SNR_dB_range = 0:2:20;  % SNR range in dB

% Initialize BER results
BER = zeros(length(M_values), length(SNR_dB_range));

for m_idx = 1:length(M_values)
    M = M_values(m_idx);
    numBitsPerSymbol = log2(M);
    
    for snr_idx = 1:length(SNR_dB_range)
        SNR_dB = SNR_dB_range(snr_idx);
        SNR_linear = 10^(SNR_dB/10);
        noisePower = 1 / SNR_linear;
        
        % Generate random bits
        bitStream = randi([0 1], numSymbols * numBitsPerSymbol * Nt, 1);
        
        % Convert bits to symbols
        symbolIndices = bi2de(reshape(bitStream, numBitsPerSymbol, []).', 'left-msb');
        txSymbols = exp(1j * (2 * pi * symbolIndices / M));
        txMatrix = reshape(txSymbols, Nt, []);
        
        % MIMO Rayleigh channel
        H = (randn(Nr, Nt, numSymbols) + 1j * randn(Nr, Nt, numSymbols)) / sqrt(2);
        
        % Add AWGN
        noise = sqrt(noisePower/2) * (randn(Nr, numSymbols) + 1j * randn(Nr, numSymbols));
        rxMatrix = zeros(Nr, numSymbols);
        for k = 1:numSymbols
            rxMatrix(:,k) = H(:,:,k) * txMatrix(:,k) + noise(:,k);
        end
        
        % ZF Equalization (assume perfect CSI)
        rxSymbols_eq = zeros(Nt, numSymbols);
        for k = 1:numSymbols
            rxSymbols_eq(:,k) = pinv(H(:,:,k)) * rxMatrix(:,k);
        end
        
%         % Demodulation: Convert symbols back to bits
%         rxAngles = mod(angle(rxSymbols_eq(:)) / (2*pi/M), M);
%         rxBits = de2bi(rxAngles, numBitsPerSymbol, 'left-msb');
%         rxBits = rxBits(:);
        % Demodulation: Correctly map received symbols to integer indices
        rxAngles = angle(rxSymbols_eq(:)); % Get phase angles [-π, π]
        rxAngles = wrapTo2Pi(rxAngles);     % Convert to [0, 2π] (requires Mapping Toolbox)
        % Or manually wrap:
        % rxAngles = mod(rxAngles, 2*pi); 
        
        % Map angles to symbol indices (0 to M-1)
        rxIndices = round(rxAngles / (2*pi/M)); % Quantize to nearest PSK symbol
        rxIndices = mod(rxIndices, M);          % Wrap around (e.g., M=8 → 0-7)
        
        % Convert indices to bits
        rxBits = de2bi(rxIndices, numBitsPerSymbol, 'left-msb')';
        rxBits = rxBits(:);
        
        % Calculate BER
        BER(m_idx, snr_idx) = sum(rxBits ~= bitStream) / length(bitStream);
    end
end

% Plot BER vs. SNR
figure;
semilogy(SNR_dB_range, BER(1,:), 'bo-', 'LineWidth', 2); hold on;
semilogy(SNR_dB_range, BER(2,:), 'rs-', 'LineWidth', 2);
semilogy(SNR_dB_range, BER(3,:), 'gx-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('Binary (M=2)','QPSK (M=4)', '8-PSK (M=8)');
title(['BER vs. SNR for MIMO (Nt=', num2str(Nt), ', Nr=', num2str(Nr), ')']);