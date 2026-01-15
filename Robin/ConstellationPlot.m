clc; clear; close all;
mod_type = 'PECS';

switch mod_type
    case 'M-PSK'
        % M-PSK Constellation Plot
        M = 8; % Modulation order (e.g., 2 for BPSK, 4 for QPSK, 8 for 8-PSK)

        % Generate M-PSK symbols
        k = 0:M-1; % Symbol indices
        symbols = exp(1j * 2 * pi * k / M); % M-PSK symbol generation

    case 'PECS'
        M = 10;
        L = 1;
        symbols = PECS_subpls(M, L);
end
% Plotting
figure;
plot(real(symbols), imag(symbols), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
grid on;
axis equal;
title([num2str(M) '-PSK Constellation Diagram']);
xlabel('In-phase');
ylabel('Quadrature');

% Annotate symbols
for i = 1:M
    text(real(symbols(i)) * 1.1, imag(symbols(i)) * 1.1, ...
        ['s_' num2str(i-1)], 'FontSize', 12);
end
phi_w = unwrap(angle(symbols));

figure; hold all; box on;
plot(phi_w);
legend('wrapped phase')
xlabel('Sample Index');
ylabel('Phase [rad]');

function [sqnce] = PECS_subpls(M, L)
Q = 3;
theta_li = [];
m = (1:M).';
a_coeff = zeros(3,L);
for l = 1:L
    a_coeff(:,l) = 1e1*rand(Q,1);
    theta_li(:,l) = a_coeff(1,l)+ a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
theta_x = theta_li(:);
x = exp(1j*theta_x);
sqnce = x;
end