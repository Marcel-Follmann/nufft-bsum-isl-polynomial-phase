% Example usage
N = 99;
n = 0:N; 
phi = mod(2 * pi * sin(0.1 * n), 2 * pi); % Simulated wrapped phase
unwrapped_phase = unwrap_L2(phi);

figure;
plot(n, phi, 'r--', 'LineWidth', 1.5); hold on;
plot(n, unwrapped_phase, 'b', 'LineWidth', 1.5);
legend('Wrapped Phase', 'Unwrapped Phase');
xlabel('Sample Index'); ylabel('Phase (rad)');
title('1D Phase Unwrapping using L2-norm');
grid on;


function psi = unwrap_L2(phi)
    N = length(phi);
    dphi = diff(phi);  % Wrapped phase differences
    dpsi = round((dphi + pi) / (2 * pi)) * 2 * pi - pi;  % Correct phase jumps
    psi = [phi(1), phi(1) + cumsum(dphi - dpsi)];  % Integrate corrected differences
end