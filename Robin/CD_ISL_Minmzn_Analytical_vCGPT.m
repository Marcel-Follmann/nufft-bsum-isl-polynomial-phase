clc; clear; close all;

%% Parameters
N = 8;
L = 1;
Q = 3;
M = N / L;

maxIter = 100;
tol = 1e-9;

a_coeff = rand(Q, L);
[x, phi] = create_sequence(a_coeff, M, L);


%% Auto-correlation equivalence test
ISL = 0;
a2 = a_coeff(3);  % Extract initial a2
a2_i = a2;
r = zeros(N-1,1);
flg_test = 't2';    %'f'; % 't1'/ 't2' / 't3' / 't4'
switch(flg_test)
    case 't1'
% r_k = \sum_{n=1}^{N-k} e^{j \psi_{nk}}
        for k = 1:N-1
            for n = 1:N-k
                phi_n = a_coeff(1) + a_coeff(2)*n + a2 * n^2;
                phi_nk = a_coeff(1) + a_coeff(2)*(n+k) + a2 * (n + k)^2;
                psi_nk = phi_n - phi_nk;
                r(k) = r(k) + exp(1j * psi_nk);
            end
            ISL = ISL + abs(r(k))^2;
        end
        rm = xcorr(x);

        lag = 1:N-1;
        figure(1); hold all; box on;
        plot(lag, abs(r), 'r', 'LineWidth',2);
        plot(lag, abs(rm(N+1:2*N-1)), '--b*', 'LineWidth',2);

        ISLm = sum(abs(rm(N+1:2*N-1)).^2);

        fprintf('ISL computed using my logic: %f \n', ISL);
        fprintf('ISL computed using Matlab: %f \n', ISLm);
    case 't2'
% \text{ISL} = \sum_{k=1}^{N-1} \Bigg[ \sum_{n=1}^{N-k} \sum_{n'=1}^{N-k} e^{j\psi_{knn'}} \Bigg]
        s1 = 0; s2 = 0;
        for k = 1:N-1
            for n = 1:N-k
                for np = 1:N-k
                    psi_knnp = 2*a2*k*(n - np);
                    s2 = s2 + exp(1j * psi_knnp);
                end
            end
        end
        ISL = s1+s2;
        rm = xcorr(x);
        ISLm = sum(abs(rm(N+1:2*N-1)).^2);

        fprintf('ISL computed using my logic: %f \n', ISL);
        fprintf('ISL computed using Matlab: %f \n', ISLm);

    case 't3'
% \text{ISL} = \sum_{k=1}^{N-1} \bigg[ \bigg(\sum_{n'=1}^{N-k}e^{-j \psi_{n'k}} \bigg). \bigg(\sum_{n=1}^{N-k}e^{j \psi_{nk}} \bigg)\bigg]
        s1 = 0; s2 = 0;
        for k = 1:N-1
            for n = 1:N-k
                for np = 1:N-k
                    % n - snippet
                    phi_n = a_coeff(1) + a_coeff(2)*n + a2 * n^2;
                    phi_nk = a_coeff(1) + a_coeff(2)*(n+k) + a2 * (n + k)^2;
                    psi_nk = phi_n - phi_nk;
                    % np - snippet
                    phi_np = a_coeff(1) + a_coeff(2)*np + a2 * np^2;
                    phi_npk = a_coeff(1) + a_coeff(2)*(np+k) + a2 * (np + k)^2;
                    psi_npk = phi_np - phi_npk;

                    s2 = s2 + exp(-1j * psi_npk)*exp(1j * psi_nk);
                end
            end

        end
        ISL = s1+s2;
        rm = xcorr(x);
        ISLm = sum(abs(rm(N+1:2*N-1)).^2);

        fprintf('ISL computed using my logic: %f \n', ISL);
        fprintf('ISL computed using Matlab: %f \n', ISLm);

    case 't4'
%\text{ISL} = \sum_{k=1}^{N-1} \bigg[ \sum_{\substack{n=1\\n=n'}}^{N-k} e^{j0} +  \sum_{\substack{n,n'=1\\ n\neq n'}}^{N-k} ( e^{j\psi_{knn'}} + e^{-j\psi_{knn'}})\bigg]
        s1 = 0; s2 = 0;
        for k = 1:N-1
            s1 = s1 + (N-k);
            for n = 1:N-k
                for np = 1:N-k
                    if n < np
                        psi_knnp = 2 * a2 * k * (n - np);
                        s2 = s2 + 2*cos(psi_knnp);
                    else
                        continue
                    end
                end
            end
        end
        ISL = s1+s2;
        rm = xcorr(x);
        ISLm = sum(abs(rm(N+1:2*N-1)).^2);

        fprintf('ISL computed using my logic: %f \n', ISL);
        fprintf('ISL computed using Matlab: %f \n', ISLm);
    otherwise
        0;
end
return
%% MM iterations
a2 = a_coeff(3);  % Extract initial a2
a2_i = a2;
normalized_ISL = zeros(maxIter, 1);
step_size = 0.0001; % or smaller if unstable
for iter = 1:maxIter
    ISL = 0;
    r = zeros(N-1,1);
    for k = 1:N-1
        for n = 1:N-k
            phi_n = a_coeff(1) + a_coeff(2)*n + a2 * n^2;
            phi_nk = a_coeff(1) + a_coeff(2)*(n+k) + a2 * (n + k)^2;
            psi_nk = phi_n - phi_nk;
            r(k) = r(k) + exp(1j * psi_nk);
        end
        ISL = ISL + abs(r(k))^2;
    end
    
    normalized_ISL(iter) = 20*log10(ISL/N^2);

    % ----------------------- %
    if false
        rm = xcorr(x);

        lag = 1:N-1;
        figure(1); hold all; box on;
        plot(lag, abs(r), 'r', 'LineWidth',2);
        plot(lag, abs(rm(N+1:2*N-1)), '--b*', 'LineWidth',2);

        ISLm = sum(abs(rm(N+1:2*N-1)).^2);

        fprintf('ISL computed using my logic: %f \n', ISL);
        fprintf('ISL computed using Matlab: %f \n', ISLm);
    end
    % ----------------------- %
    % MM Update step for each (k, n, n')
    delta = 0;
    for k = 1:N-1
        for n = 1:N-k
            for np = 1:N-k
                if n < np
                    psi_knnp = 2 * a2 * k * (n - np);
%                     delta = delta + (1/(2*k*(n - np)))*sin(psi_knnp);
                    delta = delta + (1/(2*k*(n - np)))*tan(psi_knnp);
                else
                    continue
                end
            end
        end
    end
    a2u = a2 + step_size*delta;
    a2 = a2u;
end

%% Plot normalized ISL convergence
figure;
plot(normalized_ISL, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Iteration');
ylabel('Normalized ISL');
title('Convergence of Normalized ISL via MM');

%% Final result
fprintf('Final a2: %.6f\n', a2);

%% Support functions
function [x, phi_x] = create_sequence(a_coeff, M, L)
phi_i = [];
m = (1:M).';
for l = 1:L
    phi_i(:,l) = a_coeff(1,l) + a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
phi_x = phi_i(:);
x = exp(1j*phi_x);
end


%     if iter > 1 && abs(normalized_ISL(iter) - normalized_ISL(iter-1)) < tol
%     if iter > 1 && norm(abs(a2_i) - abs(a2)) < tol
%         fprintf('Converged at iteration %d\n', iter);
%         normalized_ISL = normalized_ISL(1:iter);
%         break;
%     end