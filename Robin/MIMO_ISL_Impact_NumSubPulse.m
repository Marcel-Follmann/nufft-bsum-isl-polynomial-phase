% clc; clear; close all;
rng(1);
N = 128;
L = 8;
M = round(N/L);

lags = -N+1:N-1;

[sqnce] = PECS_subpls(M, L);

r = xcorr(sqnce);

r_abs_sq = abs(r).^2./N^2;

figure; box on;
plot(lags, 10*log10(r_abs_sq));
xlabel('lags')
ylabel('mag')

ISL = sum(r_abs_sq(1:N-1)) + sum(r_abs_sq(N+1:end));
PSL = max(max(r_abs_sq(1:N-1)), max(r_abs_sq(N+1:end)));

txt1 = sprintf('ISL of the sequence: %f', 10*log10(ISL)); disp(txt1);
txt2 = sprintf('PSL of the sequence: %f', 10*log10(PSL)); disp(txt2);


%% Support functions
function [sqnce] = PECS_subpls(M, L)
Q = 3;
theta_li = [];
m = (1:M).';
comm_sbpls = dec2bin(0:7) - '0';
comm_sbpls = comm_sbpls';
comm_sbpls = comm_sbpls([3 2 1],:);
a_coeff = zeros(3,L);
% ------------- Debug routine start ---------------- %
% for l = 1:L
%     % ----------------- %
%     if mod(l,2) == 0
% %         a_coeff(:,l) = 1e10.*comm_sbpls(:,randi(8)-1);
%         a_coeff(:,l) = 1e10.*comm_sbpls(:,l);
%         if a_coeff(3,l) == 0
%             a_coeff(3,l) = -1e10;
%         else
%             a_coeff(3,l) = +1e10;
%         end
%     else
%         a_coeff(:,l) = 1e10*rand(Q,1);
%     end
% end
% ------------- Debug routine end ---------------- %
% ============================= %
for l = 1:L
    % ----------------- %
    if mod(l,2) == 0
        a_coeff(:,l) = 1e10.*comm_sbpls(:,l);
        if a_coeff(3,l) == 0
            a_coeff(3,l) = -1e10;
        else
            a_coeff(3,l) = +1e10;
        end
    else
        a_coeff(:,l) = 1e10*rand(Q,1);
    end
    % ----------------- %
    theta_li(:,l) = a_coeff(1,l) + a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
% ============================= %
theta_x = theta_li(:);
x = exp(1j*theta_x);
sqnce = x;
end

