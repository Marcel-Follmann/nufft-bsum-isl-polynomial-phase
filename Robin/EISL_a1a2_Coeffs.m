clc; clear; close all;

N = 8;

%% Expectation of Cross-correlation of a communication subpulse with a sensing subpulse
s1 = 0;
% k=0; 
% b1 = 10; b2 = 100;
b1 = 10*rand;
b2 = 10*rand;
for k = 1-N:N-1
    for n = 1:N-k
        for n_p = 1:N-k
            n_bar1s = n_p - n;
            n_bar1a = n_p + n;
            n_bar2 = n_p^2 - n^2;
            t1 = (exp(1j*n_bar1s*pi) + exp(-1j*n_bar1s*pi))/2;
            t2 = (exp(1j*n_bar2*pi) + exp(-1j*n_bar2*pi))/2;
            C1 = exp(1j*n_bar1s*(-b1 - b2*(n_bar1a) - 2*b2*k)/pi);
            s1 = s1 + t1*t2*C1;
        end
    end
end

%% Expectation of Cross-correlation of a communication subpulse with another communication subpulse
s2 = 0;
for k = 1-N:N-1
    for n = 1:N-k
        for n_p = 1:N-k
            n_bar1s = n_p - n;
            n_bar1a = n_p + n;
            n_bar2 = n_p^2 - n^2;
            t11 = (exp(1j*n_bar1s*pi) + exp(-1j*n_bar1s*pi))/2;
            t22 = (exp(1j*n_bar2*pi) + exp(-1j*n_bar2*pi))/2;
            t33 = (exp(1j*n_bar1s*pi) + exp(-1j*n_bar1s*pi))/2;
            t44 = (exp(1j*n_bar1s*(n_bar1a+2*k)*pi) + exp(-1j*n_bar1s*(n_bar1a+2*k)*pi))/2;
            s2 = s2 + t11*t22*t33*t44;
        end
    end
end

%% Temporary snippet
% indx = -10*pi:pi:10*pi;
% figure;
% plot(indx, cos(indx))
% xlabel(indx)
% ylabel('$\cos(\theta)$')

% %% a1 coefficient Expectation
% s1 = 0;
% for k = 1-N:N-1
%     for n = 1:N-k
%         for n_p = 1:N-k
%             n_bar1 = n_p - n;
%             t1 = (exp(1j*n_bar1*pi) + exp(-1j*n_bar1*pi))/2;
%             s1 = s1 + t1;
%         end
%     end
% end
% 
% 
% %% a2 coefficient Expectation
% s2 = 0;
% for k = 1-N:N-1
%     for n = 1:N-k
%         for n_p = 1:N-k
%             n_bar2 = n_p^2 - n^2;
%             t2 = (exp(1j*n_bar2*pi) + exp(-1j*n_bar2*pi))/2;
%             s2 = s2 + t2;
%         end
%     end
% end