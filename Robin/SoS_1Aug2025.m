clc; clear; close all;
% addpath('C:\Program Files\Mosek\10.2\toolbox\r2017aom')
%% Key points to note in this simulation:
% 1) variable: "epsln" is a weight factor for "Gk_dot" and "Gk_ddot"
% (a) When [epsln = 1e-3] - the objective oscillates around the mean value with
% increasing iterations
% (b) When [epsln = 1] - the objective oscillates around the mean value with
% increasing iterations
% 2) variable: "bias" is also a weight factor.
% a) [a = 1] optimization works fine 


%% Sequence Set Design parameters
N = 32;
L = 2;
M = N/L;
% rng(1);
rng('shuffle');

% Weighting Factors for tuning optimization
% epsln = 1;
% epsln = 1e-3;  
epsln = 1e-2; 

% bias = 1; % 1e30; % bias = 1e-10;
bias = 1;

Q = 3;
max_itrns = 5;
maxitrISL = 10;
theta_li = [];
m = (1:M).';
a_coeff = zeros(3,L);
for l = 1:L
    a_coeff(:,l) = 2*pi*1e3*rand(Q,1);
    theta_li(:,l) = a_coeff(1,l)+ a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
theta_x = theta_li(:);
x = exp(1j*theta_x);

theta_y = 2*pi*1e3*rand(N,1);
xg_prime = exp(1j*theta_y);

l = 1;
a = a_coeff(:,l);
indcs = (1+(l-1)*M: l*M)';
xg = zeros(N,1);
xg(indcs) = x(indcs);

x_g = x;
x_g(indcs) = 0;

IN = eye(N);
c_dB = [];
ck = abs(xcorr(x,xg_prime)).^2;
c_dB(end+1) = 10*log10(sum(ck));



a_track = zeros(length(a_coeff(:,l)), max_itrns);
rng('shuffle');

obj_cnv = [];
%% Trigonometric (cos/sin) majorizer loop
for itr = 1:max_itrns
    coeff = zeros(35,2*N-1);
    i = 1;
    obj_sum = 0;
    acc_flc1 = []; acc_flc2 = [];
    obj_k = 0;
    %% Specific subsequence and lag design parameters
    for k = 1-N:N-1

        Uk = toeplitzC(N,k);

%         alpha_k = Uk*xg_prime;
        tmp = Uk*xg_prime;
        alpha_k = zeros(N,1);
        alpha_k(indcs) = tmp(indcs);

        gamma_k = x_g' * Uk * xg_prime;

        psi_k = sqrt(real(alpha_k).^2 + imag(alpha_k).^2);
        %     psi_k(psi_k == 0) = 1e-10;

        A_dot_mk = [];
        A_ddot_mk = [];
        A_dot_k_blkdg = [];
        A_ddot_k_blkdg = [];

        theta_l = theta_x(indcs);
        
        I = eye(N+1);
        D = [];
        zo = zeros(1,N+1); zo(end) = 1;
    
        % Subsequence elements 
        for m = 1:N

            if psi_k(m) ~= 0
                beta_mk = acos(real(alpha_k(m))./psi_k(m));

                C = theta_l(m) + beta_mk;

                theta_lk(m) = theta_l(m) - beta_mk;
                eta_1k = 0.5*psi_k(m);
                eta_2k = psi_k(m).*(-C - sin(theta_lk(m)));
                eta_3k = psi_k(m).*(0.5*C.^2 + C.*sin(theta_lk(m)) + cos(theta_lk(m))) + real(gamma_k)./M;


                eta_4k = 0.5*psi_k(m);
                eta_5k = psi_k(m).*(-C + cos(theta_lk(m)));
                eta_6k = psi_k(m).*(0.5*C.^2 - C.*cos(theta_lk(m)) + sin(theta_lk(m))) + imag(gamma_k)./M;

            else
                eta_1k = 0;
                eta_2k = 0;
                eta_3k =  real(gamma_k)./M;


                eta_4k = 0;
                eta_5k = 0;
                eta_6k = imag(gamma_k)./M;
            end

            A_dot_mk(:,:,m) = [eta_1k, 0.5*eta_2k; 0.5*eta_2k, eta_3k] ;

            A_ddot_mk(:,:,m) = [eta_4k, 0.5*eta_5k; 0.5*eta_5k, eta_6k] ;

            A_dot_k_blkdg = blkdiag(A_dot_k_blkdg, squeeze(A_dot_mk(:,:,m)));
            A_ddot_k_blkdg = blkdiag(A_ddot_k_blkdg, squeeze(A_ddot_mk(:,:,m)));

            D = [D; I(m,:); zo];

        end

        M_mtrx = [];
        pq = (1:M).';
        for q = 1:Q
            M_mtrx = [M_mtrx, pq.^(q-1)];
        end

        % 1st Term
        B_dot_k = D.'*A_dot_k_blkdg*D;
        B_1k = B_dot_k(1:M,1:M);
        b_1k = B_dot_k(1:M,M+1);
        vec_b_k_dot = b_1k.'*M_mtrx;
        sclr_b_k_dot = B_dot_k(M+1,M+1);
        Gk_dot = M_mtrx.'*B_1k*M_mtrx + epsln*eye(3);

        % 2nd Term
        B_ddot_k = D.'*A_ddot_k_blkdg*D;
        B_2k = B_ddot_k(1:M,1:M);
        b_2k = B_ddot_k(1:M,M+1);
        vec_b_k_ddot = b_2k.'*M_mtrx;
        sclr_b_k_ddot = B_ddot_k(M+1,M+1);
        Gk_ddot = M_mtrx.'*B_2k*M_mtrx + epsln*eye(3);

        [coeff1] = coeff_allcn(Gk_dot, vec_b_k_dot, sclr_b_k_dot);

        [coeff2] = coeff_allcn(Gk_ddot, vec_b_k_ddot, sclr_b_k_ddot);

        %% --------------------------------------------------- %%       
        coeff(:,i) = coeff1 + coeff2;  i = i+1;
        obj_k = obj_k + ( (a'*Gk_dot*a + 2*vec_b_k_dot*a + sclr_b_k_dot)^2     +   (a'*Gk_ddot*a + 2*vec_b_k_ddot*a + sclr_b_k_ddot)^2 );
    end

    obj_cnv(end+1) = obj_k;

s_coeff = bias.*sum(coeff,2);


[v] = func_cvx_solver(s_coeff);
a = v(2:4);

pq = (1:M).';
theta_l = a(1) + a(2).*pq + a(3).*pq.^2;

theta_x(indcs) = theta_l;
x = exp(1j*theta_x);
a_track(:,itr) = a;
end
figure; 
% plot(obj_cnv);
semilogy(obj_cnv);
xlabel('majorizer iterations')
ylabel('obj conv')

% ISL Convergence expected in this loop
ck = abs(xcorr(x,xg_prime)).^2;
c_dB(end+1) = 10*log10(sum(ck));

% disp(v)
disp('Simulation complete')

figure; plot(c_dB);
xlabel('iterations')
ylabel('ISL')

%% -----------------------------------------  %
% ---------- Support functions -------------- %
%  -----------------------------------------  %

function [v] = func_cvx_solver(coeff)

% Define the objective function
% cvx_begin sdp
cvx_begin sdp
% cvx_precision best
variable v(35,1);

V = [v(1),  v(2),   v(3),   v(4),   v(5),   v(6),    v(7),   v(8),   v(9),   v(10);
    v(2),  v(5),   v(6),   v(7),   v(11),  v(12),   v(13),  v(14),  v(15),  v(16);
    v(3),  v(6),   v(8),   v(9),   v(12),  v(14),   v(15),  v(17),  v(18),  v(19);
    v(4),  v(7),   v(9),   v(10),  v(13),  v(15),   v(16),  v(18),  v(19),  v(20);
    v(5),  v(11),  v(12),  v(13),  v(21),  v(22),   v(23),  v(24),  v(25),  v(26);
    v(6),  v(12),  v(14),  v(15),  v(22),  v(24),   v(25),  v(34),  v(28),  v(29);
    v(7),  v(13),  v(15),  v(16),  v(23),  v(25),   v(26),  v(28),  v(29),  v(27);
    v(8),  v(14),  v(17),  v(18),  v(24),  v(34),   v(28),  v(30),  v(31),  v(32);
    v(9),  v(15),  v(18),  v(19),  v(25),  v(28),   v(29),  v(31),  v(32),  v(33);
    v(10), v(16),  v(19),  v(20),  v(26),  v(29),   v(27),  v(32),  v(33),  v(35);];

obj1 = coeff(1)*v(1) + coeff(2)*v(2) + coeff(3)*v(3) + coeff(4)*v(4);
obj2 = coeff(5)*v(5) + coeff(6)*v(6) + coeff(7)*v(7) + coeff(8)*v(8) + coeff(9)*v(9) + coeff(10)*v(10);
obj3 = coeff(11)*v(11) + coeff(12)*v(12) + coeff(13)*v(13) + coeff(14)*v(14) + coeff(15)*v(15) ...
    + coeff(16)*v(16) + coeff(17)*v(17) + coeff(18)*v(18) + coeff(19)*v(19) + coeff(20)*v(20);
obj4 = coeff(21)*v(21) + coeff(22)*v(22) + coeff(23)*v(23) + coeff(24)*v(24) + coeff(25)*v(25) + coeff(26)*v(26) ...
    + coeff(27)*v(27) + coeff(28)*v(28) + coeff(29)*v(29) + coeff(30)*v(30) + coeff(31)*v(31) + coeff(32)*v(32) ...
    + coeff(33)*v(33) + coeff(34)*v(34) + coeff(35)*v(35);
obj = obj1 + obj2 + obj3 + obj4;
minimize(obj);

subject to
V >= 0;
v(1) == 1;

cvx_end

end


%  -----------------------------------------  %
function [coeff] = coeff_allcn(Gk_prime, vec_bk_prime, sclr_bk_prime)
coeff = zeros(35,1);

% Aliasing variables as per the manuscript
b0 = vec_bk_prime(1);
b1 = vec_bk_prime(2);
b2 = vec_bk_prime(3);

% Coefficient allocation
% Degree 0 & 1
coeff(1)      = sclr_bk_prime^2;
coeff(2)      = 4 * b0 * sclr_bk_prime;
coeff(3)      = 4 * b1 * sclr_bk_prime;
coeff(4)      = 4 * b2 * sclr_bk_prime;
% Degree 2
coeff(5)      = 4 * b0^2 + 2*Gk_prime(1,1)*sclr_bk_prime;
coeff(6)      = 2*Gk_prime(2,1)*sclr_bk_prime + 2*Gk_prime(1,2)*sclr_bk_prime + 8*b0*b1;
coeff(7)      = 2*Gk_prime(3,1)*sclr_bk_prime + 2*Gk_prime(1,3)*sclr_bk_prime + 8*b0*b2;
coeff(8)      = 4 * b1^2 + 2*Gk_prime(2,2)*sclr_bk_prime;
coeff(9)      = 2*Gk_prime(3,2)*sclr_bk_prime + 2*Gk_prime(2,3)*sclr_bk_prime + 8*b1*b2;
coeff(10)     = 4 * b2^2 + 2*Gk_prime(3,3)*sclr_bk_prime;
% Degree 3
coeff(11)     = 4*Gk_prime(1,1)*b0;
coeff(12)     = 4*Gk_prime(1,1)*b1 + 4*Gk_prime(2,1)*b0 + 4*Gk_prime(1,2)*b0;
coeff(13)     = 4*Gk_prime(1,1)*b2 + 4*Gk_prime(3,1)*b0 + 4*Gk_prime(1,3)*b0;
coeff(14)     = 4*Gk_prime(2,2)*b0 + 4*Gk_prime(2,1)*b1 + 4*Gk_prime(1,2)*b1;
coeff(15)     = 4*Gk_prime(2,1)*b2 + 4*Gk_prime(1,2)*b2 + 4*Gk_prime(3,1)*b1 + 4*Gk_prime(1,3)*b1 + 4*Gk_prime(3,2)*b0 + 4*Gk_prime(2,3)*b0;
coeff(16)     = 4*Gk_prime(3,3)*b0 + 4*Gk_prime(3,1)*b2 + 4*Gk_prime(1,3)*b2;
coeff(17)     = 4*Gk_prime(2,2)*b1;
coeff(18)     = 4*Gk_prime(2,2)*b2 + 4*Gk_prime(3,2)*b1 + 4*Gk_prime(2,3)*b1;
coeff(19)     = 4*Gk_prime(3,3)*b1 + 4*Gk_prime(3,2)*b2 + 4*Gk_prime(2,3)*b2;
coeff(20)     = 4*Gk_prime(3,3)*b2;

% Degree 4
coeff(21)     = Gk_prime(1,1)^2;
coeff(22)     = 2*Gk_prime(1,1)*Gk_prime(2,1) + 2*Gk_prime(1,1)*Gk_prime(1,2);
coeff(23)     = 2*Gk_prime(1,1)*Gk_prime(3,1) + 2*Gk_prime(1,1)*Gk_prime(1,3);
coeff(24)     = 2*Gk_prime(1,1)*Gk_prime(2,2) + Gk_prime(2,1)^2 + 2*Gk_prime(1,2)*Gk_prime(2,1) + Gk_prime(1,2)^2;
coeff(25)     = 2*Gk_prime(1,1)*Gk_prime(3,2) + 2*Gk_prime(2,1)*Gk_prime(3,1) + 2*Gk_prime(1,2)*Gk_prime(3,1) + 2*Gk_prime(1,1)*Gk_prime(2,3) + 2*Gk_prime(1,3)*Gk_prime(2,1) + 2*Gk_prime(1,2)*Gk_prime(1,3);
coeff(26)     = 2*Gk_prime(1,1)*Gk_prime(3,3) + Gk_prime(3,1)^2 + 2*Gk_prime(1,3)*Gk_prime(3,1) + Gk_prime(1,3)^2;
coeff(27)     = 2*Gk_prime(3,1)*Gk_prime(3,3) + 2*Gk_prime(1,3)*Gk_prime(3,3);
coeff(28)     = 2*Gk_prime(2,1)*Gk_prime(3,2) + 2*Gk_prime(1,2)*Gk_prime(3,2) + 2*Gk_prime(2,2)*Gk_prime(3,1) + 2*Gk_prime(2,1)*Gk_prime(2,3) + 2*Gk_prime(1,2)*Gk_prime(2,3) + 2*Gk_prime(1,3)*Gk_prime(2,2);
coeff(29)     = 2*Gk_prime(2,1)*Gk_prime(3,3) + 2*Gk_prime(1,2)*Gk_prime(3,3) + 2*Gk_prime(1,3)*Gk_prime(3,2) + 2*Gk_prime(2,3)*Gk_prime(3,1) + 2*Gk_prime(1,3)*Gk_prime(2,3) + 2*Gk_prime(3,1)*Gk_prime(3,2);
coeff(30)     = Gk_prime(2,2)^2;
coeff(31)     = 2*Gk_prime(2,2)*Gk_prime(2,3) + 2*Gk_prime(2,2)*Gk_prime(3,2);
coeff(32)     = 2*Gk_prime(2,2)*Gk_prime(3,3) + Gk_prime(3,2)^2 + 2*Gk_prime(2,3)*Gk_prime(3,2) + Gk_prime(2,3)^2;
coeff(33)     = 2*Gk_prime(3,2)*Gk_prime(3,3) + 2*Gk_prime(2,3)*Gk_prime(3,3);
coeff(34)     = 2*Gk_prime(2,1)*Gk_prime(2,2) + 2*Gk_prime(1,2)*Gk_prime(2,2);
coeff(35)     = Gk_prime(3,3)^2;

end

% Toeplitz function
function [Uk] = toeplitzC(N, k)
Uk = zeros(N);
for i = 1:N
    for j = 1:N
        if j-i == k
            Uk(i,j) = 1;
        else
            Uk(i,j) = 0;
        end
    end
end
end




function [flc1] = obj_check1(a, Gk_dot, vec_b_k_dot, sclr_b_k_dot, Gk_ddot, vec_b_k_ddot, sclr_b_k_ddot)
    flc1 = (a.'*Gk_dot*a + 2*a.'*vec_b_k_dot.' + sclr_b_k_dot)^2 + (a.'*Gk_ddot*a + 2*a.'*vec_b_k_ddot.' + sclr_b_k_ddot)^2;
end


function [flc2] = obj_check2(a, coeff)
v(1) = 1;
v(2) = a(1);
v(3) = a(2);
v(4) = a(3);
v(5) = a(1)^2;
v(6) = a(1)*a(2);
v(7) = a(1)*a(3);
v(8) = a(2)^2;
v(9) = a(2)*a(3);
v(10) = a(3)^2;
v(11) = a(1)^3;
v(12) = a(1)^2*a(2);
v(13) = a(1)^2*a(3);
v(14) = a(1)*a(2)^2;
v(15) = a(1)*a(2)*a(3);
v(16) = a(1)*a(3)^2;
v(17) = a(2)^3;
v(18) = a(2)^2*a(3);
v(19) = a(2)*a(3)^2;
v(20) = a(3)^3;
v(21) = a(1)^4;
v(22) = a(1)^3*a(2);
v(23) = a(1)^3*a(3);
v(24) = a(1)^2*a(2)^2;
v(25) = a(1)^2*a(2)*a(3);
v(26) = a(1)^2*a(3)^2;
v(27) = a(1)*a(3)^3;
v(28) = a(1)*a(2)^2*a(3);
v(29) = a(1)*a(2)*a(3)^2;
v(30) = a(2)^4;
v(31) = a(2)^3*a(3);
v(32) = a(2)^2*a(3)^2;
v(33) = a(2)*a(3)^3;
v(34) = a(1)*a(2)^3;
v(35) = a(3)^4;

%%%%%%%%%%%%%%
obj1 = coeff(1)*v(1) + coeff(2)*v(2) + coeff(3)*v(3) + coeff(4)*v(4);
obj2 = coeff(5)*v(5) + coeff(6)*v(6) + coeff(7)*v(7) + coeff(8)*v(8) + coeff(9)*v(9) + coeff(10)*v(10);
obj3 = coeff(11)*v(11) + coeff(12)*v(12) + coeff(13)*v(13) + coeff(14)*v(14) + coeff(15)*v(15) ...
    + coeff(16)*v(16) + coeff(17)*v(17) + coeff(18)*v(18) + coeff(19)*v(19) + coeff(20)*v(20);
obj4 = coeff(21)*v(21) + coeff(22)*v(22) + coeff(23)*v(23) + coeff(24)*v(24) + coeff(25)*v(25) + coeff(26)*v(26) ...
    + coeff(27)*v(27) + coeff(28)*v(28) + coeff(29)*v(29) + coeff(30)*v(30) + coeff(31)*v(31) + coeff(32)*v(32) ...
    + coeff(33)*v(33) + coeff(34)*v(34) + coeff(35)*v(35);
obj = obj1 + obj2 + obj3 + obj4;

flc2 = obj;
end