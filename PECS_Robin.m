clc;clear;close all;

% Parmeters
MaxIter = 1e4;
tol = 1e-4;

% Initialization
N = 50;
M = 1;
L = 1;
Q = 3;
p = 2;

[x_seed, a] = rndmSqnce_phaseEstmt(N, Q);
[x_PSL, obj_hist, PSL_hist, ISL_hist,infm] = MPSL_PECS(x_seed, N, Q, a, p);


figure("Name", "PSL Convergence");box on; hold on;
plot(0:2:length(PSL_hist)-1,PSL_hist(1:2:end),'linewidth',2);
xlabel('iteration')
ylabel('PSL [dB]')

figure("Name", "ISL Convergence");box on; hold on;
plot(0:2:length(ISL_hist)-1,ISL_hist(1:2:end),'linewidth',2);
xlabel('iteration')
ylabel('ISL [dB]')

figure("Name", "Objective Convergence");box on; hold on;
plot(0:2:length(obj_hist)-1,obj_hist(1:2:end),'linewidth',2);
xlabel('iteration')
ylabel('pth Norm - rk')

% ------------------- Correlation Plot ----------------------------
chk = xcorr(x_seed);
chkPly = 20*log10(abs(chk./chk(N)));

ac_MM_ply = xcorr(x_PSL);
corrPS = 20*log10(abs(ac_MM_ply./ac_MM_ply(N)));

figure("Name","Correlation"); hold on;
plot(1-N:N-1,chkPly, 'linewidth',1.5);
plot(1-N:N-1,corrPS,  'LineWidth',1.5);
ylim([-60, 1]);
xlabel('k')
ylabel('correlation level (dB)')
legend('Seed Sqnce','PECS - PS')

% --------------------- Phase Plot -------------------------------
argX_ply = unwrap(atan2(imag(x_seed),real(x_seed)));
argPly_PS = unwrap(atan2(imag(x_PSL),real(x_PSL)));

figure("Name","Phase"); box on; hold on; 
indx = 1:N;
plot(indx, argPly_PS, '-.o', 'linewidth', 1.5);
xlabel('Code length');ylabel('Phase');
hold off;

figure("Name", "BandWidth"); hold on;
nFFT = 256;fIndx = linspace(-0.5,0.5,nFFT);
plot(fIndx, fftshift(20*log10(abs(fft(x_PSL,nFFT))/nFFT)),'linewidth',2);
xlabel('Normalized Frequency'); ylabel('Power [dB]');

mCode = PhaseCodeType.Golomb;
x_glmb = mCode.generate(N);
[~,isl_glmb] = pslisl(x_glmb);
fprintf('ISL of Golomb sequence for length N -> %d : %f\n',N, 10*log10(isl_glmb));
[~,isl_x_PSL] = pslisl(x_PSL);
fprintf('ISL of the optimized sequence for length N -> %d : %f\n',N, 10*log10(isl_x_PSL));

% ------------------- Ambig Fn Plot ---------------------------
AmbigFn_heat(x_PSL);


%% ========================= * ==========================
% Lp Norm increment with increasing iterations
function [x_PSL,obj_PSL,PSL, ISL, infm] = MPSL_PECS(x0, N, Q, a, p, options)
    arguments
        x0 
        N
        Q
        a
        p
        options.MaxIter = 1e4
    end
    %%%%%%%% Output %%%%%%%%%%%
    % x_PSL : output sequence
    % num_iter: number of iterations
    % obj_PSL: objective values (p-norm) at each iteration
    % PSL: PSL value at each iteration

    x = x0;
    obj = [];
    PSL = [];
    ISL = [];
    f = fft(x,2*N);             % 2N point FFT of x
    r = ifft(abs(f).^2);        % r = [r0,r1,...,r_N-1,0,r_1-N,...,r_-1]
    r_abs = abs(r(2:N));
    obj(end+1) = norm(r_abs,p);
    PSL(end+1) = 20*log10(max(r_abs));
    
    for iter = 1:options.MaxIter
        xU = MPSL_update(x,p,N);
        % ------------------ PECS start ----------------------------
        [yF, a] = PECS_PerfSqr(xU, a, N, Q);
        argY = angle(yF);
        x1 = exp(1j*argY);
        % ------------------ PECS end ----------------------------
        x_fft = fft(x1,2*N); 
        rk = ifft(abs(x_fft).^2); 
        r_abs = abs(rk(2:N)); 
        obj3 = norm(r_abs,p);
        x = x1;

        obj(end+1) = obj3;
        PSL(end+1) = 20*log10(max(r_abs));
        [~,isl] = pslisl(x); 
        ISL(end+1)=20*log10(isl);
    end
    x_PSL = x;
    obj_PSL = obj;
    
    infm.actlItrns = iter;
end

%% PSL Minimization
function [x_update] = MPSL_update(x,p,N)
    % One step update of MPSL_eig
    f = fft(x,2*N);
    r = ifft(abs(f).^2); % r = [r0,r1,...,r_N-1,0,r_1-N,...,r_-1]
    r_abs = abs(r(2:N));
    t = norm(r_abs,p)+1e-5;
    % normalized version (normalize by t^p to avoid numerical issue)
    r_t_ratio = r_abs/t;
    a = (1 + (p-1)*r_t_ratio.^p - p*r_t_ratio.^(p-1))./(t-r_abs).^2;
    w = 0.5*p/t^2*r_t_ratio.^(p-2);
    
    lambda_max_L = max([N-1:-1:1]'.*a); % the maximum eigenvalue of L
    
    c = [0;w;0;w(N-1:-1:1)].*r;
    mu = real(fft(c));
    lambda_u = 0.5*(max(mu(1:2:end))+max(mu(2:2:end)));
    
    y = ifft(mu.*f); % the fator 1/2N is implicitly included
    
    x_update = (N*lambda_max_L + lambda_u)*x - y(1:N);
end

%% PECS - Perfect Square Implementation
function [z, a] = PECS_PerfSqr(y, a, N, Q)
    
    z = [];
    psi = unwrap(atan2(imag(y),real(y)));
    r = abs(y);

    indx_e = 0;
    indx_s = indx_e+1;
    indx_e = indx_s + N-1;

    psi_M = psi(indx_s:indx_e);
    r_M   = r(indx_s:indx_e);

    m = 1:N; m = m(:);

    arg_init = 0;
    for q = 1:Q
        arg_init = arg_init + a(q).*(m.^(q-1));
    end

    theta = arg_init - psi_M;

    % =========== Taylor subroutine starts ================

    am = r_M / 2 .* cos(theta);
    bm = r_M .* sin(theta);
    beta = - 2 * am .* psi_M - 2 * am .* theta + bm;
    alpha = [];
    for q = 0 : Q-1
        alpha = [alpha, 2 *  am .* ((1:N).').^q];
    end

    [a, ~] = lsqr(alpha,-beta, 1e-3);

    % =========== Taylor subroutine ends ===============

    arg_opt = 0;
    for q = 1:Q
        arg_opt =  arg_opt + a(q).*(m.^(q-1));
    end
    x_j = exp(1j*arg_opt);

    z = [z; x_j];

end



function [ psl,isl,PSLRdB,ISLRdB ] = pslisl( x , nF , pF )
    % This function calculates PSL and ISL for a radar transmit sequence 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Input: 
    %           x       : sequence
    %           nF      : normalization flag, default n = 0 
    %           pF      : pF = 1 -> calculating "Periodic" PSL and ISL,
    %                     otherwise, calculating "Aperiodic" PSL and ISL,
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Output: 
    %           psl     : calculated PSL value
    %           isl     : calculated ISL value
    %           PSLRdB  : calculated PSLR (dB) value
    %           ISLRdB  : calculated ISLR (dB) value
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % The input parameters can be chosed as follows:
    if nargin == 1
        nF      = 0;
        pF      = 0; 
    elseif nargin <= 2
       pF       = 0; 
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    x           = x(:);
    N           = length(x);   
    if pF == 0
        r       = abs(xcorr(x));
    else
        r       = abs(pcorr(x));
    end
    if nF == 1
        r       = r / max(r);
    end
    psl         = max(r(1:N-1));
    isl         = sum(r(1:N-1).^2);
    
    PSLRdB      = 10*log10(psl^2/N^2);
    ISLRdB      = 10*log10(isl/N^2);
end



function [x, a_ql] = rndmSqnce_phaseEstmt(N, Q)
    % This function is designed to keep variable sub-sequence length - M
    % only choice required is No of sub-sequences - L, minimum sub-sequence length - M_min, Q, sequence length - N   
    rng 'shuffle'
    
    z = [];
    
    mCode = PhaseCodeType.Random;
    y = mCode.generate(N);
    
    psi = angle(y);
    a = 10*rand(Q,1);
    indx_e = 0;
        
    m = 1:N; m = m(:);
    indx_s = indx_e+1;
    indx_e = indx_s + N-1;
    
    psi_M = psi(indx_s:indx_e);

    arg_init = 0;
    m_s = [ones(size(m))];
    for q = 1:Q
        arg_init = arg_init + a(q)*m.^(q-1);
    end
    
    for q = 2:Q
        m_s = [m.^(q-1), m_s];
    end
    
    A_ls = m_s;
    [vrbl,~] = lsqr(A_ls, psi_M);
    
    for q =1:Q
        a(Q+1-q) = vrbl(q);
    end
    
    arg_opt = 0;
    for q = 1:Q
        arg_opt =  arg_opt + a(q)*m.^(q-1);
    end
    x_j = exp(1j*arg_opt);
    z = [z; x_j];

    x = z; 
    a_ql = a;
end

function AmbigFn_heat(uinput)
    % Compute and plot the ambiguity function for a any code sequence
    % Compute the ambiguity function 
    % by utilizing the FFT through combining multiple range cuts 
    if isrow(uinput)
        code_seq = uinput; 
    else
        code_seq = uinput'; 
        fprintf(1, 'Note: Input should be a Row Vector!!\n');
    end
    % code_seq = uinput(:);
    N = size(code_seq,2); 
    tau = N; 
    stp = 1;                                 % Oversampling Ratio - stp = 1: Normal sampling / stp > 1: oversampling
    samp_num = size(code_seq,2) *stp;
    n = ceil(log(samp_num) / log(2)); 
    nfft = 2^n; 
    u(1:nfft) = 0; 
    j = 0; 
    
    for index = 1:stp:samp_num 
        j = j+1; 
        u(index:index+stp-1) = code_seq(j);
    end
    v = u; 
    delay = linspace(-tau, tau, nfft); 
    % Options: freq_del = 12 / tau /100;               % Why 12/tau/100    % freq_del = 2 / tau / 10;     
    freq_del = 1 / tau / 100;   
    j = 0; 
    vfft = fft(v,nfft); 
    % Options: for freq = -6/tau:freq_del:6/tau % for freq = -1/tau:freq_del:1/tau % for freq = -0.5:freq_del:0.5
    for freq = -0.1:freq_del:0.1
        j = j+1; 
        exf = exp(sqrt(-1) * 2. * pi * freq .* delay); 
        u_times_exf = u .* exf; 
        ufft = fft(u_times_exf,nfft); 
        prod = ufft .* conj(vfft); 
        ambig(:,j) = fftshift(abs(ifft(prod))');
    end
    
    norm_ambig = ambig ./ max(max(ambig));  % Normalization of the ambig function
    % Options: freq = -6/tau:freq_del:6/tau; % freq = -1/tau:freq_del:1/tau; % freq = -0.5:freq_del:0.5;
    freq = -0.1:freq_del:0.1; 
    

    delay = linspace(-1,1,nfft); 
    
    figure("Name","Ambiguity Function");
    surf(freq, delay, norm_ambig);
    xlabel('frequency') 
    ylabel('delay') 
    axis tight;
    shading interp
    % colormap hot;
    % brighten(colormap,0.7)
    colorbar;
    box on
end