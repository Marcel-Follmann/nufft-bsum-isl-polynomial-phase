clc;clear;close all;

% Parmeters
MaxIter = 1e4;

% Initialization
nSamples = 100;
polyOrder = 4;
pthNorm = 2;

tic;
[x0, params] = rndmSqnce_phaseEstmt(nSamples, polyOrder);
[x_PSL, obj_hist, PSL_hist, ISL_hist] = MPSL_PECS(x0, nSamples, polyOrder, params, pthNorm);
toc;

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
chk = xcorr(x0)/nSamples;
chkPly = 20*log10(abs(chk));

ac_MM_ply = xcorr(x_PSL);
corrPS = 20*log10(abs(ac_MM_ply./ac_MM_ply(nSamples)));

figure("Name","Correlation"); hold on;
plot(1-nSamples:nSamples-1,chkPly, 'linewidth',1.5);
plot(1-nSamples:nSamples-1,corrPS,  'LineWidth',1.5);
ylim([-60, 1]);
xlabel('k')
ylabel('correlation level (dB)')
legend('Seed Sqnce','PECS - PS')

% --------------------- Phase Plot -------------------------------
argX_ply = unwrap(atan2(imag(x0),real(x0)));
argPly_PS = unwrap(atan2(imag(x_PSL),real(x_PSL)));

figure("Name","Phase"); box on; hold on; 
indx = 1:nSamples;
plot(indx, argPly_PS, '-.o', 'linewidth', 1.5);
xlabel('Code length');ylabel('Phase');
hold off;

figure("Name", "BandWidth"); hold on;
nFFT = 256;fIndx = linspace(-0.5,0.5,nFFT);
plot(fIndx, fftshift(20*log10(abs(fft(x_PSL,nFFT))/nFFT)),'linewidth',2);
xlabel('Normalized Frequency'); ylabel('Power [dB]');

mCode = PhaseCodeType.Golomb;
x_glmb = mCode.generate(nSamples);
[~,isl_glmb] = pslisl(x_glmb);
fprintf('ISL of Golomb sequence for length N -> %d : %f\n',nSamples, 10*log10(isl_glmb));
[~,isl_x_PSL] = pslisl(x_PSL);
fprintf('ISL of the optimized sequence for length N -> %d : %f\n',nSamples, 10*log10(isl_x_PSL));

% ------------------- Ambig Fn Plot ---------------------------
AmbigFn_heat(x_PSL);


%% ========================= * ==========================
% Lp Norm increment with increasing iterations
function [x,objVec,pslVec_dB, islVec_dB] = MPSL_PECS(x0, nSamples, polyOrder, params, pthNorm, options)
    arguments
        x0 
        nSamples
        polyOrder
        params
        pthNorm
        options.MaxIter = 1e4
    end

    objVec = NaN * ones(options.MaxIter,1);
    pslVec_dB = NaN * ones(options.MaxIter,1);
    islVec_dB = NaN * ones(options.MaxIter,1);

    r = abs(xcorr(x0)/nSamples);
    objVec(1) = norm(r(1:nSamples),pthNorm);
    [pslVec_dB(1), islVec_dB(1)] = pslisl(x0, "isLogScale", true);
    
    x = x0;
    for iter = 1:options.MaxIter-1
        % Optimize
        xU = MPSL_update(x,pthNorm,nSamples);
        [x, params] = PECS_PerfSqr(xU, params, nSamples, polyOrder);

        % Eval
        r = abs(xcorr(x)/nSamples);
        objVec(iter+1) = norm(r(1:nSamples),pthNorm);
        [pslVec_dB(iter+1),islVec_dB(iter+1)] = pslisl(x, "isLogScale", true); 
    end  
end

%% PSL Minimization
function [x_update] = MPSL_update(x,pthNorm,nSamples)
    % One step update of MPSL_eig
    f = fft(x,2*nSamples);
    r = ifft(abs(f).^2); % r = [r0,r1,...,r_N-1,0,r_1-N,...,r_-1]
    r_abs = abs(r(2:nSamples));
    
    t = norm(r_abs,pthNorm)+1e-5;
    % normalized version (normalize by t^p to avoid numerical issue)
    r_t_ratio = r_abs/t;
    a = (1 + (pthNorm-1)*r_t_ratio.^pthNorm - pthNorm*r_t_ratio.^(pthNorm-1))./(t-r_abs).^2;
    w = 0.5*pthNorm/t^2*r_t_ratio.^(pthNorm-2);
    
    lambda_max_L = max((nSamples-1:-1:1)'.*a); % the maximum eigenvalue of L
    
    c = [0;w;0;w(nSamples-1:-1:1)].*r;
    mu = real(fft(c));
    lambda_u = 0.5*(max(mu(1:2:end))+max(mu(2:2:end)));
    
    y = ifft(mu.*f); % the fator 1/2N is implicitly included
    
    x_update = (nSamples*lambda_max_L + lambda_u)*x - y(1:nSamples);
end

%% PECS - Perfect Square Implementation
function [z, params] = PECS_PerfSqr(y, params, nSamples, polyOrder)
    
    z = [];
    psi = angle(y);
    r = abs(y);

    indx_e = 0;
    indx_s = indx_e+1;
    indx_e = indx_s + nSamples-1;

    psi_M = psi(indx_s:indx_e);
    r_M   = r(indx_s:indx_e);

    m = 1:nSamples; m = m(:);

    arg_init = 0;
    for q = 1:polyOrder
        arg_init = arg_init + params(q).*(m.^(q-1));
    end

    theta = arg_init - psi_M;

    % =========== Taylor subroutine starts ================

    am = r_M / 2 .* cos(theta);
    bm = r_M .* sin(theta);
    beta = - 2 * am .* psi_M - 2 * am .* theta + bm;
    alpha = [];
    for q = 0 : polyOrder-1
        alpha = [alpha, 2 *  am .* ((1:nSamples).').^q];
    end

    [params, ~] = lsqr(alpha,-beta, 1e-3);

    % =========== Taylor subroutine ends ===============

    arg_opt = 0;
    for q = 1:polyOrder
        arg_opt =  arg_opt + params(q).*(m.^(q-1));
    end
    x_j = exp(1j*arg_opt);

    z = [z; x_j];

end



function [ psl,isl ] = pslisl( x , options )
    % This function calculates PSL and ISL for a radar transmit sequence 
    arguments
        x 
        options.isNormalized = false
        options.isAperiodic = true
        options.isLogScale = false
    end

    x           = x(:);
    N           = length(x);   

    if options.isAperiodic
        x       = abs(xcorr(x));
    else
        x       = abs(pcorr(x));
    end

    if options.isNormalized
        x       = x / max(x);
    end

    psl         = max(x(1:N-1));
    isl         = sum(x(1:N-1).^2);
    
    if options.isLogScale
        psl      = 10*log10(psl^2/N^2);
        isl      = 10*log10(isl/N^2);
    end
end



function [x, a_ql] = rndmSqnce_phaseEstmt(nSamples, polyOrder)
    % This function is designed to keep variable sub-sequence length - M
    % only choice required is No of sub-sequences - L, minimum sub-sequence length - M_min, Q, sequence length - N   
    rng 'shuffle'
        
    mCode = PhaseCodeType.Random;
    y = mCode.generate(nSamples);
    
    psi = angle(y);
    a = 10*rand(polyOrder,1);
        
    m(:,1) = 1:nSamples;
    psi_M = psi(1:nSamples);

    m_s = [ones(size(m))];
    for q = 2:polyOrder
        m_s = [m.^(q-1), m_s];
    end
    
    A_ls = m_s;
    [vrbl,~] = lsqr(A_ls, psi_M);
    
    for q =1:polyOrder
        a(polyOrder+1-q) = vrbl(q);
    end
    
    arg_opt = 0;
    for q = 1:polyOrder
        arg_opt =  arg_opt + a(q)*m.^(q-1);
    end
    x = exp(1j*arg_opt);
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