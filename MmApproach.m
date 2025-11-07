classdef MmApproach
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Q, N
        D, Delta
    end

    %% Public Methods
    methods (Access = public)
        function obj = MmApproach(Q, N)
            obj.Q = Q; 
            obj.N = N;

            obj.D = obj.PrecalcD();
        end

        function x = SimulateSignal(obj, c)
            n(:,1) = 0:obj.N-1; 
            q(:,1) = 0:obj.Q-1;
            phase = c.'*Polynom(q, n.');
            x = exp(1j*phase);
        end

        function acf = SimulateAcf(obj, c) 
            x = obj.SimulateSignal(c);
            acf = abs(xcorr(x)) / obj.N; 
        end

        function isl = SimulateIsl(obj, c)
            acf = obj.SimulateAcf(c); 
            isl = sum(acf(1:obj.N)) / obj.N;
        end

        function acf = DeriveAcf(obj, c)
            phase = sum(c.*obj.D, 1);
            acfHalf = abs(sum(exp(1j*phase), 3)) / obj.N;
            acf = [acfHalf, 1, flip(acfHalf)];
        end


    end

    %% Visualization Methods
    methods (Access = public)
        function fig = VisualizeSignal(obj, c)
            x = obj.SimulateSignal(c);
            
            fig = figure("Name","Signal");
            plot(real(x));
            grid on; grid minor; 
            xlabel("Samples (-)");
            ylabel("Amplitude (-)");
        end

        function fig = VisualizeAcf(obj, c)
            acfSimu = obj.SimulateAcf(c);
            acfDerived = obj.DeriveAcf(c);
            lags = -obj.N+1:obj.N-1; 
            
            fig = figure("Name", "Acf"); hold on; 
            plot(lags, 20*log10(acfSimu), ...
                "DisplayName","Simulated");
            plot(lags, 20*log10(acfDerived), '--', ...
                "DisplayName","Theoretical");
            grid on; grid minor;
        legend("Location","northeast");
        end

        function PrintIsl(obj, c)
            isl = obj.SimulateIsl(c); 
            fprintf("ISL: %.2f dB\n", 10*log10(isl));
        end


    end

    %% Private Methods
    methods (Access = private)
        function D = PrecalcD(obj)
            q(:,1) = 0:obj.Q-1;
            D = zeros(obj.Q, obj.N-1, obj.N-1);

            for k = 1:obj.N-1
                for n = 1:obj.N-k
                    D(:, k, n) = Polynom(q, n+k) - Polynom(q, n);
                end
            end
        end

        function Delta = PrecalcDelta(obj)

        end

    
    end
end

function p = Polynom(q, n)
    p = n.^q;
end