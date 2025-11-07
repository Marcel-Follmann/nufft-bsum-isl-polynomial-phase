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
            acf = obj.SimulateAcf(c);
            lags = -obj.N+1:obj.N-1; 
            fig = figure("Name", "Acf");
            plot(lags, 20*log10(acf));
            grid on; grid minor;
        end

        function PrintIsl(obj, c)
            isl = obj.SimulateIsl(c); 
            fprintf("ISL: %.2f dB\n", 10*log10(isl));
        end


    end

    %% Private Methods
    methods (Access = private)
        function D = PrecalcD(obj)

        end

        function Delta = PrecalcDelta(obj)

        end

    
    end
end

function p = Polynom(q, n)
    p = n.^q;
end