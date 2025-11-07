classdef MmApproach < handle
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
            obj.Delta = obj.PrecalcDelta();
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
            isl = 2 * sum(acf(1:obj.N-1).^2) / obj.N;
        end

        function acf = DeriveAcf(obj, c)
            phase = sum(c .* obj.D, 1);
            acfHalf = sum(exp(1j*phase),3, 'omitnan');

            acf = abs([flip(acfHalf), obj.N, acfHalf]) / obj.N;
        end

        function isl = DeriveIsl(obj, c)
            phase = sum(c .* obj.Delta, 1);
            const = obj.N*(obj.N-1)/2;
            tmp = const + 2*sum(cos(phase), 2, 'omitnan');
            isl = 2 * tmp / obj.N^3;
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
            islSimu = obj.SimulateIsl(c); 
            islDerived = obj.DeriveIsl(c);
            fprintf("ISL (Simulative):  %.2f dB\n", 10*log10(islSimu));
            fprintf("ISL (Theoretical): %.2f dB\n", 10*log10(islDerived));
        end


    end

    %% Private Methods
    methods (Access = private)
        function D = PrecalcD(obj)
            q(:,1) = 0:obj.Q-1;
            D = NaN * ones(obj.Q, obj.N-1, obj.N-1);

            for k = 1:obj.N-1
                for n = 1:obj.N-k
                    D(:, k, n) = Polynom(q, n+k) - Polynom(q, n);
                end
            end
        end

        function Delta = PrecalcDelta(obj)
            q(:,1) = 0:obj.Q-1;
            Delta = NaN + ones(obj.Q, obj.N*(obj.N-1)*(2*obj.N-1)/3);

            i = 1;
            for k = 1:obj.N-1
                for n = 1:obj.N-k
                    for d = 1:obj.N-k-n
                        Delta(:,i) = Polynom(q, n+k) - Polynom(q,n) - ...
                            (Polynom(q, n+d+k) - Polynom(q, n+d));
                        i = i+1;
                    end
                end
            end
        end

    
    end
end

function p = Polynom(q, n)
    p = n.^q;
end