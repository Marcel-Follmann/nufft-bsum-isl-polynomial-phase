classdef PhaseCodeType
    enumeration
        Frank
        Golomb
        P1
        P2
        PX
        P4
        ZadoffChu
        Random
        Barker
        MSeq
        Gold
        Kasami
        AllOne
    end
    
    methods (Access = public)
        function signal = generate(type, nSamples)
            arguments
                type (1,1) Algorithm.internal.PhaseCodeType
                nSamples (1,1) double {mustBePositive}
            end
            switch type
                case Algorithm.internal.PhaseCodeType.Frank
                    signal = type.genFrank(nSamples);
                case Algorithm.internal.PhaseCodeType.Golomb
                    signal = type.genGolomb(nSamples);
                case Algorithm.internal.PhaseCodeType.P1
                    signal = type.genP1(nSamples);
                case Algorithm.internal.PhaseCodeType.P2
                    signal = type.genP2(nSamples);
                case Algorithm.internal.PhaseCodeType.PX
                    signal = type.genPX(nSamples);
                case Algorithm.internal.PhaseCodeType.P4
                    signal = type.genP4(nSamples);
                case Algorithm.internal.PhaseCodeType.ZadoffChu
                    signal = type.genZadoffChu(nSamples);
                case Algorithm.internal.PhaseCodeType.Random
                    signal = type.genRandom(nSamples);
                case Algorithm.internal.PhaseCodeType.Barker
                    signal = type.genBarker(nSamples);
                case Algorithm.internal.PhaseCodeType.MSeq
                    signal = type.genMSeq(nSamples);
                case Algorithm.internal.PhaseCodeType.Gold
                    signal = type.genGold(nSamples);
                case Algorithm.internal.PhaseCodeType.Kasami
                    signal = type.genKasami(nSamples);
                case Algorithm.internal.PhaseCodeType.AllOne
                    signal = type.genAllOne(nSamples);

                otherwise
                    error("No generator defined for %s\n", string(type))
            end
        end
    end

    methods (Access = private, Static)
        function x = genFrank(N)
            N       = (fix(sqrt(N)))^2;
            code    = zeros(1,N);
            k       = 1;
            for i = 1 : sqrt(N)
                for j = 1 : sqrt(N)
                    code(k) = exp(1i * 2 * pi /sqrt(N) * (i-1) * (j-1) );
                    k = k  + 1;
                end
            end
            x       = code(:);
        end

        function x = genGolomb(N)
            alpha = exp(1i * 2*pi /N);
            exponent = (0:(N-1))' .* (1:N)' / 2;
            x = alpha.^exponent;
        end

        function x = genP1(N)
            L   = sqrt(N);
            phi = zeros(L,L);
            for k = 1 : L
                for n = 1: L
                    phi(n,k) = 2*pi/L *((L+1)/2-n)*((n-1)*L+(k-1));
                end
            end
            phi = phi(:);
            x   = exp(1i * phi);
        end

        function x = genP2(N)
            L   = sqrt(N);
            phi = zeros(L,L);
            for k = 1 : L
                for n = 1: L
                    phi(n,k) = 2*pi/L *((L+1)/2-k)*((L+1)/2-n);
                end
            end
            phi = phi(:);
            x   = exp(1i * phi);
        end

        function x = genPX(N)
            L   = sqrt(N);
            phi = zeros(L,L);
            if mod(L,2) == 0
                for k = 1 : L
                    for n = 1: L
                        phi(n,k) = 2*pi/L *((L+1)/2-k)*((L+1)/2-n);
                    end
                end
            else
                for k = 1 : L
                    for n = 1: L
                        phi(n,k) = 2*pi/L *((L)/2-k)*((L+1)/2-n);
                    end
                end
            end
            phi = phi(:);
            x   = exp(1i * phi);
        end

        function x = genP4(N)
            n   = (1:N)';
            phi = 2*pi/N * (n-1) .* (n-1-N)/2;
            x   = exp(1i * phi);
        end

        function x = genZadoffChu(N)
            if mod(N,2)==0
                x = exp(1i * pi * (1:N).^2 /N);
            else
                x = exp(1i * pi * (2:(N+1)) .* (1:N) /N);
            end
        end

        function x = genRandom(N)
            phi    = (2*rand(N,1)-1)*pi;
            x      = exp(1i * phi);
        end

        function x = genBarker(N)
            arguments
                N (1,1) double {mustBeBetween(N,2,13)}
            end
            switch N
                case 2
                    x = [-1 1];
                case 3
                    x = [-1 -1 1];
                case 4
                    x = [-1 -1 1 -1];
                case 5
                    x = [-1 -1 -1 1 -1];
                case 7
                    x = [-1 -1 -1 1 1 -1 1];
                case 11
                    x = [-1 -1 -1 1 1 1 -1 1 1 -1 1];
                case 13
                    x = [-1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1] ;
            end
        end

        function x = genMSeq(N)
            x   = mseq(nextpow2(N),init);
        end

        function x = genGold(N)
            if  mod(nextpow2(N),2)~= 1 && mod(nextpow2(N),4)~=2
                disp('order of the sequence set must satisfying mod(n,2)=1 or mod(n,4)=2');
                return;
            end
            Y   = gold(nextpow2(N));
            x   = Y(:,1);
        end

        function x = genKasami(N)
            if  mod(nextpow2(N),2)~= 0
                disp('order of the sequence set, an even integer.');
                return;
            end
            
            Y   = kasami(nextpow2(N));
            x   = Y(:,1);
        end

        function x = genAllOne(N)
            x   = exp(1i * zeros(N,1));
        end
    end
    
end