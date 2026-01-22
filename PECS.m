classdef PECS < handle

    properties
        N
        Q
        M
        L
        pNorm
        tol
        nMaxIters
    end

    properties (SetAccess = protected)
        
    end

    methods
        function obj = PECS(params)
            arguments
                params.N
                params.Q
                params.M = 1
                params.pNorm = 2
                params.tol = 1e-4
                params.nMaxIters = 1e4
            end

            obj.N = params.N;
            obj.Q = params.Q;
            obj.M = params.M;
            obj.L = round(obj.N/obj.M); 
            obj.pNorm = params.pNorm;
            obj.tol = params.tol; 
            obj.nMaxIters = params.nMaxIters;

        end

        function [x, refVec, a_ql] = randomSequencePhaseEstimation(nSamples, L, M_min, Q, isRandomSubsequenceLength)
            mCode = PhaseCodeType.Random;
            signal = mCode.generate(N);
            phase = angle(signal);


        end
    end
end