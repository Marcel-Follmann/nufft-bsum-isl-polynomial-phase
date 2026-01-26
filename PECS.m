classdef PECS < handle

    properties
        nSamples
        nPolyOrder
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
                params.pNorm = 2
                params.tol = 1e-4
                params.nMaxIters = 1e4
            end

            obj.nSamples = params.N;
            obj.nPolyOrder = params.Q;
            obj.pNorm = params.pNorm;
            obj.tol = params.tol; 
            obj.nMaxIters = params.nMaxIters;

        end

        function [x0, a0] = getRandomSequence(obj)
            rng "shuffle";
            z = []

            mCode = PhaseCodeType.Random;
            signal = mCode.generate(obj.nSamples);
            phase = angle(signal);

            a = rand(obj.nPolyOrder, 1);
            m(:,1) = 1:obj.nSamples;


        end

        function mpsl(obj, x0, a0, options)
            arguments
                obj 
                x0 
                a0 
                options.logIsActive (1,1) bool = false 
            end

        end
    end
end