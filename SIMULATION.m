classdef SIMULATION < handle

	properties (GetAccess = public, SetAccess = private)
		
	end

	methods (Access = public)
		function obj = SIMULATION(params)
			arguments
				params = 0
			end

		end

		function run(obj)
			
			obj.runComplete = true; 
		end
	end

	methods (Access = private)
		
	end

	methods (Static)

	end

	% Plot Methods
	methods (Access = public)
		function fig = Plot_Function(obj)
			obj.CheckRunComplete();

			fig = figure("Name","Some Figure");
		end
	end

	% Internal parameters and functions that have nothing to do with the imulation itself
	properties (Access = private)
		runComplete = false
	end

	methods (Access = private)
		function CheckRunComplete(obj)
			if ~obj.runComplete
				error("Simulation did not run jet. Nothing to plot.");
			end
		end
	end
end