classdef TestMathematicalDerivations < matlab.unittest.TestCase

	methods (TestClassSetup)
		% Shared setup for the entire test class
	end

	methods (TestMethodSetup)
		% Setup for each test
	end

	methods (Test)
		% Test methods

		function Test_BinomialSum(tc)
			a = 1; 
			b = 2; 
			order = 2;

			expected = (a+b)^order;
			actual = BinomialSum(a, b, order);
			tc.assertEqual(actual, expected);
		end
	end

end