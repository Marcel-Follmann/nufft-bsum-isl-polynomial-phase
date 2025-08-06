function d = ComputeFrequencyIndex(n, m, k, order)
	d = ComputeShiftDifference(n,k,order) - ComputeShiftDifference(m,k,order);
end