function residual = c(n,k, params, excludeOrder)
	% Introduced in Eqn. (6)
	nOrder = numel(params);
	
	phaseSum = 0;
	for iOrder = 2:nOrder
		if iOrder == excludeOrder, continue; end

		phaseSum = phaseSum + ...
			params(iOrder) * eqn.d(n,k,iOrder);
	end

	residual = exp(1i*phaseSum);
end