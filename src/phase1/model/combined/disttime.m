function [t] = disttime(n, p, bw, psize, pcount)
	% Determine the depth of the tree.
	depth = log2(n);

	% Determine the number of bits that travel in each IKE
	% transaction.
	totalBits = psize * pcount * 8;
	
	% Determine the estimated time using the model
	etime = time(0, 0, n - 1, p, n);

	% Now compute and display the final key distribution time.
	t = 1 / ((bw / etime) / (totalBits * depth));
	t = t + (0.018 * pcount + depth); % 18ms for epoch time slot
end

