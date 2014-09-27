% Compute Pd(i,j) from the fixed probability p
function [pdval] = pd(i,j,p,d2,d3,s4)
	% If b=U (U = min(1 + d3, s4), then #sending <= #receiving
	result = 0;
	b = i+j;
	if b == min(1+d3,s4)
		% compute the result of binomial distribution
		pdval = 0;
		for k = b:s4
			pdval = pdval + (nchoosek(s4,k) * (p^k) * (1-p)^(s4-k));
		end
	% #sending > #receiving - need the ELSE case here
	else
		g = (nchoosek(d3-d2,i) * nchoosek(1 + d2, j)) / nchoosek(1 + d3, b);
		pdval = g * nchoosek(s4, b) * (p ^ b) * (1 - p)^(s4 - b);
	end
end

