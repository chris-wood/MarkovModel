function [parentList, parentCount] = readyParents(aMatrix, cMatrix, maxChildren, nNodes)
%READYNODES Summary of this function goes here
%   Detailed explanation goes here

% Assume none are ready to start
parentList = [];
tempList = zeros(1, nNodes);
parentCount = 0;
for i = 1:nNodes
    if (cMatrix(i) == 1)
        cCount = 0;
        for j = 1:nNodes
            if (aMatrix(i,j) == 1)
                cCount = cCount + 1;
            end
        end
        if (i ~= 1)
            cCount = cCount - 1; % take away the parent we are adjacent to!
        end
        if (cCount < maxChildren)
            parentCount = parentCount + 1;
            tempList(parentCount) = i;
        end
    end
end

% Copy the temp list into the parent list
parentList = zeros(1, parentCount);
for i = 1:parentCount
    parentList(i) = tempList(i);
end

end


