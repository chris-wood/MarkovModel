% readParents - Generate the list of ready parents that can accept children
% for authentication (and the length of that list for simplicity)
function [parentList, parentCount] = readyParents(aMatrix, cMatrix, authMatrix, maxChildren, nNodes, kMult)
    % Assume none are ready to start
    parentList = [];
    tempList = zeros(1, nNodes);
    parentCount = 0;
    for i = 1:nNodes
        if (cMatrix(i) == 1)
            cCount = 0;
            workingCount = 0;
            
            % Search the adjacency matrix
            for j = 1:nNodes
                if (aMatrix(i,j) == 1)
                    cCount = cCount + 1;
                end
            end
            
            % Search the authentication matrix
            for kIndex = 1:(kMult) % this was just kMult, but those in the last stage will always get the key immediately
                for cIndex = 1:nNodes
                    if (authMatrix(kIndex, i, cIndex) == 1) 
                        workingCount = workingCount + 1;
                        cCount = cCount + 1; % currently trying to authenticate some other node that isn't connected!
                    end
                end
            end
            
            % Take away the parent we are adjacent to
            if (i ~= 1)
                cCount = cCount - 1; 
            end
            
            % We're only ready if we haven't maxed out at this point
            if (cCount < maxChildren && workingCount == 0)
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


