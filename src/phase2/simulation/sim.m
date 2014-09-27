%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File: sim.m
% Author: Christopher Wood, caw4567@rit.edu
% Description: Monte carlo simulation for the key distribution times using
%   the unbounded spanning tree model for # messages per transaction and
%   the # maximum number of children.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The output file - open for overwriting with this new simulation
fName = 'sim_run_20130922.txt';
fid = fopen(fName, 'w');
if (fid == -1)
	disp('Error: could not open the file for output.');
	exit;
else
	disp('Success: output file opened successfully - proceeding to gather data.');
end

% Simulation parameters
numSamples = 100000; 
maxChildren = [2, 3, 4]; % this is k
maxMessages = [2, 3, 4]; % this is m
nodeCount = [5, 10, 25, 50, 100, 150, 250, 500];
p1Probs = [1, 0.9, 0.75, 0.5, 0.25, 0.1];
p2Probs= [1, 0.9, 0.75, 0.5, 0.25, 0.1];

% Extract sizes dynamically...
[t1, numNodes] = size(nodeCount);
[t2, numP1probs] = size(p1Probs);
[t3, numP2probs] = size(p2Probs);
[t4, numChildren] = size(maxChildren);
[t5, numMessages] = size(maxMessages);

% Result containers
times = zeros(numMessages, numChildren, numP1probs, numP2probs, numNodes, numSamples);
finalTable = zeros(numMessages, numChildren, numP1probs, numP2probs, numNodes, 4); % num nodes, avg time, stddev, error in the last coord

% Run the simulation nSamples times
disp('Starting the simulation...');
for messageIndex = 1:numMessages
    for childIndex = 1:numChildren
        disp(sprintf('k,m =  %d, %d', maxChildren(childIndex), maxMessages(messageIndex)))
        for p1Index = 1:numP1probs
            for p2Index = 1:numP2probs
                for n = 1:numNodes
                    disp(sprintf('Simulation for %d nodes with p1 = %d and p2 = %d', nodeCount(n), p1Probs(p1Index), p2Probs(p2Index)))
                    total = 0;
                    for i = 1:numSamples
                        % Initialize the adj. matrix representation for the nodes and network
                        % No one is connected at the beginning...
                        time = 0; % time = #t2 events
                        nConnected = 0;

                        % The matrix to store authentication steps in time.
                        authMatrix = zeros(maxMessages(messageIndex), nodeCount(n), nodeCount(n));

                        % The adjacency matrix stores those nodes node connections (the
                        % tree).
                        aMatrix = zeros(nodeCount(n), nodeCount(n));

                        % The connected vector that indicates whether a node has
                        % the key (it is connected).
                        cMatrix = zeros(1, nodeCount(n));

                        % Set the root node to have the key at time 0
                        cMatrix(1) = 1;

                        % Loop while we do try to establish a connection with each node
                        while (nConnected < (nodeCount(n) - 1)) % We go until connected == (n-1)
                            
                            % DEBUG
                            % fprintf('Time step: %i\n', time);

                            % Find the unconnected nodes from the connected list
                            tempList = zeros(1, nodeCount(n));
                            tempListBack = zeros(1, nodeCount(n));
                            nUnconnected = 0;
                            for j = 1:nodeCount(n)
                                tempList(j) = -1; % mark as invalid to start...
                                if (cMatrix(j) == 0)
                                    nUnconnected = nUnconnected + 1;
                                    tempList(nUnconnected) = j; % Flag as unconnected
                                end
                            end

                            % Strip out all nodes that are currently in 
                            % authentication stage.
                            for kIndex = 1:(maxMessages(messageIndex)) % all nodes in the last stage will be ready for the key...
                               for rIndex = 1:nodeCount(n)
                                  for cIndex = 1:nodeCount(n)
                                     if (authMatrix(kIndex, rIndex, cIndex) == 1)
                                        % cIndex is being authenticated by rIndex,
                                        % so take out cIndex from the list if it's
                                        % in there.
                                        for nIndex = 1:nodeCount(n)
                                            if (tempList(nIndex) == cIndex)
                                                tempList(nIndex) = -1; % set it back to invalid
                                                nUnconnected = nUnconnected - 1; % decrement since we took it out of the list
                                            end
                                        end
                                     end
                                  end
                               end
                            end

                            % Build up the unconnected list
                            unconnected = zeros(1, nUnconnected);
                            tempIndex = 1;
                            for j = 1:nUnconnected
                                % Skip over invalid entries
                                while (tempList(tempIndex) == -1)
                                    tempIndex = tempIndex + 1;
                                end

                                % Add this element to the list
                                unconnected(j) = tempList(tempIndex);
                                tempIndex = tempIndex + 1;
                            end

                            % For each node that is ready, decide with probability p
                            % if it should receive the key at this instance in time.
                            readyList = zeros(1, nUnconnected);
                            nReady = 0;
                            for j = 1:nUnconnected
                                if (rand(1) <= p1Probs(p1Index))
                                    readyList(j) = 1; % Flag it as ready for authentication
                                    nReady = nReady + 1;
                                end
                            end

                            %%% We need to work backwards through the
                            %%% authentication stages so we don't accidentally
                            %%% carry connections through the pipeline in the same
                            %%% instance in time (only one transition at a time)

                            % Handle the key distribution step now (authentication
                            % is complete at this stage in the auth matrix)
                            for rIndex = 1:nodeCount(n)
                               for cIndex = 1:nodeCount(n)
                                  if (authMatrix(maxMessages(messageIndex), rIndex, cIndex) == 1) % was maxMessages(messageIndex), not maxMessages(messageIndex) - 1
                                      % Once we're in this stage we automatically go to having the key (no questions asked)
                                      authMatrix(maxMessages(messageIndex), rIndex, cIndex) = 0; % no longer in the authentication stage...
                                      aMatrix(rIndex, cIndex) = 1;
                                      aMatrix(cIndex, rIndex) = 1;
                                      cMatrix(cIndex) = 1;
                                      nConnected = nConnected + 1;
                                  end
                               end
                            end

                            % Now check to see if the nodes doing
                            % authentication march forwards in time
                            bound = maxMessages(messageIndex) - 1;
                            for kIndex = bound:-1:1
                                for rIndex = 1:nodeCount(n)
                                   for cIndex = 1:nodeCount(n)
                                      % If a pair of nodes is attempting
                                      % authentcation, check to see if they
                                      % make progress
                                      if (authMatrix(kIndex, rIndex, cIndex) == 1)
                                          if (rand(1) <= p2Probs(p2Index))
                                              authMatrix(kIndex + 1, rIndex, cIndex) = 1;
                                              authMatrix(kIndex, rIndex, cIndex) = 0;
                                          end
                                      end
                                   end
                                end
                            end
                            
                            

                            % Compute the set of available parents at this iteration
                            [parentList, parentCount] = readyParents(aMatrix, cMatrix, authMatrix, maxChildren(childIndex), nodeCount(n), maxMessages(messageIndex));
                            %disp('parents available at this time');
                            %for j = 1:parentCount
                            %   disp(parentList(j));
                            %end

                            % Shuffle algorithm
                            for j = parentCount:-1:1
                                index = randi(j,1);
                                temp = parentList(index);
                                parentList(index) = parentList(j);
                                parentList(j) = temp;
                            end

                            % Find upper bound on connections
                            bound = min(parentCount, nReady);

                            % Start these nodes off in the authentication step
                            readyIndex = 1;
                            for j = 1:bound
                                % Skip over nodes that were deemed not ready
                                while (readyList(readyIndex) == 0)
                                    readyIndex = readyIndex + 1;
                                end

                                % Hook these guys into the auth matrix
                                child = unconnected(readyIndex);
                                parent = parentList(j);
                                authMatrix(1, parent, child) = 1; % this is a directed graph, so don't point from child->parent
                                readyIndex = readyIndex + 1;
                            end

                            time = time + 1;
                        end

                        % Take away the last step in time - since key distribution happens on the time step before
                        time = time - 1; 

                        % Record the total time for simulation
                        times(messageIndex, childIndex, p1Index, p2Index, n, i) = time;
                    end
                    avg = mean(times(messageIndex, childIndex, p1Index, p2Index, n,:));
                    stddev = std(times(messageIndex, childIndex, p1Index, p2Index, n,:));
                    stderr = 2 * (stddev / (numSamples^(1/2)));
                    fprintf('%d, %d, %d, %d, %d, %d, %d, %d\n', maxChildren(childIndex), maxMessages(messageIndex), nodeCount(n), p1Probs(p1Index), p2Probs(p2Index), avg, stddev, stderr);
					fprintf(fid, '%d, %d, %d, %d, %d, %d, %d, %d\n', maxChildren(childIndex), maxMessages(messageIndex), nodeCount(n), p1Probs(p1Index), p2Probs(p2Index), avg, stddev, stderr);
                end
            end
        end
    end
end

% Generate a plot for each one
% figureId = 1;
% for messageIndex = 1:numMessages
%     for childIndex = 1:numChildren
%         for p1Index = 1:numP1probs
%             temp = zeros(numP2probs, numNodes);
%             for p2Index = 1:numP2probs
%                for n = 1:numNodes
%                   temp(p2Index, n) = mean(times(messageIndex, childIndex, p1Index, p2Index, n,:)); % the second element is the average time
%                end
%             end
%             %figure(figureId);
%             %figureId = figureId + 1; % forward the figures... math is fun.
%             %plot(temp);
%             %set(gca,'XTickLabel',{'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7','0.8', '0.9', '1.0'});
%             %title([sprintf('Key Distribution Time for %d Children with Key Probability = %d', maxChildren(childIndex), keyProbabilities(pKeyIndex))]);
%             %xlabel('Authentication Probability');
%             %ylabel('Average Re-Key Time (epochs)');
%         end
%     end
% end

%disp('START OUTPUT');

% Calculate the average and standard deviation for each node simulation
% finalTable = zeros(numChildren, numP1probs, numP2probs, numNodes, 4);
%for messageIndex = 1:numMessages
%    for childIndex = 1:numChildren
%      for p1Index = 1:numP1probs
%        for p2Index = 1:numP2probs
%          for i = 1:numNodes
%            avg = mean(times(messageIndex, childIndex, p1Index, p2Index, i,:));
%            stddev = std(times(messageIndex, childIndex, p1Index, p2Index, i,:));
%%            stderr = 2 * (stddev / (numSamples^(1/2)));
%%            
            % Display the CSV output and build up the final table
%            fprintf('%d, %d, %d, %d, %d, %d, %d, %d\n', maxChildren(childIndex), maxMessages(messageIndex), nodeCount(i), p1Probs(p1Index), p2Probs(p2Index), avg, stddev, stderr);
%            finalTable(numMessages, childIndex, p1Index, p2Index,i,1) = nodeCount(i);
%            finalTable(numMessages, childIndex, p1Index, p2Index,i,2) = avg;
%            finalTable(numMessages, childIndex, p1Index, p2Index,i,3) = stddev;
%            finalTable(numMessages, childIndex, p1Index, p2Index,i,4) = stderr;
%          end
%        end
%      end
%    end
%end

% Display the final table
%disp(finalTable);
