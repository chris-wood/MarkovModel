%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File: runSim.m
% Author: Christopher Wood, caw4567@rit.edu
% Description: Monte carlo simulation for the key distribution times
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation parameters
NUM_SAMPLES = 10000;
MAX_CHILDREN = 2;
NUM_NODES = [5,10,15,20,25,30];
%NUM_NODES = [5];
PROBABILITIES = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
%PROBABILITIES = [0.5];
[~, numSims] = size(NUM_NODES);
[~, numProbs] = size(PROBABILITIES);

% Set up the results matrix
times = zeros(numProbs, numSims,NUM_SAMPLES); % we perform one simulation for each number of nodes
avgTimes = zeros(numProbs, numSims, 1);
finalTable = zeros(numProbs, numSims, 4);

% Run the simulation nSamples times
disp('Starting the simulation...');
totalTime = 0;
for p = 1:numProbs
    for n = 1:numSims
        for i = 1:NUM_SAMPLES
            % Initialize the adj. matrix representation for the nodes and network
            % No one is connected at the beginning.
            time = 0;
            nConnected = 0;
            aMatrix = zeros(NUM_NODES(n), NUM_NODES(n));
            cMatrix = zeros(NUM_NODES(n));
            for r = 1:NUM_NODES(n)
                cMatrix(r) = 0;
                for c = 1:NUM_NODES(n)
                    aMatrix(r,c) = 0;
                end
            end

            % Set the root node to have the key at time 0
            cMatrix(1) = 1;

            % Loop while we do try to establish a connection with each node
            while (nConnected < (NUM_NODES(n) - 1)) % We go until connected == (n-1)
                % Find the unconnected nodes
                tempList = zeros(1, NUM_NODES(n));
                nUnconnected = 0;
                for j = 1:NUM_NODES(n)
                    if (cMatrix(j) == 0)
                        nUnconnected = nUnconnected + 1;
                        tempList(nUnconnected) = j; % Flag as unconnected
                    end
                end
                unconnected = zeros(1, nUnconnected);
                for j = 1:nUnconnected
                    unconnected(j) = tempList(j);
                end
                
                % For each node that is ready, decide with probability p
                % if it should receive the key at this instance in time.
                readyList = zeros(1, nUnconnected);
                nReady = 0;
                for j = 1:nUnconnected
                    if (rand(1) < PROBABILITIES(p))
                        readyList(j) = 1; % Flag it as ready
                        nReady = nReady + 1;
                    end
                end
                
                % Compute the set of available parents at this iteration
                [parentList, parentCount] = readyParents(aMatrix, cMatrix, MAX_CHILDREN, NUM_NODES(n));
                
                % Shuffle algorithm
                unconnected;
                readyList;
                parentList;
                aMatrix;
                for j=parentCount:-1:1
                    index = randi(j,1);
                    temp = parentList(index);
                    parentList(index) = parentList(j);
                    parentList(j) = temp;
                end
                
                % Find upper bound on connections
                bound = min(parentCount, nReady);
                
                % Make the connections between ready children and available
                % parents
                readyIndex = 1;
                for j = 1:bound
                    % Skip over nodes that were deemed not ready
                    while (readyList(readyIndex) == 0)
                        readyIndex = readyIndex + 1;
                    end
                    
                    % Tie these guys together
                    %disp('connecting...');
                    child = unconnected(readyIndex);
                    parent = parentList(j);
                    aMatrix(child, parent) = 1;
                    aMatrix(parent, child) = 1;
                    cMatrix(child) = 1;
                    nConnected = nConnected + 1;
                    readyIndex = readyIndex + 1;
                end

                % Increment the time variable
                time = time + 1;
            end
            
            totalTime = totalTime + time;
            times(p,n,i) = time;
        end
    
        avgTimes(p,n,1) = totalTime / NUM_SAMPLES;
    end
end

% Display...
%times'

% Calculate the average and standard deviation for each node simulation
for p = 1:numProbs
    for i = 1:numSims
        avg = mean(times(p, i,:));
        stddev = std(times(p, i,:));
        stderr = 2 * (stddev / (NUM_SAMPLES^(1/2)));
        finalTable(p,i,1) = NUM_NODES(i);
        finalTable(p,i,2) = avg;
        finalTable(p,i,3) = stddev;
        finalTable(p,i,4) = stderr;
    end
end

% Display the final table
disp(finalTable);

%disp('Expected time');
%disp(totalTime / NUM_SAMPLES)
%xlabel('Number of Nodes');
%ylabel('Estimated Time (E(t))');
%plot(avgTimes);

