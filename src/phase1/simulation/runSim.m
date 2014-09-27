%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File: runSim.m
% Author: Christopher Wood, caw4567@rit.edu
% Description: Monte carlo simulation for the key distribution times
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation parameters
NUM_SAMPLES = 10000; %10000
%MAX_CHILDREN = 2;
MAX_CHILDREN = [2,3,4,5,6,7,8]; 
NUM_NODES = [5,10,15,20,25,30];
%NUM_NODES = [5];
PROBABILITIES = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
%PROBABILITIES = [0.5];
[~, numSims] = size(NUM_NODES);
[~, numProbs] = size(PROBABILITIES);
[~, numChildren] = size(MAX_CHILDREN);

% Set up the results matrix
times = zeros(numChildren, numProbs, numSims, NUM_SAMPLES); % we perform one simulation for each number of nodes
avgTimes = zeros(numChildren, numProbs, numSims);
finalTable = zeros(numChildren, numProbs, numSims, 4);

% Run the simulation nSamples times
disp('Starting the simulation...');
for childIndex = 1:numChildren
    disp('Children');
    disp(MAX_CHILDREN(childIndex));
    for p = 1:numProbs
        disp('Probability');
        disp(PROBABILITIES(p));
        for n = 1:numSims
            disp('Number of nodes');
            disp(NUM_NODES(n));
            totalTime = 0;
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
                    [parentList, parentCount] = readyParents(aMatrix, cMatrix, MAX_CHILDREN(childIndex), NUM_NODES(n));

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
                times(childIndex,p,n,i) = time;
            end

            totalTime;
            avgTimes(childIndex,p,n) = totalTime / NUM_SAMPLES;
        end
    end

    % Display...
    %times'

    % Calculate the average and standard deviation for each node simulation
    for p = 1:numProbs
        for i = 1:numSims
            avg = mean(times(childIndex,p,i,:));
            stddev = std(times(childIndex,p, i,:));
            stderr = 2 * (stddev / (NUM_SAMPLES^(1/2)));
            finalTable(childIndex, p,i,1) = NUM_NODES(i);
            finalTable(childIndex, p,i,2) = avg;
            finalTable(childIndex, p,i,3) = stddev;
            finalTable(childIndex, p,i,4) = stderr;
        end
    end
end

% Display the final table
%for c = 1:numChildren
    %disp(finalTable(c,:))
%end
disp(avgTimes)
%plot3(avgTimes(:,:,1))

% Generate a plot for each one
for c = 1:numChildren
    temp = zeros(numProbs, numSims);
    for p = 1:numProbs
       for n = 1:numSims
          %temp(p, n) = avgTimes(c, p, n); % the final table has the correct values
          temp(p, n) = finalTable(c, p, n, 2); % the second element is the average time
       end
    end
    figure(c);
    plot(temp);
    %axis([0.1 1.0 0 Inf])
    %set(gca, 'XTickMode', 'manual');
    %title(MAX_CHILDREN(c),'FontWeight','bold');
    set(gca,'XTickLabel',{'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'});
    title(['Key Distribution Time for ',int2str(MAX_CHILDREN(c)), ' Children']);
    xlabel('Connection Probability');
    ylabel('Average Re-Key Time (epochs)');
end

%x = -pi:.1:pi;
%y = sin(x);
%plot(x,y)
%set(gca,'XTick',0.1:0.1)
%set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})


%figure(numChildren + 1);
%scatter3(MAX_CHILDREN, PROBABILITIES, NUM_NODES, 5, avgTimes)

%surf(PROBABILITIES, NUM_NODES, avgTimes(1,:,:))

%surf(MAX_CHILDREN, PROBABILITIES, avgTimes(1,:), 1)

%p = patch(isosurface(MAX_CHILDREN,PROBABILITIES,NUM_NODES,avgTimes,60))
%isonormals(MAX_CHILDREN,PROBABILITIES,NUM_NODES,avgTimes,p)

%scatter3([1:numChildren], [1:numProbs], [1:numSims], 5, avgTimes(:))

%disp('Expected time');
%disp(totalTime / NUM_SAMPLES)
%xlabel('Number of Nodes');
%ylabel('Estimated Time (E(t))');
%plot(avgTimes);

