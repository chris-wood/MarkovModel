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
    xlabel('Probability');
    ylabel('Average Re-Key Time (epochs)');
end
