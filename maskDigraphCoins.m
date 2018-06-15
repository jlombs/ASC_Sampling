function Gm = maskDigraphCoins(G)
%% Function which uniformly masks nodes and their parents from the input digraph G, outputing the modified structure
% Masking is done with intendended uniform probability over the node space, although is this impossible to do analytically 
% for number of roots > 3. Instead, we perform a probabilistic rebalancing such that we round regions of small probability 
% mass based on the number of nodes affected upon removal. It has the properties that an isolated node will always have a coin flip
% chance of being utilized, while the complete state and the complete empty
% state have the same probabilities. See discussion in notes for details.
% supplied (as it can be precomputed).


% Uses Kahle's algorithm with uniform coin flip probabilities

%% Establish variable shortcuts
totNodes = G.numnodes;
numRoots = indegree(G,G.numnodes);
nodesPerLevel1Up = arrayfun(@(x) myNChooseK(numRoots,x),2:numRoots);

% Note that we do not remove any of the roots (probVect will be 0 in those indices).
% This is enforcing a paradigm wherein the boundary data are just the isolated roots. 
% If we could allow root removal, we would just be sampling the lesser root subspaces, since
% the existance of a root is essential for the higher dimensional pieces
% (the coupling makes the sampling factorize). So we study for the stationary 
% states on n nodes for each simulation. This also helps us from dealing
% with the empty state. More advanced boundary structures or level freezing
% can be accomplished here as well, but are not implemented yet.

%% Establish counter for level indexing and a logical vector on the removed nodes
counterL = numRoots + 1;
fullRemNodes = false(1,totNodes);

%% Loop over levels and inductively find the nodes to remove

numCoinFlips = 0;

for i = 1:(numRoots-1)
    
    % Find the total number of nodes removed so far in a given level
    removedInLevelLogical = fullRemNodes(counterL:(counterL+nodesPerLevel1Up(i)-1)); 
    totRemovedInLevel = sum(removedInLevelLogical);
    
    % If all are removed, break out as there is nothing more to do--should
    % be redunant and not used due to a break below, but left for safety.
    if nodesPerLevel1Up(i) == totRemovedInLevel
        
        break
    
    end
    
    %{
    % Compute the probability vector for removing x nodes at once. Probs
    % are uniform for all x, except x=0, which is accounted for separately
    % by normalization of the probability vector to unity through a difference method.    
    numProbVect = ones(1,nodesPerLevel1Up(i)-totRemovedInLevel)./(1+sum(nodesPerLevel1Up(i:end)) - sum(fullRemNodes(counterL:end))); 
    numProbVect = [1-sum(numProbVect),numProbVect]; %#ok<AGROW>
    
    % Draw the sample of how many nodes to remove
    numRemoveInLevel = randsample(0:(length(numProbVect)-1),1,true,numProbVect);
    
    % If no removals, continue to next loop and update level counter
    if numRemoveInLevel == 0
    
        counterL = counterL + nodesPerLevel1Up(i);
        continue
    
    end
    
    %}
    
    % Find all available nodes in the level
    availableLevelNodes = counterL:(counterL + nodesPerLevel1Up(i)-1); 
    availableLevelNodes = availableLevelNodes(~removedInLevelLogical);
    
    remNodes = zeros(1,numel(availableLevelNodes));
    
    for j = 1:numel(availableLevelNodes)
        
        numCoinFlips = numCoinFlips + 1;
        
        if rand > .5
            
            remNodes(j) = availableLevelNodes(j);
            
        end
        
    end

    remNodes = remNodes(remNodes > 0);
    
    if isempty(remNodes)
        
         counterL = counterL + nodesPerLevel1Up(i);
         continue
         
    end
    
    %{
    % If only one is available, and this point it will be removed and we
    % can do so by flagging it directly. Otherwise, if removing all nodes
    % in the level, we can also compute the effect on the logical vector
    % directly without needing to compute trees. This terminates the
    % algorithm. Lastly, random sample from the subset of nodes to remove
    % if need be.
    if numel(availableLevelNodes) == 1
        
        remNodes = availableLevelNodes;
    
    elseif numel(availableLevelNodes) == numRemoveInLevel
        
        temp = false(1,numel(fullRemNodes));
        temp(availableLevelNodes(1):end) = true;
        fullRemNodes = fullRemNodes | temp;
        break
        
    else
        
        remNodes = randsample(availableLevelNodes,numRemoveInLevel);
    
    end
    
    %}
    
    % Update the logical vector based on the maximal tree from the removed
    % nodes for continuity.
    
    for j = 1:length(remNodes)
   
        fullRemNodes = fullRemNodes | isfinite(shortestpathtree(G,remNodes(j),'OutputForm','Vector'));
    
    end
   
    counterL = counterL + nodesPerLevel1Up(i);
    
end

Gm.P = numCoinFlips;

%% If nothing is removed, we are finished.
if ~any(fullRemNodes)
    
    Gm.G = G;
    return

end

%% Remove the nodes requested and update the digraph

Gm.G = rmnode(G,find(fullRemNodes));

end