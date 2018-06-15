function [H,bndNodes,hMBFNew] = maskDigraphLatest(G,hMBF,mod,bndNodes,preprocess)
%% Function which uniformly masks nodes and their parents from the input digraph G, outputing the modified structure
% Masking is done with intendended uniform probability over the node space, although is this impossible to do analytically 
% for number of roots > 3. Instead, we perform a probabilistic rebalancing such that we round regions of small probability 
% mass based on the number of nodes affected upon removal. It has the properties that an isolated node will always have a coin flip
% chance of being utilized, while the complete state and the complete empty
% state have the same probabilities. See paper for details.

% Additionally stores the new level-node structure in H as structure for future uses
% H.G = G; H.levelMat = levelMat; H.P = probability of getting this state

% For the fixed weight simulation, we also remove any remaining nodes which would violate ptVol, since they are known in advance.
% This is done first, for speed, since the 'continue' will take care of the
% rest.

% mod.Type indicates the type of model [flag,d,flag2], with flag=true
% indicated fixed volume, and d~=inf indicating a maximal level,
% and flag2 indicates exact volume solution if d==3 also. 

% bndNodes is a list of nodes which can not be removed (boundary nodes, excludes roots).
% Note, ptVol cuts take precidence over bndNodes if mod.Type(1)=true; as
% does level cutoffs

% preprocess is a structure which accomplishes many of the tasks found in
% this first segment of code, done in the MCMC script instead to avoid
% redoing the computations for each sample. if preprocess is empty, the
% script will handle the necessary computations internally.

if ~isempty(preprocess)
    
    totNodes = preprocess.totNodes;
    numRoots = preprocess.numRoots;
    nodesPerLevel1Up = preprocess.nodesPerLevel1Up;
    bndPerLevel = preprocess.bndPerLevel;
    counterL = numRoots + 1;
    fullRemNodes = false(1,totNodes);
    
    if ~isempty(bndNodes)

        bndNodeNames = preprocess.bndNodeNames;

    end
    
    levelMax = preprocess.levelMax;
    
    
else
    
    %% Establish variable shortcuts
    totNodes = G.numnodes;
    numRoots = indegree(G,G.numnodes);
    nodesPerLevel1Up = arrayfun(@(x) myNChooseK(numRoots,x),2:numRoots);

    %% Compute boundary nodes per level (starting at level 2)
    bndPerLevel = cell(numRoots-1,1);
    if ~isempty(bndNodes)

        degrees = indegree(G,bndNodes);
        degreesU = unique(degrees);

        for i = 2:(numRoots-1)

            if i <= numel(degreesU)

                bndPerLevel{i-1} = bndNodes(degrees == degreesU(i));

            else

                break

            end

        end

    end

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

    %% Trim any higher nodes depending on ptVol and uniform edge length
    if mod.Type(1)

        for i = 2:numRoots

            vol = mod.Type(1)^(i-1)*sqrt((i-1)+1)/(2^((i-1)/2)*factorial(i-1));

            if vol < mod.ptVol || vol > mod.IR

                fullRemNodes(sum([counterL,nodesPerLevel1Up(1:(i-2))]):end) = true;

                if ~isempty(bndPerLevel{i-1})

                    tmp = bndPerLevel{i-1};
                    bndPerLevel(i-1) = {tmp(~ismember(tmp,sum([counterL,nodesPerLevel1Up(1:(i-2))]):totNodes))};
                    bndPerLevel(i:end) = cell(size(bndPerLevel(i:end)));

                end

                break

            end

        end

    end

    %% Trim any higher nodes depending on level cutoff
    if mod.Type(2) < numRoots

        fullRemNodes(sum([counterL,nodesPerLevel1Up(1:(mod.Type(2)-1))]):end) = true;

        if ~isempty(bndPerLevel{mod.Type(2)})

            bndPerLevel(mod.Type(2):end) = cell(size(bndPerLevel(mod.Type(2):end)));

        end

    end

    %% Update bndNodes depending on whether any cuts were implemented, and get node names for tracking in a cell array
    %Note, we now get rid of any roots as boundaries, since the algorithm is
    %faster if we don't have boundary pieces to worry about, and the roots are
    %always there anyway and treated separely, no our 'no boundary' shortcuts
    %are implemented if there are no boundary pieces at level 2 and higher.

    if ~isempty(bndNodes)

        bndNodes = [bndPerLevel{:}];
        bndNodeNames = G.Nodes.Name(bndNodes);

    end
   
    %% TNT trimming flag for level loops

    %If we enforce the embedding restriction that if there CAN be a structure,
    %there IS a structure, then we only need to trim at the edge level since
    %they specify the rest. This is a highly restricted space compared to the
    %full space, but what was originally coded in the classical sims before the
    %TNT move.

    if mod.TNT

        levelMax = numRoots - 1;

    else

        levelMax = 1;

    end

end

fullRemNodesStart = fullRemNodes;

%% Establish probability storage
Pf = 1;
Pb = 1;

%% Loop over levels and inductively find the nodes to remove
% Note, i doesn't index over levels, but indexes over levels-1 since we
% don't trim the roots but instead start at level 2. The indexing is
% shifted at the level of the for loop to start at 1, but the code is
% written assuming we are using the 'next level'. 

for i = 1:levelMax
    
    % Find the total number of nodes removed so far in a given level
    removedInLevelLogical = fullRemNodes(counterL:(counterL+nodesPerLevel1Up(i)-1)); 
    totRemovedInLevel = sum(removedInLevelLogical);
    
    % If all available are removed, skip to next level
    if nodesPerLevel1Up(i) - numel(bndPerLevel{i}) == totRemovedInLevel
        
        if ~any(cellfun(@numel,bndPerLevel(i:end)))
            
            break
            
        else
            
            continue
            
        end
    
    end
    
    % Compute the probability for removing x nodes at once. Probs
    % are uniform for all x, except x=0, which is accounted for separately
    % by normalization of the probability vector to unity through a
    % difference method (or in this implementation, a separate rand
    % filter).
    probRem = 1/(1+sum(nodesPerLevel1Up(i:end)) - sum(fullRemNodes(counterL:end)));
    probKeepAll = 1 - (nodesPerLevel1Up(i)-totRemovedInLevel)*probRem; 
    
    % If no removals, continue to next loop and update level counter
    if rand < probKeepAll
        
        Pf = Pf * probKeepAll;
        counterL = counterL + nodesPerLevel1Up(i);
        continue
        
    else
        
        % Draw the sample of how many nodes to remove
        avail = nodesPerLevel1Up(i)-totRemovedInLevel-numel(bndPerLevel{i});
        numRemoveInLevel = randi([1,avail]);
        Pf = Pf * probRem / myNChooseK(avail,numRemoveInLevel);
        
    end
    
    % Find all available nodes in the level
    availableLevelNodes = counterL:(counterL + nodesPerLevel1Up(i)-1); 
    availableLevelNodes = availableLevelNodes(~removedInLevelLogical);
    availableLevelNodes(ismember(availableLevelNodes,bndPerLevel{i})) = [];

    % If only one is available, at this point it will be removed and we
    % can do so by flagging it directly. Otherwise, if removing all nodes
    % in the level, we can also compute the effect on the logical vector
    % directly without needing to compute trees. This terminates the
    % algorithm. Lastly, random sample from the subset of nodes to remove
    % if need be
    if numel(availableLevelNodes) == 1
        
        remNodes = availableLevelNodes;
    
    elseif numel(availableLevelNodes) == numRemoveInLevel && isempty(bndPerLevel{i})
        
        temp = false(1,numel(fullRemNodes));
        temp(availableLevelNodes(1):end) = true;
        fullRemNodes = fullRemNodes | temp;
        break
        
    else
        
        remNodes = randsample(availableLevelNodes,numRemoveInLevel);
    
    end
    
    % Update the logical vector based on the maximal tree from the removed
    % nodes for continuity.
    
    for j = 1:numel(remNodes)
   
        fullRemNodes = fullRemNodes | isfinite(shortestpathtree(G,remNodes(j),'OutputForm','Vector'));
    
    end
   
    counterL = counterL + nodesPerLevel1Up(i);
    
end

hMBFNew = ~fullRemNodes;

%% Remove the nodes requested and update the digraph

if any(fullRemNodes)
    
    H.G = subgraph(G,find(hMBFNew));

else
    
    H.G = G;
    
end

%% Find the new boundary node numbers and return the list

if ~isempty(bndNodes)
    
    bndNodes = findnode(H.G,bndNodeNames);
    
end

%% Compute new node labels at each level, storing them in a node-level matrix for reference later

inDeg = indegree(H.G);

%Initialize matrix and handle roots
if max(inDeg) == 0
    
    levelMat = {1:numRoots};
    
else

    levelMat = cell(inDeg(end),1);
    levelMat{1} = 1:numRoots;
    
    for i = 2:max(inDeg)
    
        levelMat{i} = find(inDeg == i)';
    
    end
    
end

H.levelMat = levelMat;

%% Compute backwards probability by construction using hMBF
% Note, this could have been done by just saving the Pf for use on the next
% ratio in the MCMC script, but for uniformity of what G.P means and for
% use in hybridizing the samplers, we'll compute it explicitly and return
% the metropolis ratio instead.

counterL = numRoots + 1;
fullRemNodes = fullRemNodesStart;

for i = 1:levelMax
    
    % Find the total number of nodes removed so far in a given level
    removedInLevelLogical = fullRemNodes(counterL:(counterL+nodesPerLevel1Up(i)-1)); 
    totRemovedInLevel = sum(removedInLevelLogical);
    
    % If all available are removed, skip to next level
    if nodesPerLevel1Up(i) - numel(bndPerLevel{i}) == totRemovedInLevel
        
        if ~any(cellfun(@numel,bndPerLevel(i:end)))
            
            break
            
        else
            
            continue
            
        end
    
    end
    
    % Compute the probability for removing x nodes at once. Probs
    % are uniform for all x, except x=0, which is accounted for separately
    % by normalization of the probability vector to unity through a
    % difference method (or in this implementation, a separate rand
    % filter).
    probRem = 1/(1+sum(nodesPerLevel1Up(i:end)) - sum(fullRemNodes(counterL:end)));
    probKeepAll = 1 - (nodesPerLevel1Up(i)-totRemovedInLevel)*probRem; 
    
    % Count number of removals in level from hMBF
    remNodes = counterL + find(~hMBF(counterL:(counterL+nodesPerLevel1Up(i)-1))) - 1;
    numRemoveInLevel = numel(remNodes);
    
    if numRemoveInLevel == 0
        
        Pb = Pb * probKeepAll;
        counterL = counterL + nodesPerLevel1Up(i);
        continue
        
    else
        
        avail = nodesPerLevel1Up(i)-totRemovedInLevel-numel(bndPerLevel{i});
        Pb = Pb * probRem / myNChooseK(avail,numRemoveInLevel);
        
    end
    
    % Update the logical vector based on the maximal tree from the removed
    % nodes for continuity.
    
    for j = 1:numel(remNodes)
   
        fullRemNodes = fullRemNodes | isfinite(shortestpathtree(G,remNodes(j),'OutputForm','Vector'));
    
    end
   
    counterL = counterL + nodesPerLevel1Up(i);
    
end

H.P = Pb/Pf;

end