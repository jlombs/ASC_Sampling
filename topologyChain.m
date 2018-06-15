function [H,hMBF] = topologyChain(Gin,hMBF,mod)
% Function which performs a topology markov step on the input digraph H (also expressed as a MBF on G, hMBF),
% with respect to parameters in the structure mod, such that each step is
% drawn from an exponential distribution with respect to combinatorial
% metric distance. G is the full graph structure which
% holds the global adjacency data necessary for the nodal constraints.

%% Note, no TNT restriction is handled separately from this algorithm, and has its own notes below

global G

%% Set variable shortcuts

G = Gin;
%totNodes = numel(hMBF);
numRoots = outdegree(G,1)+1;

totChangable = myNChooseK(numRoots,floor(numRoots/2)); 

expDist = exp(-(1:totChangable)./totChangable)/sum(exp(-(1:totChangable)./totChangable));

if mod.TNT   
    
    %% Failure conditional for reversibility 
    % when jumping outside of state space (selecting too many nodes to add/rem from same distribution but now wrt current state)
    
    failure = 1;
    
    while failure
        
        %% Set add/remove action, checking for edge cases which force the situation first

        if sum(hMBF) == numRoots

            addRem = true;

        elseif all(hMBF)

            addRem = false;

        else

            addRem = true;
            if rand < .5
                addRem = false;
            end

        end

        if addRem

            %% Compute addable nodes

            uncoupledNodes = addableNodes(find(hMBF));

        else

            %% Compute removable nodes

            uncoupledNodes = removableNodes(numRoots,find(hMBF));

        end

        numUncoupledNodes = numel(uncoupledNodes);

        %% Pick nodes (number according to exponential sample distance) uniformly and update probability
        %{
        if numUncoupledNodes == 1

            dist = 1;
            bits = uncoupledNodes;

        else

            dist = datasample(1:numUncoupledNodes,1,'Replace',false,'Weights',exp(-(1:numUncoupledNodes)./numUncoupledNodes)); 
            bits = randsample(uncoupledNodes,dist,false);

        end
        %}
        
        dist = datasample(1:totChangable,1,'Replace',false,'Weights',expDist);
        
        if dist > numUncoupledNodes
            
           continue
            
        elseif numUncoupledNodes == 1
            
            bits = uncoupledNodes;
            
        else
            
            bits = randsample(uncoupledNodes,dist,false);
    
        end
        
        failure = false;

        Pf = 1/myNChooseK(numUncoupledNodes,dist);

        %% Flip the bit

        hMBF(bits) = ~hMBF(bits);

        %% Compute reverse probability from new number of uncoupledNodes, same dist

        if addRem

            %% Compute removable nodes

            uncoupledNodes = removableNodes(numRoots,find(hMBF));

        else

            %% Compute addable nodes

            uncoupledNodes = addableNodes(find(hMBF));

        end

        Pb = 1/myNChooseK(numel(uncoupledNodes),dist);  

        H.G = subgraph(G,find(hMBF));
        H.P = Pb/Pf;
        
    end
    
else
    
    %% NO TNT 
    
    %% List edge nodes
    level2Nodes = (numRoots+1):(numRoots + myNChooseK(numRoots,2));
    totEdges = numel(level2Nodes);
    
    %% Draw number of nodes to flip based on exponential distance within available edges

    dist = datasample(1:totEdges,1,'Replace',false,'Weights',exp(-(1:totEdges)./totEdges)); 
    bits = randsample(level2Nodes,dist,false);
    
    %% Re-Establish MBF based on G
    % Creates the fullly connected MBF
    
    hMBFNew = true(1,numel(hMBF));
    
    %% Flip edge bits
    % resets the edges based on hMBF and new flips
    hMBFNew(level2Nodes) = hMBF(level2Nodes);
    hMBFNew(bits) = ~hMBFNew(bits);
    
    %% Find all off edges
    
    offs = numRoots + find(~hMBFNew(level2Nodes));
    
    %% Turn off all dependent higher nodes based on edge set
    
    for i = 1:numel(offs)
        
        hMBFNew = hMBFNew & ~isfinite(shortestpathtree(G,offs(i),'OutputForm','Vector'));
        
    end
    
    %% Pull out new graph
    
    H.G = subgraph(G,find(hMBFNew));
    hMBF = hMBFNew;
    
    %% Algorithm intrinsically has symmetric transition probs
    H.P = 1;
    
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

end

function rN = removableNodes(numRoots,nodeNumH)
global G

% Get list of all removable nodes as those without children (and
% aren't roots)

tmp = subgraph(G,nodeNumH);
rN = findnode(G,tmp.Nodes.Name(outdegree(tmp) == 0));
rN(ismembc(rN,1:numRoots)) = [];

end

function aN = addableNodes(nodeNumH)
global G

 % Get list of all addable nodes as all children of H in G not already in H

aN = [];

for j = 1:numel(nodeNumH)

    aN = [aN, successors(G,nodeNumH(j))'];

end

aN = unique(aN);
aN(ismembc(aN,nodeNumH)) = [];
remFlags = false(1,numel(aN));

% Remove from list those nodes which don't have the right
% substructure: ie, not every parent is in H

for j = 1:numel(aN)

    parents = predecessors(G,aN(j));

    if ~all(ismembc(parents,nodeNumH))

        remFlags(j) = true;

    end

end

aN(remFlags) = [];

end