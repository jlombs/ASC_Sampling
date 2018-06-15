function initializeDigraph(numRoots)
%initializeDigraph a function which takes an input number of roots as a
%double valued positive integer and saves the matlab variable which stores
%the totally connected combinatorial graph K* structure associated with
%this state in Matlab's 'digraph' data structure.

% Data is saved in the 'InitialStates' folder at the same directory level
% as this function. If none exists, one is created. File naming conventions
% are based on the number of root nodes requested as initial#.mat

% If the digraph already exists, the function does not run and instead
% returns a notification.

%% Return an error for numRoots > 21, as the structure is too large to be handled
% Factorial function is only accurate to n=21 

if numRoots > 21
    
    error('numRoots too large, the digraph requires too much memory')
    
end

%% Check for proper directory structure

fileList = dir;
fileList = {fileList.name};

if ~any(cellfun(@(x) strcmp(x,'InitialStates'),fileList))
    
    mkdir('InitialStates')
    
end

%% Check if the requested state already exists

initialStateList = dir('./InitialStates/');
initialStateList = {initialStateList.name};
if length(initialStateList) > 2
    
    initialStateList = initialStateList(3:end);

    if any(cellfun(@(x) strcmp(x,strcat('initial',int2str(numRoots),'.mat')),initialStateList))

        fprintf('\nState already exists: refer to the InitialStates folder\n')
        return

    end
    
end

%% State does not exist, so we must create it

%% Shortcut algorithm/special case for totNodes = 1

if numRoots == 1
    
    G = digraph(1,{char(34)},'OmitSelfLoops'); %#ok<NASGU>
    dataName = './InitialStates/initial1';
    save(dataName,'G')
    return
    
end

%% Compute total nodes required

totNodes = sum(arrayfun(@(x) myNChooseK(numRoots,x),1:numRoots));

%% Establish source-target matrix

totEdges = sum(arrayfun(@(x) myNChooseK(numRoots,x)*x,2:numRoots));

%Typecast source-target storage matrix to the minimal integer type

if totNodes < 2^8-1
    
    storeType = 'uint8';
    
elseif totNodes < 2^16-1
    
    storeType = 'uint16';
    
elseif totNodes < 2^32-1
    
    storeType = 'uint32';
    
else
    
    error('Too many nodes for memory handling')
    
end

A = zeros(totEdges,2,storeType);

%% Establish node name list for tracking adjacency
% Coded to handle roots numbers less than 2^16-1, aka no feasible issues
% Node names are stored as unique hashes based on variable length chars
% from converting the associated root list. Can't use numerical root list
% alone since no way to distinguish [1 2 3]->123 from [1 23]->123: tested
% and causes issues.

N = cell(1,totNodes);

% Code in the initial root labels: due to cell whitespace conversion
% issues, bump up the char map by 33 characters.

N(1:numRoots) = cellstr(char(((1:numRoots) + 33)'));

%% Establish node and edge counter

counterN = numRoots;
counterE = 1;

%% Loop over levels
% Only need to include level 2 through max-1: max is handled separately
% since we know it analytically.

for i = 2:(numRoots-1)
   
    levelNodes = VChooseK(cast(1:numRoots,'uint8'),i);
    
    levelIndex =  sum(arrayfun(@(x) myNChooseK(numRoots,x),1:(i-1)));
    subLabelList = N((levelIndex-myNChooseK(numRoots,i-1)+1):levelIndex);
    
    %% Loop over nodes in each level
    
    for j = 1:size(levelNodes,1)
        
        % Update labels
        
        node = levelNodes(j,:);
        counterN  = counterN + 1;
        N{counterN} = char(node + 33);
        
        %Find Connections in level below
        
        subs = VChooseK(node,i-1);
        subLabels = char(subs + 33);

        cons = zeros(i,1,storeType);
        for k = 1:i
            
            cons(k) = find(strcmp(subLabelList,subLabels(k,:)),1);
            
        end

        cons = cons + levelIndex - myNChooseK(numRoots,i-1);

        %Connect into the s-t matrix
        
        A(counterE:(counterE+length(cons)-1),:) = [cons,ones(i,1,storeType).*counterN];
        counterE = counterE + i;
        
    end
    
end

%% Connect highest level

A(counterE:end,:) = [[(totNodes-numRoots):(totNodes-1)]',ones(numRoots,1,storeType).*totNodes]; %#ok<NBRAK>
N(end) = {char(totNodes+33)};
    
%% Save digraph 

G = digraph(A(:,1),A(:,2),[],N); %#ok<NASGU>
dataName = strcat('./InitialStates/initial',int2str(numRoots));
save(dataName,'G')

end