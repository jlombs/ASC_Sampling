% Test maskDigraph with my algorithm vs a coin flip algorithm by looking at
% a variable number of samples with a fixed number of roots
close all
clc
clear all

%% Set the number of roots to simulate on and load in the initial state
numRoots = 6;
numSamples= 10000;
numStatesExp = 5000;

newStateFound = false(1,numSamples);
newStateFoundC = false(1,numSamples);
newStateFoundT = false(1,numSamples);

modStruct.ptVol = 0;
modStruct.IR = inf;
bndNodes = [];
modStruct.TNT = true;
modStruct.Type = [2,inf,0];

GTemp = load(strcat('./InitialStates/initial',int2str(numRoots),'.mat'));

G = GTemp.G;



%% Initialize variables

totNodes = G.numnodes;

%% Perform maskings, taking the results and projecting each state as a vector in R^n, and examine the distribution for uniformity

%sampleStorage = zeros(numSamples,totNodes-numRoots);
multiplicityStorage = zeros(1,numStatesExp);
%gaugedStorage = zeros(numStatesExp,totNodes-numRoots);
graphIsoTestStorage = cell(1,numStatesExp);

multiplicityStorageC = zeros(1,numStatesExp);
%gaugedStorageC = zeros(numStatesExp,totNodes-numRoots);
graphIsoTestStorageC = cell(1,numStatesExp);

multiplicityStorageT = zeros(1,numStatesExp);
%gaugedStorageC = zeros(numStatesExp,totNodes-numRoots);
graphIsoTestStorageT = cell(1,numStatesExp);

counter = 1;
counterC = 1;
counterT = 1;

GmMBF = true(1,totNodes);

GmMBFOld = GmMBF;
GmMBFOldT = GmMBF;
GmO = G;
GmOT = G;
GmCO = G;
GmCPb = 1;

acceptedStep = 1;
acceptedStepC = 1;
acceptedStepT = 1;



for sample = 1:numSamples

    %% Perform the masking

    [Gm,~,GmMBF] = maskDigraphLatest(G,GmMBFOld,modStruct,[],[]); 
    GmC = maskDigraphCoins(G);
    [GmT,GmMBFT] = topologyChain(G,GmMBFOldT,modStruct);
    
    rejectFlag = false;
    rejectFlagC = false;
    rejectFlagT = false;

    %%{
    if sample ~= 1 && rand > Gm.P

        GmMBF = GmMBFOld;
        Gm = GmO;
        rejectFlag = true;

    else
     %}   
        GmO = Gm.G;
        Gm = Gm.G;
        GmMBFOld = GmMBF;  
        acceptedStep = acceptedStep + 1;
   %%{     
    end
    %%{
    if sample ~= 1 && rand > GmT.P

        GmMBFT = GmMBFOldT;
        GmT = GmOT;
        rejectFlagT = true;

    else
     %}   
        GmOT = GmT.G;
        GmT = GmT.G;
        GmMBFOldT = GmMBFT;  
        acceptedStepT = acceptedStepT + 1;
   %%{     
    end

    %%{
    if sample ~= 1 && rand > GmCPb/GmC.P

        GmC = GmCO;
        rejectFlagC = true;

    else
      %}  
        GmCPb = GmC.P;
        GmCO = GmC.G;
        GmC = GmC.G;
        acceptedStepC = acceptedStepC + 1;

    end


    %% Map this state (using uniform node weights) into R^(totNodes)
    % Treating each node as a coordinate, we find the projected vector.
    % We remove the roots, since they do not change

    vect = ismember(G.Nodes.Name,Gm.Nodes.Name);
    vectC = ismember(G.Nodes.Name,GmC.Nodes.Name);
    vectT = ismember(G.Nodes.Name,GmT.Nodes.Name);

    %If you want to keep track of the actual inequivalent states...
    %%{ 
    if sample == 1

        graphIsoTestStorage{counter} = Gm;
        multiplicityStorage(counter) = 1;
        newStateFound(acceptedStep) = true;
        %gaugedStorage(counter,:) = vect((numRoots+1):end);

        graphIsoTestStorageC{counterC} = GmC;
        multiplicityStorageC(counterC) = 1;
        newStateFoundC(acceptedStepC) = true;
        %gaugedStorageC(counterC,:) = vectC((numRoots+1):end);
        
        graphIsoTestStorageT{counterT} = GmT;
        multiplicityStorageT(counterT) = 1;
        newStateFoundT(acceptedStepT) = true;
        %gaugedStorageC(counterC,:) = vectC((numRoots+1):end);


    else

        %temp = ismember(sampleStorage(1:(runs-1),:),vect((numRoots+1):end)','rows');

        %if ~any(temp)
        if ~rejectFlag
            
            for i = 1:counter
                
                if  isisomorphic(graphIsoTestStorage{i},Gm)

                    multiplicityStorage(i) = multiplicityStorage(i) + 1;

                    break

                end

                if i == counter

                    counter = counter + 1;
                    graphIsoTestStorage{counter} = Gm;
                    multiplicityStorage(counter) = 1;
                    %gaugedStorage(counter,:) = vect((numRoots+1):end);
                    newStateFound(acceptedStep) = true;

                end

            end
            
        end

        
        if ~rejectFlagC
            
            for i = 1:counterC

                if isisomorphic(graphIsoTestStorageC{i},GmC)

                    multiplicityStorageC(i) = multiplicityStorageC(i) + 1;
                    break

                end

                if i == counterC

                    counterC = counterC + 1;
                    graphIsoTestStorageC{counterC} = GmC;
                    multiplicityStorageC(counterC) = 1;
                    %gaugedStorageC(counterC,:) = vectC((numRoots+1):end);
                    newStateFoundC(acceptedStepC) = true;

                end

            end
            
        end
        
        if ~rejectFlagT
            
            for i = 1:counterT

                if isisomorphic(graphIsoTestStorageT{i},GmT)

                    multiplicityStorageT(i) = multiplicityStorageT(i) + 1;
                    break

                end

                if i == counterT

                    counterT = counterT + 1;
                    graphIsoTestStorageT{counterT} = GmT;
                    multiplicityStorageT(counterT) = 1;
                    %gaugedStorageC(counterC,:) = vectC((numRoots+1):end);
                    newStateFoundT(acceptedStepT) = true;

                end

            end
            
        end

        %else

        %    find(temp,1)

    end

    %}
    %sampleStorage(runs,:) = vect((numRoots+1):end);


end

%graphIsoTestStorage((counter+1):end) = [];% = graphIsoTestStorage(~cellfun('isempty',graphIsoTestStorage));
multiplicityStorage((counter+1):end) = [];
%gaugedStorage((counter+1):end,:) = [];
%skew = skewness(sampleStorage,0);
newStateFound((acceptedStep+1):end) = [];

%graphIsoTestStorageC((counterC+1):end) = [];% = graphIsoTestStorage(~cellfun('isempty',graphIsoTestStorage));
multiplicityStorageC((counterC+1):end) = [];
%gaugedStorageC((counterC+1):end,:) = [];
newStateFoundC((acceptedStepC+1):end) = [];

%graphIsoTestStorageC((counterC+1):end) = [];% = graphIsoTestStorage(~cellfun('isempty',graphIsoTestStorage));
multiplicityStorageT((counterT+1):end) = [];
%gaugedStorageC((counterC+1):end,:) = [];
newStateFoundT((acceptedStepT+1):end) = [];

%% Compute Data

% Find unique states

%[~,ind,~] = unique(sampleStorage,'rows');

% Find how many times a unique state occurs

%mult = arrayfun(@(x) sum(ismember(sampleStorage,sampleStorage(ind(x),:),'rows')),1:length(ind))

%residuals = mult - mean(mult);

residuals = (multiplicityStorage - mean(multiplicityStorage))/numSamples; 

residualsC = (multiplicityStorageC - mean(multiplicityStorageC))/numSamples; 

residualsT = (multiplicityStorageT - mean(multiplicityStorageT))/numSamples; 

%{
figure(3)
plot(1:numel(residuals),residuals)
hold on
plot(1:numel(residualsC),residualsC)
title(sprintf('Multiplicity Residuals for Unique States with n = %d \n Normalized to the Number of Samples = %d',numRoots,numSamples))
xlabel('State Index')
ylabel('Residual')
legend(sprintf('Balanced Algorithm, %d States; std = %3.2e',numel(multiplicityStorage),std(residuals)),sprintf('Coin, %d States; std = %3.2e',numel(multiplicityStorageC),std(residualsC)))
%}


%simCounter = simCounter + 1;

%{
mom2 = moment(sampleStorage,2);
mom3 = moment(sampleStorage,3);
kurt = kurtosis(sampleStorage,0);
%}


figure(1)
plot(1:numel(newStateFound),cumsum(newStateFound))
hold on
plot(1:numel(newStateFoundC),cumsum(newStateFoundC))
plot(1:numel(newStateFoundT),cumsum(newStateFoundT))
title(sprintf('Number of Unique States Hit in a Walk on %d Roots While Varying Steps',numRoots))
title(sprintf('Number of Unique Geometric States While Sampling C_%d as a Function of Accepted Transitions, \n Linearly Interpolated',numRoots))
xlabel('Accepted Step')
ylabel('Number of Unique States')
legend('Balanced Algorithm','Coin Toss Sampling','Local Random Walk')

figure(2)
plot(1:numel(residualsC),residualsC)
hold on
plot(1:numel(residuals),residuals)
plot(1:numel(residualsT),residualsT)
title(sprintf('Multiplicity Residuals of Unique Geometric States on %d Samples Drawn from C_%d \n Linearly Interpolated',numSamples,numRoots))
xlabel('State Index')
ylabel('Residual')
legend('Coin Toss Sampling','Balanced Algorithm','Local Random Walk')
