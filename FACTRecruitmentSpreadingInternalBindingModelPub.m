%% FACT Recruitment and Spreading Model

% Full model but with the addition of FACT loading at all nucleosomes


%% Parameter sets used for figures in the publication
%If not specified, they were as in the code below


% Bare Bones with internal binding version (Figure S6A right):

%FACTBindingRate = NucNo*0.0005;
%FACTExtensionRate = 0.3;
%FACTRetractionRate = 0.3;
%FACTExtToRetSwitch = 0.3;
%FACTRetToExtSwitch = 0.3;
%FACTUnbindingRate = 0;




%% Constants of the Simulation

GeneLength = 1510;          %Length of synthetic gene (bp)
NucFootprint = 150;         %Nucleosome-DNA contact size (bp)
NucSpacing = 20;            %Length of DNA between nucleosomes (bp)
%Gene length should be an appropriate combination of Nucleosome
%footprints and spacing lengths

NucNo = floor((GeneLength+NucSpacing)/(NucFootprint+NucSpacing));
%The number of nucleosomes in the system

%Maximum and minimum DNA occupancy of the FACT complex
MaxFACTSize = 2*NucFootprint + NucSpacing;
MinFACTSize = NucFootprint;
%Assume the unextended FACT size is the same as a nucleosome

%Maximum number of FACT complexes that can coexist on the synthetic gene
MaxFACTNo = floor(GeneLength/MinFACTSize);



%Set up the nucleosome positions, boundaries and overlaps (for 
%FACT loading at all nucleosomes)
NucleosomePositions = zeros(1,GeneLength);

NucleosomeBoundaries = zeros(NucNo,2);
NucleosomeOverlaps = zeros(MaxFACTNo,1);

for ictr=1:1:NucNo
    
    temp_posn = (ictr-1)*(NucFootprint+NucSpacing) + 1;
    NucleosomePositions(temp_posn:(temp_posn+NucFootprint-1)) = 1;
    NucleosomeBoundaries(ictr,1) = temp_posn;
    NucleosomeBoundaries(ictr,2) = temp_posn+NucFootprint-1;
    
end


%Rate at which FACT complexes bind to the +1 Nucleosome (if clear):
FACTBindingRate = NucNo*0.0005;
%Rate at which a FACT complex will extend by 1bp (if possible):
FACTExtensionRate = 0.3;
%Rate at which a FACT complex will retract by 1bp (if possible)
FACTRetractionRate = 0.3;

%Rate at which a fully-extended FACT complex will switch to a retracting
%state:
FACTExtToRetSwitch = 0.03;
%Rate at which a fully-retracted FACT complex will swith to an extending
%state:
FACTRetToExtSwitch = 0.03;

%Rate at which FACT complexes will unbind:
FACTUnbindingRate = 0.00005;

%Normalising factor by which above rates are divided.
RateNorm = 10.0;

%Construct the absolute rates (and simplify the variable names)
FBRN = FACTBindingRate / RateNorm;
FERN = FACTExtensionRate / RateNorm;
FRRN = FACTRetractionRate / RateNorm;
FERSN = FACTExtToRetSwitch / RateNorm;
FRESN = FACTRetToExtSwitch / RateNorm;
FUBRN = FACTUnbindingRate / RateNorm;

%Total probability of an event occuring during a single time step must be 
%less than 1
if ( RateNorm < ...
        (FACTBindingRate + ...
        MaxFACTNo*max(FACTExtensionRate,FACTRetractionRate) + ...
        MaxFACTNo*max(FACTExtToRetSwitch,FACTRetToExtSwitch) + ...
        MaxFACTNo*FACTUnbindingRate) )

    error('InchwormScript:NormTooSmall','Warning: The Rate Normalisation Factor was too small, exiting.')
    
end


%Number of time steps used prior to taking measurements:
relaxationSteps = 1000000;

%Number of time steps used when taking measurements:
timeSteps = 1000000;

%Inteval between taking measurements:
samplingInterval = 10; 

%Number of independent simulation runs
Realisations = 10000;








%% Variables and Arrays of the Simulations

%Holder for locations and states of the FACT complexes
FACTStates = zeros(MaxFACTNo,4);
%left posn, right posn, length, extend/retract phase
%final value is 1 for extending; 0 for retracting;
% 2 for fully extended; and -1 for fully retracted;


FACTNo = 0;    %Number of FACT complexes on the synthetic gene
ExtNo = 0;      %Number that are extending
RetNo = 0;      %Number that are retracting
FullyExtNo = 0; %Number that are fully extended
FullyRetNo = 0; %Number that are fully retracted

%just look at a single realisation for now


%Holds a density heatmap of fragments that would be observed in an 
%MNase-ChIP experiment with antibodies for a FACT component
%(should not be reinitialised):
MNaseDensity = zeros(GeneLength,GeneLength); %fragment size,location


%Number of FACT complexes that initiated prior to measurements:
RelaxationInitiatedNo = 0;
%Number of FACT complexes that finished (made it to the end of the 
%synthetic gene and were removed) prior to measurements:
RelaxationFinishedNo = 0;

%Number of FACT complexes that initiated during the measurement period
SteadyStateInitiatedNo = 0;
%Number of FACT complexes that finished during the measurement period
SteadyStateFinishedNo = 0;

%Number of FACT complexes that unbound without finishing prior to 
%measurements:
RelaxationEarlyFinishedNo = 0;
%Number of FACT compleses that unbound without finishing during the 
%measurement period
SteadyStateEarlyFinishedNo = 0;


%% Initiate the RNG

rng('shuffle','twister')



%% Main Program Loop


for rctr=1:1:Realisations
    
    realtics = tic;
    
    %Reinitialise:
    FACTStates = zeros(MaxFACTNo,4);
    FACTNo = 0;
    ExtNo = 0;
    RetNo = 0;
    FullyExtNo = 0;
    FullyRetNo = 0;
    
    
    %----Relax the System to an assumed Steady State----------------------
    
    %Loop over time steps:
    for tctr=1:1:relaxationSteps
        
        %Generate a pseudo random number from a uniform distribution
        diceRoll = rand;
        
        %Use the pseudo random number to test for success of potential 
        %events
        
        if diceRoll < FBRN         
            %Attempt a new FACT binding at a randomly chosen nucleosome
            
            %First choose a potential binding site
            ChosenNuc = randi(NucNo);
            %Then bind if that site is currently unbound
            NucStart = NucleosomeBoundaries(ChosenNuc,1);
            NucEnd = NucleosomeBoundaries(ChosenNuc,2);
            NucleosomeOverlap = 0;
            for ictr=1:1:MaxFACTNo 
                if ( (FACTStates(ictr,1) >= NucStart) && (FACTStates(ictr,1) <= NucEnd) ) ...
                        || ( (FACTStates(ictr,2) >= NucStart) && (FACTStates(ictr,2) <= NucEnd) )
                    NucleosomeOverlap = 1;
                end
            end
            firstPostChosenNucFACT = 0;
            if NucleosomeOverlap == 0
                %find out which (if any) FACTs need moving in the queue
                for ictr=MaxFACTNo:-1:1
                    if(FACTStates(ictr,1) > NucleosomeBoundaries(ChosenNuc,2))
                        firstPostChosenNucFACT = ictr;
                    end
                end
                %move the post ChosenNucs FACTs and add a new one
                if firstPostChosenNucFACT > 0
                    for ictr=MaxFACTNo:-1:(firstPostChosenNucFACT+1)
                        FACTStates(ictr,:) = FACTStates(ictr-1,:);
                    end
                    FACTStates(firstPostChosenNucFACT,:) = ...
                        [NucleosomeBoundaries(ChosenNuc,:) NucFootprint 1];
                    FACTNo = FACTNo + 1;
                    ExtNo = ExtNo + 1;
                    RelaxationInitiatedNo = RelaxationInitiatedNo + 1;
                    %NB binding FACTs start off straight in the extending
                    %phase!
                else %if there are none past add in the last slot
                    FACTNo = FACTNo + 1;
                    FACTStates(FACTNo,:) = ...
                        [NucleosomeBoundaries(ChosenNuc,:) NucFootprint 1];
                    ExtNo = ExtNo + 1;
                    RelaxationInitiatedNo = RelaxationInitiatedNo + 1;
                end
                
            end
            
            
            
            
        elseif diceRoll < (FBRN+ExtNo*FERN)
            %Attempt an extension reaction
            %Choose an extending FACT
            ChosenFACT = find(FACTStates(:,4)==1);
            ChosenFACT = ChosenFACT(randi(ExtNo));
            
            %Check for ability to extend
            if (ChosenFACT < FACTNo)
                
                if (FACTStates(ChosenFACT,3) < MaxFACTSize)
                    if (FACTStates(ChosenFACT,2) < (FACTStates(ChosenFACT+1,1)-1))
                        %extend
                        FACTStates(ChosenFACT,2) = FACTStates(ChosenFACT,2) + 1;
                        %update the length
                        FACTStates(ChosenFACT,3) = FACTStates(ChosenFACT,3) + 1;
                    end
                else
                    %Become Fully Extended
                    FACTStates(ChosenFACT,4) = 2;
                    ExtNo = ExtNo - 1;
                    FullyExtNo = FullyExtNo + 1;
                end
                
            else
                %Final FACT has to check against end of gene not other FACT
                if (FACTStates(ChosenFACT,3) < MaxFACTSize) && ...
                        (FACTStates(ChosenFACT,2) < GeneLength)
                    %extend by one
                    FACTStates(ChosenFACT,2) = FACTStates(ChosenFACT,2) + 1;
                    %update length
                    FACTStates(ChosenFACT,3) = FACTStates(ChosenFACT,3) + 1;
                else
                    %Become Fully Extended
                    FACTStates(ChosenFACT,4) = 2;
                    ExtNo = ExtNo - 1;
                    FullyExtNo = FullyExtNo + 1;
                end

            end
            
            
            
        elseif diceRoll < (FBRN+(ExtNo*FERN)+(RetNo*FRRN))
            %Attempt a retraction reaction
            %Choose a retracting FACT
            ChosenFACT = find(FACTStates(:,4)==0);
            ChosenFACT = ChosenFACT(randi(RetNo));
            if (FACTStates(ChosenFACT,3) == MinFACTSize)
                %if not at end of gene, become fully retracted, otherwise
                %terminate
                if (FACTStates(ChosenFACT,2) < GeneLength)
                    %Become Fully Retracted
                    FACTStates(ChosenFACT,4) = -1;
                    RetNo = RetNo - 1;
                    FullyRetNo = FullyRetNo + 1;
                else
                    %'terminate' the FACT
                    FACTStates(ChosenFACT,:) = 0;
                    RetNo = RetNo - 1;
                    FACTNo = FACTNo - 1;
                    RelaxationFinishedNo = RelaxationFinishedNo + 1;
                end
            else
                %retract by one
                FACTStates(ChosenFACT,1) = FACTStates(ChosenFACT,1) + 1;
                %update length
                FACTStates(ChosenFACT,3) = FACTStates(ChosenFACT,3) - 1;
                
            end
            
        elseif diceRoll < (FBRN+(ExtNo*FERN)+(RetNo*FRRN)+(FullyExtNo*FERSN))
            %One of the fully extended FACTs transitions to the retracting
            %phase
            %Choose a fully extended FACT
            ChosenFACT = find(FACTStates(:,4)==2);
            ChosenFACT = ChosenFACT(randi(FullyExtNo));
            FACTStates(ChosenFACT,4) = 0;
            FullyExtNo = FullyExtNo - 1;
            RetNo = RetNo + 1;
            
            
        elseif diceRoll < (FBRN+(ExtNo*FERN)+(RetNo*FRRN)+(FullyExtNo*FERSN)+(FullyRetNo*FRESN))
            %One of the fully retracted FACTs transitions to the extending
            %phase
            %Choose a fully retracted FACT
            ChosenFACT = find(FACTStates(:,4)==-1);
            ChosenFACT = ChosenFACT(randi(FullyRetNo));
            FACTStates(ChosenFACT,4) = 1;
            FullyRetNo = FullyRetNo - 1;
            ExtNo = ExtNo + 1;
            
            
        elseif diceRoll < (FBRN+(ExtNo*FERN)+(RetNo*FRRN)+(FullyExtNo*FERSN)+(FullyRetNo*FRESN)+(FACTNo*FUBRN))
            %select a random FACT and unbind
            ChosenFACT = find(FACTStates(:,1) ~= 0);
            ChosenFACT = ChosenFACT(randi(sum(FACTStates(:,1)~=0)));
            if (FACTStates(ChosenFACT,4) == -1)
                FACTStates(ChosenFACT,:) = [0 0 0 0];
                FullyRetNo = FullyRetNo - 1;
                for ictr = ChosenFACT:1:(FACTNo-1)
                    FACTStates(ictr,:) = ...
                        FACTStates(ictr+1,:);
                    FACTStates(ictr+1,:) = [0 0 0 0];
                end
                FACTNo = FACTNo - 1;
                RelaxationEarlyFinishedNo = RelaxationEarlyFinishedNo + 1;
                
            elseif (FACTStates(ChosenFACT,4) == 0)
                FACTStates(ChosenFACT,:) = [0 0 0 0];
                RetNo = RetNo - 1;
                for ictr = ChosenFACT:1:(FACTNo-1)
                    FACTStates(ictr,:) = ...
                        FACTStates(ictr+1,:);
                    FACTStates(ictr+1,:) = [0 0 0 0];
                end
                FACTNo = FACTNo - 1;
                RelaxationEarlyFinishedNo = RelaxationEarlyFinishedNo + 1;
                
            elseif (FACTStates(ChosenFACT,4) == 1)
                FACTStates(ChosenFACT,:) = [0 0 0 0];
                ExtNo = ExtNo - 1;
                for ictr = ChosenFACT:1:(FACTNo-1)
                    FACTStates(ictr,:) = ...
                        FACTStates(ictr+1,:);
                    FACTStates(ictr+1,:) = [0 0 0 0];
                end
                FACTNo = FACTNo - 1;
                RelaxationEarlyFinishedNo = RelaxationEarlyFinishedNo + 1;
                
            elseif (FACTStates(ChosenFACT,4) == 2)
                FACTStates(ChosenFACT,:) = [0 0 0 0];
                FullyExtNo = FullyExtNo - 1;
                for ictr = ChosenFACT:1:(FACTNo-1)
                    FACTStates(ictr,:) = ...
                        FACTStates(ictr+1,:);
                    FACTStates(ictr+1,:) = [0 0 0 0];
                end
                FACTNo = FACTNo - 1;
                RelaxationEarlyFinishedNo = RelaxationEarlyFinishedNo + 1;
                
            end
            
            
        end
        
        
    end
    
    
    %----Simulate the System post Relaxation------------------------------
    %----Sample at intermediate points------------------------------------
    
    
    for tctr=1:1:timeSteps
        
        %Generate a pseudo random number from a uniform distribution
        diceRoll = rand;
        
        %Use the pseudo random number to test for success of potential
        %events
        
        if diceRoll < FBRN
            %Attempt a new FACT binding at a randomly chosen nucleosome
            
            %First choose a potential binding site
            ChosenNuc = randi(NucNo);
            %Then bind if that site is currently unbound

            NucStart = NucleosomeBoundaries(ChosenNuc,1);
            NucEnd = NucleosomeBoundaries(ChosenNuc,2);
            NucleosomeOverlap = 0;
            for ictr=1:1:MaxFACTNo
                if ( (FACTStates(ictr,1) >= NucStart) && (FACTStates(ictr,1) <= NucEnd) ) ...
                        || ( (FACTStates(ictr,2) >= NucStart) && (FACTStates(ictr,2) <= NucEnd) )
                    NucleosomeOverlap = 1;
                end
            end
            firstPostChosenNucFACT = 0;
            if NucleosomeOverlap == 0
                %find out which (if any) FACTs need moving in the queue
                for ictr=MaxFACTNo:-1:1
                    if(FACTStates(ictr,1) > NucleosomeBoundaries(ChosenNuc,2))
                        firstPostChosenNucFACT = ictr;
                    end
                end
                %move the postChosenNucs FACTs and add a new one
                if firstPostChosenNucFACT > 0
                    for ictr=MaxFACTNo:-1:(firstPostChosenNucFACT+1)
                        FACTStates(ictr,:) = FACTStates(ictr-1,:);
                    end
                    FACTStates(firstPostChosenNucFACT,:) = ...
                        [NucleosomeBoundaries(ChosenNuc,:) NucFootprint 1];
                    FACTNo = FACTNo + 1;
                    ExtNo = ExtNo + 1;
                    SteadyStateInitiatedNo = SteadyStateInitiatedNo + 1;
                    %NB binding FACTs start off straight in the extending
                    %phase!
                else %if there are none past add in the last slot
                    FACTNo = FACTNo + 1;
                    FACTStates(FACTNo,:) = ...
                        [NucleosomeBoundaries(ChosenNuc,:) NucFootprint 1];
                    ExtNo = ExtNo + 1;
                    SteadyStateInitiatedNo = SteadyStateInitiatedNo + 1;
                    
                end
                
            end
            
            
        elseif diceRoll < (FBRN+ExtNo*FERN)
            %Attempt an extension reaction
            %Choose an extending FACT
            ChosenFACT = find(FACTStates(:,4)==1);
            ChosenFACT = ChosenFACT(randi(ExtNo));
            
            %Check for ability to extend
            if (ChosenFACT < FACTNo)
                
                if (FACTStates(ChosenFACT,3) < MaxFACTSize)
                    if (FACTStates(ChosenFACT,2) < (FACTStates(ChosenFACT+1,1)-1))
                        %extend
                        FACTStates(ChosenFACT,2) = FACTStates(ChosenFACT,2) + 1;
                        %update the length
                        FACTStates(ChosenFACT,3) = FACTStates(ChosenFACT,3) + 1;
                    end
                else
                    %Become Fully Extended
                    FACTStates(ChosenFACT,4) = 2;
                    ExtNo = ExtNo - 1;
                    FullyExtNo = FullyExtNo + 1;
                end
                
            else
                %Final FACT has to check against end of gene not other FACT
                if (FACTStates(ChosenFACT,3) < MaxFACTSize) && ...
                        (FACTStates(ChosenFACT,2) < GeneLength)
                    %extend by one
                    FACTStates(ChosenFACT,2) = FACTStates(ChosenFACT,2) + 1;
                    %update length
                    FACTStates(ChosenFACT,3) = FACTStates(ChosenFACT,3) + 1;
                else
                    %Become Fully Extended
                    FACTStates(ChosenFACT,4) = 2;
                    ExtNo = ExtNo - 1;
                    FullyExtNo = FullyExtNo + 1;
                end
                
            end
            
        elseif diceRoll < (FBRN+(ExtNo*FERN)+(RetNo*FRRN))
            %Attempt a retraction reaction
            %Choose a retracting FACT
            ChosenFACT = find(FACTStates(:,4)==0);
            ChosenFACT = ChosenFACT(randi(RetNo));
            if (FACTStates(ChosenFACT,3) == MinFACTSize)
                %if not at end of gene, transition to extending, otherwise
                %terminate
                if (FACTStates(ChosenFACT,2) < GeneLength)
                    %Become Fully Retracted
                    FACTStates(ChosenFACT,4) = -1;
                    RetNo = RetNo - 1;
                    FullyRetNo = FullyRetNo + 1;
                else
                    %'terminate' the FACT
                    FACTStates(ChosenFACT,:) = 0;
                    RetNo = RetNo - 1;
                    FACTNo = FACTNo - 1;
                    SteadyStateFinishedNo = SteadyStateFinishedNo + 1;
                end
            else
                %retract by one
                FACTStates(ChosenFACT,1) = FACTStates(ChosenFACT,1) + 1;
                %update length
                FACTStates(ChosenFACT,3) = FACTStates(ChosenFACT,3) - 1;
                
            end
            
            
            
        elseif diceRoll < (FBRN+(ExtNo*FERN)+(RetNo*FRRN)+(FullyExtNo*FERSN))
            %One of the fully extended FACTs transitions to the retracting
            %phase
            %Choose a fully extended FACT
            ChosenFACT = find(FACTStates(:,4)==2);
            ChosenFACT = ChosenFACT(randi(FullyExtNo));
            FACTStates(ChosenFACT,4) = 0;
            FullyExtNo = FullyExtNo - 1;
            RetNo = RetNo + 1;
            
        elseif diceRoll < (FBRN+(ExtNo*FERN)+(RetNo*FRRN)+(FullyExtNo*FERSN)+(FullyRetNo*FRESN))
            %One of the fully retracted FACTs transitions to the extending
            %phase
            %Choose a fully retracted FACT
            ChosenFACT = find(FACTStates(:,4)==-1);
            ChosenFACT = ChosenFACT(randi(FullyRetNo));
            FACTStates(ChosenFACT,4) = 1;
            FullyRetNo = FullyRetNo - 1;
            ExtNo = ExtNo + 1;
            
            
        elseif diceRoll < (FBRN+(ExtNo*FERN)+(RetNo*FRRN)+(FullyExtNo*FERSN)+(FullyRetNo*FRESN)+(FACTNo*FUBRN))
            %select a random FACT to unbind
            ChosenFACT = find(FACTStates(:,1) ~= 0);
            ChosenFACT = ChosenFACT(randi(sum(FACTStates(:,1)~=0)));
            if (FACTStates(ChosenFACT,4) == -1)
                FACTStates(ChosenFACT,:) = [0 0 0 0];
                FullyRetNo = FullyRetNo - 1;
                for ictr = ChosenFACT:1:(FACTNo-1)
                    FACTStates(ictr,:) = ...
                        FACTStates(ictr+1,:);
                    FACTStates(ictr+1,:) = [0 0 0 0];
                end
                FACTNo = FACTNo - 1;
                SteadyStateEarlyFinishedNo = SteadyStateEarlyFinishedNo + 1;

            elseif (FACTStates(ChosenFACT,4) == 0)
                FACTStates(ChosenFACT,:) = [0 0 0 0];
                RetNo = RetNo - 1;
                for ictr = ChosenFACT:1:(FACTNo-1)
                    FACTStates(ictr,:) = ...
                        FACTStates(ictr+1,:);
                    FACTStates(ictr+1,:) = [0 0 0 0];
                end
                FACTNo = FACTNo - 1;
                SteadyStateEarlyFinishedNo = SteadyStateEarlyFinishedNo + 1;
                
            elseif (FACTStates(ChosenFACT,4) == 1)
                FACTStates(ChosenFACT,:) = [0 0 0 0];
                ExtNo = ExtNo - 1;
                for ictr = ChosenFACT:1:(FACTNo-1)
                    FACTStates(ictr,:) = ...
                        FACTStates(ictr+1,:);
                    FACTStates(ictr+1,:) = [0 0 0 0];
                end
                FACTNo = FACTNo - 1;
                SteadyStateEarlyFinishedNo = SteadyStateEarlyFinishedNo + 1;
                
            elseif (FACTStates(ChosenFACT,4) == 2)
                FACTStates(ChosenFACT,:) = [0 0 0 0];
                FullyExtNo = FullyExtNo - 1;
                for ictr = ChosenFACT:1:(FACTNo-1)
                    FACTStates(ictr,:) = ...
                        FACTStates(ictr+1,:);
                    FACTStates(ictr+1,:) = [0 0 0 0];
                end
                FACTNo = FACTNo - 1;
                SteadyStateEarlyFinishedNo = SteadyStateEarlyFinishedNo + 1;
                
            end
            
            
        end
        
        
        %Sample the protected regions of DNA:
        
        abuttingFACT = [];
        edgePoints = reshape(FACTStates(1:FACTNo,1:2)',[1 2*FACTNo]);
        for sctr=1:1:FACTNo-1
            if (edgePoints(2*sctr+1)-edgePoints(2*sctr) == 1)
                abuttingFACT = [abuttingFACT sctr];
            end
        end
        if ~isempty(abuttingFACT)
            for sctr=1:1:length(abuttingFACT)
                edgePoints(2*abuttingFACT(sctr):(2*abuttingFACT(sctr)+1)) = [];
                abuttingFACT = abuttingFACT - 1;
            end
        end

        
        %Process Location and Length from edgepoints
        
        densityLocs = round(mean([edgePoints(2:2:length(edgePoints)); ...
            edgePoints(1:2:length(edgePoints))]));
        
        densityLengths = edgePoints(2:2:length(edgePoints)) - ...
            edgePoints(1:2:length(edgePoints));
        
        for sctr=1:1:length(densityLocs)
            
            MNaseDensity(densityLengths(sctr),densityLocs(sctr)) = ...
                MNaseDensity(densityLengths(sctr),densityLocs(sctr)) + 1;
            
        end
        
        
        
        
        
    end
    
    disp(['Completed Realisation ' num2str(rctr) ...
        ' of ' num2str(Realisations) ...
        ' in ' num2str(toc(realtics)) ' secs.'])
    
    
    
end





%% Block Sum Plot of Density Heatmap

figpos = [200  200  560  420];

mfig = figure('position',figpos);


avgingWindow = 20;
maxFragSizeToShow = 400;

tl = tiledlayout(1,1);

nt = nexttile;

p = pcolor(conv2(ones(avgingWindow),MNaseDensity(:,:))/(avgingWindow^2));
set(p,'EdgeColor','None')


set(nt,'ydir','normal')
set(nt,'ylim',[50 maxFragSizeToShow])
set(nt,'xlim',[1 1000+1+MinFACTSize])

set(nt,'ytick',[50 100 150 200 250 300 350 400]);

set(nt,'xtick',[0 500 1000]+MinFACTSize/2);
set(nt,'xticklabels',{'0','0.5','1'});
 
 
 
colorbar





















