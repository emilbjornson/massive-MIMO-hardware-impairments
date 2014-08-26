%This Matlab script can be used to generate Figure 14 in the article:
%
%Emil Björnson, Jakob Hoydis, Marios Kountouris, Mérouane Debbah, “Massive
%MIMO Systems with Non-Ideal Hardware: Energy Efficiency, Estimation, and
%Capacity Limits,” To appear in IEEE Transactions on Information Theory.
%
%Download article: http://arxiv.org/pdf/1307.2584
%
%This is version 1.0 (Last edited: 2014-08-26)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.

%Initialization
close all;
clear all;


%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));

%Range of number of BS antennas
Nantennas = [1:5 10:5:100 150:50:1000];

maxN = max(Nantennas); %Maximal number of antennas in simulation

BSs = 16; %Number of sites (one cell under study surrounded by 15 interfering cells)
Ksite = 6; %Number of users per site

%Define the level of hardware impairments at the UE and BS
kappaUE = 0.1^2;
kappaBS = 0.1^2;

%Portion of all symbols used for either DL or UL data transmission.
trafficPortion = 0.45;

%Bandwidth of the flat-fading subcarrier
subcarrierBandwidth = 15e3; %Hz in the LTE-like system

%Transmit power (in mW) per UE over the subcarrier bandwidth (this
%corresponds to 200 mW over a 9 MHz bandwidth or 0.0222 microJ/channel use
pUE = 0.0222e-3*subcarrierBandwidth;

%Receiver noise properties
noiseFigure = 5; %Noise amplification in the receiver hardware in dB
noiseFloordBm = -174+10*log10(subcarrierBandwidth)+noiseFigure; %Noise floor in dBm

%Define the size of the square cells in the simulation
intersiteDistance = 0.4; %Distance between the middle of two adjacent cells (in vertical or horizontal direction)
intersiteDistanceHalf = intersiteDistance/2; %Half the inter-cell distance
lengthOfRegion = intersiteDistance*sqrt(BSs); %This is the length of the region where all the BSs are placed

%Generate grid locations for the cells
locationHorizontal = repmat(-3*intersiteDistance/2 : intersiteDistance : 3*intersiteDistance/2,[4 1]);
locationVertical = locationHorizontal';
BSlocations = locationHorizontal(:) + 1i*locationVertical(:);  %Real part is horizontal coordinate and imaginary part is vertical coordinate

%Distance between UEs and their serving BS
userDistance = 0.1;

%Put the users equally spaced on a circle
rotations = exp(1i*(0:2*pi/Ksite:2*pi-2*pi/Ksite)); %Angles of the equally spaced users
UElocations = (userDistance*ones(BSs,Ksite)) .* repmat(rotations,[BSs 1]); %Put out the users on the circle
UElocations = UElocations + repmat(BSlocations,[1 Ksite]); %Finalize the user location by putting the users around the serving BS

%Compute distances between each UE and BS in a system with wrap-around
distances = zeros(BSs,BSs,Ksite); %(j,l,k): Channel from BS_j to UE_k in Cell l.

for j = 1:BSs
    
    %Compute the shortest distance between a UE and each of the users. In a
    %system with wrap-around, there are four different locations that each
    %BS can have and these are tested here.
    move = reshape([0 sign(real(BSlocations(j)))*lengthOfRegion sign(imag(BSlocations(j)))*1i*lengthOfRegion sign(real(BSlocations(j)))*lengthOfRegion+sign(imag(BSlocations(j)))*1i*lengthOfRegion],[1 1 4]);
    distances(j,:,:)= min(abs(repmat(UElocations - BSlocations(j),[1 1 4]) + repmat(move,[BSs,Ksite 1])),[],3);
    
end

%Compute the ratio between pathloss and noise for every user location
pathlossesOverNoise = 10.^( -(128.1+37.6*log10(distances) + noiseFloordBm)/10);
pathlossesOverNoiseSqrt = sqrt(pathlossesOverNoise);


%%Initialize Monte Carlo simulations

%Number realizations in Monte Carlo simulations
nbrOfMonteCarloRealizations = 10000;

%Index of the cell that is studied in the simulations (it is
%computationally efficient to only study one of the cells, and we can
%choose any cell due to symmetry.
jBSstudied = 1;


%Placeholders for storing simulation results, with either hardware
%impairments or ideal hardware, and with either unique pilots or reuse

firstMoment_impair_unique = zeros(length(Nantennas),Ksite);
secondMoments_impair_unique = zeros(length(Nantennas),Ksite);
distortionTerm_impair_unique = zeros(length(Nantennas),Ksite);

firstMoment_impair_reuse = zeros(length(Nantennas),Ksite);
secondMoments_impair_reuse = zeros(length(Nantennas),Ksite);
distortionTerm_impair_reuse = zeros(length(Nantennas),Ksite);

firstMoment_ideal_unique = zeros(length(Nantennas),Ksite);
secondMoments_ideal_unique = zeros(length(Nantennas),Ksite);

firstMoment_ideal_reuse = zeros(length(Nantennas),Ksite);
secondMoments_ideal_reuse = zeros(length(Nantennas),Ksite);


%Go through all Monte Carlo realizations
for iter = 1:nbrOfMonteCarloRealizations
    
    disp(['Iteration ' num2str(iter) '/' num2str(nbrOfMonteCarloRealizations)]);
    
    %Generate random realizations
    H = (randn(maxN,1,BSs,Ksite)+1i*randn(maxN,1,BSs,Ksite))/sqrt(2); %Generate channel realizations
    etaUEs = sqrt(kappaUE*pUE) * (randn(BSs,1,Ksite)+1i*randn(BSs,1,Ksite))/sqrt(2); %Generate distortion noise at UEs
    etaBS = sqrt(kappaBS*pUE) * (randn(maxN,1,Ksite)+1i*randn(maxN,1,Ksite))/sqrt(2); %Generate distortion noise at BS
    noise = (randn(maxN,1,Ksite)+1i*randn(maxN,1,Ksite))/sqrt(2);
    
    %Prepare to store channel estimates
    Hhat_impair_unique = zeros(maxN,Ksite,BSs);
    Hhat_impair_reuse = zeros(maxN,Ksite,BSs);
    Hhat_ideal_unique = zeros(maxN,Ksite,BSs);
    Hhat_ideal_reuse = zeros(maxN,Ksite,BSs);

    
    %Go through all UEs in the cell that is studied
    for k = 1:Ksite
        
        %Extract the gain of the channel from the served UE and the sum of
        %the gains of the channels from all active UEs with pilot reuse
        ownChannelGain = pathlossesOverNoise(jBSstudied,jBSstudied,k);
        sumAllChannelGains = sum(pathlossesOverNoise(jBSstudied,:,k));

        %Compute matrix A in the LMMSE estimator (see Eq. (9)), including
        %the interference from pilot reuse when this is applicable
        A_impair_reuse = sqrt(pUE)* ownChannelGain / (    pUE*(1 + kappaUE)*sumAllChannelGains + pUE*kappaBS*sumAllChannelGains + 1);
        A_impair_unique = sqrt(pUE) * ownChannelGain / (    pUE*(1 + kappaUE)*ownChannelGain + pUE*kappaBS*ownChannelGain + 1);
        A_ideal_reuse = sqrt(pUE) * ownChannelGain / (    pUE*sumAllChannelGains +1);
        A_ideal_unique = sqrt(pUE) * ownChannelGain / (    pUE*ownChannelGain + 1);
        
        %Calculate channel estimates from pilot signaling
        Hhat_impair_unique(:,k,jBSstudied) = A_impair_unique * ( H(:,jBSstudied,jBSstudied,k) * ((sqrt(pUE) + etaUEs(jBSstudied,1,k) )*pathlossesOverNoiseSqrt(jBSstudied,jBSstudied,k)) + abs(H(:,jBSstudied,jBSstudied,k)) .* etaBS(:,1,k) * pathlossesOverNoiseSqrt(jBSstudied,jBSstudied,k) + noise(:,1,k) );
        Hhat_impair_reuse(:,k,jBSstudied) = A_impair_reuse * ( reshape(H(:,jBSstudied,:,k),[maxN BSs]) * (sqrt(pUE)*pathlossesOverNoiseSqrt(jBSstudied,:,k)' + etaUEs(:,1,k).*(pathlossesOverNoiseSqrt(jBSstudied,:,k)') ) + sqrt((abs(reshape(H(:,jBSstudied,:,k),[maxN BSs])).^2)*pathlossesOverNoiseSqrt(jBSstudied,:,k)') .* etaBS(:,1,k) + noise(:,1,k) );
        Hhat_ideal_unique(:,k,jBSstudied) = A_ideal_unique * ( H(:,jBSstudied,jBSstudied,k) * (sqrt(pUE)*pathlossesOverNoiseSqrt(jBSstudied,jBSstudied,k))  + noise(:,1,k) );
        Hhat_ideal_reuse(:,k,jBSstudied) = A_ideal_reuse * ( (reshape(H(:,jBSstudied,:,k),[maxN BSs])*pathlossesOverNoiseSqrt(jBSstudied,:,k)') * sqrt(pUE) + noise(:,1,k) );
        
    end
    
    
    %Go through different number of antennas
    for n = 1:length(Nantennas)
        
        %Extract current number of antennas
        N = Nantennas(n);
        
        %Extract the useful channels between the BS and its UEs
        usefulChannels = reshape(H(1:N,jBSstudied,jBSstudied,:),[N Ksite])*sparse(diag(reshape(pathlossesOverNoiseSqrt(jBSstudied,jBSstudied,:),[Ksite 1])));
        
        %Compute the average strength of the intercell interference
        intercellInterference = pUE*sum(sum(pathlossesOverNoise(jBSstudied,[1:jBSstudied-1 jBSstudied+1:BSs],:),2),3);
        
        
        
        %%Compute relizations of the expectations in Eq. (37) for different
        %%setups (ideal hardware or hardware impairments, unique pilots or
        %%pilot reuse across cells).
        
        
        %Case: Hardware impairments, unique pilots
        estimatedChannel = Hhat_impair_unique(1:N,:,jBSstudied);
        
        %Compute MMSE receive beamforming
        V_MMSE = ((1+intercellInterference)*eye(N) + pUE*(estimatedChannel*estimatedChannel'))\estimatedChannel;
        V_MMSE = V_MMSE ./ repmat(sqrt(sum(abs(V_MMSE).^2,1)),[N 1]);
        
        %Compute the received useful channel gains after MMSE receive
        %beamforming filters have been applied
        filteredChannelGains = usefulChannels'*V_MMSE;
        
        %Store a realization of the first moment of the inner product
        %between beamforming and useful channel
        firstMoment_impair_unique(n,:) = firstMoment_impair_unique(n,:) + diag(filteredChannelGains)';
        
        %Compute and store a realization of the second moment of the inner
        %product between beamforming vectors and channels, both for useful
        %and interfering channels
        for m = 1:BSs
            secondMoments_impair_unique(n,:) = secondMoments_impair_unique(n,:) + sum(sparse(diag(reshape(pathlossesOverNoise(jBSstudied,m,:),[Ksite 1])))*abs(reshape(H(1:N,jBSstudied,m,:),[N,Ksite])'*V_MMSE).^2,1);
        end
        
        %The elementwise product between channel and beamforming
        %vectors (and sum over these elements) appear in the
        %distortion noise term. This computes a realization.
        distortionTerm_impair_unique(n,:) = distortionTerm_impair_unique(n,:) + (reshape(abs(H(1:N,jBSstudied,:,:)).^2,[N BSs*Ksite])*reshape(pathlossesOverNoise(jBSstudied,:,:),[BSs*Ksite 1]))' * abs(V_MMSE).^2;
        
        
        
        %Case: Hardware impairments, pilot reuse
        estimatedChannel = Hhat_impair_reuse(1:N,:,jBSstudied);
        
        %Compute MMSE receiver filter
        V_MMSE = ((1+intercellInterference)*eye(N) + pUE*(estimatedChannel*estimatedChannel'))\estimatedChannel;
        V_MMSE = V_MMSE ./ repmat(sqrt(sum(abs(V_MMSE).^2,1)),[N 1]);
        
        %Compute the received useful channel gains after MMSE receive
        %beamforming filters have been applied
        filteredChannelGains = usefulChannels'*V_MMSE;
        
        %Store a realization of the first moment of the inner product
        %between beamforming and useful channel
        firstMoment_impair_reuse(n,:) = firstMoment_impair_reuse(n,:) + diag(filteredChannelGains)';
        
        %Compute and store a realization of the second moment of the inner
        %product between beamforming vectors and channels, both for useful
        %and interfering channels
        for m = 1:BSs
            secondMoments_impair_reuse(n,:) = secondMoments_impair_reuse(n,:) + sum(sparse(diag(reshape(pathlossesOverNoise(jBSstudied,m,:),[Ksite 1])))*abs(reshape(H(1:N,jBSstudied,m,:),[N,Ksite])'*V_MMSE).^2,1);
        end
        
        %The elementwise product between channel and beamforming
        %vectors (and sum over these elements) appear in the
        %distortion noise term. This computes and stores a realization.
        distortionTerm_impair_reuse(n,:) = distortionTerm_impair_reuse(n,:) + (reshape(abs(H(1:N,jBSstudied,:,:)).^2,[N BSs*Ksite])*reshape(pathlossesOverNoise(jBSstudied,:,:),[BSs*Ksite 1]))' * abs(V_MMSE).^2;
        
        
        
        %Case: Ideal hardware, unique pilots
        estimatedChannel = Hhat_ideal_unique(1:N,:,jBSstudied);
        
        %Compute MMSE receiver filter
        V_MMSE = ((1+intercellInterference)*eye(N) + pUE*(estimatedChannel*estimatedChannel'))\estimatedChannel;
        V_MMSE = V_MMSE ./ repmat(sqrt(sum(abs(V_MMSE).^2,1)),[N 1]);
        
        %Compute the received useful channel gains after MMSE receive 
        %beamforming filters have been applied
        filteredChannelGains = usefulChannels'*V_MMSE;
        
        %Store a realization of the first moment of the inner product
        %between beamforming and useful channel
        firstMoment_ideal_unique(n,:) = firstMoment_ideal_unique(n,:) + diag(filteredChannelGains)';
        
        %Compute and store a realization of the second moment of the inner
        %product between beamforming vectors and channels, both for useful
        %and interfering channels
        for m = 1:BSs
            secondMoments_ideal_unique(n,:) = secondMoments_ideal_unique(n,:) + sum(sparse(diag(reshape(pathlossesOverNoise(jBSstudied,m,:),[Ksite 1])))*abs(reshape(H(1:N,jBSstudied,m,:),[N,Ksite])'*V_MMSE).^2,1);
        end
        

        
        %Case: Ideal hardware, pilot reuse
        estimatedChannel = Hhat_ideal_reuse(1:N,:,jBSstudied);
        
        %Compute MMSE receiver filter
        V_MMSE = ((1+intercellInterference)*eye(N) + pUE*(estimatedChannel*estimatedChannel'))\estimatedChannel;
        V_MMSE = V_MMSE ./ repmat(sqrt(sum(abs(V_MMSE).^2,1)),[N 1]);
        
        %Compute the received useful channel gains after MMSE receive
        %beamforming filters have been applied
        filteredChannelGains = usefulChannels'*V_MMSE;
        
        %Store a realization of the first moment of the inner product
        %between beamforming and useful channel
        firstMoment_ideal_reuse(n,:) = firstMoment_ideal_reuse(n,:) + diag(filteredChannelGains)';
        
        %Compute and store a realization of the second moment of the inner
        %product between beamforming vectors and channels, both for useful
        %and interfering channels
        for m = 1:BSs
            secondMoments_ideal_reuse(n,:) = secondMoments_ideal_reuse(n,:) + sum(sparse(diag(reshape(pathlossesOverNoise(jBSstudied,m,:),[Ksite 1])))*abs(reshape(H(1:N,jBSstudied,m,:),[N,Ksite])'*V_MMSE).^2,1);
        end

    end
    
end


%Compute lower bounds on the achievable rates of all the UEs in the cell
rates_impair_unique = real(log2(1+ pUE*abs(firstMoment_impair_unique/nbrOfMonteCarloRealizations).^2 ./ ( pUE*(1+kappaUE)*secondMoments_impair_unique/nbrOfMonteCarloRealizations + pUE*kappaBS*distortionTerm_impair_unique/nbrOfMonteCarloRealizations - pUE*abs(firstMoment_impair_unique/nbrOfMonteCarloRealizations).^2 + 1) ));
rates_impair_reuse = real(log2(1+ pUE*abs(firstMoment_impair_reuse/nbrOfMonteCarloRealizations).^2 ./ ( pUE*(1+kappaUE)*secondMoments_impair_reuse/nbrOfMonteCarloRealizations + pUE*kappaBS*distortionTerm_impair_reuse/nbrOfMonteCarloRealizations - pUE*abs(firstMoment_impair_reuse/nbrOfMonteCarloRealizations).^2 + 1) ));
rates_ideal_unique = real(log2(1+ pUE*abs(firstMoment_ideal_unique/nbrOfMonteCarloRealizations).^2 ./ ( pUE*secondMoments_ideal_unique/nbrOfMonteCarloRealizations - pUE*abs(firstMoment_ideal_unique/nbrOfMonteCarloRealizations).^2 + 1) ));
rates_ideal_reuse = real(log2(1+ pUE*abs(firstMoment_ideal_reuse/nbrOfMonteCarloRealizations).^2 ./ ( pUE*secondMoments_ideal_reuse/nbrOfMonteCarloRealizations - pUE*abs(firstMoment_ideal_reuse/nbrOfMonteCarloRealizations).^2 + 1) ));


%Plot Figure 14 from the paper

figure; hold on; box on;

plot(Nantennas,trafficPortion*mean(rates_ideal_unique,2),'b-.');
plot(Nantennas,trafficPortion*mean(rates_ideal_reuse,2),'b--');

plot(Nantennas,trafficPortion*mean(rates_impair_unique,2),'k-');
plot(Nantennas,trafficPortion*mean(rates_impair_reuse,2),'k:');

legend('Ideal: Unique Pilots','Ideal: Inter-Cell Pilot Reuse','Non-ideal: Unique Pilots','Non-ideal: Inter-Cell Pilot Reuse','Location','SouthEast')
xlabel('Number of Antennas per BS (N)')
ylabel('Average Spectral Efficiency per User [bit/channel use]')

axis([0 maxN 0 5]);
