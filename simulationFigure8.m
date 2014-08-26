%This Matlab script can be used to generate Figure 8 in the article:
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
N = [1:50 60:10:150 175:25:1000];

Nmax = max(N); %Maximal number of antennas in simulation

%Define the normalized channel covariance matrix R as an identity matrix
R = eye(Nmax);

%Define the range of level of hardware impairments at the transmitter and
%receiver. The values are the same for BS and UE, but can be different.
kappaBS = 0.05^2;
kappaUE = 0.05^2;

%Portion of all symbols used for either DL or UL data transmission.
trafficPortion = 0.45;

%Range of SNRs in simulation (we have normalized sigma2 to 1)
SNRdB = 20; %In decibel scale
SNR = 10.^(SNRdB/10); %In linear scale



%%Initialize Monte Carlo simulations

%Number realizations in Monte Carlo simulations
nbrOfMonteCarloRealizations = 10000;

%Number of different channel covariance models
nbrOfCovarianceModels = 4;

%Generate random realizations
huncorrelated = sqrtm(R)*(randn(Nmax,nbrOfMonteCarloRealizations)+1i*randn(Nmax,nbrOfMonteCarloRealizations))/sqrt(2); %Generate channel realizations
etaUE = (randn(1,nbrOfMonteCarloRealizations)+1i*randn(1,nbrOfMonteCarloRealizations))/sqrt(2); %Generate distortion noise at UE (to be scaled by kappa and SNR)
etaBSuncorrelated = (randn(Nmax,nbrOfMonteCarloRealizations)+1i*randn(Nmax,nbrOfMonteCarloRealizations))/sqrt(2); %Generate distortion noise at BS (to be scaled by kappa, SNR and channel)
noise = (randn(Nmax,nbrOfMonteCarloRealizations)+1i*randn(Nmax,nbrOfMonteCarloRealizations))/sqrt(2);


%Initialize the matrices for asymptotic capacity limits
asymptoticLimitLower = zeros(length(kappaBS),length(SNR));
asymptoticLimitUpper = log2(1+1./kappaUE); %This is only well-defined for kappa>0

%Placeholders for storing simulation results
capacityLowerBound = zeros(length(N),nbrOfCovarianceModels,length(SNR));
capacityUpperBound = zeros(length(N),nbrOfCovarianceModels,length(SNR));



%Go through the different channel covariance models
for c = 1:nbrOfCovarianceModels
    
    %Output the progress of the simulation
    disp(['Covariance models: ' num2str(c) '/' num2str(nbrOfCovarianceModels)]);
    
    %Generate channel covariance matrix for largest number of antennas
    if c == 1 %Uncorrelated covariance matrix
        R = eye(Nmax);
    elseif c == 2 %Exponential correlation model
        correlationFactor = 0.7;
        R = toeplitz(correlationFactor.^(0:Nmax-1));
    elseif c == 3 %One-ring model with 20 degrees angular spread
        R = functionOneRingModel(Nmax,20);
    elseif c == 4 %One-ring model with 10 degrees angular spread
        R = functionOneRingModel(Nmax,10);
    end
    
    %Modify channel and distortion noise realizations based on the current 
    %channel covariance matrix
    h = sqrtm(R)*huncorrelated;
    etaBS = abs(h) .* etaBSuncorrelated;
    
    %Go through different number of antennas
    for n = 1:length(N)
        
        %Go through SNR values
        for m = 1:length(SNR)
            
            %A typical pilot signal
            d = sqrt(SNR(m));
            
            %Compute matrix A in the LMMSE estimator (see Eq. (9))
            A_LMMSE = conj(d) * R(1:N(n),1:N(n)) / (abs(d)^2*(1+kappaUE)*R(1:N(n),1:N(n)) + abs(d)^2*kappaBS*diag(diag(R(1:N(n),1:N(n))))+eye(N(n)));
            
            %Compute the MSE
            MSE = trace(d*  A_LMMSE *R(1:N(n),1:N(n)));
            
            %Compute the lower capacity limit in Corollaries 4 and 5 for
            %the largest of the N-values
            if n == length(N)
                
                %Compute the phi in Eq. (41) by Monte-Carlo simulations
                realizationsPhi = 1000000;
                
                randomness1 = 1+sqrt(kappaUE)*(randn(1000000,1)+1i*randn(1000000,1))/sqrt(2);
                randomness2 = abs(randomness1).^2;
                
                phiFirstMomentSquared = abs(mean( randomness1 * MSE ./ sqrt(  MSE - trace(A_LMMSE*R(1:N(n),1:N(n))*A_LMMSE')*abs(d)^2*((1+kappaUE)-randomness2))))^2;
                phiSecondMoment  = real(mean( randomness2 * MSE^2 ./ (  MSE - trace(A_LMMSE*R(1:N(n),1:N(n))*A_LMMSE')*abs(d)^2*((1+kappaUE)-randomness2))));
                
                %Copmute the lower capacity limit
                asymptoticLimitLower(c,m) = log2(1+ phiFirstMomentSquared ./ ((1+kappaUE)*phiSecondMoment-phiFirstMomentSquared ) );
                
            end
            
            %Placeholders for storing Monte Carlo simulation results
            firstMoment = zeros(nbrOfMonteCarloRealizations,1);
            secondMoment = zeros(nbrOfMonteCarloRealizations,1);
            distortionTerm = zeros(nbrOfMonteCarloRealizations,1);
            boundBeamforming = zeros(nbrOfMonteCarloRealizations,1);
            
            %Go through all Monte Carlo realizations
            for k = 1:nbrOfMonteCarloRealizations
                
                %Compute received signal
                z = h(1:N(n),k) * ( d + sqrt(kappaUE*SNR(m))*etaUE(k) ) + sqrt(kappaBS*SNR(m))*etaBS(1:N(n),k) + noise(1:N(n),k);
                
                %Compute channel estimates
                hhat = A_LMMSE*z;
                
                %Compute the beamforming vector (MRT/MRC)
                beamforming = sqrt(SNR(m))*hhat/norm(hhat);
                
                %Compute a realization of the first and second moments of
                %the inner product between beamforming and channel
                firstMoment(k) = h(1:N(n),k)'*beamforming;
                secondMoment(k) = abs(firstMoment(k)).^2;
                
                %The elementwise product between channel and beamforming
                %vectors (and sum over these elements) appear in the
                %distortion noise term. This computes a realization.
                distortionTerm(k) = sum( (abs(h(1:N(n),k)).^2) .* (abs(beamforming).^2));
                
                %Compute the unnormalized beamforming vector according to
                %Eq. (24) and Eq. (25) in the upper bound
                boundBeamforming(k) = sum(abs(h(1:N(n),k)).^2 ./ ( kappaBS*abs(h(1:N(n),k)).^2+ 1/SNR(m)));
                
            end
            
            %Finalize the Monte Carlo simulations by computing lower and
            %upper bounds on the capacity
            capacityLowerBound(n,c,m) = log2(1+ abs(mean(firstMoment,1)).^2 ./ ( (1+kappaUE) * var(firstMoment) + kappaUE* abs(mean(firstMoment,1)).^2 + kappaBS*mean(distortionTerm,1) + 1) );
            capacityUpperBound(n,c,m) = mean(log2(1+ boundBeamforming ./ ( 1 + kappaUE*boundBeamforming ) )); %This is based on Lemma 1, but uses Lemma 2 to make the computation more efficient
            
        end
        
    end
    
    
end


%Plot Figure 8 from the paper

for m = 1:length(SNR)
    
    figure; hold on; box on;
    
    plot(N,trafficPortion*capacityUpperBound(:,1,m),'k--');
    

    plot(N,trafficPortion*capacityLowerBound(:,1,m),'k');
    plot(N,trafficPortion*capacityLowerBound(:,2,m),'r-.');
    plot(N,trafficPortion*capacityLowerBound(:,3,m),'b-.');
    plot(N,trafficPortion*capacityLowerBound(:,4,m),'b--');
        
    plot(1:Nmax,trafficPortion*asymptoticLimitUpper*ones(1,Nmax),'r:');
    
    for c = 1:nbrOfCovarianceModels
        
        plot(1:Nmax,trafficPortion*asymptoticLimitLower(c,m) *ones(1,Nmax),'r:');
        
    end
    
    legend('Upper Bound','Case 1: Uncorrelated','Case 2: Expontial Mod. r=0.7','Case 3: One-Ring, 20 degrees','Case 4: One-Ring, 10 degrees','Location','SouthEast');
    xlabel('Number of Base Station Antennas (N)');
    ylabel('Spectral Efficiency [bit/channel use]');
    
    axis([0 1000 0 4]);
    
end
