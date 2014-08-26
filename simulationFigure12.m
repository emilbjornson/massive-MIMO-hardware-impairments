%This Matlab script can be used to generate Figure 12 in the article:
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

N = 200; %Number of BS antennas

%Define the normalized channel covariance matrix R as an identity matrix
R = eye(N);

%Number of co-users that interfere during data transmission
nbrOfInterferingUsers = 2;

%Number of co-users that interfere during pilot transmission
nbrOfContaminatingUsers = 1;

%Define the range of level of hardware impairments at the BS and
%UE. The values are the same for BS and UE, but can be different.
kappaValuesBS = [0 0.05^2 0.1^2];
kappaValuesUE = [0 0.05^2 0.1^2];

%Portion of all symbols used for either DL or UL data transmission.
trafficPortion = 0.45;

%Range of SNRs in simulation (we have normalized sigma2 to 1)
SNRdB = 20; %In decibel scale
SNR = 10.^(SNRdB/10); %In linear scale

%Range of how much weaker the interfering signals are as compared to the
%useful signal (on average).
weakerIndB = -50:2:0; %Range of weaker SNRs for pilot contamination
weakerInterferenceCont = SNR*10.^([-Inf weakerIndB]/10); %Add minus infinity to model the interference-free case
weakerInterferenceRegular = SNR*10.^([-Inf -10*ones(size(weakerIndB))]/10); %Regular interferer is -10 dB weaker


%%Initialize Monte Carlo simulations

%Number realizations in Monte Carlo simulations
nbrOfMonteCarloRealizations = 10000;

%Generate random realizations
h = sqrtm(R)*(randn(N,nbrOfMonteCarloRealizations)+1i*randn(N,nbrOfMonteCarloRealizations))/sqrt(2); %Generate channel realizations for useful channel
h_cousers = (randn(N,nbrOfInterferingUsers,nbrOfMonteCarloRealizations)+1i*randn(N,nbrOfInterferingUsers,nbrOfMonteCarloRealizations))/sqrt(2); %Generate channel realizations for interfering users
etaUE = (randn(1,nbrOfMonteCarloRealizations)+1i*randn(1,nbrOfMonteCarloRealizations))/sqrt(2);  %Generate distortion noise at UE (to be scaled by kappa and SNR)
etaBS = ( h .* (randn(N,nbrOfMonteCarloRealizations)+1i*randn(N,nbrOfMonteCarloRealizations)))/sqrt(2);  %Generate distortion noise at BS (to be scaled by kappa and SNR)
noise = (randn(N,nbrOfMonteCarloRealizations)+1i*randn(N,nbrOfMonteCarloRealizations))/sqrt(2);


%Placeholders for storing simulation results
rate_onlycont = zeros(length(kappaValuesBS),length(weakerInterferenceCont));
rate_interf = zeros(length(kappaValuesBS),length(weakerInterferenceCont));


for m = 1:length(weakerInterferenceCont)
    
    %Output the progress of the simulation
    disp(['Interference: ' num2str(m) '/' num2str(length(weakerInterferenceCont))]);
    
    %Go through all level of impairments values
    for impairIndex = 1:length(kappaValuesBS)
        
        %Extract current level of impairments at BS and UE
        kappaUE = kappaValuesUE(impairIndex);
        kappaBS = kappaValuesBS(impairIndex);
        
        %A typical pilot signal
        d = sqrt(SNR);
        
        %Compute matrix A in the LMMSE estimator (see Eq. (9)) by extending
        %them to also take interference into account
        A_LMMSE_interf = conj(d) * R / (abs(d)^2*(1+kappaUE)*R + abs(d)^2*kappaBS*diag(diag(R)) + (1+ nbrOfContaminatingUsers*weakerInterferenceCont(m)) * eye(N) );
        
        %Placeholders for storing Monte Carlo simulation results
        firstMoment_interf = zeros(nbrOfMonteCarloRealizations,1);
        distortionTerm_interf = zeros(nbrOfMonteCarloRealizations,1);
        interference_all = zeros(nbrOfMonteCarloRealizations,1);
        interference_onlycont = zeros(nbrOfMonteCarloRealizations,1);
        
        for k = 1:nbrOfMonteCarloRealizations
            
            %Compute received signal
            z_interf = h(:,k) * ( d + d*sqrt(kappaUE)*etaUE(k) ) + d*sqrt(kappaBS)*etaBS(:,k) + noise(:,k) + sqrt(weakerInterferenceCont(m)) * sum(h_cousers(:,1:nbrOfContaminatingUsers,k),2);
            
            %Compute channel estimates
            hhat_interf = A_LMMSE_interf*z_interf;
            
            %Compute the beamforming vector (MRT/MRC)
            beamforming_interf = sqrt(SNR)*hhat_interf/norm(hhat_interf);
            
            %Compute a realization of the first and second moments of
            %the inner product between beamforming and channel
            firstMoment_interf(k) = h(:,k)'*beamforming_interf;
            
            %The elementwise product between channel and beamforming
            %vectors (and sum over these elements) appear in the
            %distortion noise term. This computes a realization.
            distortionTerm_interf(k) = sum( (abs(h(:,k)).^2) .* (abs(beamforming_interf).^2));
            
            %Compute the interference power during data transmission, when
            %there is only pilot contamination or also regular interference
            interference_onlycont(k) = (weakerInterferenceCont(m)/SNR) * norm( h_cousers(:,1:nbrOfContaminatingUsers,k)'*beamforming_interf)^2;
            interference_all(k) = interference_onlycont(k) + (weakerInterferenceRegular(m)/SNR) * norm( h_cousers(:,nbrOfContaminatingUsers+1:end,k)'*beamforming_interf)^2;
            
        end
        
        %Compute averages and variances by the Monte Carlo simulation
        meanvalue_interf = mean(firstMoment_interf);
        variance_interf = var(firstMoment_interf);
        meanDistortion_interf = mean(distortionTerm_interf);
        
        %Finalize the Monte Carlo simulations by computing the achievable
        %uplink spectral efficiency
        rate_onlycont(impairIndex,m) = log2(1+ abs(meanvalue_interf).^2 ./ ( (1+kappaUE) * variance_interf + kappaUE*abs(meanvalue_interf).^2 + kappaBS*meanDistortion_interf + mean(interference_onlycont) + 1/SNR) );
        rate_interf(impairIndex,m) = log2(1+ abs(meanvalue_interf).^2 ./ ( (1+kappaUE) * variance_interf + kappaUE*abs(meanvalue_interf).^2 + kappaBS*meanDistortion_interf + mean(interference_all) + 1/SNR) );
        
    end
    
end


%Plot Figure 12 from the paper

figure; hold on; box on;

for impairIndex = 1:length(kappaValuesBS)
    
    plot(weakerIndB,trafficPortion*rate_onlycont(impairIndex,1)*ones(size(weakerIndB)),'r--');
    plot(weakerIndB,trafficPortion*rate_onlycont(impairIndex,2:end),'b-.');
    plot(weakerIndB,trafficPortion*rate_interf(impairIndex,2:end),'k-');
    
end

xlabel('Relative Channel Gain of Pilot Contamination [dB]');
ylabel('Uplink Spectral Efficiency [bit/channel use]');
axis([min(weakerIndB) max(weakerIndB) 0 4.5]);
legend('No Interference','Only Pilot Cont.','Pilot Con. & Regular Interf. (Relative Gain -10 dB)','Location','SouthWest');
