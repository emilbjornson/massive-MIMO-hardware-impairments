%This Matlab script can be used to generate Figure 11 in the article:
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
N = [1:50 60:10:150 175:25:500];

Nmax = max(N); %Maximal number of antennas in simulation

%Compute normalized channel covariance matrix R according to the exponential
%correlaton model in Eq. (17).
correlationFactor = 0.7;
R = toeplitz(correlationFactor.^(0:Nmax-1));

%Define the original level of hardware impairments at the transmitter and
%receiver, before the hardware scaling law is applied. The values are the
%same for BS and UE, but can be different.

%Level of hardware impairments at the BS, before the hardware scaling law
%has been applied.
kappaBSOriginal = 0.05^2;

%Fixed level of of hardware impairments at the UE
kappaUE = 0.05^2;

%Portion of all symbols used for either DL or UL data transmission.
trafficPortion = 0.45;

%Range of tau-values in the hardware scaling law of Corollary 7
tauValues = [0 1/4 1/2 1 2];

%Range of SNRs in simulation (we have normalized sigma2 to 1)
SNRdB = 20; %In decibel scale
SNR = 10.^(SNRdB/10); %In linear scale


%%Initialize Monte Carlo simulations

%Number realizations in Monte Carlo simulations
nbrOfMonteCarloRealizations = 10000;

%Generate random realizations
h = sqrtm(R)*(randn(Nmax,nbrOfMonteCarloRealizations)+1i*randn(Nmax,nbrOfMonteCarloRealizations))/sqrt(2); %Generate channel realizations
etaUE = (randn(1,nbrOfMonteCarloRealizations)+1i*randn(1,nbrOfMonteCarloRealizations))/sqrt(2);  %Generate distortion noise at UE
etaBS = ( abs(h) .* (randn(Nmax,nbrOfMonteCarloRealizations)+1i*randn(Nmax,nbrOfMonteCarloRealizations)))/sqrt(2); %Generate distortion noise at BS
noise = (randn(Nmax,nbrOfMonteCarloRealizations)+1i*randn(Nmax,nbrOfMonteCarloRealizations))/sqrt(2);


%Placeholders for storing simulation results
capacityLowerBound = zeros(length(N),length(tauValues));


%Go through different number of antennas
for n = 1:length(N)
    
    %Output the progress of the simulation
    disp(['Antennas: ' num2str(n) '/' num2str(length(N))]);
    
    %Go through different tau-values in the hardware scaling law
    for m = 1:length(tauValues)
        
        %Compute the level of hardware impairments at the BS, according to
        %the scaling law in Corollary 7
        kappaBS=kappaBSOriginal*(N(n))^(tauValues(m));
        
        %A typical pilot signal
        d = sqrt(SNR);
        
        %Compute matrix A in the LMMSE estimator (see Eq. (9))
        A_LMMSE = conj(d) * R(1:N(n),1:N(n)) / (abs(d)^2*(1+kappaUE)*R(1:N(n),1:N(n)) + abs(d)^2*kappaBS*diag(diag(R(1:N(n),1:N(n))))+eye(N(n)));
        
        %Placeholders for storing Monte Carlo simulation results
        firstMoment = zeros(nbrOfMonteCarloRealizations,1);
        distortionTerm = zeros(nbrOfMonteCarloRealizations,1);
        
        
        %Go through all Monte Carlo realizations
        for k = 1:nbrOfMonteCarloRealizations
            
            %Compute received signal
            z = h(1:N(n),k) * ( d + d*sqrt(kappaUE)*etaUE(k) ) + sqrt(kappaBS)*d*etaBS(1:N(n),k) + noise(1:N(n),k);
            
            %Compute channel estimates
            hhat = A_LMMSE*z;
            
            %Compute the beamforming vector (MRT/MRC)
            beamforming = sqrt(SNR)*hhat/norm(hhat);
            
            %Compute a realization of the first moment of the inner product
            %between beamforming and channel
            firstMoment(k) = h(1:N(n),k)'*beamforming;

            %The elementwise product between channel and beamforming
            %vectors (and sum over these elements) appear in the
            %distortion noise term. This computes a realization.
            distortionTerm(k) = sum( (abs(h(1:N(n),k)).^2) .* (abs(beamforming).^2));
            
        end
        
        %Finalize the Monte Carlo simulations by computing the lower bound
        %on the capacity
        capacityLowerBound(n,m) = log2(1+ abs(mean(firstMoment,1)).^2 ./ ( (1+kappaUE) * var(firstMoment) + kappaUE* abs(mean(firstMoment,1)).^2 + kappaBS*mean(distortionTerm,1) + 1) );
        
    end
end



%Plot Figure 11 from the paper

figure; hold on; box on;

plot(N,trafficPortion*capacityLowerBound(:,1),'b-.');

for m=2:length(tauValues)
    
    plot(N,trafficPortion*capacityLowerBound(:,m),'k');
    
end

legend('Fixed Hardware (\tau=0)','Hardware Degradation','Location','SouthEast')
xlabel('Number of Base Station Antennas (N)')
ylabel('Spectral Efficiency [bit/channel use]')

axis([0 500 0 3.5]);
