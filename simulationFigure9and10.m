%This Matlab script can be used to generate Figures 9 and 10 in the article:
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
N = [1:10 15:5:40 50:10:100 125:25:500];

Nmax = max(N); %Maximal number of antennas in simulation

%Compute normalized channel covariance matrix R according to the exponential
%correlaton model in Eq. (17).
correlationFactor = 0.7;
R = toeplitz(correlationFactor.^(0:Nmax-1));

%Define the level of hardware impairments at the transmitter and receiver
kappaBS = 0.05^2;
kappaUE = 0.05^2;

%Range of tau parameters (tUE and tBS) in scaling law of Corollary 5
tScalings = [0 1/2];

%Maximal transmit power
powerMax = 2.2222e-8; %in J per channel user
maxSNRdB = 20; %SNR at full output power (in decibel scale)
maxSNR = 10.^(maxSNRdB/10); %SNR at full output power (in linear scale)

%Efficiency of power amplifiers at BS and UE, divided by maximal power
omegaScaled = 0.3 / powerMax;

%Different sums of zeta and rho, divided by maximal power
sumZetaRho = [2e-6 2e-8]; 

%Different splittings between zeta and rho
quotients = [0 0.01 0.1]; %This is rho/(zeta+rho)

%Compute the different combinations of rho and zeta
rho = kron(sumZetaRho,quotients);
zeta = kron(sumZetaRho,1-quotients);

%Equal portion of uplink and downlink transmission
ULDLratio = 0.5;

%Portion of the coherence block used for uplink tranmission (and equally
%for the downlink transmission)
trafficPortion = 0.45;



%%Initialize Monte Carlo simulations

%Number realizations in Monte Carlo simulations
nbrOfMonteCarloRealizations = 10000;

%Generate random realizations
h = sqrt(maxSNR)*sqrtm(R)*(randn(Nmax,nbrOfMonteCarloRealizations)+1i*randn(Nmax,nbrOfMonteCarloRealizations))/sqrt(2); %Generate channel realizations (normalized towards maximal SNR)
etaUE = sqrt(kappaUE)*(randn(1,nbrOfMonteCarloRealizations)+1i*randn(1,nbrOfMonteCarloRealizations))/sqrt(2); %Generate distortion noise at UE
etaBS = sqrt(kappaBS)*( abs(h) .* (randn(Nmax,nbrOfMonteCarloRealizations)+1i*randn(Nmax,nbrOfMonteCarloRealizations)))/sqrt(2); %Generate distortion noise at BS
noise = (randn(Nmax,nbrOfMonteCarloRealizations)+1i*randn(Nmax,nbrOfMonteCarloRealizations))/sqrt(2);


%Placeholders for storing simulation results

%EE and power with hardware impairments, using optimized power allocation
%or the scaling law in Corollary 5
EEImpairmentsOptimized = zeros(length(N),length(rho));
powerImpairmentsOptimized = zeros(length(N),length(rho));

EEImpairmentsScalingLaw = zeros(length(N),length(rho),length(tScalings));
powerImpairmentsScalingLaw = zeros(length(N),length(rho),length(tScalings));

%EE and power with ideal hardware, using optimized power allocation or
%the scaling law in Corollary 5
EEIdealOptimized = zeros(length(N),length(rho));
powerIdealOptimized = zeros(length(N),length(rho));

EEIdealScalingLaw = zeros(length(N),length(rho),length(tScalings));
powerIdealScalingLaw = zeros(length(N),length(rho),length(tScalings));


%Go through different number of antennas
for n = 1:length(N)
    
    %Output the progress of the simulation
    disp(['Antennas: ' num2str(n) '/' num2str(length(N))]);
    
    
    %Go through all combinations of circuit power consumptions
    for r = 1:length(rho)
        
        %Optimize the power allocation (with powerMax as the maximal value,
        %which represents x=1 due to normalization) and store the results
        
        %Hardware impairments
        optimalPowerFraction = fminbnd(@(x) -functionEnergyEfficiency(x,h(1:N(n),:),etaUE,etaBS(1:N(n),:),noise(1:N(n),:),R(1:N(n),1:N(n)),kappaUE,kappaBS,omegaScaled,rho(r),zeta(r),trafficPortion,ULDLratio),0,1);
        EEImpairmentsOptimized(n,r) = functionEnergyEfficiency(optimalPowerFraction,h(1:N(n),:),etaUE,etaBS(1:N(n),:),noise(1:N(n),:),R(1:N(n),1:N(n)),kappaUE,kappaBS,omegaScaled,rho(r),zeta(r),trafficPortion,ULDLratio);
        powerImpairmentsOptimized(n,r) = powerMax*optimalPowerFraction;
        
        %Ideal hardware
        optimalPowerFraction = fminbnd(@(x) -functionEnergyEfficiency(x,h(1:N(n),:),etaUE,etaBS(1:N(n),:),noise(1:N(n),:),R(1:N(n),1:N(n)),0,0,omegaScaled,rho(r),zeta(r),trafficPortion,ULDLratio),0,1);
        EEIdealOptimized(n,r) = functionEnergyEfficiency(optimalPowerFraction,h(1:N(n),:),etaUE,etaBS(1:N(n),:),noise(1:N(n),:),R(1:N(n),1:N(n)),0,0,omegaScaled,rho(r),zeta(r),trafficPortion,ULDLratio);
        powerIdealOptimized(n,r) = powerMax*optimalPowerFraction;
        
        
        %Go through all scaling exponents t of the transmit power
        for m = 1:length(tScalings)
            
            %Use the optimal power allocation for N=1 and then follow the
            %scaling law in Corollary 5 for the given t.

            %Hardware impairments
            powerFractionScaled = powerImpairmentsOptimized(1,r)/powerMax/N(n)^tScalings(m);
            EEImpairmentsScalingLaw(n,r,m) = functionEnergyEfficiency(powerFractionScaled,h(1:N(n),:),etaUE,etaBS(1:N(n),:),noise(1:N(n),:),R(1:N(n),1:N(n)),kappaUE,kappaBS,omegaScaled,rho(r),zeta(r),trafficPortion,ULDLratio);
            powerImpairmentsScalingLaw(n,r,m) = powerMax*powerFractionScaled;
            
            %Ideal hardware
            powerFractionScaledIdeal = powerIdealOptimized(1,r)/powerMax/N(n)^tScalings(m);
            EEIdealScalingLaw(n,r,m) = functionEnergyEfficiency(powerFractionScaledIdeal,h(1:N(n),:),etaUE,etaBS(1:N(n),:),noise(1:N(n),:),R(1:N(n),1:N(n)),0,0,omegaScaled,rho(r),zeta(r),trafficPortion,ULDLratio);
            powerIdealScalingLaw(n,r,m) = powerMax*powerFractionScaledIdeal;
            
        end
        
    end
    
end



%Plot Figure 9 from the paper
figure; hold on; box on;

for r = 1:length(rho)
    
    plot(N,EEIdealOptimized(:,r),'k','LineWidth',1);
    plot(N,EEImpairmentsOptimized(:,r),'r:','LineWidth',1);
    plot(N,EEImpairmentsScalingLaw(:,r,1),'b-.','LineWidth',1);
    plot(N,EEImpairmentsScalingLaw(:,r,2),'k--','LineWidth',1);
    
end

legend('Ideal: Optimized','Non-Ideal: Optimized','Non-Ideal: t=0','Non-Ideal: t=1/2','Location','SouthWest');
xlabel('Number of Base Station Antennas (N)');
ylabel('Energy Efficiency [bit/Joule]');
set(gca,'YScale','log');


%Plot Figure 10 from the paper
figure; hold on; box on;

for r = 1:length(rho)
    
    %These are scaled to get micro Joule
    plot(N,powerIdealOptimized(:,r)*1e6,'k','LineWidth',1);
    plot(N,powerImpairmentsOptimized(:,r)*1e6,'r:','LineWidth',1);
    plot(N,powerImpairmentsScalingLaw(:,r,1)*1e6,'b-.','LineWidth',1);
    plot(N,powerImpairmentsScalingLaw(:,r,2)*1e6,'k--','LineWidth',1);
    
end

legend('Ideal: Optimized','Non-Ideal: Optimized','Non-Ideal: t=0','Non-Ideal: t=1/2','Location','NorthWest');
xlabel('Number of Base Station Antennas (N)');
ylabel('Transmit Power [\mu J/channel use]');
set(gca,'YScale','log');
