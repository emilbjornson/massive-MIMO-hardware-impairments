%This Matlab script can be used to generate Figure 4, in the article:
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

N = 50; %Number of BS antennas

%Compute normalized channel covariance matrix R according to the exponential
%correlaton model in Eq. (17).
correlationFactor = 0.7;
R = toeplitz(correlationFactor.^(0:N-1));

%Maximal pilot length
Bmax = 10;

%Define the level of hardware impairments at the UE and BS
kappatUE = 0.05^2;
kapparBS = 0.05^2;


%Range of SNRs in simulation (we have normalized sigma2 to 1)
SNRdB = [5 30]; %In decibel scale
SNR = 10.^(SNRdB/10); %In linear scale


%%Initialize Monte Carlo simulations

%Number realizations in Monte Carlo simulations
nbrOfMonteCarloRealizations = 100000;

%Generate random realizations
h = sqrtm(R)*(randn(N,nbrOfMonteCarloRealizations)+1i*randn(N,nbrOfMonteCarloRealizations))/sqrt(2); %Generate channel realizations
etatBS = (randn(1,nbrOfMonteCarloRealizations,Bmax)+1i*randn(1,nbrOfMonteCarloRealizations,Bmax))/sqrt(2); %Generate distortion noise at transmitter (to be scaled by kappatUE)
etarUE = ( repmat(abs(h),[1 1 Bmax]) .* (randn(N,nbrOfMonteCarloRealizations,Bmax)+1i*randn(N,nbrOfMonteCarloRealizations,Bmax)))/sqrt(2); %Generate distortion noise at receiver (to be scaled by kapparBS)
nu = (randn(N,nbrOfMonteCarloRealizations,Bmax)+1i*randn(N,nbrOfMonteCarloRealizations,Bmax))/sqrt(2); %Generate receiver realizations


%Placeholders for storing simulation results
normalizedMSE_corr_distortion = zeros(length(SNR),Bmax); %Normalized MSE for estimator in Eq. (15) for fully correlated distortion noise
normalizedMSE_uncorr_distortion = zeros(length(SNR),Bmax); %Normalized MSE for estimator in Eq. (15) for uncorrelated distortion noise
normalizedMSE_ideal = zeros(length(SNR),Bmax); %Normalized MSE for LMMSE estimator with ideal hardware


%Go through SNR values
for m = 1:length(SNR)
    
    %Output the progress of the simulation
    disp(['SNR: ' num2str(m) '/' num2str(length(SNR))]);
    
    %Compute the pilot signal for the given SNR value
    d = sqrt(SNR(m));
    
    %Go through different pilot lengths
    for B = 1:Bmax
        
        %Compute matrix A in the LMMSE estimator (see Eq. (9))
        A_LMMSE = conj(d) * R(1:N,1:N) / (abs(d)^2*(1+kappatUE)*R(1:N,1:N) + abs(d)^2*kapparBS*diag(diag(R(1:N,1:N)))+eye(N));
        
        %Compute matrix A in the LMMSE estimator (see Eq. (9)) for ideal hardware
        A_ideal = conj(d) * R(1:N,1:N) / (abs(d)^2*R(1:N,1:N) +eye(N));
        
        
        %Placeholders for storing squared estimation errors at current SNR
        errors_corr_distortion = zeros(nbrOfMonteCarloRealizations,1);
        errors_uncorr_distortion = zeros(nbrOfMonteCarloRealizations,1);
        errors_ideal = zeros(nbrOfMonteCarloRealizations,1);
        
        %Go through all Monte Carlo realizations
        for k = 1:nbrOfMonteCarloRealizations

            %Compute received signals 
            z_corr = h(1:N,k) * ( d + abs(d)*sqrt(kappatUE)*etatBS(1,k,1) ) + abs(d)*sqrt(kapparBS)*etarUE(1:N,k,1) + sum(nu(1:N,k,1:B),3)/B;
            z_uncorr = h(1:N,k) * ( d + abs(d)*sqrt(kappatUE)*sum(etatBS(1,k,1:B),3)/B ) + abs(d)*sqrt(kapparBS)*sum(etarUE(1:N,k,1:B),3)/B + sum(nu(1:N,k,1:B),3)/B;
            z_ideal = h(1:N,k) * d  + sum(nu(1:N,k,1:B),3)/B;
            
            %Compute channel estimates
            hhat_corr = A_LMMSE*z_corr;
            hhat_uncorr = A_LMMSE*z_uncorr;
            hhat_ideal = A_ideal*z_ideal;
            
            %Compute the squared norms of the channel estimation errors
            errors_corr_distortion(k) = norm(hhat_corr -  h(1:N,k)).^2/N;
            errors_uncorr_distortion(k) = norm(hhat_uncorr -  h(1:N,k)).^2/N;
            errors_ideal(k) = norm(hhat_ideal -  h(1:N,k)).^2/N;
            
        end
        
        %Compute normalized MSEs as the average of the squared norms of the 
        %estimation errors over the Monte Carlo realizations
        normalizedMSE_corr_distortion(m,B) = mean(errors_corr_distortion);
        normalizedMSE_uncorr_distortion(m,B) = mean(errors_uncorr_distortion);
        normalizedMSE_ideal(m,B) = mean(errors_ideal);
        
    end
    
end



%Plot Figure 4 from the paper
figure; hold on; box on;

for m = 1:length(SNR)
    
    plot(1:Bmax,normalizedMSE_corr_distortion(m,:),'ro-','LineWidth',1);
    plot(1:Bmax,normalizedMSE_uncorr_distortion(m,:),'bd-.','LineWidth',1);
    plot(1:Bmax,normalizedMSE_ideal(m,:),'k*--','LineWidth',1);
    
end

set(gca,'YScale','log');
xlabel('Pilot Length (B)');
ylabel('Relative Estimation Error per Antenna');
legend('Fully-Correlated Distortion Noise','Uncorrelated Distortion Noise','Ideal Hardware','Location','NorthEast');
