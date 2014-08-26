function [EE,rate,power] = functionEnergyEfficiency(power,h,etaUE,etaBS,noise,R,kappaUE,kappaBS,omega,rho,zeta,trafficPortion,ULDLratio)
%This is an implementation of energy efficiency metric used in the article:
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
%INPUT:
%power          = Transmit power
%h              = Channel realizations (dimension: N x nbrOfMonteCarloRealizations)
%etaUE          = Distortion noise realization at UE 
%                 (dimension 1 x nbrOfMonteCarloRealizations)
%etaBS          = Distortion noise realization at UE
%                 (dimension N x nbrOfMonteCarloRealizations)
%noise          = Matrix with noise realizations in UL 
%                (dimension: N x nbrOfMonteCarloRealizations)
%R              = N x N channel covariance matrix
%kappaUE        = Level of hardware impairments at UE
%kappaBS        = Level of hardware impairments at BS
%omega          = Efficiency of power amplifiers at UE and BS
%rho            = Circuit power that scales with N
%zeta           = Circuit power independent of N
%trafficPortion = Portion of the total resources available for data
%                 transmission in the direction under study (smaller than
%                 ULDLratio)
%ULDLratio      = Fraction of total resources allocated to the transmission
%                 direction under study (set it to 0.5 to study both UL/DL)
%
%OUTPUT:
%EE             = Energy efficiency according to model in Definition 1


%A typical pilot signal
d = sqrt(power);

%Extract number of antennas
N = size(h,1);

%Extract number of Monte Carlo simulations
nbrOfMonteCarloRealizations = size(h,2);
            
%Compute matrix A in the LMMSE estimator (see Eq. (9))
A_LMMSE = conj(d) * R / (abs(d)^2*(1+kappaUE)*R + abs(d)^2*kappaBS*diag(diag(R))+eye(N));


%Placeholders for storing Monte Carlo simulation results
firstMoment = zeros(nbrOfMonteCarloRealizations,1);
distortionTerm = zeros(nbrOfMonteCarloRealizations,1);

%Go through all Monte Carlo realizations
for k = 1:nbrOfMonteCarloRealizations
    
    %Compute received signal
    z = h(:,k) * ( d + d*etaUE(k) ) + d*etaBS(:,k) + noise(:,k);
    
    %Compute channel estimates
    hhat = A_LMMSE*z;
    
    %Compute the beamforming vector (MRT/MRC)
    beamforming = sqrt(power)*hhat/norm(hhat);

    %Compute a realization of the first moment of the inner product between
    %beamforming and channel 
    firstMoment(k) = h(:,k)'*beamforming;
    
    %The elementwise product between channel and beamforming vectors (and
    %sum over these elements) appear in the distortion noise term. This
    %computes a realization.
    distortionTerm(k) = sum( (abs(h(:,k)).^2) .* (abs(beamforming).^2));
    
end

%Finalize the Monte Carlo simulations by computing lower bounds on the capacity
rate = trafficPortion * log2(1+ abs(mean(firstMoment,1)).^2 ./ ( (1+kappaUE) * var(firstMoment) + kappaUE* abs(mean(firstMoment,1)).^2 + kappaBS*mean(distortionTerm,1) + 1) );

%Compute the energy efficiency according to Definition 1
EE = rate/((power/omega + N*rho + zeta)*ULDLratio);
