function [P] = piezocoefficient_pi44(T,Na)
%=========================================================================%
% PURPOSE:
%           the normalized piezocoefficient pi44 
%           as a function of temperature T and acceptor density Na
%
% INPUT:    T  in K
%           Na in cm^-3
%           
% OUTPUT:   P == normalized piezocoefficient pi44
%           
% REFERENCE: 
%           Piezoresistance in p-type silicon revisited
%           J. Richter, J. Pedersen, M. Brandbyge, E. V. Thomsen, O. Hansen
%           Journal of Applied Physics 104, 023715 (2008); 
%           doi: 10.1063/1.2960335
%-------------------------------------------------------------------------%
% Kolchuzhin V.A., MMT, TU Chemnitz
% <vladimir.kolchuzhin@etit.tu-chemnitz.de>
% 28.11.2014 test passed in R2009a
%=========================================================================%
if nargin==0 % selftest

    T=300; Na=1e20;
    
    Na_set=[...
        [1.0; 2.0; 4.0; 6.0; 8.0;].*1e17;
        [1.0; 2.0; 4.0; 6.0; 8.0;].*1e18;
        [1.0; 2.0; 4.0; 6.0; 8.0;].*1e19;
        1.0e20;];
    
    T_set=[200; 225; 250; 275; 300; 325; 350; 375; 400; 425; 450;];
    
    color_set=[...
    0.0 0.0 1.0
    0.1 0.2 0.9
    0.2 0.4 0.8
    0.3 0.6 0.7
    0.4 0.8 0.6
    0.5 1.0 0.5
    0.6 0.8 0.4
    0.7 0.6 0.3
    0.8 0.4 0.2
    0.9 0.2 0.1
    1.0 0.0 0.0
    ];
%-------------------------------------------------------------------------%
        for i=1:numel(Na_set)
            Pset(i)=piezocoefficient_pi44(T_set(1),Na_set(i));    
        end
        semilogx(Na_set,Pset,'LineWidth',2,'Color',color_set(1,:));
    
    
    for j=1:numel(T_set)
        for i=1:numel(Na_set)
           Pset(i)=piezocoefficient_pi44(T_set(j),Na_set(i));    
        end
        hold on;
        semilogx(Na_set,Pset,'LineWidth',2,'Color',color_set(j,:));
    end
    
    box on; grid on;
    xlabel('N_A, cm^-^3'); ylabel('P(N_A,T)'); 
    legend('200 K','225 K','250 K','275 K','300 K','325 K','350 K','375 K','400 K','425 K','450 K');
    %axis([1e17 1e20 0 1.5]);
    
end
%=========================================================================%
To=300; % K
theta=T/To;
%-------------------------------------------------------------------------%
% fitting parameters
Nb=6e19;    % cm^-3
Nc=7e20;    % cm^-3
% power coefficients
alpfa=0.43;
beta=0.1;
gamma=1.6;
eta=3;
upsilon=0.9;
%-------------------------------------------------------------------------%
P = theta^-upsilon / ( 1 + ((Na/Nb)^alpfa)*theta^-beta  + ((Na/Nc)^gamma)*theta^-eta);
%=========================================================================%
return
