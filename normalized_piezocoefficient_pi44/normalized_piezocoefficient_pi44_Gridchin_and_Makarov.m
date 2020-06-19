function [P] = normalized_piezocoefficient_pi44_Gridchin_and_Makarov(T,Na)
%=========================================================================%
% PURPOSE:
%           the normalized piezocoefficient pi44 for Si p-type
%           as a function of temperature T and acceptor density Na
%
% INPUT:    T  in K
%           Na in cm^-3
%           
% OUTPUT:   P == normalized piezocoefficient pi44
%           
% REFERENCE: 
%           Gridchin V.A., Makarov E.A. "Raschet temperaturnoj i koncentracionnoj
%           zavisimosti p'ezosoprotivlenija diffuzionnyh sloev kremnija" 
%           Izv. Severo-Kavkazkogo nauchn. centra Vysshej shkoly. 1976, Nr3. 
%           (in Russian)
%
%           Na[1e15-1e20] & T[100-400]: p44 ~ 1/T
%-------------------------------------------------------------------------%
% <vladimir.kolchuzhin@ieee.org>
% 19.06.2020 Hohenschaeftlarn
%=========================================================================%
if nargin==0 % selftest

    T=300;
    Na=1e17;
    
    %[P] = normalized_piezocoefficient_pi44_Gridchin_and_Makarov(T,Na); % normalization
    
    Na_set=[...
        [1.0; 2.0; 4.0; 6.0; 8.0;].*1e17;
        [1.0; 2.0; 4.0; 6.0; 8.0;].*1e18;
        [1.0; 2.0; 4.0; 6.0; 8.0;].*1e19;
        1.0e20;];
    
    T_set=[200; 225; 250; 275; 300; 325; 350; 375; 400; 425; 450;];
    
    %ind_color=['r';'g';'b';'c';'m';'y';'k';];
    %ind_marker=['o';'^';'d';'*';'<';'>';'k';];
    
%-------------------------------------------------------------------------%
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
            %Pset(i)=piezocoefficient_pi44(T_set(1),Na_set(i));
            Pset(i)=normalized_piezocoefficient_pi44_Gridchin_and_Makarov(T_set(1),Na_set(i));
        end
        semilogx(Na_set,Pset,'LineWidth',2,'Color',color_set(1,:));
    
    
    for j=1:numel(T_set)
        for i=1:numel(Na_set)
           %Pset(i)=piezocoefficient_pi44(T_set(j),Na_set(i));
           Pset(i)=normalized_piezocoefficient_pi44_Gridchin_and_Makarov(T_set(j),Na_set(i));
        end
        hold on;
        %plot(Na_set,Pset);
        semilogx(Na_set,Pset,'LineWidth',2,'Color',color_set(j,:));
    end
    
    box on; grid on;
    xlabel('N_A, cm^-^3'); ylabel('P(N_A,T)'); 
    legend('200 K','225 K','250 K','275 K','300 K','325 K','350 K','375 K','400 K','425 K','450 K');
    %axis([1e17 1e20 0 1.5]);
    
end
%=========================================================================%
To=300; % K
V=To./T;
%-------------------------------------------------------------------------%
k=1e-11;    % 1/Pa
norm=[116.480267095506];

% fitting parameters
D=0.0; C=0.0;
if (Na < 1e15)              D= 91+47.7*V; C= 0.0;       end
if (Na > 1e15)&&(Na < 3e18) D=398-85.8*V; C= 8.9-3.9*V; end
if (Na > 3e18)              D=670*V-202;  C=13.9*V-5.2; end
%-------------------------------------------------------------------------%
P = D - C.*log(Na);
P=P./norm;
%=========================================================================%
return