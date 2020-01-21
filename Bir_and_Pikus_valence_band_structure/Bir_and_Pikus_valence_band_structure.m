function [E1,E2]=Bir_and_Pikus_valence_band_structure(x,y,z,teps,eps)
%=========================================================================%
% <vladimir.kolchuzhin@ieee.org>
% 2017-07-01 Munich 
% 2020-01-21 realise
%-------------------------------------------------------------------------%
%		FILE: Bir_and_Pikus_valence_band_structure.m
%
%       PURPOSE: energy calculation for p-type silicon valence bands
%		
%		INPUT:
%           x,y,z - wave vector k, 1/nm
%            eps - starin value
%           teps - starin tensor
%
%		OUTPUT:
%           E1 band1 energy, eV
%           E2 band2 energy, eV
%
%       REFERENCE: 
%           G. L. Bir and G. E. Pikus,
%           "Symmetry and Strain-Induced Effects in Semiconductors",
%           Wiley, New York,1974
%=========================================================================%
%=========================================================================%
Sp_eps=teps(1,1)+teps(2,2)+teps(3,3);

x2=x.*x;
y2=y.*y;
z2=z.*z;
r2=x2+y2+z2;

% scale factor: Ei => eV; k => 1/nm
k=0.038;    

% energy band parameters for p-Si
A=4.1*k;
B=1.6*k;
C=3.3*k;
D=4.3*k;

% deformation potentials for Si
Eb=2.30;
Ed=5.10;

% energy
Ee1=(Eb*Eb/2.)*((teps(1,1)-teps(2,2)).^2+(teps(2,2)-teps(3,3)).^2+(teps(3,3)-teps(1,1)).^2);
Ee2=Ed*Ed*(teps(1,2).^2+teps(1,3).^2+teps(2,3).^2);

Ee=eps*eps*(Ee1+Ee2);
Eek=eps.*(B*Eb.*(3.*(x2.*teps(1,1)+y2.*teps(2,2)+z2.*teps(3,3))-r2.*Sp_eps)+2.0*D*Ed.*(x.*y.*teps(1,2)+x.*z.*teps(1,3)+y.*z.*teps(2,3)));
Ek=B.*B.*r2.*r2+C.*C.*(x2.*y2+x2.*z2+y2.*z2);

root=sqrt(abs(Ee+Eek+Ek));

E1=A*r2+root;
E2=A*r2-root;
%=========================================================================%
return