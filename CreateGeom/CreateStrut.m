function S=CreateStrut(NElem)

% Create empty strut with specified number of elements, NElem. Strut element
% locations defined at zero turbine rotation phase.

% NElem: Number of strut elements.
% TtoC: Strut thickness to chord ratio.
% MC: Mid chord location (over ref. radius) for each element end (size NElem + 1).
% CtoR: Chord to ref. radius for each element end (size NElem + 1).
% PE: Location (over ref. radius) for each element center (size NElem).
% sE: Element spanwise vector (size NElem).
% ECtoR: Chord to ref. radius for each element (size NElem).
% EAreaR: Area over ref. radius squared for each element (size NElem).
% BIndS: Index of the blade to which the first strut element is attached.
% EIndS: Index of the element on blade BIndS where the first strut element is attached.
% BIndE: Index of the blade to which the last strut element is attached.
% EIndE: Index of the element on blade BInd where the last strut element is attached.
%   Blade and element indicies used in CACTUS to identify the relevant
%   blade data to use in calculating strut-blade interference drag.
%   For struts that are attached to the rotor shaft at one end (not to
%   another blade), set the appropriate BInd and EInd values to zero.

S.NElem=NElem;
S.TtoC=0;
% Element end geometry
S.MCx=zeros(1,NElem+1);
S.MCy=zeros(1,NElem+1);
S.MCz=zeros(1,NElem+1);
S.CtoR=zeros(1,NElem+1);
% Element geometry
S.PEx=zeros(1,NElem);
S.PEy=zeros(1,NElem);
S.PEz=zeros(1,NElem);
S.sEx=zeros(1,NElem);
S.sEy=zeros(1,NElem);
S.sEz=zeros(1,NElem);
S.ECtoR=zeros(1,NElem);
S.EAreaR=zeros(1,NElem);
% Blade interference parameters
S.BIndS=0;
S.EIndS=0;
S.BIndE=1;
S.EIndE=1;