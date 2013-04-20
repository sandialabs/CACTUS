function S=CreateStrut(NElem)

% Create empty strut with specified number of elements, NElem. Strut element
% locations defined at zero turbine rotation phase.

% SE: Location (over ref. radius) for each element end (size NElem + 1).
% CtoR: Chord to ref. radius for each element end (size NElem + 1).
% AreaR: Element area over ref. radius squared for each element (size NElem).
% TtoC: Strut thickness to chord ratio.
% BIndS: Index of the blade to which the first strut element is attached.
% EIndS: Index of the element on blade BIndS where the first strut element is attached.
% BIndE: Index of the blade to which the last strut element is attached.
% EIndE: Index of the element on blade BInd where the last strut element is attached.
%   Blade and element indicies used in CACTUS to identify the relevant
%   blade data to use in calculating strut-blade interference drag.
%   For struts that are attached to the rotor shaft at one end (not to
%   another blade), set the appropriate BInd and EInd values to zero.

S.NElem=NElem;
S.SEx=zeros(1,NElem+1);
S.SEy=zeros(1,NElem+1);
S.SEz=zeros(1,NElem+1);
S.CtoR=zeros(1,NElem+1);
S.AreaR=zeros(1,NElem);
S.TtoC=0;
S.BIndS=0;
S.EIndS=0;
S.BIndE=1;
S.EIndE=1;