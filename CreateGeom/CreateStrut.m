function S=CreateStrut(NElem)

% Create empty strut with specified number of elements, NElem. Strut element
% locations defined at zero turbine rotation phase.

% SE: Location (over ref. radius) for each element end (size NElem + 1).
% CtoR: Chord to ref. radius for each element end (size NElem + 1).
% AreaR: Element area over ref. radius squared for each element (size NElem).
% TtoC: Strut thickness to chord ratio.
% BInd: Index of the blade to which the strut is attached.
% EInd: Index of the element on blade BInd where the strut is attached.
%   Blade and element indicies used in CACTUS to identify the relevant
%   blade data to use in calculating strut-blade interference drag.

S.NElem=NElem;
S.SEx=zeros(1,NElem+1);
S.SEy=zeros(1,NElem+1);
S.SEz=zeros(1,NElem+1);
S.CtoR=zeros(1,NElem+1);
S.AreaR=zeros(1,NElem);
S.TtoC=0;
S.BInd=1;
S.EInd=1;