function SU=UpdateSElemArea(S)

% Updates AreaR component of given strut structure (element areas). Must
% set element locations, and chord to radius ratio components in S before 
% running this function.

for i=1:length(S.AreaR)
    P1=[S.SEx(i),S.SEy(i),S.SEz(i)];
    P2=[S.SEx(i+1),S.SEy(i+1),S.SEz(i+1)];
    dSpR=sqrt(sum((P2-P1).^2));
    S.AreaR(i)=dSpR*(S.CtoR(i)+S.CtoR(i+1))/2;
end

SU=S;