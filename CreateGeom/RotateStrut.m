function SR=RotateStrut(S,Theta,nR,Origin)

% Rotates strut structure around normal vector nR (size 1 x 3) through
% angle Theta (rad), using specified origin point (size 1 x 3).

% Rotate element locations
P=[S.SEx;S.SEy;S.SEz];
PR=QuatRot(P',Theta,nR,Origin);
S.SEx=PR(:,1)';
S.SEy=PR(:,2)';
S.SEz=PR(:,3)';

SR=S;