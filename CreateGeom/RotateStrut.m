function SR=RotateStrut(S,Theta,nR,Origin)

% Rotates strut structure around normal vector nR (size 1 x 3) through
% angle Theta (rad), using specified origin point (size 1 x 3).

% Rotate element end geometry and recalculate element geom for consistency.

% Rotate element locations
P=[S.MCx;S.MCy;S.MCz];
PR=QuatRot(P',Theta,nR,Origin);
S.MCx=PR(:,1)';
S.MCy=PR(:,2)';
S.MCz=PR(:,3)';

% Calc element geometry
S=CalcSEGeom(S);

SR=S;