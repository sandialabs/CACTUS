function TR=RotateTurbine(T,Theta,nR,Origin)

% Rotates turbine structure around normal vector nR (size 1 x 3) through
% angle Theta (rad), using specified origin point (size 1 x 3).

% Rotate turbine rotation axis vector
T.RotN=QuatRot(T.RotN,Theta,nR,[0,0,0]);

% Rotate turbine rotation origin point
T.RotP=QuatRot(T.RotP,Theta,nR,Origin);

% Rotate blades
for i=1:T.NBlade
    T.B(i)=RotateBlade(T.B(i),Theta,nR,Origin);
end

% Rotate struts
for i=1:T.NStrut
    T.S(i)=RotateStrut(T.S(i),Theta,nR,Origin);
end

TR=T;