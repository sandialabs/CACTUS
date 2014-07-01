function WriteTurbineGeom(FN,T)

% Write a CACTUS turbine geometry input file.
%
% Format example:
% NBlade: 3
% NStrut: 3
% RotN: 0 0 1
% RotP: 0 0 0
% RefAR: 2.0
% RefR: 10.0
% Type: VAWT
% Blade 1:
%     NElem: 5
%     FlipN: 0
%     QCx: 0 0 0 0 0 0
%     QCy: 1 2 3 4 5 6
%     QCz: 1 1 1 1 1 1
%     tx: 1 1 1 1 1 1
%     ty: 0 0 0 0 0 0
%     tz: 0 0 0 0 0 0
%     CtoR: .1 .1 .1 .1 .1 .1
%     PEx: 0 0 0 0 0 
%     PEy: 1 2 3 4 5 
%     PEz: 1 1 1 1 1 
%     tEx: 0 0 0 0 0 
%     tEy: 1 2 3 4 5 
%     tEz: 1 1 1 1 1 
%     nEx: 0 0 0 0 0
%     nEy: 1 2 3 4 5
%     nEz: 1 1 1 1 1
%     sEx: 1 1 1 1 1
%     sEy: 0 0 0 0 0
%     sEz: 1 2 3 4 5
%     ECtoR: .1 .1 .1 .1 .1
%     EAreaR: .1 .1 .1 .1 .1
%     iSect: 1 1 1 1 1
% Blade 2:
%     ...
% Strut 1:
%     NElem: 5
%     TtoC: .15
%     MCx: 0 0 0 0 0 0
%     MCy: 1 2 3 4 5 6
%     MCz: 1 1 1 1 1 1
%     CtoR: .1 .1 .1 .1 .1 .1
%     PEx: 0 0 0 0 0
%     PEy: 1 2 3 4 5
%     PEz: 1 1 1 1 1
%     ECtoR: .1 .1 .1 .1 .1
%     EAreaR: .1 .1 .1 .1 .1
%     BIndS: 0
%     EIndS: 0
%     BIndE: 1
%     EIndE: 3
% Strut 2:
%     ...

fid=fopen(FN,'w');

fprintf(fid,'NBlade: ');
fprintf(fid,'%3d ',T.NBlade);
fprintf(fid,'\n');

fprintf(fid,'NStrut: ');
fprintf(fid,'%3d ',T.NStrut);
fprintf(fid,'\n');

fprintf(fid,'RotN: ');
fprintf(fid,'%13.5e ',T.RotN);
fprintf(fid,'\n');

fprintf(fid,'RotP: ');
fprintf(fid,'%13.5e ',T.RotP);
fprintf(fid,'\n');

fprintf(fid,'RefAR: ');
fprintf(fid,'%13.5e ',T.RefAR);
fprintf(fid,'\n');

fprintf(fid,'RefR: ');
fprintf(fid,'%13.5e ',T.RefR);
fprintf(fid,'\n');

fprintf(fid,'Type: ');
if ischar(T.Type)
    fprintf(fid,T.Type);
end
fprintf(fid,'\n');

for i=1:T.NBlade
    fprintf(fid,['Blade ',num2str(i),': ']);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tNElem: ');
    fprintf(fid,'%3d ',T.B(i).NElem);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tFlipN: ');
    fprintf(fid,'%3d ',T.B(i).FlipN);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tQCx: ');
    fprintf(fid,'%13.5e ',T.B(i).QCx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tQCy: ');
    fprintf(fid,'%13.5e ',T.B(i).QCy);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tQCz: ');
    fprintf(fid,'%13.5e ',T.B(i).QCz);
    fprintf(fid,'\n');
    
    fprintf(fid,'\ttx: ');
    fprintf(fid,'%13.5e ',T.B(i).tx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tty: ');
    fprintf(fid,'%13.5e ',T.B(i).ty);
    fprintf(fid,'\n');
    
    fprintf(fid,'\ttz: ');
    fprintf(fid,'%13.5e ',T.B(i).tz);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tCtoR: ');
    fprintf(fid,'%13.5e ',T.B(i).CtoR);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tPEx: ');
    fprintf(fid,'%13.5e ',T.B(i).PEx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tPEy: ');
    fprintf(fid,'%13.5e ',T.B(i).PEy);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tPEz: ');
    fprintf(fid,'%13.5e ',T.B(i).PEz);
    fprintf(fid,'\n');

    fprintf(fid,'\ttEx: ');
    fprintf(fid,'%13.5e ',T.B(i).tEx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\ttEy: ');
    fprintf(fid,'%13.5e ',T.B(i).tEy);
    fprintf(fid,'\n');
    
    fprintf(fid,'\ttEz: ');
    fprintf(fid,'%13.5e ',T.B(i).tEz);
    fprintf(fid,'\n');
  
    fprintf(fid,'\tnEx: ');
    fprintf(fid,'%13.5e ',T.B(i).nEx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tnEy: ');
    fprintf(fid,'%13.5e ',T.B(i).nEy);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tnEz: ');
    fprintf(fid,'%13.5e ',T.B(i).nEz);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tsEx: ');
    fprintf(fid,'%13.5e ',T.B(i).sEx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tsEy: ');
    fprintf(fid,'%13.5e ',T.B(i).sEy);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tsEz: ');
    fprintf(fid,'%13.5e ',T.B(i).sEz);
    fprintf(fid,'\n');

    fprintf(fid,'\tECtoR: ');
    fprintf(fid,'%13.5e ',T.B(i).ECtoR);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tEAreaR: ');
    fprintf(fid,'%13.5e ',T.B(i).EAreaR);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tiSect: ');
    fprintf(fid,'%3d ',T.B(i).iSect);
    fprintf(fid,'\n');
    
end

for i=1:T.NStrut
    fprintf(fid,['Strut ',num2str(i),': ']);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tNElem: ');
    fprintf(fid,'%3d ',T.S(i).NElem);
    fprintf(fid,'\n');

    fprintf(fid,'\tTtoC: ');
    fprintf(fid,'%13.5e ',T.S(i).TtoC);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tMCx: ');
    fprintf(fid,'%13.5e ',T.S(i).MCx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tMCy: ');
    fprintf(fid,'%13.5e ',T.S(i).MCy);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tMCz: ');
    fprintf(fid,'%13.5e ',T.S(i).MCz);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tCtoR: ');
    fprintf(fid,'%13.5e ',T.S(i).CtoR);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tPEx: ');
    fprintf(fid,'%13.5e ',T.S(i).PEx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tPEy: ');
    fprintf(fid,'%13.5e ',T.S(i).PEy);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tPEz: ');
    fprintf(fid,'%13.5e ',T.S(i).PEz);
    fprintf(fid,'\n');

    fprintf(fid,'\tsEx: ');
    fprintf(fid,'%13.5e ',T.S(i).sEx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tsEy: ');
    fprintf(fid,'%13.5e ',T.S(i).sEy);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tsEz: ');
    fprintf(fid,'%13.5e ',T.S(i).sEz);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tECtoR: ');
    fprintf(fid,'%13.5e ',T.S(i).ECtoR);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tEAreaR: ');
    fprintf(fid,'%13.5e ',T.S(i).EAreaR);
    fprintf(fid,'\n');
   
    fprintf(fid,'\tBIndS: ');
    fprintf(fid,'%3d ',T.S(i).BIndS);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tEIndS: ');
    fprintf(fid,'%3d ',T.S(i).EIndS);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tBIndE: ');
    fprintf(fid,'%3d ',T.S(i).BIndE);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tEIndE: ');
    fprintf(fid,'%3d ',T.S(i).EIndE);
    fprintf(fid,'\n');
    
end

fclose(fid);