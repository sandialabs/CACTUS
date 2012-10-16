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
%     QCx: 0 0 0 0 0 0
%     QCy: 1 2 3 4 5 6
%     QCz: 1 1 1 1 1 1
%     nx: 0 0 0 0 0 0
%     ny: 0 0 0 0 0 0
%     nz: -1 -1 -1 -1 -1 -1
%     tx: 1 1 1 1 1 1
%     ty: 0 0 0 0 0 0
%     tz: 0 0 0 0 0 0
%     CtoR: .1 .1 .1 .1 .1 .1
%     AreaR: .1 .1 .1 .1 .1
%     iSect: 1 1 1 1 1
% Blade 2:
%     ...
% Strut 1:
%     NElem: 5
%     SEx: 0 0 0 0 0 0
%     SEy: 1 2 3 4 5 6
%     SEz: 1 1 1 1 1 1
%     CtoR: .1 .1 .1 .1 .1 .1
%     AreaR: .1 .1 .1 .1 .1
%     TtoC: .15
%     BInd: 1
%     EInd: 3
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
    
    fprintf(fid,'\tQCx: ');
    fprintf(fid,'%13.5e ',T.B(i).QCx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tQCy: ');
    fprintf(fid,'%13.5e ',T.B(i).QCy);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tQCz: ');
    fprintf(fid,'%13.5e ',T.B(i).QCz);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tnx: ');
    fprintf(fid,'%13.5e ',T.B(i).nx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tny: ');
    fprintf(fid,'%13.5e ',T.B(i).ny);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tnz: ');
    fprintf(fid,'%13.5e ',T.B(i).nz);
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
    
    fprintf(fid,'\tAreaR: ');
    fprintf(fid,'%13.5e ',T.B(i).AreaR);
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
    
    fprintf(fid,'\tSEx: ');
    fprintf(fid,'%13.5e ',T.S(i).SEx);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tSEy: ');
    fprintf(fid,'%13.5e ',T.S(i).SEy);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tSEz: ');
    fprintf(fid,'%13.5e ',T.S(i).SEz);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tCtoR: ');
    fprintf(fid,'%13.5e ',T.S(i).CtoR);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tAreaR: ');
    fprintf(fid,'%13.5e ',T.S(i).AreaR);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tTtoC: ');
    fprintf(fid,'%13.5e ',T.S(i).TtoC);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tBInd: ');
    fprintf(fid,'%3d ',T.S(i).BInd);
    fprintf(fid,'\n');
    
    fprintf(fid,'\tEInd: ');
    fprintf(fid,'%3d ',T.S(i).EInd);
    fprintf(fid,'\n');
    
end

fclose(fid);