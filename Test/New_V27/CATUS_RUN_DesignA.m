close all
clear all
tic
if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 11
end
workingDir = '/Users/clkell/Documents/Modeling_Tools/CACTUS_gfortran/Test/New_V27';
cd(workingDir);
global beta_off  ix R CRr rB bTwist HubRR Tilt eta bCone bi NBlade af NElem;

%%%%%%%%%% Define wind turbine blade geometry
Im_vortex = importdata('Design_A_dec.dat','\t',1);
% Params
R = 44.2913; % radius (ft)
HubRR = Im_vortex.data(1,1); % hub radius ratio
Tilt = 0; % Rotor tilt angle (deg, positive windward axis tilted up)
eta = 0; % Blade mount point ratio ((distance behind leading edge of the blade mount point) / (root chord))
rB = Im_vortex.data(:,1)'; % r/R blade stations
CRr = Im_vortex.data(:,2)'; % c/R blade chord stations
bCone=0; % Blade coning angle (deg, positive tip into the wind (-x))
bi=0; % Blade planform incidence (deg, w.r.t. rotor disk plane, positive LE into the wind (-x))
bTwist = Im_vortex.data(:,3)'; % Blade station twist in degrees relative to rotation direction
af = Im_vortex.data(:,4)'; % airfoil station type
NBlade=3;
NElem = length(CRr)-1;
wake_cutoff = 5;

Rm = R.*0.3048;
lambda = [4:1:13];
nr  = 30.*ones(1,length(lambda));
U = 8.*ones(1,length(lambda));
RPM = U.*lambda./Rm .*60./2./pi;
beta_off = 0.*ones(1,length(lambda));
% ELEMENT geometry
for i = 1:length(rB)-1
    rBe(i) = (rB(i)+rB(i+1))./2;
    cRre(i) = (CRr(i)+CRr(i+1))./2;
    bTwiste(i) = (bTwist(i)+bTwist(i+1))./2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA = cell(1,length(lambda));
DATA_ALL = cell(1,length(lambda));
% DATA_WAKE = cell(1,length(lambda));
DATA_ELEMENT = cell(1,length(lambda));

for ix = 1:length(lambda)
    Create_New_V27_A;
    fid  = fopen(['New_V27' num2str(ix) '.geom']  ,'r');
    f = char(fread(fid,'char*1'))';
    fclose(fid);
    f = strrep(f,[repmat('1   ',1,length(rB)-2) '1'],num2str(af(2:end)));
    geomfile{ix} = char(['New_V27' num2str(ix) '.geom']);
    fid  = fopen(geomfile{ix},'w');
    fprintf(fid,'%s',f);
    fclose(fid);
        
    
fid = fopen('INPUT.in','r');
f = char(fread(fid,'char*1'))';
fclose(fid);
% Replace number revolutions, TSR, xstop with variables from MATLAB
f = strrep(f,char('nr = x'),char(strcat({'nr = '},char(num2str(nr(ix))))));
f = strrep(f,char('RPM = x'),char(strcat({'RPM = '},char(num2str(RPM(ix))))));
f = strrep(f,char('Ut = x'),char(strcat({'Ut = '},char(num2str(lambda(ix))))));
f = strrep(f,char('xstop = x'),char(strcat({'xstop = '},char(num2str(wake_cutoff)))));
% f = strrep(f,char('xstop = x'),char(strcat({'xstop = '},char(num2str(20*2*pi/lambda(ix))))));
f = strrep(f,char('New_V27.geom'),char(['New_V27' num2str(ix) '.geom']));

% Write to new INPUT1.in file
inputfile{ix} = ['INPUT',num2str(ix),'.in'];
fid  = fopen(inputfile{ix},'w');
fprintf(fid,'%s',f);
fclose(fid);
end



parfor i = 1:length(lambda)
% Run CACTUS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cmd = char(strcat({'cactus '},inputfile{i}));
system(cmd);
display(char(strcat({'TSR = '},num2str(lambda(i)),{' complete'})))
end


for i = 1:length(lambda)
% Load data from csv file generated 'INPUT1_RevData.csv'
filename_rev = strcat('INPUT',num2str(i),'_RevData.csv');
filename_param = strcat('INPUT',num2str(i),'_Param.csv');
filename_time = strcat('INPUT',num2str(i),'_TimeData.csv');
filename_el = strcat('INPUT',num2str(i),'_ElementData.csv');
% filename_wake = strcat('INPUT',num2str(ix),'_WakeDefData.csv');
% filename_wake_el = strcat('INPUT',num2str(ix),'_WakeData.csv');

delimiterIn = ',';
headerlinesIn = 1;
import = importdata(filename_rev,delimiterIn,headerlinesIn);
importel = importdata(filename_el,delimiterIn,headerlinesIn);
% import_wake = importdata(filename_wake,delimiterIn);
DATA{1,i} = import.data(end,:);
DATA_ALL{1,i} = import.data;
DATA_ELEMENT{1,i} = importel.data;
% DATA_WAKE{ix,iy} = import_wake.data;
% Delete CACTUS generated files
% delete(inputfile,filename_rev,filename_param,filename_time,filename_el,filename_wake,filename_wake_el);
delete(inputfile{i},geomfile{i},filename_rev,filename_param,filename_time,filename_el);


end

    
% Write to new DATA_ALL to .mat file
save('CACTUS_OUT_REV_CLEAN_A','DATA_ALL')

% LARGE FILE
save('CACTUS_OUT_ELEM_CLEAN_A','DATA_ELEMENT','-v7.3') 

% 
% save('CACTUS_OUT_WAKE_CLEAN','DATA_WAKE')


toc



