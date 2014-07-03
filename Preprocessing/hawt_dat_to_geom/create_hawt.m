% create_hawt.m : Reads in parameters from './hawt_params.m', builds a HAWT geometry structure, and writes out 
%				  to a CACTUS-formatted .geom file.
%
%                 The aerodynamic schedule file (tab-delimited: r/R, c/R, twist, airfoil_id) defining each blade
%                 is specified in dat_filename.

clc
clear all
close all

%% Parameters %%
dat_filename  = 'DESIGN_A_dec.dat'
geom_filename = 'DESIGN_A_dec.geom'

%% Import hawt_params.m
% hawt_params holds three structs:
%		rotor_params
%		blade_params
%		grid_params
hawt_params

%% Generate turbine geometry and write to file
T = hawt_dat_to_geom(dat_filename, geom_filename, rotor_params, blade_params, grid_params)

% % Plot controls
% % Plot animated turbine rotation
% XLim=[-1,1];
% YLim=[-1,1];
% ZLim=[-1,1];

% % Plot controls
% PlotVec=1;
% SFVec=.5;
% Trans=.5;

% hf=figure(1);
% set(hf,'Position',[303   124   956   610]) 
% set(gca,'Position',[5.2743e-002  5.1245e-002  8.9979e-001  8.8141e-001])
% set(gca,'CameraPosition',[-52.1999   30.4749   62.2119])
% set(gca,'CameraUpVector',[1.8643e-001  9.7433e-001 -1.2615e-001])
% set(gca,'CameraViewAngle',6.3060e+000)
% grid on
% set(gcf,'Color','white');
% % hl=light('Position',[0,-1,0]);
% set(gca,'Color','white');
% set(gca,'DataAspectRatio',[1,1,1])
% set(gca,'XLim',XLim,'YLim',YLim,'ZLim',ZLim)

% HIn=[];
% PhasePlot=linspace(0,20*pi,1500);
% for i=1:length(PhasePlot)
%    H=PlotTurbineGeom(T,hf,PhasePlot(i),HIn,Trans,PlotVec,SFVec);
%    HIn=H;
   
% end