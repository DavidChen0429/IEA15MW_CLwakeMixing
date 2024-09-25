clear
close all
addpath('.\Functions');

%% Load model
sys = load('ModelOrder4.mat.mat');
A = sys.OLi.A;
B = sys.OLi.B;
C = sys.OLi.B;
D = sys.OLi.B;

% decouple matrix 
% delay factor 
% controller tuning