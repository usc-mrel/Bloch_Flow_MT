%%% function to generate Ernst signal curve for given flip (deg), TR, R1
%
% This function is copied from the ihMT_steadystate package published by 
% Shaihan Malik (King's College London), available here:
% https://github.com/mriphysics/ihMT_steadystate

function sig = ernst(flip,TR,R1)
e1=exp(-TR*R1);
sig = sind(flip).*(1-e1)./(1-e1.*cosd(flip));
end