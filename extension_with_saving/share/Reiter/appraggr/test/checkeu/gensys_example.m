% Neoclassical growth model solved using gensys
% Solve the following social planner problem:
%   utility function:  (c^(1 - eta) - 1)/(1 - eta)
%   resouce constraint: c + k' = z*k^alpha + (1 - delta)*k
%   technology: ln(z') = (1-psi)ln(zbar) + psi ln(z) + epsilon
% To run this program, you will need to download Sims' gensys.m, qzdiv.m,
% and qzswitch.m from http://sims.princeton.edu/yftp/gensys/
% Last Updated: Dec. 2, 2008 by Nora Traum

clear all; clc;
rand('state', 87);
randn('state', 87);

%------------------------
% Calibrated Parameters 
%------------------------
beta = 0.99;
psi = 0.85;
alpha = 0.36;
eta = 1;
delta = 0.025;
zbar = 1;

%-------------------------
% Steady state values
%-------------------------
rss = 1/beta;
kss = ((alpha*zbar)/(rss - 1 + delta))^(1/(1 - alpha));
yss = zbar*(kss^alpha);
css = yss - delta*kss;


% Specify size of system of equations 
nvar = 5;           % number of variables
nshocks = 1;        % number of shocks
nforerrors = 2;     % number of forecast errors

%--------------------------------------------------
% Cast the system in Gensys form:
% G0*X(t) = G1*X(t-1) + Psi*e(t) + Pi*eta(t) + CC
%--------------------------------------------------
% Initializations
G0 = zeros(nvar,nvar);
G1 = zeros(nvar,nvar);
CC = zeros(nvar,1);
Psi = zeros(nvar,nshocks);
Pi = zeros(nvar,nforerrors);

% Define variables
nC = 1;
nK = 2;
nZ = 3;
nY = 4;
nR = 5;

% Define Log-Linearized Equations
%--------------------------------------------
% (1)	Aggregate Resource Constraint
%--------------------------------------------
G0(1,nC) = 1;
G0(1,nK) = kss/css;
G0(1,nZ) = -yss/css;
G1(1,nK) = (kss*rss)/css;
%--------------------------------------------
% (2)	Marginal cost of capital
%--------------------------------------------
G0(2,nR) = 1;
G0(2,nZ) = -(1 - beta*(1 - delta));
G1(2,nK) = -(1 - beta*(1 - delta))*(1 - alpha);
%--------------------------------------------
% (3)	Euler Equation
%--------------------------------------------
G0(3,nR) = 1;
G0(3,nC) = -eta;
G1(3,nC) = -eta;
Pi(3,1) = 1;
Pi(3,2) = -eta;
%--------------------------------------------
% (4)	Process for Z
%--------------------------------------------
G0(4,nZ) = 1;
G1(4,nZ) = psi;
Psi(4,1) = 1;
%--------------------------------------------
% (5)	Production function
%--------------------------------------------
G0(5,nY) = 1;
G0(5,nZ) = -1;
G1(5,nK) = alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution methods

% Call to Sims' gensys
addpath /home/amckay/Dropbox/computation/gensys;
[G,C0,M,fmat,fwt,ywt,gev,eu] = gensys(G0,G1,CC,Psi,Pi);
disp(['eu:    ' mat2str(eu)])
% M is often called 'impact'
%G is sometimes called G1


%call to checkeu (Reiter's implementation)
div = 1+1e-10;
[eu_reit,eigvals_reit,G1_reit,impact_reit] = checkeu(G0,G1,Psi,Pi, div);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%-----------------------------------
% Compute impulse responses
%-----------------------------------
% These are the responses following a 1% increase in the shocks.
steps = 40;
imf = zeros(nvar,nshocks,steps);	
imf(:,:,1) = M;
for t = 2:steps
   imf(:,:,t) = G*imf(:,:,t-1);
end

%eliminate very small responses
for i = 1:nvar
   for j = 1:nshocks
      for k = 1:steps
         if abs(imf(i,j,k))<10^-10
            imf(i,j,k) = 0;
         end
      end
   end
end

Time = (0:(steps-1))/4;  % converts time dimension to years
vnames = {'consumption';'capital';'technology'; 'output'}; 

for k = 1:nshocks
    figure(k);
    for j = 1:(nvar-1)
        subplot(2,2,j);
        imf1(1,:)=imf(j,k,1:steps);
        plot(Time,imf1,'b-','MarkerSize',3);
        grid on;
        title([char(vnames(j))], 'fontsize', 8);
    end
end


%------------------------
% Simulate data 
%------------------------
simlength = 220;        % number of simulated draws
burn = 20;              % burn-in period
cov = 0.712;            % variance-covariance matrix of shocks, units: %

sig = chol(cov);
epsvec = sig*randn(nshocks,simlength);

gendata = zeros(nvar,simlength);
for i = 2:simlength
    gendata(:,i) = G*gendata(:,i-1) + M*epsvec(:,i);
end

gendata = gendata(:,burn+1:simlength);

% Plot Data
figure(2)
for j = 1:(nvar-1)
    subplot(2,2,j);
    plot(gendata(j,:));
    grid on;
    title([char(vnames(j))], 'fontsize', 8);
end






%% show that the two methods give the same answer

disp 'Are the two methods the same?'
if all(all(abs([(G - G1_reit) (M - impact_reit)])<1e-15)) == 1
    disp 'true'
else
    disp 'false'
end












