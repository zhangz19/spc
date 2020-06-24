function [] = spc_demo()
% Demo for SPC: Spatial clustering for method described in: 
% Feng, W., Lim, C. Y., Maiti, T., & Zhang, Z. (2016). Spatial regression and estimation of disease risks: 
% A clustering?based approach. Statistical Analysis and Data Mining: The ASA Data Science Journal, 9(6), 417-434.

rng('default');  rng(8);  %set random seeds for reproducibility of the example result

% Observe/Expect: N by 1 observed/expected counts (Expect=1 if not available) 
% X: design matrix: [1, poverty/1e4, urban/1e4, abovehigh/1e4] for demo
% W: N by N spatial adjacency matrix
% Lat/Long: latitude/longitude of the centrods for computing distance
% ID: unique ID code for each region, = 1:numel(Observe) if not available
load('spc_demo_data.mat','Expect','Observe','X','W','Lat','Long','ID') %exmaple Michigan demo data

% set environmenal variables (ev)
ev.Expected = Expect;  ev.Observed = Observe;
ev.X = X;   ev.W = W;
ev.Latitude = Lat;  ev.Longitude = Long;
ev.DistIndex = ID;  ev.n=length(ID);
ev.nrmin = 5; %number of minimal size of a spatial cluster
% set prior
ev.A = 2.000122;   ev.B = 0.03500429;  %for inverse gamma prior for variance
ev.a = 2.000122;   ev.b = 0.03500429;
ev.mean=zeros(size(ev.X,2),1);   ev.var=10*eye(size(ev.X,2));  %prior for regression coefficients
ev.invvar=inv(ev.var); %precision
ev.Sigma2Range = 0.001:0.001:2;

% general MCMC setting
ev.verbose = true;
ev.tot = 100; %50000;  % used a small number of runs for demo only
ev.burnin = ceil(ev.tot/2);  %use second half
ev.nchain = 1; %for this demo. For multiple chains, need to modify spc_summary.m accordingly 

% run the model, store results
result = cell(1, ev.nchain);   for ch=1:ev.nchain;  result{ch} = spc(ev);  end

% check clustering results
out = spc_summary(result, ev); 
tabulate(out.labs)  %frequency of clustering labels

% clsuter-wise summary of regression coefficients (row: cluster, column: p-covariates)
out.BetaR_mean

end



