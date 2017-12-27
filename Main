% main
% Main function to call other functions that are application specific
% 
% Usage: m = main
%
% Where: 
%   There are no return values
%
% Promblem: General perpouse nested sampling algorithm
%
% Inputs to the Promblem: 
%   N - Number of Objects
%   Max - Number of iterates
%   Obj - Collection of n objects
%   Samples - Objects stored for posterior results
%   
% Outputs to the Promblem:
%   logwidth - width in prior mass
%   logLstar - likelyhood constant
%   H_ent - prior information
%   logz - Evidence
%   logznew - New Evidence
%
% Originally written in C
% Modified: 
%           William Dahl
%           12 December 2017 
%           Converted to Matlab

function main

n = 100; % Num Objects
MAX = 1000;

Obj = struct([]);
Samples = struct([]);
for i = 1:n
  Obj(i).u = 0.;
  Obj(i).v = 0.;
  Obj(i).x = 0.;
  Obj(i).y = 0.;
  Obj(i).logL = 0.;
  Obj(i).logWt = 0.;
end

for i = 1:MAX
  Samples(i).u = 0.;
  Samples(i).v = 0.;
  Samples(i).x = 0.;
  Samples(i).y = 0.;
  Samples(i).logL = 0.;
  Samples(i).logWt = 0.;
end
H = 0.0;         % Information, init 0
logZ = log(1e-100);  % ln(Evidence Z, init 0);

% Set prior objects
for i = 1:n  
  Obj(i) = Prior(Obj(i));
end
% Outermost interval of prior mass
logwidth = log(1.0 - exp(-1.0 / n));

%% Begin nested sampling loop
nest = 1;
while (nest <= MAX) 
    % (find) Worst object in collection, with Weight = width * likelihood
    worst = 1;
    i = 2;
    while(i <= n)
        if (Obj(i).logL < Obj(worst).logL)  
            worst = i;
        end
        i = i + 1;
    end
    Obj(worst).logWt = logwidth + Obj(worst).logL;

    logZnew = log(exp(logZ) + exp(Obj(worst).logWt)); % Updating Evindence Z and Information H_ent
    H = exp(Obj(worst).logWt - logZnew) * Obj(worst).logL + exp(logZ - logZnew) * (H + logZ) - logZnew;
    logZ = logZnew;
    
    Samples(nest) = Obj(worst);   % Posterior samples (optional, care with storage overflow;
    
    copy = worst;
    while (copy == worst && n > 1) % Kill worst object in favor of copy of different survivor
        copy = int32((n - 1) * rand()) + 1; 
    end   
    logLstar = Obj(worst).logL;         % new likelihood constraint
    Obj(worst) = Obj(copy);             % overwrite worst object
    Obj(worst) = Explore(Obj(worst), logLstar);     % Evolve copied object within constraint
    logwidth = logwidth - (1.0 / n);                     % Shrink Interval
    nest = nest + 1;
end %% end nested sampling loop ============================
nest = nest - 1;
fprintf("# iterates = %d\n", nest);
fprintf("Evidence: ln(Z) = %g +- %g\n", logZ, sqrt(H/n));
fprintf("information: H = %g nats = %g bits\n", H, H/log(2.0));
Results(Samples, nest, logZ);
return
