function data = metaModelFit_generate_data(ntrials, coh, theta1_index, theta2_index, conf_sigma, model, varargin)
% Generate per-condition predictions from different models in change of mind
% study
% Takes as inputs parameters from STAN fits, and a model type
% (see plot_simpleModel_fits.m)
%
% Stephen Fleming stephen.fleming@ucl.ac.uk 

% parse inputs
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch lower(arg)
            case 'k1', k1 = varargin{i+1};
            case 'm', m = varargin{i+1};
            case 'w', w = varargin{i+1};
            case 'gamma', gamma = varargin{i+1};
        end
    end
end

theta1 = coh(theta1_index);
theta2 = coh(theta2_index);
muT = sum(k1.*coh)/3;
varT = (sum((k1.*coh - muT).^2)./3) + 1;

for i = 1:ntrials
   
   flip_d = rand;
   if flip_d > 0.5
       d = 1;
   else
       d = -1;
   end
   
   x1 = normrnd(d*k1*theta1, 1);
   x2 = normrnd(d*k1*theta2, 1);
   
   prob_a = 1/(1+exp(-(100*(x1-m))));
   if prob_a > 0.5
       a = 1;
   else
       a = 0;
   end
   
   if (d == 1 & a == 1) | (d == -1 & a == 0)
       acc = 1;
   else
       acc = 0;
   end   
   
   % Compute log-odds of d=1 marginalizing over thetas
   loglikdir_pre = (2*muT*x1)/varT;
   loglikdir_post = (2*muT*x2)/varT;

   if a == 1
       loglikC_pre = loglikdir_pre;
       loglikC_post = loglikdir_post;
   else
       loglikC_pre = -loglikdir_pre;
       loglikC_post = -loglikdir_post;
   end
   
   switch model
       case 'ideal'
            loglikC = loglikC_pre + loglikC_post;
            conf = 1/(1+exp(-loglikC));
       case 'weighted'
            loglikC = (2*w*loglikC_pre) + ((2-(2*w))*loglikC_post);
            conf = 1/(1+exp(-loglikC));
       case 'accweighted'
            if acc == 1
                loglikC = loglikC_pre + loglikC_post;
            else
                loglikC = (2*w*loglikC_pre) + ((2-(2*w))*loglikC_post);
            end
            conf = 1/(1+exp(-loglikC));
       case 'mapping'
            loglikC = loglikC_pre + loglikC_post;
            temp_conf = 1/(1+exp(-loglikC));
            LO_pi_conf = gamma.*(log(temp_conf./(1-temp_conf)));
            conf = 1/(1+exp(-LO_pi_conf));
       case 'weightedmapping'
            loglikC = (2*w*loglikC_pre) + ((2-(2*w))*loglikC_post);
            temp_conf = 1/(1+exp(-loglikC));
            LO_pi_conf = gamma.*(log(temp_conf./(1-temp_conf)));
            conf = 1/(1+exp(-LO_pi_conf));
       case 'accweightedmapping'
            if acc == 1
                loglikC = loglikC_pre + loglikC_post;
            else
                loglikC = (2*w*loglikC_pre) + ((2-(2*w))*loglikC_post);
            end
            temp_conf = 1/(1+exp(-loglikC));
            LO_pi_conf = gamma.*(log(temp_conf./(1-temp_conf)));
            conf = 1/(1+exp(-LO_pi_conf));
   end
   
   % Store data
   r = normrnd(conf, conf_sigma);
   data.conf(i) = r;
   if data.conf(i) > 1
       data.conf(i) = 1;
   elseif data.conf(i) < 0
       data.conf(i) = 0;
   end
   data.loglik_pre(i) = loglikC_pre;
   data.loglik_post(i) = loglikC_post;
   data.acc(i) = acc;
   data.d(i) = d;
   data.a(i) = a;
end

