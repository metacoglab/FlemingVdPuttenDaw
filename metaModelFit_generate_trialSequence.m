function data = metaModelFit_generate_trialSequence(d, action, conf, coh, theta1_index, theta2_index, conf_sigma, model, varargin)
% Generate trial series of predictions from model parameters conditional on
% known stimulus, action and confidence
% Returns log-probability of each datapoint (both decision and confidence)
%
% SF 2017

% parse inputs
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch lower(arg)
            case 'k1', k1 = varargin{i+1};
            case 'm', m = varargin{i+1};
            case 'b', b = varargin{i+1};
            case 'w1', w1 = varargin{i+1};
            case 'w2', w2 = varargin{i+1};
            case 'gamma', gamma = varargin{i+1};
            case 'brt',  brt = varargin{i+1};
            case 'logrt', logrt = varargin{i+1};
        end
    end
end

theta1 = coh(theta1_index);
theta2 = coh(theta2_index);
muT = sum(k1.*coh)/3;
varT = (sum((k1.*coh - muT).^2)./3) + 1;

for i = 1:length(d)
    
    x1 = normrnd(d(i)*k1*theta1(i), 1);
    x2 = normrnd(d(i)*k1*theta2(i), 1);
    
    % Define model action and accuracy for this trial
    prob_a = 1/(1+exp(-(100*(x1-m))));
    if prob_a > 0.5
        a(i) = 1;
    else
        a(i) = 0;
    end
    
    if (d(i) == 1 & a(i) == 1) | (d(i) == -1 & a(i) == 0)
        acc = 1;
    else
        acc = 0;
    end
    
    % Compute log-odds of d=1 marginalizing over thetas
    loglikdir_pre = (2*muT*x1)/varT;
    loglikdir_post = (2*muT*x2)/varT;
    
    if a(i) == 1
        loglikC_pre = loglikdir_pre;
        loglikC_post = loglikdir_post;
    else
        loglikC_pre = -loglikdir_pre;
        loglikC_post = -loglikdir_post;
    end
    
    switch model
        case 'ideal'
            loglikC = loglikC_pre + loglikC_post;
            modelConf = 1/(1+exp(-loglikC));
            model_loglikC_pre = loglikC_pre;
            model_loglikC_post = loglikC_post;
        case 'weighted'
            loglikC = (w1*loglikC_pre) + (w2*loglikC_post);
            modelConf = 1/(1+exp(-loglikC));
            model_loglikC_pre = w1*loglikC_pre;
            model_loglikC_post = w2*loglikC_post;
        case 'choiceweighted'
            if ((x2 > 0) && a(i) == 1) || ((x2 <= 0) && a(i) == 0)
                loglikC = loglikC_pre + w1*loglikC_post;
                model_loglikC_post = w1*loglikC_post;
            else
                loglikC = loglikC_pre + w2*loglikC_post;
                model_loglikC_post = w2*loglikC_post;
                
            end
            modelConf = 1/(1+exp(-loglikC));
            model_loglikC_pre = loglikC_pre;
        case 'choicebias'
            if a(i) == 1
                loglikdir_bias = log(b./(1-b)); % positive towards rightward choice
            else
                loglikdir_bias = log((1-b)./b);  % negative towards leftward choice
            end
            loglikdir_total = loglikdir_pre + loglikdir_post + loglikdir_bias;
            if a(i) == 1
                loglikC = loglikdir_total;
            else
                loglikC = -loglikdir_total;
            end
            modelConf = 1/(1+exp(-loglikC));
            model_loglikC_pre = loglikC_pre;
            model_loglikC_post = loglikC_post;
        case 'mapping'
            loglikC = loglikC_pre + loglikC_post;
            temp_conf = 1/(1+exp(-loglikC));
            LO_pi_conf = gamma.*(log(temp_conf./(1-temp_conf)));
            modelConf = 1/(1+exp(-LO_pi_conf));
            model_loglikC_pre = loglikC_pre;
            model_loglikC_post = loglikC_post;
        case 'ideal_dt'
            loglikC = loglikC_pre + loglikC_post + brt.*(logrt(i));
            modelConf = 1/(1+exp(-loglikC));
            model_loglikC_pre = loglikC_pre;
            model_loglikC_post = loglikC_post;
        case 'weighted_dt'
            loglikC = (w1*loglikC_pre) + (w2*loglikC_post) + brt.*(logrt(i));
            modelConf = 1/(1+exp(-loglikC));
            model_loglikC_pre = w1*loglikC_pre;
            model_loglikC_post = w2*loglikC_post;
        case 'choiceweighted_dt'
            if ((x2 > 0) && a(i) == 1) || ((x2 <= 0) && a(i) == 0)
                loglikC = loglikC_pre + w1*loglikC_post + brt.*(logrt(i));
                model_loglikC_post = w1*loglikC_post;
            else
                loglikC = loglikC_pre + w2*loglikC_post  + brt.*(logrt(i));
                model_loglikC_post = w2*loglikC_post;
                
            end
            modelConf = 1/(1+exp(-loglikC));
            model_loglikC_pre = loglikC_pre;
        case 'choicebias_dt'
            if a(i) == 1
                loglikdir_bias = log(b./(1-b)); % positive towards rightward choice
            else
                loglikdir_bias = log((1-b)./b);  % negative towards leftward choice
            end
            loglikdir_total = loglikdir_pre + loglikdir_post + loglikdir_bias;
            if a(i) == 1
                loglikC = loglikdir_total + brt.*(logrt(i));
            else
                loglikC = -loglikdir_total + brt.*(logrt(i));
            end
            modelConf = 1/(1+exp(-loglikC));
            model_loglikC_pre = loglikC_pre;
            model_loglikC_post = loglikC_post;
        case 'mapping_dt'
            loglikC = loglikC_pre + loglikC_post + brt.*(logrt(i));
            temp_conf = 1/(1+exp(-loglikC));
            LO_pi_conf = gamma.*(log(temp_conf./(1-temp_conf)));
            modelConf = 1/(1+exp(-LO_pi_conf));
            model_loglikC_pre = loglikC_pre;
            model_loglikC_post = loglikC_post;
        case 'random'
            modelConf = rand;
            prob_a = rand;
            model_loglikC_pre = 0;
            model_loglikC_post = 0;
    end
    
    % Store output as vectors
    r = normrnd(modelConf, conf_sigma);
    if r > 1
        r = 1;
    elseif r < 0
        r = 0;
    end
    data.modelConf(i) = r;
    data.modelProbA(i) = prob_a;
    data.modelAcc(i) = acc;
    data.loglik_pre(i) = model_loglikC_pre;
    data.loglik_post(i) = model_loglikC_post;
    
    % Store log-probability of model prediction with reference to actual data
    if action(i) == 1
        data.logprobAction(i) = max(log(prob_a), log(eps(0)));
    else
        data.logprobAction(i) = max(log(1-prob_a), log(eps(0)));
    end
    data.logprobConf(i) = max(log(normpdf(conf(i),modelConf,conf_sigma)), log(eps(0)));
end

