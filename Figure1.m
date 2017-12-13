% Simulation of Bayesian change of mind model
% Reproduces Figure 1 in Fleming, van der Putten & Daw
%
% Steve Fleming stephen.fleming@ucl.ac.uk 2017

clear all
close all

% ntrials = 10000;  % N trials per condition
ntrials = 100;
k1 = 4; % mapping from coherence to evidence 
m = 0; % choice bias
coh = linspace(0.1,0.5,3); % coherence values 

% Likelihood distributions marginalizing coherence
muT = sum(k1.*coh)/3;
varT = (sum((k1.*coh - muT).^2)./3) + 1;

% Simulate trials
for pre = 1:length(coh)
    for post = 1:length(coh)
        
        theta1 = coh(pre);
        theta2 = coh(post);
        
        for i = 1:ntrials
            
            s = rand;
            if s < 0.5 
                d(i) = -1;
            else
                d(i) = 1;
            end
            
            x1 = normrnd(d(i)*k1*theta1, 1);
            x2 = normrnd(d(i)*k1*theta2, 1);
            
            if x1 > m
                a = 1;
            else
                a = 0;
            end
            
            if (d(i) == 1 & a == 1) | (d(i) == -1 & a == 0)
                acc(i) = 1;
            else
                acc(i) = 0;
            end
            
            % Compute log-odds of d=1 (rightward motion) marginalizing over thetas
            LORpre(i) = (2*muT*x1)/varT;
            LORpost(i) = (2*muT*x2)/varT;
            LORtotal(i) = LORpre(i) + LORpost(i);
            
            % Change in coordinate frame for LO correct
            if a == 1
                LOCpre(i) = LORpre(i);
                LOCpost(i) = LORpost(i);
                LOCtotal(i) = LORtotal(i);
            else
                LOCpre(i) = -LORpre(i);
                LOCpost(i) = -LORpost(i);
                LOCtotal(i) = -LORtotal(i);
            end
            
        end
        % Mean change in log-odds induced by post-decision evidence split
        % by correct/error, and right/left motion direction
        mean_LORpost_cor_R(pre,post) = mean(LORpost(acc == 1 & d == 1));
        mean_LORpost_err_R(pre,post) = mean(LORpost(acc == 0 & d == 1));
        mean_LORpost_cor_L(pre,post) = mean(LORpost(acc == 1 & d == -1));
        mean_LORpost_err_L(pre,post) = mean(LORpost(acc == 0 & d == -1));
        mean_LOCpost_cor_R(pre,post) = mean(LOCpost(acc == 1 & d == 1));
        mean_LOCpost_err_R(pre,post) = mean(LOCpost(acc == 0 & d == 1));
        mean_LOCpost_cor_L(pre,post) = mean(LOCpost(acc == 1 & d == -1));
        mean_LOCpost_err_L(pre,post) = mean(LOCpost(acc == 0 & d == -1));

        mean_LOCcollapse_cor(pre,post) = mean(LOCpost(acc == 1));
        mean_LOCcollapse_err(pre,post) = mean(LOCpost(acc == 0));
    end
end

% % get means of columns to average over pre
post_LORacc = [mean(mean_LORpost_cor_L) mean(mean_LORpost_cor_R) mean(mean_LORpost_err_L) mean(mean_LORpost_err_R)];
post_LOCacc = [mean(mean_LOCpost_cor_L) mean(mean_LOCpost_cor_R) mean(mean_LOCpost_err_L) mean(mean_LOCpost_err_R)];

%% Summary figures collapsing over levels of pre coherence
h1 = figure;
set(gcf, 'Position', [400 400 1200 300])
subplot(1,3,1)
plot([1 2 3], post_LORacc(1:3), 'g--', 'LineWidth', 3)
hold on
plot([1 2 3], post_LORacc(4:6), 'g', 'LineWidth', 3)
plot([1 2 3]+0.2, post_LORacc(7:9), 'r--', 'LineWidth', 3)
plot([1 2 3]+0.2, post_LORacc(10:12), 'r', 'LineWidth', 3)
set(gca, 'XLim', [0.5 3.5], 'XTick', [1 2 3], 'XTickLabel', {'Low', 'Med', 'High'}, 'FontSize', 16);
xlabel('Postdecision motion strength','FontSize',20);
ylabel('\Delta Log-odds rightward','FontSize',20)
legend('Left, correct', 'Right, correct', 'Left, error', 'Right, error', 'Location', 'East')
legend boxoff
box off

subplot(1,3,2)
plot([1 2 3], post_LOCacc(1:3), 'g--', 'LineWidth', 3)
hold on
plot([1 2 3]+0.2, post_LOCacc(4:6), 'g', 'LineWidth', 3)
plot([1 2 3], post_LOCacc(7:9), 'r--', 'LineWidth', 3)
plot([1 2 3]+0.2, post_LOCacc(10:12), 'r', 'LineWidth', 3)
set(gca, 'XLim', [0.5 3.5], 'XTick', [1 2 3], 'XTickLabel', {'Low', 'Med', 'High'}, 'FontSize', 16);
xlabel('Postdecision motion strength','FontSize',20);
ylabel('\Delta Log-odds correct','FontSize',20)
legend('Left, correct', 'Right, correct', 'Left, error', 'Right, error', 'Location', 'East')
legend boxoff
box off

% Confidence and value, show nonlinear transform from log-odds correct
baseLOC = linspace(-4,4,500);
conf = 1./(1+exp(-baseLOC));
acc = conf > 0.5;
qsr = 1-((acc - conf).^2);
subplot(1,3,3)
% line([0 0], [0 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
plot(baseLOC, conf, 'k', 'LineWidth', 3)
hold on
plot(baseLOC, qsr, 'k--', 'LineWidth', 3)
set(gca, 'FontSize', 16, 'XLim', [min(baseLOC) max(baseLOC)]);
xlabel('Final log-odds correct','FontSize',20);
ylabel('Confidence / value','FontSize',20)
fill([min(baseLOC) 0 0 min(baseLOC)], [0 0 1 1],[0.6 0.6 1]);
fill([0 max(baseLOC) max(baseLOC) 0], [0 0 1 1],[0.8 0.8 0.8]);
plot(baseLOC, conf, 'k', 'LineWidth', 3)
plot(baseLOC, qsr, 'k--', 'LineWidth', 3)
box off
legend('Confidence', 'Value', 'Location', 'West')
legend boxoff