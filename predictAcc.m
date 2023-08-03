function [predAcc,predAccCI,BER,EbNo] = predictAcc(corrs,method,CIflag)
% PREDICTACC Predict the accuracy of a correlation-based stimulus
% decoder for AAD using BER estimation in a BPSK system with AWGN.
%
%   Input parameters:
%       corrs [DOUBLE]: matrix of correlations (nb of decision windows x 2)
%       method [STRING]: 'mle' (maximum likelihood estimation) or 'mom'
%           (method of moments)
%       CIflag [BOOLEAN]: true when to compute confidence intervals, false
%           when not
%
%   Output:
%       predAcc [DOUBLE]: the predicted accuracy
%       predAccCI [DOUBLE]: the 95%-confidence interval (1 x 2)
%       BER [DOUBLE]: the bit-error rate in log-scale
%       EbNo [DOUBLE] : the computed Eb/No-ratio in dB

% Author: Simon Geirnaert, KU Leuven, ESAT & Dept. of Neurosciences
% Correspondence: simon.geirnaert@esat.kuleuven.be

options = optimoptions('fsolve','Display','none');
B = 1000; % number of resamplings for bootstrapping


%% Unsupervised estimation
% estimate standard deviation of difference via sum
sigmaDiff = std(sum(corrs,2));

% estimate mean of differences via folded normal distribution
X = abs(diff(corrs,[],2));
switch method
    case 'mle' % maximum likelihood estimation
        mDiff = fsolve(@(x)mMLE(x,sigmaDiff,X),1,options);
    case 'mom' % method of moments
        mDiff = fsolve(@(x)mMoM(x,sigmaDiff)-mean(X),1,options);
end

% predict accuracy via BER for BPSK and AWGN
predAcc = 1-1/2*erfc(mDiff/(sqrt(2)*sigmaDiff));

if CIflag
    %% 95%-confidence interval via bootstrapping
    switch method
        case 'mle' % maximum likelihood estimation
            predAccCI = bootci(B,{@(C)1-1/2*erfc(fsolve(@(x)mMLE(x,std(sum(C,2)),abs(diff(C,[],2))),1,options)/(sqrt(2)*std(sum(C,2)))),corrs},'Type','bca'); % bias-corrected and accelerated
        case 'mom' % method of moments
            predAccCI = bootci(B,{@(C)1-1/2*erfc(fsolve(@(x)mMoM(x,std(sum(C,2)))-mean(abs(diff(C,[],2))),1,options)/(sqrt(2)*std(sum(C,2)))),corrs},'Type','bca'); % bias-corrected and accelerated
    end
else
    predAccCI = [nan,nan];
end

BER = log10(1-predAcc); % BER
EbNo = 10*log10(mDiff^2/(2*sigmaDiff^2)); % Eb/No-ratio

end