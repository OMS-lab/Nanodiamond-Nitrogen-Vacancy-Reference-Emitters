%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs ===
% correlations = array of times within precalculated coherence window [s]
% N1 = Countrate from detector 1
% N2 = Countrate from detector 2
% tbin = Bin duration for coincidence counts
% coherence_w = coherence window length [s]
% ID = ID of object (e.g. 29 for ND no. 29)
% 
% Outputs ===
% tau_c = delay times [s], adjusted for t0 - i.e. systematic delay between 2 detectors.
% g2y = g2(t) y-values
% N = binned coincidences
% par = g2(t) fitting parameters
% RMSE = Root-mean-square error between 1 µs < |tau| < 5 µs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitresult,X,Y,N,g2t_error,RMSE,GOF] = g2_calc(correlations,N1,N2,tbin,coherence_w,dwell_time,filtering)
    
    if nargin < 4
        tbin = 1;
        coherence_w = 100;
    end

    if nargin < 5
        coherence_w = 100;
    end

    if nargin < 7
        filtering = 0;
    end

    %[ns] -> [s]
    tbin = 1e-9*tbin;
    correlations = 1e-9*correlations;
    coherence_w = 1e-9*coherence_w;


    [N,tau] = histcounts(correlations,-coherence_w:tbin:coherence_w);
    tau = tau(2:end)-tbin/2; %adjust x values to midpoint of bins
    
    % J. Opt. Soc. Am. B / Vol. 24, No. 12 / December 2007 (M. Beck,
    % Optica Publishing)


%### EQN: g2y = t_bin*c(t)/(T*P1*P2) where;
%     g2y = g2(t) values normalised to a poissonian source (var = tbin)
%     t_bin = bin duration (var = tbin)
%     c(t) = coincidences at each delay time (var = tbin)
%     T = total acquisition change (var = tbin)
%     P(1,2) = Poissonian Distribution Probability of detecting a photon from detector (1, 2) (var = tbin)
     g2y = tbin*N./(dwell_time*poisspdf(1,tbin*N1)*poisspdf(1,tbin*N2));

     if filtering == true
      %Remove "devil horn" peaks due to reflection of photons at SPAD
      %collection spot interface (these peaks occur at a specific delay
      %time linked to the length of the optical fibre at collection
          g0_mask = and(tau < 50e-9,tau > -50e-9);
          subset_g2y = g2y(g0_mask);
          h_window = 3; %window length to median around maximum values
          for i = 1:4 %since there are always two peaks (with up to 4 data points) due to this issue
              [~,I] = max(subset_g2y);
              if (I <= 3 || I > numel(subset_g2y)-h_window)
                  warning("g(2)(t) outlier is at edge of ±50 ns window")
              else
                subset_g2y(I) = median(subset_g2y(I-h_window:I+h_window)); %replace outlier
              end
          end
          g2y(g0_mask) = subset_g2y;
     end
    %   g^{2}_{Fit}(\tau) = 1-a[bexp(-|tau-t0|/c + (1-b)exp(-|\tau-t0|/d)]
    [fitresult, GOF] = g2fit(tau, g2y);
    g2t_error = confint(fitresult,0.68);
    g2t_error(2,1) = min([1,g2t_error(2,1)]); %ensure max "a" upper value is =<1
    tau_c = tau-fitresult.e; %Adjust g2(0) point at t0
    yfit = fitresult(tau_c);
    
    b_tau = 1e-9*abs(1e9*tau_c) > 1e-6;
    RMSE_yval = sqrt((g2y - yfit).^2);
    RMSE = mean(RMSE_yval(b_tau));

    X = tau;
    Y = g2y;
end
