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

function fitresult = g2_plot(X,Y,inset,fitresult,plot_components)
    tau = X;
    g2y = Y;

    %   g^{(2)}_{Fit}(\tau) = 1-a[bexp(-|tau-t0|/c + (1-b)exp(-|\tau-t0|/d)]
    if nargin < 3
        inset = 0;
    end
    if nargin < 4
        [fitresult, ~] = g2fit(tau, g2y);
    end
    if nargin < 5
        plot_components = false;
    end
    g2t_error = confint(fitresult,0.68);
    g2t_error(2,1) = min([1,g2t_error(2,1)]); %ensure max "a" upper value is =<1
    tau_c = tau-fitresult.e; %Adjust g2(0) point at t0
    yfit = fitresult(tau);
    

    %PLOT
    figure(1)
    if plot_components == false
        plot(1e9.*tau_c,g2y,'kx-')
        hold on
        plot(1e9.*tau_c,yfit,'r-','LineWidth',2)
        yline(1 - fitresult.a,'r--')
        xlabel('Delay Time, \tau (ns)');ylabel('g^{(2)}(\tau)');
    else
        plot(tau_c,g2y,'kx-')
        hold on
        plot(tau_c,yfit,'r-','LineWidth',2)
        yline(1 - fitresult.a,'r--')
        %components
        f_antibunching = @(tau) 1-fitresult.a*(fitresult.b*exp(-abs(tau)/fitresult.c) - fitresult.a); %subtract "a" at end, so dip isn't below 0
        f_bunching = @(tau) 1-fitresult.a*(1-fitresult.b)*exp(-abs(tau)/fitresult.d);
        fplot(f_antibunching,'b-','LineWidth',2)
        fplot(f_bunching,'g-')
        xlabel('Delay Time, \tau (s)');ylabel('g^{(2)}(\tau)');
    end
    hold off
    title(strcat("g^{(2)}(0) = ",num2str(1-fitresult.a),", \tau = ",num2str(round(fitresult.c*1e9,2,"decimal"))," ns"))
    legend(["g^{(2)}_{data}(\tau)","g^{(2)}_{fit}(\tau)","g^{(2)}(0)"],"Location","southwest")
    figform
    if plot_components == false
        xlim([-100 100])
    else
        xlim([-100*1e-9 100*1e-9])
    end
    ylim([0 1.5*max(yfit)])
        
    if inset == true
        axes('Position',[.63 .20 .25 .25],'TickDir','in','YTick',[0,0.5,1],'YTickLabel',{'0','0.5','1'})
        hold on
        plot(1e9*tau,g2y,'kx','MarkerSize',2)
        plot(1e9*tau_c,yfit,'r-','LineWidth',2)
        ylim([0 1.2*max(yfit)])
        hold off
        axis tight
        box on
        figform
        xlim([-1000 1000])
        ylim([0 1.6])
        set(gca,'FontSize',8,'TickDir','in','linewidth',2,'FontWeight','bold')
    end
    %set(gca,'FontSize',8,'TickDir','in','linewidth',2,'FontWeight','bold')
    
    b_tau = 1e-9*abs(1e9*tau_c) > 1e-6;
    RMSE_yval = sqrt((g2y - yfit).^2);
    RMSE = mean(RMSE_yval(b_tau));
    
    %Output values to console
    disp('___________________________________________________')
    disp(['Bunching lifetime = ',num2str(fitresult.d*1e9),' ns'])
    disp(['Offset = ',num2str(fitresult.e*1e9),' ns'])
    disp(['NV Lifetime = ',num2str(fitresult.c*1e9),' ns'])
    disp(['g2(0) = ',num2str(1 - fitresult.a),'(+',num2str(fitresult.a - g2t_error(1,1)),')','(',num2str(fitresult.a - g2t_error(2,1)),')'])
    disp(['RMSE = ',num2str(RMSE)])
end

function figform
    axis tight
    box on
    set(gca,'FontSize',12,'TickDir','in','linewidth',2,'FontWeight','bold')
end

function output = multmovmean(input,n,win)
    for i = 1:n
        input = movmedian(input,win);
    end
    output = input;
end