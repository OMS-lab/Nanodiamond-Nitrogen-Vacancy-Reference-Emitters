%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to plot power dependence measurements for an NV center along with
its components (linear laser offset, NV emission)
error is calculated from the 68% confidence interval
function fitresult = plot_pdep(X,Y,fitresult)
    if nargin < 3
        [fitresult, ~] = psat_fit(X,Y);
    end
    coeffs = [coeffvalues(fitresult); confint(fitresult,0.68)]';
    coeffs = [coeffs(:,1)-coeffs(:,2),coeffs(:,1),coeffs(:,3)-coeffs(:,1)]; % sorting in to lower,val,upper
    k_inf = coeffs(1,:);
    Psat = coeffs(2,:);

    psat_full = @(P) k_inf(2)*P/(P+Psat(2)) + coeffs(3,2)*P;
    psat_ND = @(P) k_inf(2)*P/(P+Psat(2));
    psat_laser = @(P) coeffs(3,2)*P;
        %Fitting pdep curve
        fplot(psat_full,"LineWidth",2,LineStyle="-")
        hold on
        fplot(psat_ND,LineStyle="-")
        fplot(psat_laser,LineStyle="-")
        plot(X,Y,'kx',MarkerSize=5)
        xline(coeffs(2,2),'r--')
        hold off
        xlabel("Power (µW)")
        ylabel("Countrate, I (cps)")
        axis tight
        title(strcat("P_{sat} = ",num2str(round(fitresult.b,3,"significant"))," µW"))
        legend(["I","I_{ND}","I_{laser}","","P_{sat}"],Location="best")
        figform
        ylim([0 1.1*max(Y)])
        xlim([0 max(X)])
end