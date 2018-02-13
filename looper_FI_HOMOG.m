%Compute the Fisher information for many different population sizes
%And strengths of correlations

%Makes Fig. 7BC of 
%Zylberberg, Cafaro, Turner, et al. Neuron 2016
%Direction Selective Circuits Shape Noise to Ensure a Precise Population Code
%contact: joel.zylberberg@ucdenver.edu to report bugs or issues

%Homogeneous tuning curves

%pick the parameter range to consider
rhomaxvals = [0.8 0.6 0.4 0.2 0];
Ncellvals = [5 10 20 50 100 500];

%set booleans specifying homogeneous tuning
homog = 1;
heterog = 0;

for rrr = 1:length(rhomaxvals) %for many different values of rhomax
    for NNN = 1: length(Ncellvals) %and many different population sizes
        rhomax = rhomaxvals(rrr);
        Ncells = Ncellvals(NNN);
        do_FI_calc_STIMDEP; %compute the Fisher info
        FI_stimdep_saver(rrr,NNN) = meanFI_stimdep; %save the calculated values
        FI_const_saver(rrr,NNN) = meanFI_const;
        close all
        meancorrval(rrr,NNN) = mean(meancorrel);  %save the mean correlations -- to use in Fig. legend

    end
end

%make the plots
figure()
set(gca,'fontsize',16)
plot(Ncellvals,FI_stimdep_saver','linewidth',2)
legend('SD, \rho_{max}=0.8, <\rho> = 0.11', 'SD, \rho_{max}=0.6, <\rho> = 0.081','SD, \rho_{max}=0.4, <\rho> = 0.054','SD, \rho_{max}=0.2, <\rho> = 0.027','Uncorrelated')
xlabel('Number of Cells')
ylabel('Fisher Information')


figure()
set(gca,'fontsize',16)
plot(Ncellvals,transpose((FI_const_saver')),'linewidth',2)
legend('CONST, \rho_{max}=0.8','CONST, \rho_{max}=0.6','CONST, \rho_{max}=0.4','CONST, \rho_{max}=0.2','Uncorrelated')
xlabel('Number of Cells')
ylabel('Fisher Information')