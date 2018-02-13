%Compute the Fisher information for many different population sizes
%And strengths of correlations

%Makes Fig. 7EF of 
%Zylberberg, Cafaro, Turner, et al. Neuron 2016
%Direction Selective Circuits Shape Noise to Ensure a Precise Population Code
%contact: joel.zylberberg@ucdenver.edu to report bugs or issues

%Heterogeneous tuning curves

%pick the parameter range to consider
rhomaxvals = [0.8 0.6 0.4 0.2 0];
Ncellvals = [5 10 50 100 500];

%set booleans specifying heterogenous tuning
homog = 0;
heterog = 1;

%repeat over several random sets of tuning curves
for repeater = 1:20
repeater
for rrr = 1:length(rhomaxvals) %for many different values of rhomax
    for NNN = 1: length(Ncellvals) %and many different population sizes
        rhomax = rhomaxvals(rrr);
        Ncells = Ncellvals(NNN);
        do_FI_calc_STIMDEP; %compute the Fisher info
        FI_stimdep_saver(rrr,NNN) = meanFI_stimdep; %save the calculated values
        FI_const_saver(rrr,NNN) = meanFI_const;
        
        meancorrval(rrr,NNN) = mean(meancorrel); %save the mean correlations -- to use in Fig. legend
        close all
    end
end

%record the results for this "run" of random tuning curves
mFI_stimdep(repeater,:,:) = (FI_stimdep_saver); 
mFI_const(repeater,:,:)= (FI_const_saver);
meanrho(repeater,:,:) = meancorrval;

end

%get average and SEM over the runs
FI_stimdep_saver = reshape(mean(mFI_stimdep),rrr,NNN);
FI_const_saver = reshape(mean(mFI_const),rrr,NNN);

FI_stimdep_ERR = reshape(std(mFI_stimdep),rrr,NNN)/sqrt(repeater);
FI_const_ERR = reshape(std(mFI_const),rrr,NNN)/sqrt(repeater);


%make the plots
figure()
set(gca,'fontsize',16)
errorbar(transpose(repmat(Ncellvals,rrr,1)),FI_stimdep_saver',FI_stimdep_ERR','linewidth',2)
hold on
legend('SD, \rho_{max}=0.8, <\rho> ~ 0.4','SD, \rho_{max}=0.6, <\rho> ~ 0.3','SD, \rho_{max}=0.4, <\rho> ~ 0.2','SD, \rho_{max}=0.2, <\rho> ~ 0.1','Uncorrelated','location','NorthWest')
xlabel('Number of Cells')
ylabel('Fisher Information (rad.^{-2})')


figure()
set(gca,'fontsize',16)
errorbar(transpose(repmat(Ncellvals,rrr,1)),FI_const_saver',FI_const_ERR','-.','linewidth',2)
hold on
legend('CONST, \rho_{max}=0.8','CONST, \rho_{max}=0.6','CONST, \rho_{max}=0.4','CONST, \rho_{max}=0.2', 'Uncorrelated','location','NorthWest')
xlabel('Number of Cells')
ylabel('Fisher Information (rad.^{-2})')




