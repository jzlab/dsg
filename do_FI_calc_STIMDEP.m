%Compute the Fisher information
%For stimulus dependent correlations, and for "matched" constant correlations
%Tuning curves can be either homogeneous or heterogeneous

%Used in Zylberberg, Cafaro, Turner, et al. Neuron 2016
%Direction Selective Circuits Shape Noise to Ensure a Precise Population Code
%contact: joel.zylberberg@ucdenver.edu to report bugs or issues

%note: Ncells needs to be set -- it is in the "looper" scripts

nangles = 50; %at which to eval (then avg) FI
dtheta = 0.00001;  %range for doing local linear est
plotting = 0; %bolean to suppress plot output of tuning curves

%%%%%%%parameters for the TCs
if(homog)
alphas = 0*ones(Ncells,1);  
gammas = 3*ones(Ncells,1);
betas = 25*ones(Ncells,1);
phis = linspace(0,2*pi,Ncells + 1);
phis = transpose(phis(1:end-1));
end

if(heterog)
alphas = 5*rand(Ncells,1);  
gammas = 3*rand(Ncells,1);
betas = 25*rand(Ncells,1);
phis = 2*pi*rand(Ncells,1);
end


%make the tuning curve plots
angles = linspace(0,2*pi,nangles);
clear TCcalc
for stimvals = 1:nangles
    sval = angles(stimvals);
    TCcalc(stimvals,:) = TC(sval,alphas,betas,gammas,phis);
end

%show the tuning curves...
if(plotting)
figure()
set(gca,'fontsize',16)
plot(angles*180/pi, TCcalc,'linewidth',2)
xlabel('Stimulus')
ylabel('Mean Neural Response')
end

%make the max product matrix
maxprodmatrix = transpose(max(TCcalc))*(max(TCcalc));
clear corrmatrices

%loop over all stimuli, getting FI each time
stimangles = angles;
for stimvals = 1:nangles
    
    %stimulus value
    sval = stimangles(stimvals);
    
    %compute mean inputs -- scale them to max 25
    mean_inputs = TC(sval,alphas,betas,gammas,phis)';
    mean_inputs2 = TC(sval + dtheta,alphas,betas,gammas,phis)'; 

    %put in the correlation matrix -- the off diagonal elements are Eq. 2
    %diagonals of the correlation matrix are letter set to 1
    correlmat_offdiag = real(rhomax*sqrt(mean_inputs'*mean_inputs)./sqrt(maxprodmatrix));
    hollowcorrel = correlmat_offdiag - diag(diag(correlmat_offdiag));
    
	%covariance matrix
	covmat = sqrt(diag(mean_inputs))*(eye(Ncells,Ncells) + hollowcorrel)*sqrt(diag(mean_inputs));
    
    %save the correlations for later
    corrmatrices(stimvals,:,:) = corrcov(covmat);
    cmmm = corrcov(covmat);
    
    %compute the FI
    input_derivs = (mean_inputs - mean_inputs2)/dtheta;
    FI_value(stimvals) = input_derivs * inv(covmat) * input_derivs';
      
    %save the mean correlation over all cells, to use in figure legends
    meancorrel(stimvals) = mean(cmmm(cmmm<1)); %the <1 flag means you don't take the diagonals, which are always 1, and aren't cross-correlations
    
end

%compute the stim averaged correlation matrices
avcorr = reshape(mean((corrmatrices)),Ncells,Ncells);

%can display some plots here if desired.
if(plotting)
figure()
set(gca,'fontsize',16)
TC_diffs_mat = bsxfun(@minus,phis',phis);
TC_diffs_mat(TC_diffs_mat > pi) = 2*pi - TC_diffs_mat(TC_diffs_mat > pi);


plot(TC_diffs_mat(1,2:end),avcorr(1,2:end),'linewidth',2)
xlabel('Difference in Preferred Stimulus')
ylabel('Mean Correlation Coefficient')
end

%compute the mean FI for the stim dep correlations
meanFI_stimdep = mean(FI_value);


%%%%%%%% now repeat for constant correlations
for stimvals = 1:nangles
    
    %stimulus value
    sval = stimangles(stimvals);
    
    %compute mean inputs -- scale them to max 25
    mean_inputs = TC(sval,alphas,betas,gammas,phis)';
    mean_inputs2 = TC(sval + dtheta,alphas,betas,gammas,phis)';  

    %put in the correlation matrix!
    %the diagonal is later put in as ones.
    correlmat_off_diag = avcorr;
    hollowcorrel = correlmat_off_diag - diag(diag(correlmat_off_diag));
    
	%covariance matrix of inputs -- put in the stim dep means
	covmat = sqrt(diag(mean_inputs))*(eye(Ncells,Ncells) + hollowcorrel)*sqrt(diag(mean_inputs)); 

    %compute the FI
    input_derivs = (mean_inputs - mean_inputs2)/dtheta;
    FI_value(stimvals) = input_derivs * inv(covmat) * input_derivs';
   

    
end

%save the FI
meanFI_const = mean(FI_value);