%ooDS Fisher code
%calculates Fisher information from the ooDS cell responses, as in Zylberberg, Cafaro, Turner, et al., Neuron 2016
%this code calls the ooDS_model

%this is the number of times to repeat the analysis (the models differ a
%bit each time due to random cell parameters)
numruns = 10;


for run= 1:numruns
    %first, run the simulation of the model
    ooDS_model
    
    %then compute the Fisher information with the model outputs
    %Fisher info is defined locally, for each stim value
    %loop over all the stim values and average the Fisher Info
    for jjj=1:length(anglelist)-1
        
        %get mean neural for this stim, and for the next one "up" -- use
        %these for finite different derivative
        a_bunch_of_means = means(:,jjj);
        a_bunch_of_means2 = means(:,jjj+1);
        derivs = (a_bunch_of_means2 - a_bunch_of_means)/(anglelist(jjj + 1) - anglelist(jjj));
        
        %extract the covariance matrix for these responses
        covmat = reshape(covariances(:,:,jjj),Ncells,Ncells);
        
        %and also get the covariance that has the "matched" constant
        %correlations, but with the stim-dependent variances
        meancor = squeeze(mean(correlations,3));
        cov_const_corr = sqrt(diag(diag(covmat)))*meancor*sqrt(diag(diag(covmat)));
        
        %remove any non-firing cells that giver zero information
        %this step can probably be removed without any problems
        slicer = find(abs(derivs)>1E-6);
        
        %compute the fisher info
        fisher_save(jjj) = transpose(derivs(slicer))*inv(covmat(slicer,slicer))*derivs(slicer);
        
        %and for the shuffle data (no correlations)
        fisher_save_shuff(jjj) = transpose(derivs(slicer))*inv(diag(diag(covmat(slicer,slicer))))*derivs(slicer);
    
        %and for the "matched" constant correlations
        fisher_save_constcorr(jjj) = transpose(derivs(slicer))*inv(cov_const_corr(slicer,slicer))*derivs(slicer);
    end
    
    %save the (mean) Fisher info for this model population
    fisher_info(run) = mean(fisher_save);
    fisher_info_shuff(run) = mean(fisher_save_shuff);
    fisher_info_constcorr(run) = mean(fisher_save_constcorr);
end


%average the result over the set of model populations
FISH_INFO = mean(fisher_info)
FISH_INFO_SHUFF = mean(fisher_info_shuff) %for the trial shuffled
FISH_INFO_CONSTCORR = mean(fisher_info_constcorr)

beep