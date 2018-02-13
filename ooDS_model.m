%simulates the ooDS cell model from Zylberberg, Cafaro, Turner, et al., Neuron 2016
%model population contains 8 cells: 2 each of the 4 sub-types
%this code loads params from ooDS_params

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin model set-up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first, load up the model parameters -- fitted by the Genetic Algo
%get the params from the GA
load ooDS_params
Ncells = 8;
 
%baseline value for threshold in the threshold-linear neuron model
basethresh = GAparams(1); 
%cell-by-cell variance in the thresholds
threshstd = 0.1*basethresh;
thresh = basethresh + threshstd*randn(Ncells,1);

%magnitudes of constant upstream E and I "signals"
Ein = GAparams(2);
Iin = GAparams(3);

%mean and std of the gE value -- untuned excitatory gain
gebase = GAparams(4);
gestd = GAparams(5);

gE = (gebase + gestd*randn(Ncells,1))*0.67; 
gE(gE<0) = 0;


%inhibition is tuned: several parameters. gI = gI_max*(((1-sin(angle - TC_centers(acell) - pi/2))/2)^alpha) + GI_offset;
%mean and std of gImax: amplitude of tuned part
gimaxbase = GAparams(6);
gimaxstd = GAparams(7);

gI_max = gimaxbase + gimaxstd*randn(Ncells,1); %5.8
gI_max(gI_max<0) = 0;

    
%mean and std of gioff: baseline "offset" of gI
gioffbase = GAparams(8);
gioffstd = GAparams(9);

GI_offset = gioffbase + gioffstd*randn(Ncells,1); %1.8


%gamma: multiplier on E varianc relative to I
gamma = GAparams(10);

%sig2c: variance of common "upstream" noise: shared btw. cells
sig2c = GAparams(11);

%variance of E and I noise that is downstream of multiplicative gain stage
%this noise is independent between cels and channels
post_noise_var_E = GAparams(12);
post_noise_var_I = GAparams(13);
    
%this is actually the slope of the thresh-lin model.
TCamp = GAparams(14);
slope = TCamp; %for thresh-linear..s
 
%minimum "beta": overlap between cells in the upstream noise sources.
minbeta = GAparams(15);
    
%alpha is sharpness of the stimulus tuning of inhibition
%baseline value is the mean.
base_alpha = GAparams(16);
alpha = base_alpha + 0.5*randn(Ncells,1);
alpha(alpha<0) = 0;

%K is the different multiplier on E vs I inputs. This comes about because
%of different driving force through E vs I channels in the cells.
k = 3; 

%These are the locations (angles) of the tuning curve peaks
TC_centers = [0 0 pi/2 pi/2 pi pi 3*pi/2 3*pi/2] ; 
TC_centers = TC_centers + 0.2*randn(size(TC_centers));

%RF overlap (as a fraction...). This sets the max fraction of shared
%upstream noise
RF_overlap = 0.95;

while(1)
    BetaArray = RF_overlap*ones(Ncells,Ncells) - minbeta*rand(Ncells,Ncells);
    BetaArray(BetaArray>1) = 1; %no beta value can exceed 1
    
    BetaArray = triu(BetaArray) + transpose(triu(BetaArray)); %beta must be symmetric;
    BetaArray = eye(size(BetaArray)) + BetaArray - diag(diag(BetaArray)); %diags are 1

    tester = sqrtm(BetaArray);

    %check that this matrix has a square root: if not, nonphysical overlaps...
    if(isreal(tester))
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of model set-up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simulate the model  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for repeats of each stim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of repeats of each stim
ntrials = 1000;

%number of stimulus values to get responses for
nangles = 100;

%list of stimulus to consider
anglelist = linspace(0,2*pi,nangles);

%correlation matrix of the upstream noise going to the cells
in_corr_mat = BetaArray;

%"mixing matrix" to generate correlated inputs from independent random
%numbers
input_mix = sqrtm(in_corr_mat);


%simulate the model for all stim values
for i=1:nangles %each index is a new stimulus
  
    %define the stim
    angle = anglelist(i);
    
    %generate the vectors of (correlated) random numbers
    input_noise = input_mix*randn(Ncells,ntrials); %this is the "upstream" noise
    
    %compute the Excitatory inputs to the cells. 
    %includes (signal + upstream noise)*gain
    %plus the downstream noise (independent between cells)
    E = diag(gE)*(Ein + gamma*sqrt(sig2c)*input_noise) + post_noise_var_E*randn(Ncells,ntrials);
    E(E<0) = 0; %truncate E values so non-negative
 
  
    %make a matrix of stim-dep gains for the Inhibitory channel
    for acell = 1:Ncells
        stim_dep_gain(acell) = gI_max(acell)*(((1-sin(angle - TC_centers(acell) - pi/2))/2)^alpha(acell)) + GI_offset(acell);
    end
    
    %then generate the sets of I inputs to the cells
    %includes (signal + upstream noise)*gain
    %plus the downstream noise (independent between cells)
    I = (transpose(repmat(stim_dep_gain,ntrials,1))).*(Iin + sqrt(sig2c)*input_noise) + post_noise_var_I*randn(Ncells,ntrials);
    I(I<0) = 0; %truncate I values so non-negative
        
    %inputs to the cells are KE-I on each trial
    theinput = k*E - I;
    
    %now use some awkward code to implement the treshold-linear neural
    %model
    aresp = theinput;
    threshmat = repmat(thresh,1,ntrials);
    %anything below threshold is set to zero    
    aresp(aresp<= threshmat) = 0;
    %anything above threshold is multiplied by the slop of the
    %thresh-linear function
    aresp(theinput>threshmat) = slope*(aresp(theinput>threshmat) - threshmat(theinput>threshmat));    
    
    %add a small amount of random noise to the outputs
    aresp = aresp + randn(size(aresp)); 
    aresp(aresp<=0) = 0; %truncate so non-negative
    
    %now compute some statistics for the responses
    cov_for_this_stim = cov(aresp'); %cov matrix
    corr_for_this_stim = corrcoef(aresp'); %correlation matrix
    mean_for_this_stim = mean(aresp'); %mean responses
    
    %save the responses and associated statistics
    response(:,:,i) = aresp; %save all the responses
    stimangle(:,:,i) = repmat(angle,size(aresp)); %save the stim variables (for decoders)
   
    means(:,i) = mean_for_this_stim;
    covariances(:,:,i) = cov_for_this_stim;
    correlations(:,:,i) = corr_for_this_stim;
    
    correlation_saver(i,:,:) = corr_for_this_stim;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end simulate the model  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
      