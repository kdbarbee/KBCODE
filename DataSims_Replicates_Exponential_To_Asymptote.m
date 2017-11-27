clc
clear all
close all

rep_num = 2; % INPUT: Define the number of replicates to be generated per unique inputs defined in the following xvals vector
xvals = [5,15,25,35]; % INPUT: Defines the number of unique inputs as well as the range
sz_xvals = size(xvals); % Get the size of the xvals vector defined above

rep_vect = ones(1,rep_num); % Establish vector to represent the number of replicates for each unique input value
xrep = []; % Create empty vector that will be populated with all x-inputs and their replicates, if applicable
% Loop to generate 
for rv = 1:sz_xvals(1,2) % Loop index shall go up to number of unique inputs in xvals vector
    xrepv = xvals(1,rv)*rep_vect; % Temporary xrepv vector is the product of the scalar element of the xvals input vector and the replicates vector 
    xrep = [xrep xrepv]; % Append the temporary xrepv vector generated for each xvals element to the xrep vector. 
    % This xrep vector will grow by rep_num elements for each iteration in this loop
end

upr_range = max(xrep); % Get upper limit of xrep vector to help create the 'singles' data set
lwr_range = min(xrep); % Get lower limit of xrep vector to help create the 'singles' data set
sz_xrep = size(xrep); % Get total number of data points from xrep data set in matrix format
data_points = sz_xrep(1,2); % Get total number of data points from xrep data set in scalar format
scale_factor = 0.5; % INPUT: Set scale factor for y-axis bounds on plots. (Multiply '1+/-scale_factor' by ysin_min and ysin_max from model)

xsin_spacing = (upr_range - lwr_range)/(data_points-1); % Calculate the spacing to use between elements in the singles data set
xsin_raw = lwr_range:xsin_spacing:upr_range; % Define x values for non-replicate data set (singles). These may not be integer values
xsin = round(xsin_raw); % Define x values for non-replicate data set (singles), integer values (normal rounding)


A = 10000; % INPUT: Define initial value for model exponential fnc
k = -0.10; % INPUT: Define exp constant for model exponential fnc
params_fnc = [A;k]; % Store in column vector for later use as starting points in fminsearch
params_fnc_noise_bandwidth = 0.9; % INPUT: Define noise band to be applied to noise on starting values of fminsearch
upr2 = 1 + params_fnc_noise_bandwidth; % Upper limit on noise band to be applied to model parameters in creation of fminsearch starting values
lwr2 = 1 - params_fnc_noise_bandwidth; % Lower limit on noise band to be applied to model parameters in creation of fminsearch starting values

ysin = A*(1-exp(k*xsin)); % Define model exponential fnc for non-replicate data set
yrep = A*(1-exp(k*xrep)); % Define model exponential fnc for replicate data set

ysin_min = min(ysin);
ysin_max = max(ysin);

xfin = round(linspace(lwr_range,upr_range,100)); % Create fine x-data vector of integers for smooth curve of model exponential
yfin = A*(1-exp(k*xfin)); % Define model exponential fnc for fine x data set

opt_param1 = 10000; % INPUT: optimset param1
opt_param2 = 10000; % INPUT: optimset param2

% % % % % % % % % % % % % % % % % % % % 
sims = 25;  % INPUT: Set number of iterations of secondary loop
nb_int = 5; % INPUT: Set number of intervals for noise band
nb_div = 10; % INPUT: Set divisor for noise band calculations
trials=100; % INPUT: Define number of trials, where each trial generates a new noisy data set for subsequent curve fitting
% % % % % % % % % % % % % % % % % % % % 

% Set up matrices for use in main loops:
fits_error_summary = zeros(sims, 4); % Set up matrix for compiling errors from each sim 
ratio_errors_summary = zeros(sims,2); % Set up matrix for compiling ratios of errors from each sim 
% Ratios less than 1 imply that error for singles sets is greater than that from replicates set
ratio_means = zeros(nb_int, 4); % Set up matrix for compiling means of ratios generated from each sim and for each noise band
ratio_sd = zeros(nb_int, 4); % Set up matrix for compiling std dev of ratios generated from each sim and for each noise band
% ratios_CV = zeros(nb_int, 6);


for nb=1:nb_int;
    figure % Create new figure for each noise band iteration
    for p=1:sims;
        noise_band = nb/nb_div; % Define noise band
        lwr = 1 - noise_band; % Define lower bound for noise vector
        upr = 1 + noise_band; % Define upper bound for noise vector
        fits=zeros(trials,4); % Define matrix for storing fit parameters generated in the following loop
        ysin_noisy_matrix = zeros(trials,numel(ysin));
        yrep_noisy_matrix = zeros(trials,numel(yrep));

        for i=1:trials; % Begin loop for generating noisy data and curve fitting to each noisy data set
            
            noise_fin = (upr - lwr)*rand(size(xfin)) + lwr; % Define noise vector using bounds set above
            noise = (upr - lwr)*rand(size(xsin)) + lwr; % Define noise vector using bounds set above
            noise_range = [min(noise),max(noise)]; % Check that noise vector is within expected range
            
            ysin_noisy=ysin.*noise; % Apply noise to non-replicate data (use same noise vector on singles and replicates)
            ysin_noisy_matrix(i,:) = ysin_noisy; % Add noisy data to main noisy data matrix for later calculation of min, max and mean

            yrep_noisy=yrep.*noise; % Apply noise to replicate data (use same noise vector on singles and replicates)
            yrep_noisy_matrix(i,:) = yrep_noisy; % Add noisy data to main noisy data matrix for later calculation of min, max and mean

            noise_fmin = (upr2 - lwr2)*rand(size(params_fnc)) + lwr2; % Create noise column vector with bounds set by lwr2 and upr2
            noisy_fmin_start = params_fnc.*noise_fmin; % Create noisy starting points for fminsearch, but base them on model parameters

            % Fitting data to an exponential:
            yrep_obj = @(b,x) b(1)*(1-exp(b(2)*x));             % Objective function for replicate data
            OLS = @(b) sum((yrep_obj(b,xrep) - yrep_noisy).^2);          % Ordinary Least Squares cost function
            opts = optimset('MaxFunEvals',opt_param1, 'MaxIter',opt_param2);
            B = fminsearch(OLS, noisy_fmin_start, opts);       % Use ‘fminsearch’ to minimise the ‘OLS’ function

            fits(i,1) = B(1,1); % Store replicates fit param 'A' in col 1 of fits matrix
            fits(i,2) = B(2,1); % Store replicates fit param 'k' in clo 2 of fits matrix

            ysin_obj = @(c,x) c(1)*(1-exp(c(2)*x));             % Objective function for single data points
            OLS2 = @(c) sum((ysin_obj(c,xsin) - ysin_noisy).^2);          % Ordinary Least Squares cost function
            opts = optimset('MaxFunEvals',opt_param1, 'MaxIter',opt_param2);
            B2 = fminsearch(OLS2, noisy_fmin_start, opts);       % Use ‘fminsearch’ to minimise the ‘OLS’ function
            fits(i,3) = B2(1,1); % Store singles fit param 'A' in col 3 of fits matrix
            fits(i,4) = B2(2,1); % Store singles fit param 'k' in col 4 of fits matrix
            
        end
        ysin_noise_max = max(ysin_noisy_matrix);   % Calculate max of noisy data from all trials in each sim
	    ysin_noise_min = min(ysin_noisy_matrix);   % Calculate min of noisy data from all trials in each sim
        ysin_noise_mean = mean(ysin_noisy_matrix); % Calculate mean of noisy data from all trials in each sim
        
        yrep_noise_max = max(yrep_noisy_matrix);   % Calculate max of noisy data from all trials in each sim
	    yrep_noise_min = min(yrep_noisy_matrix);   % Calculate min of noisy data from all trials in each sim
        yrep_noise_mean = mean(yrep_noisy_matrix); % Calculate mean of noisy data from all trials in each sim

        plot1 = plot(xrep,yrep_noise_min,  '*b');  % Plot min of noisy data from all trials in each sim
        hold on
        plot(xrep, yrep_noise_max, '*b');          % Plot max of noisy data from all trials in each sim
        hold on
        plot(xrep, yrep_noise_mean, '*b');         % Plot mean of noisy data from all trials in each sim
        hold on
        plot2 = plot(xsin, ysin_noise_min,'+m'); 
        hold on
        plot(xsin, ysin_noise_max, '+m'); 
        hold on
        plot(xsin, ysin_noise_mean, '+m'); 
        hold on
        
        fits_avg = mean(fits); % Calculate averages of each fit parameter for both replicates and singles data sets for the above trial
        Arep_fit = fits_avg(1,1); % Extract average of fit parameter 'A' for replicates data set for the above trial
        krep_fit = fits_avg(1,2); % Extract average of fit parameter 'k' for replicates data set for the above trial
        Asin_fit = fits_avg(1,3); % Extract average of fit parameter 'A' for singles data set for the above trial
        ksin_fit = fits_avg(1,4);% Extract average of fit parameter 'k' for singles data set for the above trial
        
        yrep_fit = Arep_fit.*(1-exp(krep_fit.*xfin)); % Generate y values for exponential fnc from fit parameters (replicate)
        ysin_fit = Asin_fit.*(1-exp(ksin_fit.*xfin)); % Generate y values for exponential fnc from fit parameters (singles)
               
        plot3 = plot(xfin, yrep_fit, '--b', 'LineWidth',2); % Add a fitted curve using avg fit params from all trials in a sim
        hold on
        plot4 = plot(xfin, ysin_fit, ':m','LineWidth',2);
        hold on
        
        fits_error = zeros(size(fits_avg)); % Define vector for error of averages of fit parameters
        fits_error(1,1) = abs(((fits_avg(1,1)-A)/A)*100); % Calculate error of fit for A parameter for replicate data
        fits_error(1,2) = abs(((fits_avg(1,2)-k)/k)*100); % Calculate error of fit for k parameter for replicate data
        fits_error(1,3) = abs(((fits_avg(1,3)-A)/A)*100); % Calculate error of fit for A parameter for singles data
        fits_error(1,4) = abs(((fits_avg(1,4)-k)/k)*100); % Calculate error of fit for A parameter for singles data
        ratio_errors(1,1) = fits_error(1,1)/fits_error(1,3); % Calculate ratio of error of fit for replicate(A) to single(A)
        ratio_errors(1,2) = fits_error(1,2)/fits_error(1,4); % Calculate ratio of error of fit for replicate(k) to single(k)
        fits_error_summary(p,:)= fits_error;
        ratio_errors_summary(p,:)= ratio_errors;
      
    end
    
    plot5 = plot(xfin,yfin,...
            '-r',...
            'LineWidth',2); % Plot fine model data
    hold on
	plot(xfin,yfin,...
            'or',...
            'LineWidth',2); % Plot fine model data

    hold on
    lgd = legend([plot1 plot2 plot3 plot4 plot5],...
            {'Noisy replicate data (min, mean, max)', 'Noisy singles data (min, mean, max)',...
            'Fitted curve (replicates)','Fitted curve (singles)', 'Model curve'},'FontSize', 11, 'Location','northwest');
    
    axis([lwr_range upr_range (1-scale_factor)*ysin_min (1+scale_factor)*ysin_max])
	hold off
    xlabel('Concentration','FontSize', 12)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 11)
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 11)
    ylabel('AFU','FontSize', 12)
    lbl_lwr = num2str(round(100*lwr)/100);
    lbl_upr = num2str(round(100*upr)/100);
    plot_title = title({'Impact of Replicates on Curve Fitting'; ['Noise Band = ', lbl_lwr, ' to ', lbl_upr, ', Replicates = ', num2str(rep_num)]},'FontSize', 12);
    grid
    
    
    
    ratio_means(nb,:) = [lwr, upr, mean(ratio_errors_summary)];
    ratio_sd(nb,:) = [lwr, upr, std(ratio_errors_summary)];
%   ratios_CV(nb,:) = [lwr, upr, 100*ratio_sd./ratio_means];
end

