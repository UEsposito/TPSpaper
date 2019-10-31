
function [predicted_time,predicted_time_STD,gen_distance_closest_ref_pop] = Run_TPS(exp_txt_GEN,exp_GEN,exp_TEM,exp_TEM_STD,ad_coef_test,displaySwitch)

% Predict dates

%% Parameters
%Parameters for the predictive algorithm
PredMethod = 'InvGenDist';
n_used_ref_pop = 3;

%% Variables and parameters
n_test_ID = size(ad_coef_test,1);
closest_ref_pop = cell(n_test_ID,1);
ind_used_ref_pop = zeros(n_test_ID,n_used_ref_pop);
predicted_time = zeros(n_test_ID, 1);
predicted_time_STD = zeros(n_test_ID, 1);
gen_distance_closest_ref_pop = zeros(n_test_ID, 1);
W1_all = zeros(n_test_ID, 1);
W2_all = zeros(n_test_ID, 1);
W3_all = zeros(n_test_ID, 1);

%% Algorithm
for i=1:n_test_ID
    
    %Print progress on screen
    if strcmp(displaySwitch,'On')
        fprintf('\t\tAnalysing individual #%d/%d\n',i,n_test_ID);
    end
    
    %Calculate genetic distances
    dist_matrix = dist([ad_coef_test(i,:);exp_GEN]');
    genetic_distances = dist_matrix(1,2:end);
    
    %Sort genetic distances
    [~,sorting_ind] = sort(genetic_distances,'ascend');

    %Store indices of closest reference populations
    n_used_ref_pop_temp = min(n_used_ref_pop,size(sorting_ind,2));
    ind_used_ref_pop = sorting_ind(1:n_used_ref_pop_temp);
    
    %Extract their ID
    IDs_used_ref_pop = exp_txt_GEN(ind_used_ref_pop,:);
    
    %Extract their genetic distance with the current test individual
    gen_distance_used_ref_pop = genetic_distances(ind_used_ref_pop);
    gen_distance_closest_ref_pop(i) = gen_distance_used_ref_pop(1);
    
    %Extract their time period coordinates
    time_period_used_ref_pop = exp_TEM(ind_used_ref_pop,:);
    time_period_STD_used_ref_pop = exp_TEM_STD(ind_used_ref_pop,:);
    time_period_used_ref_pop_unique = unique(time_period_used_ref_pop,'rows');
    
    %Extract closest reference population
    closest_ref_pop{i} = IDs_used_ref_pop{1};
    
    %Calculate predictions
    if gen_distance_used_ref_pop(1) == 0
        predicted_time(i,:) = time_period_used_ref_pop(1,:);
        predicted_time_STD(i,:) = time_period_STD_used_ref_pop(1,:);
    else
        if size(time_period_used_ref_pop_unique, 1) == 1
            predicted_time(i,:) = time_period_used_ref_pop_unique;
            predicted_time_STD(i,:) = time_period_STD_used_ref_pop(1,:);
        else
            
            %---Weighting predicting algorithm
            if strcmp(PredMethod,'80w') %Give 80% to the closest reference population and 20% to the second
                %Define the weight W to assign to the closest population            
                W = 0.8;
                predicted_time(i,:) = time_period_used_ref_pop(1,:) .* W + time_period_used_ref_pop(2,:) .* (1-W);
                predicted_time_STD(i,:) = sqrt( (W * time_period_STD_used_ref_pop(1,:))^2 + ((1-W) * time_period_STD_used_ref_pop(2,:))^2 );
                
            elseif strcmp(PredMethod,'80w10w10w') %Give 50% to the closest reference population, 30% to the second and 20% to the third
                %Define the weight W to assign to the closest population
                W1 = 0.8;
                W2 = 0.1;
                W3 = 0.1;
                predicted_time(i,:) = time_period_used_ref_pop(1,:) .* W1 + time_period_used_ref_pop(2,:) .* W2 + time_period_used_ref_pop(3,:) .* W3;
                predicted_time_STD(i,:) = sqrt( (W1 * time_period_STD_used_ref_pop(1,:))^2 + (W2 * time_period_STD_used_ref_pop(2,:))^2 +(W3 * time_period_STD_used_ref_pop(3,:))^2 );
            
            elseif strcmp(PredMethod,'InvGenDist')
                W1 = 1/gen_distance_used_ref_pop(1);
                W2 = 1/gen_distance_used_ref_pop(2);
                W3 = 1/gen_distance_used_ref_pop(3);
                norm = W1 + W2 + W3;
                W1_all(i) = W1/norm;
                W2_all(i) = W2/norm;
                W3_all(i) = W3/norm;
                predicted_time(i,:) = time_period_used_ref_pop(1,:) .* W1_all(i) + time_period_used_ref_pop(2,:) .* W2_all(i) + time_period_used_ref_pop(3,:) .* W3_all(i);
                predicted_time_STD(i,:) = sqrt( (W1_all(i) * time_period_STD_used_ref_pop(1,:))^2 + (W2_all(i) * time_period_STD_used_ref_pop(2,:))^2 +(W3_all(i) * time_period_STD_used_ref_pop(3,:))^2 );
            end
        end
    end
end

