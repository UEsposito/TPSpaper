
function [testID_subPop_name,exp_txt_GEN,exp_txt_TEM,exp_GEN,exp_TEM,exp_TEM_STD] = Build_reference_panel(Dataset_training,testID_adCoef,displaySwitch)

% Generate a reference panel for TPS

%% Parameters
%Number of repetitions for kmeans algorithm
n_rep_kmeans = 100;

%Method for calculating reference panel populations
populations_method_setting = 'Average time';

%Parameters controlling the sizes of the populations and sub-populations allowed in the reference panel
pop_min_size = 2;
pop_min_size_for_subpop = 3;
subpop_min_size = 2;
subpop_min_number = 1;
f_max_k = @(x) round(x/2);

%Parameters used to split populations into sub-populations within the reference panel
clusteringAlg_iterNum = 2;
f_k_to_use = @(x) x;

%% Variables
testID_supPop_name = Dataset_training.Population;
testID_time = Dataset_training.DateBP;
testID_timeSTD = Dataset_training.DeltaT95CI;

counter1 = 1;
alphabet = {'a';'b';'c';'d';'e';'f';'g';'h';'i';'j';'k';'l';'m';'n';'o';'p';'q';'r';'s';'t';'u';'v';'w';'x';'y';'z'};

%% Algorithm Part 1: Identify the reference populations - using name and time
if strcmp(displaySwitch,'On')
    fprintf('\n\tRunning algorithm part1: Identifying reference populations...\n');
end
[unique_supPop_names, ~, ~] = unique(testID_supPop_name,'stable');
n_supPops = size(unique_supPop_names,1);
testID_subPop_name = cell(length(testID_supPop_name),1);

exp_txt_GEN_pop = string(repmat({' '},n_supPops,1));
exp_txt_TEM_pop = string(repmat({' '},n_supPops,1));
test_pop_membership = string(repmat({' '},size(testID_supPop_name,1),1));

for i=1:n_supPops

    %Print progress on screen
    if strcmp(displaySwitch,'On')
        fprintf('\t\tAnalysing superpopulation''s name #%d/%d\n',i,n_supPops);
    end

    matching_test_supPop_logical = ismember(testID_supPop_name,unique_supPop_names(i));
    
    %Each unique population name is a population in the reference panel
    exp_txt_GEN_pop(counter1,1) = unique_supPop_names(i);
    exp_txt_TEM_pop(counter1,1) = unique_supPop_names(i);
    test_pop_membership(matching_test_supPop_logical) = unique_supPop_names(i);
    counter1 = counter1 + 1;
end
pop_panel_size = counter1 - 1;

%% Algorithm Part 2: Divide each population in several subpopulations whenever possible - using admixture components
for ni = 1:clusteringAlg_iterNum
    
    if ni==1
        panel_size = pop_panel_size;
    elseif ni==2
        if strcmp(displaySwitch,'On')
            fprintf('\n\tStarting nested subpopulations analysis...\n');
        end
        panel_size = subpop_panel_size;
        test_pop_membership = testID_subPop_name;   %contains for each individual the name of the reference population to which he belongs
        exp_txt_GEN_pop = exp_txt_GEN;
    end
    
    exp_txt_GEN = string(repmat({' '},n_supPops,1));
    exp_txt_TEM = string(repmat({' '},n_supPops,1));
    exp_GEN = [];
    exp_TEM = [];
    exp_TEM_STD = [];
    size_subpop = [];
    
    counter2 = 1;
    
    %testID_subPop_name = cell(length(testID_supPop_name),1);
    testID_subPop_name = string(repmat({' '},size(testID_subPop_name)));
    
    for i=1:panel_size

        t_iter = tic;

        %Print progress on screen
        if strcmp(displaySwitch,'On')
            if ni==1
                fprintf('\t\tEvaluating reference population #%d/%d\n',i,panel_size);
            elseif ni==2
                fprintf('\t\tEvaluating reference subpopulation #%d/%d\n',i,panel_size);
            end
        end

        %Extract information about current reference population
        matching_test_pop_logical = ismember(test_pop_membership, exp_txt_GEN_pop(i));
        ad_coef_matching_test = testID_adCoef(matching_test_pop_logical,:);
        time_matching_test = testID_time(matching_test_pop_logical,:);
        timeSTD_matching_test = testID_timeSTD(matching_test_pop_logical,:);
        ad_coeff_subpop = [];

        switch populations_method_setting
            case 'Unique time'
                unique_time_matching_test = unique(time_matching_test,'rows');    
                if size(unique_time_matching_test,1) ~= 1
                    error('More than one time for the current reference population')
                end
        end
        
        %Determine the reference subpopulations
        if sum(matching_test_pop_logical) >= pop_min_size
            %A subpopulation or population can be part of the reference panel only if it has a minimum number of individuals
            if sum(matching_test_pop_logical) < pop_min_size_for_subpop
                %---Population size does not allow for sub-divisions ---> There is only one reference subpopulation
                if ni==1
                    ID_subpop = strcat(exp_txt_GEN_pop{i},'_0');
                elseif ni==2
                    ID_subpop = strcat(exp_txt_GEN_pop{i},'a');
                end
                ad_coeff_subpop = mean(ad_coef_matching_test,1);
                time_subpop = mean(time_matching_test,1);
                timeSTD_subpop = sqrt(sum(timeSTD_matching_test.^2)) ./ size(timeSTD_matching_test,1);
                exp_txt_GEN(counter2,1) = {ID_subpop};
                exp_txt_TEM(counter2,1) = {ID_subpop};
                exp_GEN(counter2,:) = ad_coeff_subpop;
                exp_TEM(counter2,:) = time_subpop;
                exp_TEM_STD(counter2,:) = timeSTD_subpop;
                size_subpop(counter2,:) = sum(matching_test_pop_logical);
                testID_subPop_name(matching_test_pop_logical,:) = {ID_subpop};
                counter2 = counter2 + 1;
            else
                %---Population size allows for sub-divisions ---> There could be several reference subpopulations
                %Determine the maximum possible number of subpopulations
                max_k = f_max_k(sum(matching_test_pop_logical));
                
                %Apply the selected method to extract the most likely value of k
                %Method1: multiple pairwise F-test. It performs all the possible pairwise F-test between the groups (the more robust criterion)
                p_pair = cell(1,max_k);
                k_ind = 3;
                if k_ind > max_k
                    [C_ind, ~, ~, ~] = kmeans(ad_coef_matching_test,2,'start','sample','replicates',n_rep_kmeans);
                    p_pair{2} = mpFtest(ad_coef_matching_test,C_ind);
                    binary_Ftest = p_pair{2}>0.05;
                    if any(binary_Ftest==1)
                        k_ind = 2;
                    end
                else
                    while k_ind <= max_k
                        [C_ind, ~, ~, ~] = kmeans(ad_coef_matching_test,k_ind,'start','sample','replicates',n_rep_kmeans);
                        p_pair{k_ind} = mpFtest(ad_coef_matching_test,C_ind);
                        binary_Ftest = p_pair{k_ind}>0.05;
                        if any(binary_Ftest==1)
                            if k_ind == 3
                                [C_ind, ~, ~, ~] = kmeans(ad_coef_matching_test,2,'start','sample','replicates',n_rep_kmeans);
                                p_pair{2} = mpFtest(ad_coef_matching_test,C_ind);
                                binary_Ftest = p_pair{2}>0.05;
                                if any(binary_Ftest==1)
                                    k_ind = 2;
                                end
                                break;
                            else
                                break;
                            end
                        else
                            k_ind = k_ind + 1;
                        end
                    end
                end
                optimal_k = k_ind - 1;                                                 

                %Apply the selected criterion to establish the value of k to be used, in relation with the above optmial_k
                final_k = f_k_to_use(optimal_k);
                %---In doing this an error may occur when final_k exceeds the possible values [1,size of the population]. The code below fixes these issues:
                final_k = subpop_min_number .* (final_k<subpop_min_number) + final_k .* (final_k>=subpop_min_number);   %Lower boundary
                final_k = optimal_k .* (final_k>sum(matching_test_pop_logical)) + final_k .* (final_k<=sum(matching_test_pop_logical));   %Upper boundary
                %Apply the value of k to enforce the subpopulations division
                if final_k == 1
                    if size(ad_coef_matching_test,1) >= subpop_min_size
                        %There is only one reference subpopulation
                        if ni==1
                            ID_subpop = strcat(exp_txt_GEN_pop{i},'_0');
                        elseif ni==2
                            ID_subpop = strcat(exp_txt_GEN_pop{i},'a');
                        end
                        ad_coeff_subpop = mean(ad_coef_matching_test,1);
                        time_subpop = mean(time_matching_test,1);
                        timeSTD_subpop = sqrt(sum(timeSTD_matching_test.^2)) ./ size(timeSTD_matching_test,1);
                        exp_txt_GEN(counter2,1) = {ID_subpop};
                        exp_txt_TEM(counter2,1) = {ID_subpop};
                        exp_GEN(counter2,:) = ad_coeff_subpop;
                        exp_TEM(counter2,:) = time_subpop;
                        exp_TEM_STD(counter2,:) = timeSTD_subpop;
                        size_subpop(counter2,:) = sum(matching_test_pop_logical);
                        testID_subPop_name(matching_test_pop_logical,:) = {ID_subpop};
                        counter2 = counter2 + 1;
                    end
                else
                    %There are several reference subpopulations
                    [C_ind, ad_coeff_subpop, ~, ~] = kmeans(ad_coef_matching_test,final_k,'start','sample','replicates',n_rep_kmeans);
                    test_used_logical = zeros(size(C_ind,1),1);
                    for j = 1:size(ad_coeff_subpop,1)
                        if sum(C_ind==j) >= subpop_min_size
                            if ni==1
                                ID_subpop = strcat(exp_txt_GEN_pop{i},'_',num2str(j));
                            elseif ni==2
                                ID_subpop = strcat(exp_txt_GEN_pop{i},alphabet{j});
                            end
                            time_subpop = mean(time_matching_test(C_ind==j,:),1);
                            timeSTD_subpop = sqrt(sum(timeSTD_matching_test(C_ind==j,:).^2)) ./ size(timeSTD_matching_test(C_ind==j,:),1);
                            exp_txt_GEN(counter2,1) = {ID_subpop};
                            exp_txt_TEM(counter2,1) = {ID_subpop};
                            exp_GEN(counter2,:) = ad_coeff_subpop(j,:);
                            exp_TEM(counter2,:) = time_subpop;
                            exp_TEM_STD(counter2,:) = timeSTD_subpop;
                            size_subpop(counter2,:) = sum(C_ind==j);
                            test_used_logical = test_used_logical + (C_ind==j);
                            counter2 = counter2 + 1;
                        end
                    end
                    if ni==1
                        test_current_subPop_membership = strcat(repmat(strcat(exp_txt_GEN_pop(i),'_'),length(C_ind),1),num2str(C_ind));
                    elseif ni==2
                        test_current_subPop_membership = strcat(repmat(exp_txt_GEN_pop{i},length(C_ind),1),alphabet(C_ind));
                    end
                    test_current_subPop_membership(~test_used_logical) = {''};
                    testID_subPop_name(matching_test_pop_logical,:) = regexprep(test_current_subPop_membership, '\W', '');
                end
            end
        else
            %testID_subPop_name(matching_test_pop_logical,:) = {''};
        end
        if strcmp(displaySwitch,'On')
            toc(t_iter);
        end
    end
    subpop_panel_size = counter2 - 1;
    logical_index_testID_to_delete = cell2mat(cellfun(@isempty,testID_subPop_name,'UniformOutput',false));
    n_testID_to_delete = sum(logical_index_testID_to_delete);
end

%% Algorithm Part 3: Sort dataset rows according to subpopulations
if strcmp(displaySwitch,'On')
    fprintf('\n\n\tRunning algorithm part4: Sort dataset rows by subpopulations...\n');
end
%Exclude testId to delete
testID_adCoef = testID_adCoef(~logical_index_testID_to_delete,:);
testID_time = testID_time(~logical_index_testID_to_delete,:);
testID_timeSTD = testID_timeSTD(~logical_index_testID_to_delete,:);
testID_subPop_name = testID_subPop_name(~logical_index_testID_to_delete,:);
%Extract list of populations name for test individuals
testID_pop_name = cellfun(@(x) x(1:find(ismember(x,'_'),1,'last')-1),testID_subPop_name,'UniformOutput',false);
temp2 = cellfun(@strsplit,testID_pop_name,repmat({'_pop'},size(testID_pop_name)),'UniformOutput',false);
testID_supPop_NameList = cellfun(@ (x) x{1},temp2,'UniformOutput',false);
%Find indicies to sort subpopulations within each population
supPop_NameList = unique(testID_supPop_NameList,'stable');
test_subpop_ID_sort_ind = zeros(size(testID_subPop_name));
k = 1;
for i=1:length(supPop_NameList)
    current_membership = ismember(testID_supPop_NameList,supPop_NameList(i));
    current_indices = find(current_membership==1);
    [~, sort_ind] = sort(testID_subPop_name(current_membership));
    current_indices_sorted = current_indices(sort_ind);
    test_subpop_ID_sort_ind(k:(k+sum(current_membership)-1),1) = current_indices_sorted;
    k = k + sum(current_membership);
end
