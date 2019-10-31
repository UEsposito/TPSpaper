
function PopulationLabelsStats = Make_pop_stats(testID_pop,testID_yearsDist,displaySwitch)

% Analyse TPS results

%% Algorithm
% Remove individuals for which there is no prediction
testID_noPrediction_index = isnan(testID_yearsDist);
testID_pop(testID_noPrediction_index) = [];
testID_yearsDist(testID_noPrediction_index) = [];

%%%---*** 1. Extract statistics for single populations
%---Average distance per population
%Extract the date in numerical format for the sorting function to properly sort
testID_pop_unique = unique(testID_pop,'stable');
temp1 = cellfun(@strsplit,testID_pop_unique,repmat("BP",size(testID_pop_unique)),'UniformOutput',false);
temp2 = cellfun(@(x) x{1}(6:end),temp1,'UniformOutput',false);
temp3 = cellfun(@str2double,temp2,'UniformOutput',false);

[~,sort_indx] = sort([temp3{:}],'d');
testID_pop_unique_sorted = testID_pop_unique(sort_indx);
n_pop = length(testID_pop_unique_sorted);
n_testIDs_per_pop = nan(n_pop,1);
bin_edges = [0,200:200:1000,inf];
distance_statistics_pop = nan(n_pop,length(bin_edges)-1);
for i = 1:n_pop

    %Print progress on screen
    if strcmp(displaySwitch,'On')
        fprintf('\t\tExaminng superpopulation #%d/%d\n',i,n_pop)
    end

    %Find matching individuals for the current population
    testID_pop_match = ismember(testID_pop,testID_pop_unique_sorted(i));
    if testID_pop_match == 0
        error('No match found for the current superpopulation')
    end

    %Calculate number of matching individuals in the current population
    n_testIDs_per_pop(i) = sum(testID_pop_match);

    %Calculate current population average distances
    distance_statistics_pop(i,:) = 100 .* histcounts(testID_yearsDist(testID_pop_match),bin_edges) ./ n_testIDs_per_pop(i);
end

%Extract the average distance
t1 = distance_statistics_pop .* repmat(n_testIDs_per_pop,1,size(distance_statistics_pop,2));


%%%---*** 2. Extract global quantities
%---Global average statistics of the distance across all individuals
distance_statistics_average = sum(t1,1) ./ size(testID_pop,1);


%%%---*** 3. Assign quantities to final variables
testID_pop_unique_sorted_plot = regexprep(testID_pop_unique_sorted,'_',' ');
distance_statistics = [distance_statistics_average;distance_statistics_pop];
PopulationLabelsStats = table;
PopulationLabelsStats.Labels = ['Average';testID_pop_unique_sorted_plot];
PopulationLabelsStats.Dist_0_200 = distance_statistics(:,1);
PopulationLabelsStats.Dist_200_400 = distance_statistics(:,2);
PopulationLabelsStats.Dist_400_600 = distance_statistics(:,3);
PopulationLabelsStats.Dist_600_800 = distance_statistics(:,4);
PopulationLabelsStats.Dist_800_1000 = distance_statistics(:,5);
PopulationLabelsStats.Dist_1000more = distance_statistics(:,6);
PopulationLabelsStats.nIDs = [sum(n_testIDs_per_pop);n_testIDs_per_pop];

