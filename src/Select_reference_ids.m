
function [Dataset_training, Dataset_training_indexIDs] = Select_reference_ids(Dataset,Dataset_training_nIDs)

% Construct training (reference) and testing (generalisation) datasets

%% Algorithm    
%Define the unique groups
IDs_groupLabel = Dataset.Population;
IDs_groupLabel_unique = unique(IDs_groupLabel);
n_groups = length(IDs_groupLabel_unique);

%For each group, assign 2 random IDs to the reference dataset (this is in order to guarantee that all groups will be represented in the reference panel)
IDs_ref_index = [];
for i=1:n_groups
    match_index = find(ismember(IDs_groupLabel,IDs_groupLabel_unique{i}));
    if isempty(match_index)
        error('No match')
    elseif length(match_index) <= 2
        IDs_ref_index = [IDs_ref_index; match_index];
    else
        p = randperm(length(match_index),2);
        IDs_ref_index = [IDs_ref_index; match_index(p)];
    end
end
refPanel_p1 = Dataset(IDs_ref_index,:);

%Assign the remaining individuals
Dataset_temp = Dataset;
Dataset_temp(IDs_ref_index,:) = [];
refPanel_p2 = datasample(Dataset_temp,Dataset_training_nIDs - length(IDs_ref_index),'Replace',false);
Dataset_training = [refPanel_p1;refPanel_p2];

%Extract their indices in the main dataset
[~,Dataset_training_indexIDs] = ismember(Dataset_training.ID,Dataset.ID);
Dataset_training_indexIDs = Dataset_training_indexIDs';

