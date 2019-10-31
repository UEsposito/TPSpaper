
function [] = Paper_results_main(varargin)

%This function reproduces the results in Fig.3 of the manuscript Esposito et al. "A genomic dating tool for ancient genomes resolves the origins of hundreds of Eurasian genomes"

% Inputs (all optional):
%    input_file - xlsx file with the aDNA dataset (it needs to follow a specific format)
%               Default value: "../data/TPSpaper.xlsx"
%    n_rand_repetitions - Number of random training/testing repetitions
%               Default value: 500
%    training_size - Number of individuals used to build the reference panel for each random iteration
%               Default value: 400
%
% Outputs (stored in the folder ../results/):
%    Output1 - xlsx file with TPS results
%    Output2 - multiple figures
%
% How to call this script:
%   - All arguments are optional
%   - Recommended call (which uses default values): Paper_results_main
%
% Subfunctions called: 
%   Select_reference_ids.m
%   Build_reference_panel.m (which calls mpFtest.m)
%   Run_TPS.m
%   Make_pop_stats.m
%
% Author: Umberto Esposito
% email: umberto.espos@gmail.com
% Created: October 2019
% Last edited: 30 October 2019


%% Assign input arguments
%Default values
inputs_defaults = {'../data/TPSpaper.xlsx', 500, 400};

%Assign values to use
inputs = inputs_defaults;
idx = ~cellfun('isempty', varargin);
inputs(idx) = varargin(idx);

%Define input variables
input_file = inputs{1};
n_rand_repetitions = inputs{2};
training_size = inputs{3};


%% Model parameters - General architecture
N_bootstrap_runs = 2;
%Calibration discards individual predicted more than 1,000 years from their radiocarbon date
max_av_dist = 1000;
%Outliers criterion for calibration
linFitStd = 1.5;
%Consider only individuals younger than 14,000 years
apriori_thresholding_value = 14000; %years


%% Read inputs
%Read data
opts = detectImportOptions(input_file);
IDsTable_all = readtable(input_file,opts);


%% Variables and data
%Select data to use according to the parameters'value chosen above
IDs_toDiscard_logicalInd = IDsTable_all.DateBP > apriori_thresholding_value;
IDsTable_all(IDs_toDiscard_logicalInd,:) = [];

%Names of the distinct datasets in the input spreadsheet (Column E "Dataset")
Datasets_LabelsList = {'Radiocarbon dating';'Other dating';'Uncertain'}; %sorted according to the order I want them to appear in the structure
Datasets_NamesList = {'RadioDated';'OtherDating';'UncertainDate'};
Datasets_LabelsList_fromTable = unique(IDsTable_all.Dataset);
if any( strcmp(sort(Datasets_LabelsList),Datasets_LabelsList_fromTable)==0 )
    error('Datasets names need to be updated');
end
N_testDatasets = length(Datasets_LabelsList);

%Admixture components in the input spreadsheet (Columns F-M)
adComponents_names = {'AncientComponent1','AncientComponent2','AncientComponent3','AncientComponent4','AncientComponent5','ModernComponent1','ModernComponent2','ModernComponent3'};

%Group individuals into populations by rounding to the nearest 500
Population = strcat('Date_',string(round(IDsTable_all.DateBP*2,-3)/2),'BP');
Population(ismissing(Population)) = 'No_Date';
IDsTable_all = addvars(IDsTable_all,Population,'After','ID');
strcat('Date_',num2str(round(IDsTable_all.DateBP*2,-3)));

%Define the data structure
Datasets = struct();
%---Initialise dataset #1 (radiocarbon dates, both training and test dataset #1), 
Datasets(1).Name = Datasets_NamesList{1};
Datasets(1).Content = IDsTable_all(strcmp(IDsTable_all.Dataset,Datasets_LabelsList{1}),:);
Datasets(1).nIDs = size(Datasets(1).Content,1);
Datasets(1).RefPanelSinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
Datasets(1).RefPanelStdSinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
Datasets(1).RefPanelMinGenDistSinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
Datasets(1).SinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
Datasets(1).StdSinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
Datasets(1).MinGenDistSinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
%---Initialise the other datasets
for i_d = 2:N_testDatasets
    %---Datasets name
    Datasets(i_d).Name = Datasets_NamesList{i_d};
    %---Datasets content
    Datasets(i_d).Content = IDsTable_all(strcmp(IDsTable_all.Dataset,Datasets_LabelsList{i_d}),:);
    %---Datasets size
    Datasets(i_d).nIDs = size(Datasets(i_d).Content,1);
    %---Initialise Datasets single predictions fields
    Datasets(i_d).SinglePredictions = nan(Datasets(i_d).nIDs,n_rand_repetitions);
    Datasets(i_d).StdSinglePredictions = nan(Datasets(i_d).nIDs,n_rand_repetitions);
    Datasets(i_d).MinGenDistSinglePredictions = nan(Datasets(i_d).nIDs,n_rand_repetitions);
end

%Results folder
output_folder = '../results/';


%% TPS bootstrapping routine
fprintf('\nRoutine started:\n\n');

%Start routine
for i_run = 1:N_bootstrap_runs
    
    %Update datasets composition
    if i_run == 2
        %---Determine the individuals to rmeove
        IDsToRemove_index = TPS_RefPanelResults.AccuracyDistMeanTime > max_av_dist;

        %---Last dataset (noDating) moves one position in the Datasets structure, to make room for the new one of IDs remove by the calibration
        %(this is to keep using foor loops below, as the last dataset (noDating) is the only one that does not enter some for loops)
        Datasets(N_testDatasets+1) = Datasets(N_testDatasets);
        N_testDatasets = N_testDatasets + 1;
        %---Define the new dataset, that goes in the (end-1) position  of the Datasets structure
        Datasets(N_testDatasets-1).Name = 'CalibrationRemovals';
        Datasets(N_testDatasets-1).Content = Datasets(1).Content(IDsToRemove_index,:);
        Datasets(N_testDatasets-1).Content.Dataset(:) = {'CalibRemov'};
        TPS_CalibDiscarded = Datasets(N_testDatasets-1).Content.RefPanel;
        Datasets(N_testDatasets-1).Content.RefPanel = [];
        Datasets(N_testDatasets-1).nIDs = size(Datasets(N_testDatasets-1).Content,1);
        Datasets(N_testDatasets-1).SinglePredictions = nan(Datasets(N_testDatasets-1).nIDs,n_rand_repetitions);
        Datasets(N_testDatasets-1).StdSinglePredictions = nan(Datasets(N_testDatasets-1).nIDs,n_rand_repetitions);
        Datasets(N_testDatasets-1).MinGenDistSinglePredictions = nan(Datasets(N_testDatasets-1).nIDs,n_rand_repetitions);
        %---Clean the fields of the radiocarbon dated individuals for the new iteration (test dataset #1)
        Datasets(1).Content(IDsToRemove_index,:) = [];
        Datasets(1).Content.RefPanel = [];
        Datasets(1).nIDs = size(Datasets(1).Content,1);
        Datasets(1).RefPanelSinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
        Datasets(1).RefPanelStdSinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
        Datasets(1).RefPanelMinGenDistSinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
        Datasets(1).SinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
        Datasets(1).StdSinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
        Datasets(1).MinGenDistSinglePredictions = nan(Datasets(1).nIDs,n_rand_repetitions);
        Datasets(1).RefPanelLinFit1 = [];
        Datasets(1).RefPanelLinFit2 = [];
    end
    
    
    %% TPS single runs
    for i_rep = 1:n_rand_repetitions

        %Print progress on screen
        t_iter = tic;
        fprintf('\tRun #%d iteration #%d\n',i_run,i_rep);
        
        %Randomly choose 400 individuals among the ones with radiocarbon date, according to the chosen randomisation mode (training dataset/reference panel)
        [Dataset_training, Dataset_training_indexIDs] = Select_reference_ids(Datasets(1).Content,training_size);

        %Create the TPS reference panel
        adCoef_training = table2array(Dataset_training(:,adComponents_names));
        [~,exp_txt_GEN,exp_txt_TEM,exp_GEN,exp_TEM,exp_TEM_STD] = Build_reference_panel(Dataset_training, adCoef_training, 'Off');

        %Control statement
        if any(strcmp(exp_txt_GEN,exp_txt_TEM) == 0)
            error('Reference panel populations names are not the same');
        end

        %Run TPS predictive routine
        if i_run ~= N_bootstrap_runs %Calibration
            %adCoef_refPanel = table2array(Datasets(1).Content(:,adComponents_names));
            [predicted_time,predicted_time_STD,min_gen_dist] = Run_TPS(exp_txt_GEN, exp_GEN, exp_TEM, exp_TEM_STD, adCoef_training, 'Off');
            %Store the results
            Datasets(1).RefPanelSinglePredictions(Dataset_training_indexIDs,i_rep) = predicted_time;
            Datasets(1).RefPanelStdSinglePredictions(Dataset_training_indexIDs,i_rep) = predicted_time_STD;
            Datasets(1).RefPanelMinGenDistSinglePredictions(Dataset_training_indexIDs,i_rep) = min_gen_dist;

        elseif i_run == N_bootstrap_runs
            %Build a unique dataset
            adCoef_all = [];
            for i_d = 1:N_testDatasets
                adCoef_dataset = table2array(Datasets(i_d).Content(:,adComponents_names));
                adCoef_all = [adCoef_all; adCoef_dataset];
            end
            %Run TPS
            [predicted_time,predicted_time_STD,min_gen_dist] = Run_TPS(exp_txt_GEN, exp_GEN, exp_TEM, exp_TEM_STD, adCoef_all, 'Off');
            %Store the results
            %---Reference panel
            Datasets(1).RefPanelSinglePredictions(Dataset_training_indexIDs,i_rep) = predicted_time(Dataset_training_indexIDs);
            Datasets(1).RefPanelStdSinglePredictions(Dataset_training_indexIDs,i_rep) = predicted_time_STD(Dataset_training_indexIDs);
            Datasets(1).RefPanelMinGenDistSinglePredictions(Dataset_training_indexIDs,i_rep) = min_gen_dist(Dataset_training_indexIDs);
            %---Test dataset #1
            Dataset_NonTraining_indexIDs = setdiff(1:Datasets(1).nIDs,Dataset_training_indexIDs);
            Datasets(1).SinglePredictions(Dataset_NonTraining_indexIDs,i_rep) = predicted_time(Dataset_NonTraining_indexIDs);
            Datasets(1).StdSinglePredictions(Dataset_NonTraining_indexIDs,i_rep) = predicted_time_STD(Dataset_NonTraining_indexIDs);
            Datasets(1).MinGenDistSinglePredictions(Dataset_NonTraining_indexIDs,i_rep) = min_gen_dist(Dataset_NonTraining_indexIDs);
            %---Other test datasets
            ind_first = 1;
            ind_last = Datasets(1).nIDs;
            for i_d = 2:N_testDatasets
                ind_first = ind_first + Datasets(i_d-1).nIDs;
                ind_last = ind_last + Datasets(i_d).nIDs;
                Datasets(i_d).SinglePredictions(:,i_rep) = predicted_time(ind_first:ind_last);
                Datasets(i_d).StdSinglePredictions(:,i_rep) = predicted_time_STD(ind_first:ind_last);
                Datasets(i_d).MinGenDistSinglePredictions(:,i_rep) = min_gen_dist(ind_first:ind_last);
            end
        end
        toc(t_iter)
    end

    %Control statements
    if any( sum( ~isnan(Datasets(1).RefPanelSinglePredictions) ) ~= training_size)
       error('Predictions not correctly stored');
    end
    if any( sum( ~isnan(Datasets(1).RefPanelStdSinglePredictions) ) ~= training_size)
       error('Predictions errors not correctly stored');
    end
    
    
    %% TPS individuals performance (averages across single runs)
    %Assess TPS performance on the reference panel's individuals
    %---Define empty table
    TPS_RefPanelResults = table;
    %---TPS average quantities. Note that even if StdSinglePredictions correspond to 2sigma, both mean value (TPSpredictions) and 1std (TPSstdPredictions) do not change
    TPS_RefPanelResults.TPSpredictions = nansum(Datasets(1).RefPanelSinglePredictions ./ (Datasets(1).RefPanelStdSinglePredictions.^2),2) ./ nansum(1 ./ (Datasets(1).RefPanelStdSinglePredictions.^2),2);
    TPS_RefPanelResults.TPSstdPredictions = sqrt(nansum( (1 ./ (Datasets(1).RefPanelStdSinglePredictions.^2)) .* ((Datasets(1).RefPanelSinglePredictions - TPS_RefPanelResults.TPSpredictions).^2), 2) ./ (nansum(1 ./ (Datasets(1).RefPanelStdSinglePredictions.^2),2)));
    for i=1:size(Datasets(1).RefPanelSinglePredictions,1)
        SinglePredictions = sort(Datasets(1).RefPanelSinglePredictions(i,~isnan(Datasets(1).RefPanelSinglePredictions(i,:))));
        margin = (length(SinglePredictions) * 5 / 100);
        TPS_RefPanelResults.TPSCIupper(i,:) = SinglePredictions(end-ceil(margin)+1);
        TPS_RefPanelResults.TPSCIlower(i,:) = SinglePredictions(floor(margin));
    end
    TPS_RefPanelResults.TPSminGenDist = nanmean(Datasets(1).RefPanelMinGenDistSinglePredictions,2);
    %---TPS final accuracy estimators
    TPS_RefPanelResults.AccuracyDistMeanTime = abs(Datasets(1).Content.DateBP - TPS_RefPanelResults.TPSpredictions);
    TPS_RefPanelResults.AccuracyOverlap2StdMeanTime = Datasets(1).Content.DateBP <= (TPS_RefPanelResults.TPSpredictions + 2*TPS_RefPanelResults.TPSstdPredictions) & Datasets(1).Content.DateBP >= (TPS_RefPanelResults.TPSpredictions - 2*TPS_RefPanelResults.TPSstdPredictions);
    TPS_RefPanelResults.AccuracyOverlap3StdMeanTime = Datasets(1).Content.DateBP <= (TPS_RefPanelResults.TPSpredictions + 3*TPS_RefPanelResults.TPSstdPredictions) & Datasets(1).Content.DateBP >= (TPS_RefPanelResults.TPSpredictions - 3*TPS_RefPanelResults.TPSstdPredictions);
    TPS_RefPanelResults.AccuracyOverlap2Std = (TPS_RefPanelResults.TPSpredictions + 2*TPS_RefPanelResults.TPSstdPredictions) >= (Datasets(1).Content.DateBP - Datasets(1).Content.DeltaT95CI) & (Datasets(1).Content.DateBP + Datasets(1).Content.DeltaT95CI) >= (TPS_RefPanelResults.TPSpredictions - 2*TPS_RefPanelResults.TPSstdPredictions);
    TPS_RefPanelResults.AccuracyOverlap3Std = (TPS_RefPanelResults.TPSpredictions + 3*TPS_RefPanelResults.TPSstdPredictions) >= (Datasets(1).Content.DateBP - 3*Datasets(1).Content.DeltaT95CI/2) & (Datasets(1).Content.DateBP + 3*Datasets(1).Content.DeltaT95CI/2) >= (TPS_RefPanelResults.TPSpredictions - 3*TPS_RefPanelResults.TPSstdPredictions);
    %---Merge back into the Dataset main table
    Datasets(1).Content.RefPanel = TPS_RefPanelResults;

    if i_run == N_bootstrap_runs
        %Assess TPS performance on the test datasets' individuals
        for i_d = 1:N_testDatasets
            %---Define empty table
            TPS_TestDatasetResults = table;
            %---TPS average quantities
            TPS_TestDatasetResults.TPSpredictions = nansum(Datasets(i_d).SinglePredictions ./ (Datasets(i_d).StdSinglePredictions.^2),2) ./ nansum(1 ./ (Datasets(i_d).StdSinglePredictions.^2),2);
            TPS_TestDatasetResults.TPSstdPredictions = sqrt(nansum( (1 ./ (Datasets(i_d).StdSinglePredictions.^2)) .* ((Datasets(i_d).SinglePredictions - TPS_TestDatasetResults.TPSpredictions).^2), 2) ./ (nansum(1 ./ (Datasets(i_d).StdSinglePredictions.^2),2)));
            for i=1:size(Datasets(i_d).SinglePredictions,1)
                SinglePredictions = sort(Datasets(i_d).SinglePredictions(i,~isnan(Datasets(i_d).SinglePredictions(i,:))));
                if isempty(SinglePredictions)
                    TPS_TestDatasetResults.TPSCIupper(i,:) = nan;
                    TPS_TestDatasetResults.TPSCIlower(i,:) = nan;
                else
                    margin = (length(SinglePredictions) * 5 / 100);
                    upper_element = length(SinglePredictions)-ceil(margin)+1;
                    upper_element = (upper_element>length(SinglePredictions)) .* length(SinglePredictions) + (upper_element<=length(SinglePredictions)) .* upper_element;
                    lower_element = floor(margin);
                    lower_element = (lower_element<1) .* 1 + (lower_element>=1) .* lower_element;
                    TPS_TestDatasetResults.TPSCIupper(i,:) = SinglePredictions(upper_element);
                    TPS_TestDatasetResults.TPSCIlower(i,:) = SinglePredictions(lower_element);
                end
            end
            TPS_TestDatasetResults.TPSminGenDist = nanmean(Datasets(i_d).MinGenDistSinglePredictions,2);
            %---TPS final accuracy estimators, except for last test dataset (noDating)
            if i_d < N_testDatasets
                TPS_TestDatasetResults.AccuracyDistMeanTime = abs(Datasets(i_d).Content.DateBP - TPS_TestDatasetResults.TPSpredictions);
                TPS_TestDatasetResults.AccuracyOverlap2StdMeanTime = Datasets(i_d).Content.DateBP <= (TPS_TestDatasetResults.TPSpredictions + 2*TPS_TestDatasetResults.TPSstdPredictions) & Datasets(i_d).Content.DateBP >= (TPS_TestDatasetResults.TPSpredictions - 2*TPS_TestDatasetResults.TPSstdPredictions);
                TPS_TestDatasetResults.AccuracyOverlap3StdMeanTime = Datasets(i_d).Content.DateBP <= (TPS_TestDatasetResults.TPSpredictions + 3*TPS_TestDatasetResults.TPSstdPredictions) & Datasets(i_d).Content.DateBP >= (TPS_TestDatasetResults.TPSpredictions - 3*TPS_TestDatasetResults.TPSstdPredictions);
                TPS_TestDatasetResults.AccuracyOverlap2Std = (TPS_TestDatasetResults.TPSpredictions + 2*TPS_TestDatasetResults.TPSstdPredictions) >= (Datasets(i_d).Content.DateBP - 2*Datasets(i_d).Content.DeltaT95CI) & (Datasets(i_d).Content.DateBP + 2*Datasets(i_d).Content.DeltaT95CI) >= (TPS_TestDatasetResults.TPSpredictions - 2*TPS_TestDatasetResults.TPSstdPredictions);
                TPS_TestDatasetResults.AccuracyOverlap3Std = (TPS_TestDatasetResults.TPSpredictions + 3*TPS_TestDatasetResults.TPSstdPredictions) >= (Datasets(i_d).Content.DateBP - 3*Datasets(i_d).Content.DeltaT95CI) & (Datasets(i_d).Content.DateBP + 3*Datasets(i_d).Content.DeltaT95CI) >= (TPS_TestDatasetResults.TPSpredictions - 3*TPS_TestDatasetResults.TPSstdPredictions);
            end
            %---Merge back into the Dataset main table
            if i_d == 1
                Datasets(i_d).Content.TestDataset = TPS_TestDatasetResults;
            else
                Datasets(i_d).Content = [Datasets(i_d).Content,TPS_TestDatasetResults];
            end
        end
    end
    

    %% TPS global performance (Linear fit)
    %Linear fit on the reference panel's individuals
    %---Fit 1
    Datasets(1).RefPanelLinFit1 = fitlm(Datasets(1).Content.DateBP, Datasets(1).Content.RefPanel.TPSpredictions);
    %---Outliers
    Datasets(1).Content.RefPanel.ResidFit1 = Datasets(1).RefPanelLinFit1.Residuals.Standardized;
    outliers_index = Datasets(1).RefPanelLinFit1.Residuals.Standardized>=linFitStd | Datasets(1).RefPanelLinFit1.Residuals.Standardized<=-linFitStd;
    %---Fit 2
    Datasets(1).RefPanelLinFit2 = fitlm(Datasets(1).Content.DateBP, Datasets(1).Content.RefPanel.TPSpredictions,'Exclude',outliers_index);

    if i_run == N_bootstrap_runs
        %Linear fit on the test #1 datasets' individuals
        %---Fit 1
        Datasets(1).LinFit1 = fitlm(Datasets(1).Content.DateBP, Datasets(1).Content.TestDataset.TPSpredictions);
        %---Outliers
        Datasets(1).Content.TestDataset.ResidFit1 = Datasets(1).LinFit1.Residuals.Standardized;
        outliers_index = Datasets(1).LinFit1.Residuals.Standardized>=linFitStd | Datasets(1).LinFit1.Residuals.Standardized<=-linFitStd;
        %---Fit 2
        Datasets(1).LinFit2 = fitlm(Datasets(1).Content.DateBP, Datasets(1).Content.TestDataset.TPSpredictions,'Exclude',outliers_index);

        %Linear fit on the other test datsets, except for last test dataset (noDating)
        for i_d = 2:N_testDatasets-1
            %---Fit 1
            Datasets(i_d).LinFit1 = fitlm(Datasets(i_d).Content.DateBP, Datasets(i_d).Content.TPSpredictions);
            %---Outliers
            Datasets(i_d).Content.ResidFit1 = Datasets(i_d).LinFit1.Residuals.Standardized;
            outliers_index = Datasets(i_d).LinFit1.Residuals.Standardized>=linFitStd | Datasets(i_d).LinFit1.Residuals.Standardized<=-linFitStd;
            %---Fit 2
            Datasets(i_d).LinFit2 = fitlm(Datasets(i_d).Content.DateBP, Datasets(i_d).Content.TPSpredictions,'Exclude',outliers_index);
        end
    end
    
    
    %% TPS population performance (averages across individuals)
    if i_run == N_bootstrap_runs
        Datasets(1).RefPanelPopulationLabelsStats = Make_pop_stats(Datasets(1).Content.Population, Datasets(1).Content.RefPanel.AccuracyDistMeanTime, 'Off');
        Datasets(1).PopulationLabelsStats = Make_pop_stats(Datasets(1).Content.Population, Datasets(1).Content.TestDataset.AccuracyDistMeanTime, 'Off');
        for i_d = 2:N_testDatasets-1
            Datasets(i_d).PopulationLabelsStats = Make_pop_stats(Datasets(i_d).Content.Population, Datasets(i_d).Content.AccuracyDistMeanTime, 'Off');
        end
    end
    
end

%% Save & print on file
%Create folder if needed
if exist(output_folder, 'dir') ~= 7
    mkdir(output_folder);
end

%Save mat file
save([output_folder,'TPS_results.mat']);

%Delete existing file
outputFilename = ['TPS_results_',num2str(n_rand_repetitions),'repetitions.xlsx'];
if exist([output_folder,outputFilename], 'file') == 2
    delete([output_folder,outputFilename]);
end

%Write xlsx results file
writetable(Datasets(1).Content,[output_folder,outputFilename],'Sheet','RefPanel','Range','A1:N1000');
writetable(Datasets(1).Content.RefPanel,[output_folder,outputFilename],'Sheet','RefPanel','Range','O1:Y1000');
writetable(Datasets(1).Content,[output_folder,outputFilename],'Sheet',Datasets(1).Name,'Range','A1:N1000');
writetable(Datasets(1).Content.TestDataset,[output_folder,outputFilename],'Sheet',Datasets(1).Name,'Range','O1:Y1000');
for i_d = 2:N_testDatasets
    writetable(Datasets(i_d).Content,[output_folder,outputFilename],'Sheet',Datasets(i_d).Name,'Range','A1');
end


%% Plot figures - Bar plot with accuracies of predictions
for i_d = 1:2
    
    if i_d==1
        Results_Struc = Datasets(1).RefPanelPopulationLabelsStats;
        Results_avDist = mean(Datasets(1).Content.RefPanel.AccuracyDistMeanTime);
        figFilename = 'RadioDated_Training_Accuracy';
    elseif i_d==2
        Results_Struc = Datasets(1).PopulationLabelsStats;
        Results_avDist = nanmean(Datasets(1).Content.TestDataset.AccuracyDistMeanTime);
        figFilename = 'RadioDated_Testing_Accuracy';
    end

    fi = figure;
    %---Bar plot
    bar_handle = bar(table2array(Results_Struc(:,2:end-1)),'stacked');
    bar_handle(1).FaceColor = [.9 1 1];
    bar_handle(2).FaceColor = [0 1 1];
    bar_handle(3).FaceColor = [0. 0.8 .7];
    bar_handle(4).FaceColor = [0. 0.6 .7];
    bar_handle(5).FaceColor = [0.5 0 .5];
    bar_handle(6).FaceColor = [0. 0. 1];
    
    %---Plot appearance
    box off;
    grid off;
    ax = gca;
    set(ax,'xLim',[0,size(Results_Struc.Labels,1)+1],'XTickLabel',[]);
    ax.YLim = ([0, 100]);
    ylabel('Prediction accuracy (%)');
    hold on
    plot([0,size(Results_Struc.Labels,1)+1],[25,25],'--k');
    plot([0,size(Results_Struc.Labels,1)+1],[50,50],'--k');
    plot([0,size(Results_Struc.Labels,1)+1],[75,75],'--k');
    set(ax,'xLim',[0,size(Results_Struc.Labels,1)+1],'XTick',1:size(Results_Struc.Labels,1),'XTickLabel',Results_Struc.Labels,'fontsize',4);
    ax.XTickLabelRotation = 90;
    leg = legend('0-200 Years','200-400 Years','400-600 Years','600-800 Years','800-1000 Years','1000+ Years');
    set(leg,'location','northeastoutside');

    %---Add text
    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[425,200,110,15]);
    set(mTextBox,'String','Number of predictions');
    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[425,187,110,15]);
    set(mTextBox,'String',[num2str(Results_Struc.nIDs(1)),' individuals']);

    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[425,160,110,15]);
    set(mTextBox,'String','Average distance');
    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[425,147,110,15]);
    set(mTextBox,'String',[num2str(round(Results_avDist,1)),' Years']);

    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[425,120,110,15]);
    set(mTextBox,'String','Within 200 years');
    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[425,107,110,15]);
    set(mTextBox,'String',[num2str(round(Results_Struc.Dist_0_200(1),1)),'%']);

    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[425,80,110,15]);
    set(mTextBox,'String','Within 400 years');
    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[425,67,110,15]);
    set(mTextBox,'String',[num2str(round(Results_Struc.Dist_0_200(1) + Results_Struc.Dist_200_400(1),1)),'%']);

    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[425,40,110,15]);
    set(mTextBox,'String','Above 1000 years');
    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[425,27,110,15]);
    set(mTextBox,'String',[num2str(round(Results_Struc.Dist_1000more(1),1)),'%']);

    %---Print figure
    set(fi,'color','w');
    set(fi,'units','centimeters');
    op = get(fi,'OuterPosition');
    %create a page the same size as the figure
    set(fi, 'PaperUnits','centimeters');
    set(fi, 'PaperSize', [op(3) op(4)]);
    %locate the figure in the page
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperPosition',[0 0 op(3) op(4)]);
    print(fi,[output_folder,figFilename],'-dpng','-r600');
end


%% Plot figures - Correlation between radiocarbon dates and predictions
for i_d = 1:2
    
    if i_d==1
        linearFit1 = Datasets(1).RefPanelLinFit1;
        %linearFit2 = Datasets(1).RefPanelLinFit2;
        GivenTime = Datasets(1).Content.DateBP;
        PredictedTime = Datasets(1).Content.RefPanel.TPSpredictions;
        PredictedTimeStd = Datasets(1).Content.RefPanel.TPSstdPredictions;
        CI_upper = abs(Datasets(1).Content.RefPanel.TPSCIupper-PredictedTime);
        CI_lower = abs(Datasets(1).Content.RefPanel.TPSCIlower-PredictedTime);
        figFilename = 'RadioDated_Training_Correlation';
    elseif i_d==2
        linearFit1 = Datasets(1).LinFit1;
        %linearFit2 = Datasets(1).LinFit2;
        GivenTime = Datasets(1).Content.DateBP;
        PredictedTime = Datasets(1).Content.TestDataset.TPSpredictions;
        PredictedTimeStd = Datasets(1).Content.TestDataset.TPSstdPredictions;
        CI_upper = abs(Datasets(1).Content.TestDataset.TPSCIupper-PredictedTime);
        CI_lower = abs(Datasets(1).Content.TestDataset.TPSCIlower-PredictedTime);
        figFilename = 'RadioDated_Testing_Correlation';
    end

    fi = figure;
    %---Scatter plot
    scatter(GivenTime,PredictedTime,15,'filled');
    
    %---Error bars
    hold on
    errorbar(GivenTime,PredictedTime,CI_upper,CI_lower,'LineStyle','none','Color','b');
    
    %---Plot appearance
    axis_UpperLim = max([GivenTime;PredictedTime]) + 1000;
    plot([0,axis_UpperLim],[0,axis_UpperLim]);
    axis square;
    xlabel('Radiocarbon date');
    ylabel('TPS predicted date');
    
    %---Linear fits
    x_plot = 0:10:axis_UpperLim;
    plot(x_plot, linearFit1.Coefficients.Estimate(2) .* x_plot + linearFit1.Coefficients.Estimate(1),'r')
    %plot(x_plot, linearFit2.Coefficients.Estimate(2) .* x_plot + linearFit2.Coefficients.Estimate(1),'b')
    %outliers_index = linearFit1.Residuals.Standardized>=linFitStd | linearFit1.Residuals.Standardized<=-linFitStd;
    %scatter(GivenTime(outliers_index), PredictedTime(outliers_index),15,'filled','r');
    
    %---Add text
    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[300,109,200,15]);
    set(mTextBox,'String','Fit 1');
    mTextBox = uicontrol('style','text');
    set(mTextBox,'BackgroundColor',[1 1 1]);
    set(mTextBox,'Position',[300,96,200,15]);
    set(mTextBox,'String',strcat('Rsquared = ',num2str(round(linearFit1.Rsquared.Ordinary,2)),' b = ',num2str(round(linearFit1.Coefficients.Estimate(1),2)),' a = ',num2str(round(linearFit1.Coefficients.Estimate(2),2))));

    %mTextBox = uicontrol('style','text');
    %set(mTextBox,'BackgroundColor',[1 1 1]);
    %set(mTextBox,'Position',[300,70,200,15]);
    %set(mTextBox,'String','Fit 2');
    %mTextBox = uicontrol('style','text');
    %set(mTextBox,'BackgroundColor',[1 1 1]);
    %set(mTextBox,'Position',[300,57,200,15]);
    %set(mTextBox,'String',strcat('Rsquared = ',num2str(round(linearFit2.Rsquared.Ordinary,2)),' b = ',num2str(round(linearFit2.Coefficients.Estimate(1),2)),' a = ',num2str(round(linearFit2.Coefficients.Estimate(2),2))));
    
    %---Print figure
    set(fi,'color','w');
    set(fi,'units','centimeters');
    op = get(fi,'OuterPosition');
    %create a page the same size as the figure
    set(fi, 'PaperUnits','centimeters');
    set(fi, 'PaperSize', [op(3) op(4)]);
    %locate the figure in the page
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperPosition',[0 0 op(3) op(4)]);
    print(fi,[output_folder,figFilename],'-dpng','-r600');
end

