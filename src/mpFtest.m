
function p = mpFtest(ad_coef_matching_test, C_ind)

% Multiple, pairwise F-tests

%% Algorithm
k = length(unique(C_ind));
clusters_combinations = combnk(1:k,2);
n_combinations = size(clusters_combinations,1);
p = zeros(1,n_combinations);

for i=1:n_combinations
    c1_ind = clusters_combinations(i,1);
    c2_ind = clusters_combinations(i,2);
    
    c1_points = ad_coef_matching_test(C_ind==c1_ind,:);
    c2_points = ad_coef_matching_test(C_ind==c2_ind,:);
    
    c1_n = size(c1_points,1);
    c2_n = size(c2_points,1);
    
    c1_centroid = mean(c1_points,1);
    c2_centroid = mean(c2_points,1);
    
    points_tot = [c1_points;c2_points];
    n_tot = c1_n + c2_n;
    mean_tot = mean(points_tot,1);
    
    ss_total_pair = sumsqr(points_tot - repmat(mean_tot,n_tot,1));
    ss_between_pair = c1_n * sumsqr(c1_centroid - mean_tot) + c2_n * sumsqr(c2_centroid - mean_tot);
    ratio_bt_pair = ss_between_pair / ss_total_pair;
    ss_within_pair = sumsqr(c1_points-repmat(c1_centroid,c1_n,1)) + sumsqr(c2_points-repmat(c2_centroid,c2_n,1));
    
    %Control statement
    if (abs(ss_within_pair+ss_between_pair) - ss_total_pair) > 0.0001
        error('Sums of squares do not sums up correctly');
    end
    
    if ( (ss_between_pair/(c1_n+c2_n)) < 1E-11 ) && ( (ss_within_pair/(c1_n+c2_n)) < 1E-11 )
       p(i) = 1;
       continue;
    end
    
    ratio_bw_pair = ss_between_pair / ss_within_pair;
    ratio_dof_pair = (n_tot-2)/(2-1);
    F_pair = ratio_bw_pair * ratio_dof_pair;
    p(i) = fcdf(F_pair,1,n_tot-2,'upper');
end
