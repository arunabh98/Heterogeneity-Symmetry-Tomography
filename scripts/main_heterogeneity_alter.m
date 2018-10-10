% Increase the number of parpool workers.
parpool('local', 14)
% warning('off', 'MATLAB:rankDeficientMatrix');

% Include the moment based estimation scripts and noise scripts.
addpath(genpath('../data'));
addpath(genpath('single_axis_symmetry'));
addpath(genpath('horizontal_symmetry_utilities'));
addpath(genpath('moment_based_estimation'));
addpath(genpath('noise_scripts'));
addpath(genpath('utilities'));

% Get the images
image_size = 100;
P1 = read_process_image('refs_011.png', image_size);
P2 = read_process_image('refs_015.png', image_size);
P3 = read_process_image('refs_018.jpg', image_size);

% Constants.
non_uniform_distribution = 0;
sigmaNoiseFraction = 0.25;
if non_uniform_distribution == 0
    filename = ...
        strcat('../results/heterogeneity/', num2str(sigmaNoiseFraction*100), '_percent_noise/');
else
    filename = ...
        strcat('../results/heterogeneity/', num2str(sigmaNoiseFraction*100), '_percent_noise/non_uniform_distribution/');
end
output_size = max(size(P1));

% Experimentatal conditions.
symmetry_prior = 1;
noisy_orientations = 0;
symmetry_method = 4;
include_clustering = 1;
num_clusters = 540;
num_theta = 40000;
max_angle_error = 0;
max_shift_amplitude = 0;

% Create the folder to hold the results of the experiment.
mkdir(strcat(filename, num2str(num_theta), '/all_variables/'));
if include_clustering ~= 1
    theta_to_write = zeros(10, num_theta);
end

% The file which contains all the errors.
fileID = fopen(strcat(filename,...
    num2str(num_theta), '/result.txt'),'w');

% Write the original images.
imwrite(P1, strcat(filename,...
    num2str(num_theta), '/original_image_1.png'));
imwrite(P2, strcat(filename,...
    num2str(num_theta), '/original_image_2.png'));
imwrite(P3, strcat(filename,...
    num2str(num_theta), '/original_image_3.png'));

% Define ground truth angles and take the tomographic projection.
if non_uniform_distribution == 0
    theta = datasample(0:0.005:179, num_theta);
elseif non_uniform_distribution == 1
    theta = [datasample(0:0.005:20, num_theta/5)...
    datasample(40:0.005:60, num_theta/5)...
    datasample(80:0.005:120, 2*num_theta/5)...
    datasample(140:0.005:160, num_theta/5)];
end

% Generate the projections including outliers of class 1 and 2.
[projections, svector, original_class_of_projections] = ...
    get_projections(theta, P1, P2, P3);

% [projections, svector] = radon(P,theta);
original_projections = projections;

% Normalize s to a unit circle
smax = max(abs(svector));
svector = svector / smax;

% Define the length of the projection.
projection_length = size(projections, 1);

% Add noise to projections.
[measured_projections, sigmaNoise] = add_noise(projections, sigmaNoiseFraction);

% Initial error between the projections.
disp('**** L2-norm error between the original projections and measured projections ****');
disp(norm(measured_projections - original_projections, 'fro'));
disp('');

disp('**** Initial classification of projections ****');
[projection_classes, ~, ~, ~] =...
    classify_projections(measured_projections, num_clusters, theta, original_class_of_projections,...
        sigmaNoise);

disp('**** Initial classification performance ****')
fprintf('Number of projections classified correctly: %d',...
        sum(projection_classes ~= original_class_of_projections));

% No. of clusters to create while estimating the structure.
num_clusters = 150;

%% Estimate the first image %%
class_measured_projections = measured_projections(:, projection_classes == 1);
class_original_projections = original_projections(:, projection_classes == 1);
actual_classes = original_class_of_projections(projection_classes == 1);
class_theta = theta(projection_classes == 1);
class_num_theta = size(class_theta, 2);
% idx = detect_outliers(projections_1, theta_1, num_clusters,...
%                       sigmaNoise, actual_classes, noisy_orientations,...
%                       svector, max_angle_error, output_size);

disp('**** Reconstruct the first image ****')
reconstructed_image_1 = reconstruct_image_symmetry(...
    class_measured_projections, class_original_projections,...
    num_clusters, class_theta, sigmaNoise, class_num_theta, noisy_orientations,...
    max_shift_amplitude, svector, output_size);

%% Estimate the second image %%
class_measured_projections = measured_projections(:, projection_classes == 2);
class_original_projections = original_projections(:, projection_classes == 2);
actual_classes = original_class_of_projections(projection_classes == 2);
class_theta = theta(projection_classes == 2);
class_num_theta = size(class_theta, 2);
% idx = detect_outliers(projections_1, theta_1, num_clusters,...
%                       sigmaNoise, actual_classes, noisy_orientations,...
%                       svector, max_angle_error, output_size);

disp('**** Reconstruct the second image ****')
reconstructed_image_2 = reconstruct_image_symmetry(...
    class_measured_projections, class_original_projections,...
    num_clusters, class_theta, sigmaNoise, class_num_theta, noisy_orientations,...
    max_shift_amplitude, svector, output_size);

%% Estimate the third image %%
class_measured_projections = measured_projections(:, projection_classes == 3);
class_original_projections = original_projections(:, projection_classes == 3);
actual_classes = original_class_of_projections(projection_classes == 3);
class_theta = theta(projection_classes == 3);
class_num_theta = size(class_theta, 2);

% idx = detect_outliers(projections_1, theta_1, num_clusters,...
%                       sigmaNoise, actual_classes, noisy_orientations,...
%                       svector, max_angle_error, output_size);

disp('**** Reconstruct the third image ****')
reconstructed_image_3 = reconstruct_image_symmetry(...
    class_measured_projections, class_original_projections,...
    num_clusters, class_theta, sigmaNoise, class_num_theta, noisy_orientations,...
    max_shift_amplitude, svector, output_size);

%% Store the results %%
formatSpec = 'Final image rmse error: %4.2f \r\n';
fprintf(fileID, formatSpec, calculate_rmse_error(reconstructed_image_1, P1));
fprintf(fileID, formatSpec, calculate_rmse_error(reconstructed_image_2, P2));
fprintf(fileID, formatSpec, calculate_rmse_error(reconstructed_image_3, P3));

% Save the images.
imwrite(reconstructed_image_1, ...
    strcat(filename, num2str(num_theta), '/reconstructed_image_1.png'));
imwrite(reconstructed_image_2, ...
    strcat(filename, num2str(num_theta), '/reconstructed_image_2.png'));
imwrite(reconstructed_image_3, ...
    strcat(filename, num2str(num_theta), '/reconstructed_image_3.png'));

% Save the important variables.
save(strcat(filename, num2str(num_theta), '/all_variables/all_variables.mat'),...
    '-regexp',...
    '^(?!(measured_projections|original_projections|projections)$).');

