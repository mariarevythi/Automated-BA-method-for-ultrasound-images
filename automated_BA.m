close all;
clear;
clc;

path = 'C:\Users\maria\Documents\MARIA\BioMed\Thesis\PART2\week6\'; % input path
Paths.OrigImgPath = fullfile(path, 'SRI17\'); % input frames folder
ImgPara.ImagePrefix = 'IMG'; % image file name prefix
ImgPara.ImageFormat =  '.tif'; % image format
ImgPara.FirstFrameNum = 1; % the No. of the first frame
ImgPara.LastFrameNum = 2046; % the No. of the end frame
ImgPara.Digits4Enum = 4; % number of digits used for frame enumeration (1-4), now the filenames are arranged from

numFrames = ImgPara.LastFrameNum - ImgPara.FirstFrameNum + 1; % Number of frames in the sequence
filename1 = strcat(ImgPara.ImagePrefix, '%0', num2str(ImgPara.Digits4Enum), 'd', ImgPara.ImageFormat);
totalElementsOutsideBorder = 0;
new=zeros(117,224);
%% average intensity
for i = ImgPara.FirstFrameNum:ImgPara.LastFrameNum
filename = sprintf(filename1, i);
    image = imread(fullfile(Paths.OrigImgPath, filename));
    image = double(image);
    new=image+new;
end
image=new/2070;
image=uint8(image);
%% remove the outside of the prostate
frame = image;

% Step 2: Preprocess the frame (if needed)
% You can apply any necessary preprocessing steps such as filtering, enhancing contrast, or denoising.

% Step 3: Threshold the frame to obtain a binary mask
threshold = 0; % Set an appropriate threshold value
binaryMask = imbinarize(frame, threshold);

% Step 4: Label connected components in the binary mask
labeledImage = bwlabel(binaryMask);

% Step 5: Find the largest connected component (assumed to be the prostate region)
prostateLabel = mode(labeledImage(:)); % Assuming the largest component is the most frequent label

% Step 6: Create a binary mask for the prostate region
prostateMask = (labeledImage == prostateLabel);

% Step 7: Dilate the prostate mask to include a border around the region
borderSize = 5; % Set the desired size of the border
borderMask = imdilate(prostateMask, strel('disk', borderSize));
% Count the number of elements outside the border
    elementsOutsideBorder = sum(borderMask(:) == 0);

    totalElementsOutsideBorder = totalElementsOutsideBorder + elementsOutsideBorder;
disp(['Total elements outside the border for the image sequence: ' num2str(totalElementsOutsideBorder)]);
%% remove the elements 
zero_indices_1 = find(image == 0, totalElementsOutsideBorder);
image(zero_indices_1) = [];
%% ecdf calculation
% Calculate the ECDF
image=double(image);
ecdf_values= ecdf(image);

% Plot the ECDF
figure;
plot(ecdf_values);
grid;
title('Empirical Cumulative Distribution Function (ECDF) Original Image');
xlabel('Values');
ylabel('Probability');

%% Define the optimal BA
% Define the target value
target_value = 0.37;

% Calculate the absolute differences between each element and the target value
differences = abs(ecdf_values - target_value);

% Find the index of the element with the minimum difference
[~, index] = min(differences);

% Retrieve the value and index
closest_value = ecdf_values(index);

fprintf('The optimal BA value is %d\n', index);
