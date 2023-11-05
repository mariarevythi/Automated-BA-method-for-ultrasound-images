function Optimal_BA = computeOptimalBA(filename1, Paths, ImgPara,I)
totalElementsOutsideBorder = 0;
numFrames = ImgPara.LastFrameNum - ImgPara.FirstFrameNum + 1; % Number of frames in the sequence
[y, x] = size(I(:,:,1));
new=zeros(y,x);
%% average intensity
for i = ImgPara.FirstFrameNum:ImgPara.LastFrameNum
filename = sprintf(filename1, i);
    image = imread(fullfile(Paths.OrigImgPath, filename));
    image = double(image);
    new=image+new;
end
image=new/numFrames;
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
%disp(['Total elements outside the border for the image sequence: ' num2str(totalElementsOutsideBorder)]);
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
xlabel('BA Values');
ylabel('ECDF values');

%% Define the optimal BA
% Define the target value
target_value = 0.37;

% Calculate the absolute differences between each element and the target value
differences = abs(ecdf_values - target_value);

% Find the index of the element with the minimum difference
[~, Optimal_BA] = min(differences);

% Retrieve the value and index
closest_value = ecdf_values(Optimal_BA);

fprintf('The optimal BA value is %d\n', Optimal_BA);
end 