% Clear workspace, close figures, and clear command window
clear all     % Remove all variables from the workspace
close all     % Close all open figures
clc           % Clear the command window
% ================== PARAMETERS ==================

% Define motion estimation parameters
opts.BlockSize   = 8;  % Size of blocks for motion estimation
opts.SearchLimit = 10; % Maximum displacement allowed for motion search

% ================== LOAD IMAGES ==================

% Read and convert the first image to double precision
img0 = im2double(imread('./imgs/foreman001.png')); 

% Read and convert the second image to double precision
img1 = im2double(imread('./imgs/foreman002.png')); 

% ================== MOTION ESTIMATION ==================

tic % Start measuring time

% Perform bidirectional motion estimation 
% This function estimates motion vectors MVx and MVy between img0 and img1
[MVx, MVy] = Bidirectional_ME(img0, img1, opts);

toc % Stop measuring time and display elapsed time

% ================== MOTION COMPENSATION ==================

tic % Start measuring time

% Perform motion compensation to reconstruct the image using estimated motion vectors
% 0.5 is typically a weighting factor for bidirectional estimation
imgMC = reconstruct(img0, MVx, MVy, 0.5);

toc % Stop measuring time and display elapsed time

% ================== EVALUATION ==================

% Get image dimensions (M: rows, N: columns, C: color channels)
[M, N, C] = size(imgMC);

% Compute residual image (difference between motion-compensated image and actual second image)
Res  = imgMC - img1(1:M, 1:N, 1:C);

% Compute Mean Squared Error (MSE)
MSE  = norm(Res(:), 'fro')^2 / numel(imgMC);

% Compute Peak Signal-to-Noise Ratio (PSNR)
PSNR = 10 * log10(max(imgMC(:))^2 / MSE);

% ================== DISPLAY RESULTS ==================

% Plot motion vector field
figure(1);
quiver(MVx(end:-1:1,:), MVy(end:-1:1,:)); % Visualize motion vectors
title('Motion Vector Field');

% Display original, target, and reconstructed images
figure(2);

subplot(221);
imshow(img0); title('img_0'); % Display first image

subplot(222);
imshow(img1); title('img_1'); % Display second image

subplot(223);
imshow(imgMC); title('img_M compensated'); % Display motion-compensated image

subplot(224); 
T = sprintf('img_M - img_1, PSNR %3g dB, MSE %f', PSNR, MSE); % Format title with PSNR value

% Show the difference between reconstructed and actual images (grayscale)
imshow(rgb2gray(imgMC) - rgb2gray(img1(1:M, 1:N, :)), []); 
title(T);
