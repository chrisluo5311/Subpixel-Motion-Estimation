clc; clear; close all;

% 1: Load image 
img0 = im2double(imread('./imgs/foreman001.png')); 
img1 = im2double(imread('./imgs/foreman002.png')); 


% 2:  Motion Blur PSF
motion_length = 10; %  the blur extends 20 pixels
motion_angle = 30;  % the blur is diagonal (30°)
% generates a motion blur kernel.
PSF = fspecial('motion', motion_length, motion_angle); 

% Step 3: Apply motion blur using convolution
% convolve the image with the PSF.
img1 = imfilter(img1, PSF, 'conv', 'circular'); 
% figure, imshow(img1), title('Blurred Image with Motion PSF');

% ================== PARAMETERS ==================

% Define motion estimation parameters
opts.BlockSize   = 8;  % Size of blocks for motion estimation
opts.SearchLimit = 10; % Maximum displacement allowed for motion search

% ================== MOTION ESTIMATION ==================

tic % Start measuring time

% Perform bidirectional motion estimation 
% This function estimates motion vectors MVx and MVy between img0 and img1
[MVx, MVy] = Bidirectional_ME(img0, img1, opts);

toc % Stop measuring time and display elapsed time

% **Fix 1: Apply Gaussian Smoothing to Motion Vectors**
MVx = imgaussfilt(MVx, 2.5);
MVy = imgaussfilt(MVy, 2.5);
% ================== ESTIMATE POINT SPREAD FUNCTION (PSF) ==================
PSF_est = estimate_psf(MVx, MVy, 15);

% ================== APPLY MOTION DEBLURRING (Using deconvreg) ==================
deblurred_img = deconvreg(img1, PSF_est,0.03);
%deblurred_img = apply_deblurring(img1, PSF_est);
%deblurred_img = deconvlucy(img1, PSF_est, 10);

% ================== DISPLAY RESULTS ==================
img2 = im2double(imread('./imgs/foreman002.png')); 
figure;
subplot(2,3,1), imshow(img0),title("img0 original")
subplot(2,3,2), imshow(img2), title("img1")
subplot(2,3,4), imshow(img1), title('Blurred Image');
subplot(2,3,5), quiver(MVx, MVy), title('Estimated Motion Vectors');
subplot(2,3,6), imshow(deblurred_img), title('Deblurred Image');

%% Function: Estimate Point Spread Function (PSF)
function psf = estimate_psf(MVx, MVy, psf_size)
    % Initialize the PSF matrix
    psf = zeros(psf_size, psf_size);
    center = floor(psf_size / 2) + 1;

    % Iterate over each motion vector
    for i = 1:size(MVx, 1)
        for j = 1:size(MVx, 2)
            dx = MVx(i, j);
            dy = MVy(i, j);

            % Compute the subpixel positions
            x = center + dx;
            y = center + dy;

            % Determine the integer and fractional parts
            x_int = floor(x);
            y_int = floor(y);
            x_frac = x - x_int;
            y_frac = y - y_int;

            % Distribute the PSF energy to the neighboring pixels
            if x_int >= 1 && x_int < psf_size && y_int >= 1 && y_int < psf_size
                psf(x_int, y_int) = psf(x_int, y_int) + (1 - x_frac) * (1 - y_frac);
                if x_int + 1 <= psf_size
                    psf(x_int + 1, y_int) = psf(x_int + 1, y_int) + x_frac * (1 - y_frac);
                end
                if y_int + 1 <= psf_size
                    psf(x_int, y_int + 1) = psf(x_int, y_int + 1) + (1 - x_frac) * y_frac;
                end
                if x_int + 1 <= psf_size && y_int + 1 <= psf_size
                    psf(x_int + 1, y_int + 1) = psf(x_int + 1, y_int + 1) + x_frac * y_frac;
                end
            end
        end
    end

    % Normalize the PSF to sum to 1
    psf = psf / sum(psf(:));
end



%% Function: Apply Wiener Deconvolution for Deblurring
function deblurred_img = apply_deblurring(blurred_img, psf)
    blurred_fft = fft2(blurred_img);
    psf_fft = fft2(psf, size(blurred_img, 1), size(blurred_img, 2));
    deblurred_fft = blurred_fft ./ (psf_fft + 0.01); % Avoid division by zero
    deblurred_img = real(ifft2(deblurred_fft));
    deblurred_img = imadjust(deblurred_img, stretchlim(deblurred_img), []);
end