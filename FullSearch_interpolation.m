function [MVy, MVx] = FullSearch_interpolation(Block, img_ref, xc, yc, SearchLimit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full Search Motion Estimation
%
% [MVy, MVx] = FullSearch(Block, img_ref, xc, yc, SearchLimit)
% Finds the motion vector of the (yc, xc)-th Block in the reference image.
%
% Input:
%   Block       - Current block being searched (from test image)
%   img_ref     - Reference image
%   xc, yc      - Center coordinates of the current block
%   SearchLimit - Maximum search range for motion estimation
%
% Output:
%   MVy, MVx    - Estimated motion vectors (vertical and horizontal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ================== PARAMETERS ==================

% Get image dimensions (M: rows, N: columns, C: color channels)
[M, N, C] = size(img_ref); 

% Get block size from the first dimension of the Block
BlockSize = size(Block, 1); 

% Define half block size for indexing
L = floor(BlockSize / 2); 

% Define block index range relative to the center
BlockRange = -L:L-1; 

% Set the search range from -SearchLimit to +SearchLimit
SearchRange = SearchLimit; 

% Initialize minimum Sum of Absolute Differences (SAD) value
SADmin = 1e6; 

%% ================== BOUNDARY CHECK ==================

% Check if the block is within the valid search area vertically
if (yc <= SearchRange + L) || (yc >= M - (SearchRange + L))
    error('Please set yc > %3g pixels from the boundary.\n', SearchRange + L);
end

% Check if the block is within the valid search area horizontally
if (xc <= SearchRange + L) || (xc >= N - (SearchRange + L))
    error('Please set xc > %3g pixels from the boundary.\n', SearchRange + L);
end

%% ================== FULL SEARCH LOOP ==================

% Loop through all possible displacements within the search range
for i = -SearchRange:SearchRange
    for j = -SearchRange:SearchRange
        
        % Calculate candidate block center coordinates
        xt = xc + j;
        yt = yc + i;
        
        % Extract candidate block from reference image
        Block_ref = img_ref(yt + BlockRange, xt + BlockRange, :);
        
        % Calculate Sum of Absolute Differences (SAD)
        SAD = sum(abs(Block(:) - Block_ref(:))) / (BlockSize^2);
        
        % Update minimum SAD and best matching coordinates
        if SAD < SADmin
            SADmin = SAD;
            x_min = xt;
            y_min = yt;
        end
        
        % ================== INTEGER MOTION VECTOR ==================
        
        % Calculate the integer motion vector
        MVx_int = xc - x_min;
        MVy_int = yc - y_min;
    end
end

%% ================== LINEAR INTERPOLATION REFINEMENT (Subpixel Motion) ==================
[MVx_frac, MVy_frac] = BilinearInterpolation(img_ref, x_min, y_min);

%% ================== OVERALL MOTION VECTOR ==================
MVx = MVx_int + MVx_frac;
MVy = MVy_int + MVy_frac;
end


%% ================== BILINEAR INTERPOLATION FUNCTION ==================
function [MVx_frac, MVy_frac] = BilinearInterpolation(img_ref, x_min, y_min)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilinear Interpolation for Subpixel Motion Estimation
% Computes the fractional motion vector using interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure valid indices for interpolation
[M, N] = size(img_ref);
x_min = max(2, min(N-1, x_min)); % Prevents out-of-bounds errors
y_min = max(2, min(M-1, y_min));

% Get pixel intensities around (x_min, y_min)
I00 = double(img_ref(y_min, x_min)); % Top-left pixel
I10 = double(img_ref(y_min, x_min+1)); % Top-right pixel
I01 = double(img_ref(y_min+1, x_min)); % Bottom-left pixel
I11 = double(img_ref(y_min+1, x_min+1)); % Bottom-right pixel

% Compute interpolated subpixel motion
MVx_frac = ((I10 - I00) + (I11 - I01)) / (2 * (abs(I10 - I00) + abs(I11 - I01) + 1e-6));
MVy_frac = ((I01 - I00) + (I11 - I10)) / (2 * (abs(I01 - I00) + abs(I11 - I10) + 1e-6));

end