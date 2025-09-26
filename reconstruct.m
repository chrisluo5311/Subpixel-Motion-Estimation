function g = reconstruct(img0, MVx, MVy, pel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integer Pixel Motion Compensation
%
% g = reconstruct(img0, MVx, MVy, pel)
% Constructs a motion-compensated frame of img0 according to the motion
% vectors specified by MVx and MVy.
%
% Inputs:
%   img0 - Reference image (previous frame)
%   MVx  - Horizontal motion vectors
%   MVy  - Vertical motion vectors
%   pel  - Precision level (e.g., 0.5 for half-pixel accuracy)
%
% Output:
%   g - Motion-compensated image
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ================== IMAGE RESCALING ==================

% Rescale the image by the inverse of the precision level
% For example, if pel = 0.5, the image is doubled in size
img = imresize(img0, 1/pel, 'bilinear'); 

% ================== BLOCK SIZE CALCULATION ==================

% Determine block size by dividing image height by number of motion vectors
BlockSize = floor(size(img, 1) / size(MVx, 1)); 

% ================== IMAGE CROPPING ==================

% Get image dimensions (m: rows, n: columns, C: color channels)
[m, n, C] = size(img); 

% Adjust dimensions to be multiples of BlockSize
M = floor(m / BlockSize) * BlockSize; % Rows
N = floor(n / BlockSize) * BlockSize; % Columns

% Crop image to adjusted dimensions to avoid boundary issues
f = img(1:M, 1:N, 1:C); 

% Initialize the output image with zeros
g = zeros(M, N, C); 

% ================== MOTION VECTOR RESIZING ==================

% Resize motion vectors to match the pixel grid of the scaled image
MVxmap = imresize(MVx, BlockSize); 
MVymap = imresize(MVy, BlockSize); 

% Scale and round motion vectors for integer pixel shifts
Dx = round(MVxmap * (1/pel)); % Horizontal displacement
Dy = round(MVymap * (1/pel)); % Vertical displacement

% ================== GRID COORDINATES ==================

% Create a grid of pixel coordinates for the image
[xgrid, ygrid] = meshgrid(1:N, 1:M); 

% Compute the shifted grid coordinates using motion vectors
X = min(max(xgrid + Dx, 1), N); % Ensure X is within image width
Y = min(max(ygrid + Dy, 1), M); % Ensure Y is within image height

% Convert 2D grid coordinates to 1D linear indices
% (X-1)*M + Y converts (x, y) to a linear index for column-major order
idx = (X(:) - 1) * M + Y(:);

% ================== PIXEL MAPPING ==================

% Loop through each color channel (for RGB images)
for coloridx = 1:C
    % Extract the current color channel
    fc = f(:, :, coloridx); 
    
    % Map pixels according to motion vectors and reshape to original size
    g(:, :, coloridx) = reshape(fc(idx), M, N); 
end

% ================== IMAGE RESCALING ==================

% Resize the motion-compensated image back to the original resolution
g = imresize(g, pel); 
