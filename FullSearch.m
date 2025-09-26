function [MVy, MVx] = FullSearch(Block, img_ref, xc, yc, SearchLimit)
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

%% ================== TAYLOR SERIES REFINEMENT ==================

% Extract the best matching block for subpixel refinement
Block_ref = img_ref(y_min + BlockRange, x_min + BlockRange, :);

% Apply Taylor approximation for subpixel accuracy
Taylor_sol = Taylor_App(Block, Block_ref);

%

% Calculate fractional motion vectors
MVx_frac = Taylor_sol(1);
MVy_frac = Taylor_sol(2);

%% ================== OVERALL MOTION VECTOR ==================

% Combine integer and fractional parts
MVx = MVx_int + MVx_frac;
MVy = MVy_int + MVy_frac;
end

% Taylor Refinement
function x = Taylor_App(f, g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taylor Refinement
% 
% This function computes the motion vector using Taylor series
% approximation.
% f(x + dx, y + dy) ~= f(x,y) + dx df/dx + dy df/dy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dfx dfy] = gradient(f);

a = sum(dfx(:).^2);
b = sum(dfx(:).*dfy(:));
d = sum(dfy(:).^2);

z = g-f;
p = sum(z(:).*dfx(:));
q = sum(z(:).*dfy(:));

A = [a b; b d];
rhs = [p;q];

if cond(A)>1e6
    x = [0 0]';
else
    x = A\rhs;
end
end