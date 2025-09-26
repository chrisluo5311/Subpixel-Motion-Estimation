function [MVy, MVx] = LogSearch(Block, img_ref, xc, yc, SearchLimit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log Search Algorithm for Motion Estimation
%
% [MVy, MVx] = LogSearch(Block, img_ref, xc, yc, SearchLimit)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ================== INITIALIZATION ==================

% Get image dimensions (M: rows, N: columns, C: color channels)
[M, N, C] = size(img_ref); 

% Get block size from the first dimension of the Block
BlockSize = size(Block, 1); 

% Define half block size for indexing
L = floor(BlockSize/2); 

% Define block index range relative to the center
BlockRange = -L:L-1; 

% Initialize minimum Sum of Absolute Differences (SAD) value
SADmin = 1e6; 

% Initialize the best matching block's coordinates
y_min = yc; 
x_min = xc;

% ================== BOUNDARY CHECK ==================

% Check if the block is within the valid search area
if (yc < SearchLimit + L) || (yc > M - (SearchLimit + L))
    error('Please set yc > %3g pixels from the boundary.\n', SearchLimit + L);
end

if (xc < SearchLimit + L) || (xc > N - (SearchLimit + L))
    error('Please set xc > %3g pixels from the boundary.\n', SearchLimit + L);
end

% ================== INITIALIZATION FOR SEARCH ==================

% Starting point for the search
x0 = xc;
y0 = yc;

% Number of levels for logarithmic search
LevelMax = 2; 

% Initialize level-specific search limits
LevelLimit = zeros(1, LevelMax+1);

% ================== LOGARITHMIC SEARCH ==================

% Loop through each level in the logarithmic search
for k = 1:LevelMax
    
    % Define the search step size for the current level
    LevelLimit(k+1) = max(2^(floor(log2(SearchLimit)) - k + 1), 1);
    
    % Generate the search pattern for the current level
    c = 2.^(0:log2(LevelLimit(k+1))); 
    
    % Remove search points that exceed the SearchLimit
    c(c + sum(LevelLimit(1:k)) > SearchLimit) = []; 
    
    % Define the search range in x and y directions
    x_range = zeros(1, 2*length(c)+1); 
    x_range(1) = 0;
    x_range(2:2:2*length(c)) = c;
    x_range(3:2:2*length(c)+1) = -c;
    y_range = x_range;
    
    % ================== SEARCH OVER CANDIDATE BLOCKS ==================
    
    % Loop through each candidate block position
    for i = 1:length(y_range)
        for j = 1:length(x_range)
            
            % Check if a better match is still possible
            if SADmin > 1e-3 
                
                % Calculate candidate block center coordinates
                xt = x0 + x_range(j);
                yt = y0 + y_range(i);
                
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
            else
                % Stop searching if SAD is very small
                SADmin = 0;
                x_min = xc;
                y_min = yc;
            end
        end
    end
    
    % Update the search center to the best match found at this level
    x0 = x_min;
    y0 = y_min;
end

% ================== INTEGER MOTION VECTOR ==================

% Calculate the integer motion vector
MVx_int = xc - x_min;
MVy_int = yc - y_min;

% ================== TAYLOR SERIES REFINEMENT ==================

% Extract the best matching block for subpixel refinement
Block_ref = img_ref(y_min + BlockRange, x_min + BlockRange, :);

% Apply Taylor approximation for subpixel accuracy
Taylor_sol = Taylor_App(Block, Block_ref);

% Calculate fractional motion vectors
MVx_frac = Taylor_sol(1);
MVy_frac = Taylor_sol(2);

% ================== OVERALL MOTION VECTOR ==================

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
