function [MVx, MVy] = Motion_Est(img_test, img_ref, opts)
% ================== SET DEFAULT PARAMETERS ==================

% Set default block size if not specified
if ~isfield(opts,'BlockSize')
    opts.BlockSize = 10; 
end

% Set default search limit if not specified
if ~isfield(opts,'SearchLimit')
    opts.SearchLimit = 20; 
end

% Ensure that BlockSize and SearchLimit are at least 6
BlockSize   = max(opts.BlockSize,   6);
SearchLimit = max(opts.SearchLimit, 6);

% Ensure BlockSize is even (motion estimation works better with even blocks)
if mod(BlockSize,2)~=0
    error('It is better to have BlockSize as an even number.');
end

% Ensure SearchLimit is even
if mod(SearchLimit,2)~=0
    error('It is better to have SearchLimit as an even number.');
end

% ================== PREPROCESSING ==================

% Crop the images so that their dimensions are multiples of BlockSize
M        = floor(size(img_ref, 1)/BlockSize) * BlockSize; % Adjust rows
N        = floor(size(img_ref, 2)/BlockSize) * BlockSize; % Adjust columns
img_ref  = img_ref(1:M, 1:N, :); % Crop reference image (row from 1 to M, column from 1 to N)
img_test = img_test(1:M, 1:N, :); % Crop test image

% Enlarge image boundaries by BlockSize/2 pixels (to prevent boundary issues)
img_ref  = padarray(img_ref,  [BlockSize/2, BlockSize/2], 'replicate'); 
img_test = padarray(img_test, [BlockSize/2, BlockSize/2], 'replicate');

% Pad images with zeros to accommodate search limit range
img_ref  = padarray(img_ref,  [SearchLimit, SearchLimit]); 
img_test = padarray(img_test, [SearchLimit, SearchLimit]); 

% ================== INITIALIZE PARAMETERS ==================

% Get new image dimensions after padding
[M, N, C] = size(img_ref); 

% Define half-block size for indexing
L = floor(BlockSize/2); 

% Define block range (indices relative to block center)
BlockRange = -L:L-1; 

% Define coordinate ranges for block centers (ensuring blocks fit within image)
xc_range = SearchLimit+L+1 : BlockSize : N-(SearchLimit+L); % X-coordinates
yc_range = SearchLimit+L+1 : BlockSize : M-(SearchLimit+L); % Y-coordinates

% Initialize motion vector matrices
MVx = zeros(length(yc_range), length(xc_range)); % Horizontal motion vectors
MVy = zeros(length(yc_range), length(xc_range)); % Vertical motion vectors

% ================== MAIN LOOP FOR BLOCK MATCHING ==================

%disp("length of yc_range")
%disp(length(yc_range))
%disp("length of xc_range")
%disp(length(xc_range))

% Loop over each block's center location
for i = 1:length(yc_range)
    for j = 1:length(xc_range)
        % Get current block's center coordinates
        xc = xc_range(j); 
        yc = yc_range(i);
       
        
        % Extract the block from the test image at the current location
        Block = img_test(yc + BlockRange, xc + BlockRange, :); 
        
        % Perform motion estimation using either Full Search or Log Search
        % Uncomment one of the following lines to choose the method:
        
        % [MVy1, MVx1] = FullSearch(Block, img_ref, xc, yc, SearchLimit); % Full Search method

        % [MVy1, MVx1] = FullSearch_interpolation(Block, img_ref, xc, yc, SearchLimit); % Full Search method with linear interpolation

        % [MVy1, MVx1] = LogSearch(Block, img_ref, xc, yc, SearchLimit); % Logarithmic Search method

        [MVy1, MVx1] = LogSearch_interpolation(Block, img_ref, xc, yc, SearchLimit);

        %[MVy1, MVx1] = ThreeStepSearch(Block, img_ref, xc, yc, SearchLimit);
        
        % Store motion vectors for this block
        MVx(i,j) = MVx1;
        MVy(i,j) = MVy1;
    end
end

% ================== POST-PROCESSING ==================

% Limit motion vectors to be within the search range
MVx(MVx >  SearchLimit) =  SearchLimit;
MVx(MVx < -SearchLimit) = -SearchLimit;
MVy(MVy >  SearchLimit) =  SearchLimit;
MVy(MVy < -SearchLimit) = -SearchLimit;

% Apply a median filter to smooth the motion vectors
MVx = medfilt2(MVx, [3 3]); 
MVy = medfilt2(MVy, [3 3]);
