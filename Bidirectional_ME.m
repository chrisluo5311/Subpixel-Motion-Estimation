function [MVx, MVy] = Bidirectional_ME(img0, img1, opts)
% Bidirectional Motion Estimation
% Estimates motion vectors in both forward and backward directions
%
% Inputs:
%   img0 - Reference image (previous frame)
%   img1 - Target image (current frame)
%   opts - Options structure containing motion estimation parameters
%
% Outputs:
%   MVx - Estimated horizontal motion vectors
%   MVy - Estimated vertical motion vectors

% ================== FORWARD MOTION ESTIMATION ==================

% Compute motion estimation from img0 to img1 (forward motion)
[MVx1, MVy1] = Motion_Est(img0, img1, opts);

% ================== BACKWARD MOTION ESTIMATION ==================

% Compute motion estimation from img1 to img0 (backward motion)
[MVx2, MVy2] = Motion_Est(img1, img0, opts);

% ================== MOTION REFINEMENT ==================

% Store forward motion vectors in the first channel of MVx
MVx(:,:,1) =  MVx1;  

% Store the negation of backward motion vectors in the second channel
MVx(:,:,2) = -MVx2;  

% Store forward motion vectors in the first channel of MVy
MVy(:,:,1) =  MVy1;  

% Store the negation of backward motion vectors in the second channel
MVy(:,:,2) = -MVy2;  

% Take the maximum value across the two motion estimations for final motion vectors
MVx = max(MVx, [], 3); % Select the dominant horizontal motion vector
MVy = max(MVy, [], 3); % Select the dominant vertical motion vector
