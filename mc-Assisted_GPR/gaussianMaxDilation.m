function B = gaussianMaxDilation(A, sigma, filterSize)
% GAUSSIANMAXDILATION Spreads non-zero points using a Gaussian profile.
%   B = gaussianMaxDilation(A, sigma, filterSize) takes matrix A and 
%   for each non-zero element, creates a Gaussian spread. The final 
%   value at each pixel is the maximum of all generated Gaussian spreads.
%
%   INPUTS:
%       A          - Input sparse-valued matrix
%       sigma      - Standard deviation of the Gaussian
%       filterSize - Size of the neighborhood (should be odd, e.g., 6*sigma)

    % 1. Ensure filterSize is odd for a centered peak
    if mod(filterSize, 2) == 0
        filterSize = filterSize + 1;
    end

    % 2. Create the Gaussian Kernel
    % We use meshgrid to define coordinates relative to the center
    lim = (filterSize - 1) / 2;
    [X, Y] = meshgrid(-lim:lim);
    
    % Calculate Gaussian values
    % Peak is 1.0 at the center (0,0)
    kernel = exp(-(X.^2 + Y.^2) / (2 * sigma^2));
    
    % 3. Convert to an Offset Structuring Element
    % We subtract 1 so the peak is exactly 0. 
    % In dilation: NewValue = max(OriginalValue + Offset)
    % Peak center: Value + 0 = Value (Retention)
    se = offsetstrel(kernel - 1);
    
    % 4. Apply Grayscale Dilation
    % This handles the 'maximum value at each index' automatically
    B = imdilate(A, se);
    
end