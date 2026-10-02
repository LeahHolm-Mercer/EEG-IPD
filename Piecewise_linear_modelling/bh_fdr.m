function q = bh_fdr(p)
%BH_FDR Benjamini-Hochberg false discovery rate adjustment.
%
%   q = BH_FDR(p) returns Benjamini-Hochberg FDR-adjusted p-values
%   corresponding to the input p-values.
%
%   NaN p-values are excluded from the multiple-testing correction and
%   remain NaN in the output. The number of tests is therefore the number
%   of non-NaN p-values.
%
%   Input:
%       p - vector of raw p-values
%
%   Output:
%       q - vector of BH FDR-adjusted p-values, in the same order and
%           shape as p
%
%   Reference:
%       Benjamini Y, Hochberg Y. Controlling the false discovery rate:
%       a practical and powerful approach to multiple testing.
%       Journal of the Royal Statistical Society: Series B.
%       1995;57(1):289-300.

    % Preserve the original shape
    originalSize = size(p);

    % Convert to a column vector for processing
    p = p(:);

    % Identify valid p-values
    valid = ~isnan(p);
    pValid = p(valid);

    % Initialise output
    q = NaN(size(p));

    % Nothing to correct
    m = numel(pValid);
    if m == 0
        q = reshape(q, originalSize);
        return;
    end

    % Sort valid p-values
    [pSorted, order] = sort(pValid);

    % Benjamini-Hochberg adjustment
    qSorted = pSorted .* m ./ (1:m)';

    % Enforce monotonicity
    for i = m-1:-1:1
        qSorted(i) = min(qSorted(i), qSorted(i+1));
    end

    % Adjusted p-values cannot exceed 1
    qSorted = min(qSorted, 1);

    % Restore original ordering
    qValid = NaN(m,1);
    qValid(order) = qSorted;

    q(valid) = qValid;

    % Restore original input shape
    q = reshape(q, originalSize);
end