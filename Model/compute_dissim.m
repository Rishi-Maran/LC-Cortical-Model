function [dissim_time] = compute_dissim(ts1, ts2)
%COMPUTE_DISSIM Summary of this function goes here
%   Detailed explanation goes here
    if ~all(size(ts1) == size(ts2))
        error('time series must be same dimension');
    end
    dissim_time = zeros(1, size(ts1, 2));
    for n = 1:size(ts1, 2)
        vec1 = reshape(ts1(:, n), [], 1)';
        vec2 = reshape(ts2(:, n), [], 1)';
        
        % Check if either vector is all zeros
        if all(vec1 == 0) || all(vec2 == 0)
            dissim_time(n) = 0;
        else
            dissim_time(n) = pdist2(vec1, vec2, 'cosine');
        end
    end
end

