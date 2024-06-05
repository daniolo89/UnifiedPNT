% Inputs:
%   Y - Received signal (matrix of size N x L), where N is the size of
%   subcarriers and L is the nr. of transmissions
%   w - Vector of size 2
% Output:
%   Z - Processed signal (matrix of size N x L/2)

function Z = match(Y, w)

% Get the dimensions of Y
[N, L] = size(Y);

% Initialize the output matrix Z with zeros
Z = zeros(N, L/2);

% Loop over k from 0 to L/2 - 1
for k = 0:(L/2 - 1)
    % Compute the k-th column of Z
    Z(:, k+1) = w(1) * Y(:, 2*k + 1) + w(2) * Y(:, 2*k + 2);
end
end
