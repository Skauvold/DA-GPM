function Hx = gr_obs_model(x,varargin)
% GR_OBS_MODEL computes thickness and synthetic/expected gamma ray
% observations based on the state variables x.

% Assume x is 4-by-ne matrix where top row is z and the bottom 3 rows are
% s1, s2 and s3

% Check whether input has correct dimensions
if size(x,1) ~= 4
    error('x must be a 4-by-ne matrix.');
end

% Use default values if none are specified
if nargin==2
    grvals = varargin{1};
else
    grvals = [0.05, 0.15, 0.65, 0.95];
end

% Separate out elevations
z = x(1,:);

% Compute proportions
props = eta2p(x(2:4,:));
props = [props;ones(1,size(props,2))-sum(props,1)];

% Compute synthetic gamma ray measurements
gr = grvals*props;

% Combine elevations and gamma ray measurements
Hx = [z;gr];

end



