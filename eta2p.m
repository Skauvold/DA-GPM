function p = eta2p(eta)
% ETA2P: convert real-valued eta vector
% to proportion vector p

% Unpack input
eta1 = eta(1,:);
eta2 = eta(2,:);
eta3 = eta(3,:);

% Apply inverse logistic transformation
sv = ones(1,size(eta,2)) + sum(exp(eta));

% Normalise
p1 = exp(eta1)./sv;
p2 = exp(eta2)./sv;
p3 = exp(eta3)./sv;

% Assemble output
p = [p1;p2;p3];

end