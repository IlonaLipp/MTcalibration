% Inputs of form:
%   MTw.data
%   MTw.fa
% N.B.: there must only be one echo time for T1 and A!

function MTsat=computeMTsatApprox_using_exact_FAs(A,R1,MTw,TR,percB1)

MTsat=MTw; % preserve fields

% singular dimension expansion as T1.data and MTW.fa should have only one "echo".
% Solution from Helms, et al., MRM (2008), Equation (8).
% use actual flip angles
ft = percB1 / 100;
MTsat.data=bsxfun(@minus,bsxfun(@rdivide,bsxfun(@minus,A.data.*(MTw.fa*ft),MTw.data),getfield(R1toT1(R1),'data')/TR)./MTw.data,(MTw.fa*ft).^2/2);

MTsat.data=MTsat.data*100; % convert to percentage units

end
