% Inputs of form:
%   weighted.data
%   weighted.fa

function A=PDwT1w2A(PDw,T1w,method,relativeB1)

if ~exist('relativeB1','var')
    relativeB1=1;
end

A=PDw; % preserve fields
A.data=0*A.data;
A.fa=[]; % no well-defined flip angle 

switch method
    case {'smallFlipAngle','Helms2008'}
        PDw.t=relativeB1.*PDw.fa;
        T1w.t=relativeB1.*T1w.fa;
    case {'exact','smallTR*R1'}
        PDw.t=2*tan(relativeB1.*PDw.fa/2);
        T1w.t=2*tan(relativeB1.*T1w.fa/2);
    otherwise 
        error('unknown method')
end

% Zero elements with no signal
mask=(T1w.data>0)&(PDw.data>0)&(relativeB1>1e-1);

if length(T1w.t)>1
    T1w.t=T1w.t(mask);
    PDw.t=PDw.t(mask);
end

% Equation (18) in Dathe and Helms, Physics in Medicine and Biology (2010)
A.data(mask)=bsxfun(@times,PDw.data(mask).*T1w.data(mask),T1w.t./PDw.t-PDw.t./T1w.t)...
    ./(bsxfun(@times,T1w.data(mask),T1w.t)-bsxfun(@times,PDw.data(mask),PDw.t));

end