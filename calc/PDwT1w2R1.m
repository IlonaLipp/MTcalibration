% Inputs of form:
%   weighted.data
%   weighted.fa

function R1=PDwT1w2R1(PDw,T1w,TR,method,relativeB1)

if ~exist('relativeB1','var')
    relativeB1=1;
end

R1=T1w; % preserve fields
R1.data=0*R1.data;
R1.fa=[]; % no well-defined flip angle

switch method
    case 'Helms2008'
        PDw.t=PDw.fa;
        T1w.t=T1w.fa;
    case 'smallFlipAngle'
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


% Equation (17) in Dathe and Helms, Physics in Medicine and Biology (2010)
rho=0.5*(bsxfun(@times,T1w.data(mask),T1w.t)-bsxfun(@times,PDw.data(mask),PDw.t))...
    ./(bsxfun(@rdivide,PDw.data(mask),PDw.t)-bsxfun(@rdivide,T1w.data(mask),T1w.t));

% Remove values which can't give sensible results
rho(isnan(rho))=0;
rho(rho<0)=0;

switch method
    case 'Helms2008' 
        if length(relativeB1)>1, relativeB1=relativeB1(mask); end
        R1.data(mask)=relativeB1.^2.*rho/TR;
    case 'smallTR*R1'
        R1.data(mask)=rho/TR;
    case {'exact','smallFlipAngle'}
        % Equation (11) in Dathe and Helms, Physics in Medicine and Biology (2010)
        rho(rho>=2)=0; % Remove values which can't give sensible results
        R1.data(mask)=(2/TR)*atanh(rho/2);
    otherwise 
        error('unknown method')
end

end