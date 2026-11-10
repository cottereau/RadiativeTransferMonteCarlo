function prediction = fonteModalRelaxation(bands,materialMetadata,time)
%FONTEMODALRELAXATION Exact all-space two-mode ODE, band by band.
% Sigma(1,2) is P->S; Sigma(2,1) is S->P, both in inverse seconds.
% See Fonte/main.tex eq:ode. No attenuation, boundaries or receiver window.
required={'energyP','energyS','eventTime'};
if ~isstruct(bands) || ~isscalar(bands) || ~all(isfield(bands,required))
    error('Fonte:RelaxationInput','Physical modal band energies are required.');
end
p=double(bands.energyP(:)); s=double(bands.energyS(:)); n=numel(p);
if numel(s)~=n || numel(materialMetadata)~=n || ...
        any(~isfinite([p;s])) || any([p;s]<0)
    error('Fonte:RelaxationInput','Band energies or material metadata are invalid.');
end
validateattributes(time,{'numeric'},{'real','finite','vector','nonempty'});
t=double(time(:)); if any(diff(t)<=0), error('Fonte:RelaxationTime','Times must increase.'); end
validateattributes(bands.eventTime,{'numeric'},{'real','finite','scalar'});
age=max(0,t-bands.eventTime); emitted=double(t>=bands.eventTime);
energy=zeros(numel(t),2); individual=zeros(numel(t),2,n);
for j=1:n
    if p(j)+s(j)==0, continue; end
    if ~isfield(materialMetadata(j),'Sigma') || ...
            ~isequal(size(materialMetadata(j).Sigma),[2 2])
        error('Fonte:RelaxationRates','Each active band needs a prepared 2x2 Sigma.');
    end
    sig=double(materialMetadata(j).Sigma);
    if any(~isfinite(sig(:))) || any(sig(:)<0)
        error('Fonte:RelaxationRates','Scattering rates must be finite and nonnegative.');
    end
    ps=sig(1,2); sp=sig(2,1); total=p(j)+s(j);
    if ps+sp==0
        ep=p(j)*ones(size(t)); es=s(j)*ones(size(t));
    else
        epEq=total*sp/(ps+sp);
        ep=epEq+(p(j)-epEq)*exp(-(ps+sp)*age);
        es=total-ep;
    end
    individual(:,:,j)=[ep.*emitted es.*emitted];
    energy=energy+individual(:,:,j);
end
prediction=struct('t',t,'eventTime',bands.eventTime,'energy',energy, ...
    'bandEnergy',individual,'totalEnergy',sum(energy,2), ...
    'ratioSP',energy(:,2)./max(energy(:,1),realmin), ...
    'energyUnits','J','model','all-space two-state conversion ODE');
end
