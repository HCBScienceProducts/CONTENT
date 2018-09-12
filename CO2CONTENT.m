function out=CO2CONTENT(gtime,lat,lon,pres,temp,psal,doxy,epres,etemp,epsal,edoxy)
% function out=CO2CONTENT(gtime,lat,lon,pres,temp,psal,doxy,epres,etemp,epsal,edoxy)
% 
% CONTENT calculation according to Bittig et al. (2018):
% Combination of CANYON-B Bayesian neural network mappings of AT, CT, pH,
% and pCO2 made consistent with carbonate chemistry constraints 
% for any set of water column P, T, S, O2, location data as an alternative
% to (spatial) climatological interpolation. Based on GLODAPv2/GO-SHIP
% bottle data.
%
% Please cite the paper when using CONTENT/CANYON-B.
%
% input:
% gtime - date (UTC) as matlab time (fractional days; 01-Jan-0000 00:00:00 = 1)
% lat   - latitude / °N  [-90 90]
% lon   - longitude / °E [-180 180] or [0 360]
% pres  - pressure / dbar
% temp  - in-situ temperature / °C
% psal  - salinity
% doxy  - dissolved oxygen / umol kg-1 (!)
% e(var)- input error for pres, temp, psal, doxy
%
% default values for input errors:
% epres:   0.5 dbar
% etemp:   0.005 °C
% epsal:   0.005 psu
% edoxy:   1 % of doxy (<- This is a rather optimistic default, meant for
% GO-SHIP quality bottle data. A reasonable default sensor doxy error would
% be 3 % of doxy (or higher!).)
%
% output:
% out.(var)           - estimate of variable (same size as 'pres' input)
% out.(var)_sigma     - 1 sigma of variable (same size as 'pres' input)
% out.(var)_sigma_min - 1 sigma_min of variable (same size as 'pres' input)
% out.(var)_raw       - raw CANYON and 3 indirect estimates of variable (size numel(pres)x4)
% out.ncruise         - number of cruises within...[see 1..6 below] (size numel(pres)x6)
% out.nstation        - number of stations within...[see 1..6 below] (size numel(pres)x6)
%                       1: +/-10° and +/-15 d window (doy-based)
%                       2: +/-10° and +/-30 d window (doy-based)
%                       3: +/-10° and +/-45 d window (doy-based)
%                       4: +/-20° and +/-15 d window (doy-based)
%                       5: +/-20° and +/-30 d window (doy-based)
%                       6: +/-20° and +/-45 d window (doy-based)
%
% check values for 
% 09-Dec-2014 08:45, 17.6° N, -24.3° E, 180 dbar, 16 °C, 36.1 psu, 104 umol O2 kg-1:
%                (var)    (var)_sigma  (var)_sigma_min      unit
%         AT:    2357.817     10.215        7.464      umol kg-1
%         CT:    2199.472      9.811        7.281      umol kg-1
%         pH:    7.870137      0.021367     0.016767   (insitu total scale, comparable to "calculated pH values" as described in Carter et al., 2018)
%         pCO2:  639.8477     34.1077      26.8143          uatm
%
% examples:
%out=CO2CONTENT(datenum([2014 12 09 08 45 00]),17.6,-24.3,180,16,36.1,104);  % default input errors
%out=CO2CONTENT(datenum([2014 12 09 08 45 00]),17.6,-24.3,180,16,36.1,104,[],[],[],0.03*104); % 3 % input error on doxy, otherwise default input errors
%out=CO2CONTENT(datenum([2014 12 09 08 45 00])*[1 1],17.6*[1 1],-24.3*[1 1],180*[1 1],16*[1 1],36.1*[1 1],104*[1 1]); % more than one single input; default input errors
%
%
% references:
%
% - CONTENT and CANYON-B method:
% Bittig et al. (2018). An alternative to static climatologies: Robust
% estimation of open ocean CO2 variables and nutrient concentrations from
% T, S and O2 data using Bayesian neural networks. Front. Mar. Sci. 5:328. doi:
% 10.3389/fmars.2018.00328 
%
% requires CO2SYS-matlab:
% van Heuven et al. (2011). MATLAB Program Developed for CO2 System
% Calculations. ORNL/CDIAC-105b. Carbon Dioxide Information Analysis
% Center, Oak Ridge National Laboratory, US Department of Energy,  Oak
% Ridge, Tennessee. http://dx.doi.org/10.3334/CDIAC/otg.CO2SYS_MATLAB_v1.1
%
% requires CO2SYS-matlab v2 extensions derivnum and errors:
% Orr et al. (2018). Routine uncertainty propagation for the marine carbon
% dioxide system. Mar. Chem. subm.
% https://github.com/jamesorr/CO2SYS-MATLAB
%
% Henry Bittig, LOV
% v0.8, 02.03.2018, included correlations between carbon parameters;  uses CANYON-B
% v1.0, 11.09.2018, initial publication

%  The software is provided "as is", without warranty of any kind, 
%  express or implied, including but not limited to the warranties of
%  merchantability, fitness for a particular purpose and noninfringement.
%  In no event shall the authors or copyright holders be liable for any
%  claim, damages or other liability, whether in an action of contract,
%  tort or otherwise, arising from, out of or in connection with the
%  software or the use or other dealings in the software.      

% No input checks! Assumes informed use, e.g., same dimensions for all
% inputs, ...

% GLODAPv2 training data density 'GLODAPtraining_coverage_CONTENT_v4.mat'
% can be created from GLODAPv2 files and
% map_coverage_of_GLODAP_training_v4.m (takes a while..);

%% Nothing below here should need to be changed by the user %%

if nargin<1, gtime=datenum([2014 12 09 08 45 00])*[1;1];lat=17.6*[1;1];lon=-24.3*[1;1];pres=180*[1;1];temp=16*[1;1];psal=36.1*[1;1];doxy=104*[1;1]; end % test values
if nargin< 8 || isempty(epres), epres=.5; end
if nargin< 9 || isempty(etemp), etemp=0.005; end
if nargin<10 || isempty(epsal), epsal=0.005; end
if nargin<11 || isempty(edoxy), edoxy=0.01*doxy(:); end

% No input checks! Assumes informed use, e.g., same dimensions for all
% inputs, ...; only eX can either be scalar or same dimension as inputs

% get number of elements
nol=numel(pres);
% and expand input errors if needed
if length(epres)==1, epres=repmat(epres,nol,1); end
if length(etemp)==1, etemp=repmat(etemp,nol,1); end
if length(epsal)==1, epsal=repmat(epsal,nol,1); end
if length(edoxy)==1, edoxy=repmat(edoxy,nol,1); end
% calculate (pure) carbonate system estimates (#1) and include estimates of
% nutrients from CANYON-B

% add CANYON-B carbonate estimates as well as phosphate and silicate
paramnames={'AT';'CT';'pH';'pCO2';'SiOH4';'PO4'};
cy=CANYONB_private(gtime(:),lat(:),lon(:),pres(:),temp(:),psal(:),doxy(:),paramnames,epres(:),etemp(:),epsal(:),edoxy(:));

% copy to raw output and start mixing calculations
paramnames={'AT';'CT';'pH';'pCO2'};
for i=1:4, 
    raw.(paramnames{i})=ones(nol,4)*NaN; 
    raw.(paramnames{i})(:,1)=cy.(paramnames{i});
end

% define arrangements of calculations
inpar=[4 1; 3 1; 4 3; 4 2; 3 2; 1 2]; % (1): 4,1; (2): 3,1; (3): 4,3; (4): 4,2; (5): 3,2; (6): 1,2
outpar=ones(6,2)*NaN;for p=1:6, outpar(p,:)=setdiff(1:4,inpar(p,:)); end
svi=[2 4; 3 3; 4 4; 2 3; 3 4; 2 2]; % save indices for outpar

% deal with weights / uncertainties
for i=1:4, sigma.(paramnames{i})=ones(nol,4)*NaN; end
% nutrient uncertainty (var(inputs) + var(training data) + var(MLP))
cy.eSiOH4=cy.SiOH4_ci;cy.ePO4=cy.PO4_ci;
% covariance matrix of CANYON-B AT, CT, pH, pCO2 due to common inputs
cycov=ones(4,4,nol)*NaN;
for i=1:4
    for j=i:4
        cycov(i,j,:)=sum(cy.inx.(paramnames{i}).*cy.inx.(paramnames{j}).*[etemp(:) epsal(:) edoxy(:) epres(:)].^2,2);cycov(j,i,:)=cycov(i,j,:);
    end
end

% correlation matrix of CANYON-B AT, CT, pH, pCO2 due to common inputs and
% CANYON-B estimation uncertainty
cyr=zeros(size(cycov));
% and full variance on diagonal: 
% var(inputs)[local] + var(training data)[global] + var(MLP)[local]
sigma.AT(:,1)=cy.AT_ci; cycov(1,1,:)=cy.AT_ci.^2;
sigma.CT(:,1)=cy.CT_ci; cycov(2,2,:)=cy.CT_ci.^2;
sigma.pH(:,1)=cy.pH_ci; cycov(3,3,:)=cy.pH_ci.^2;
sigma.pCO2(:,1)=cy.pCO2_ci; cycov(4,4,:)=cy.pCO2_ci.^2;
% correlation matrix with full variance 
for i=1:nol, cyr(:,:,i)=cycov(:,:,i)./(ones(4,1)*sqrt(diag(cycov(:,:,i)))')./(sqrt(diag(cycov(:,:,i)))*ones(1,4)); end

% calculate derivatives of carbonate system calculations (#2-#4) and store
% dx/dy's separately for each calculation (1)-(6)
dcout=ones(4,4,2,nol)*NaN;dctit=cell(4,4,2);
for p=1:6
[deriv, dheaders]=derivnum('par1',cy.(paramnames{inpar(p,1)}),cy.(paramnames{inpar(p,2)}),inpar(p,1),inpar(p,2),psal(:),temp(:),NaN,pres(:),NaN,cy.SiOH4,cy.PO4,1,10,1); % (1): 4,1
dcout(inpar(p,1),inpar(p,2),1,:)=deriv(:,1); dcout(inpar(p,1),inpar(p,2),2,:)=deriv(:,2); % dTCO2/dpCO2, dHin/dpCO2
dctit(inpar(p,1),inpar(p,2),:)= dheaders(1:2);
[deriv, dheaders]=derivnum('par2',cy.(paramnames{inpar(p,1)}),cy.(paramnames{inpar(p,2)}),inpar(p,1),inpar(p,2),psal(:),temp(:),NaN,pres(:),NaN,cy.SiOH4,cy.PO4,1,10,1); % (1): 4,1
dcout(inpar(p,2),inpar(p,1),1,:)=deriv(:,1); dcout(inpar(p,2),inpar(p,1),2,:)=deriv(:,2); % dTCO2/dALK, dHin/dALK
dctit(inpar(p,2),inpar(p,1),:)= dheaders(1:2);
clear deriv dheaders dunits
end


% do carbonate system calculations (#2-#4)
% using the K1K2 constants of Lueker et al, 2000, KSO4 of Dickson 1990 & TB of Uppstrom 1979
% add localized error calculations incl. parameter correlation due to common inputs
for p=1:6
outcalc=CO2SYS(cy.(paramnames{inpar(p,1)}),cy.(paramnames{inpar(p,2)}),inpar(p,1),inpar(p,2),psal(:),temp(:),NaN,pres(:),NaN,cy.SiOH4,cy.PO4,1,10,1); % (1): 4,1
outerr= real( errors(cy.(paramnames{inpar(p,1)}),cy.(paramnames{inpar(p,2)}),inpar(p,1),inpar(p,2),psal(:),temp(:),NaN,pres(:),NaN,cy.SiOH4,cy.PO4,...
  sqrt(cycov(inpar(p,1),inpar(p,1),:)),sqrt(cycov(inpar(p,2),inpar(p,2),:)),epsal,etemp,cy.eSiOH4,cy.ePO4,'',0,cyr(inpar(p,1),inpar(p,2),:),1,10,1) ) ; % (1): 4,1
raw.(paramnames{outpar(p,1)})(:,svi(p,1))=outcalc(:,outpar(p,1)); % outpar 1
raw.(paramnames{outpar(p,2)})(:,svi(p,2))=outcalc(:,outpar(p,2)); % outpar 2
sigma.(paramnames{outpar(p,1)})(:,svi(p,1))=outerr(:,1); % eoutpar 1
sigma.(paramnames{outpar(p,2)})(:,svi(p,2))=outerr(:,2); % eoutpar 2
for i=1:2 % check if one of the two calculated variables is pH
    if outpar(p,i)==3, 
        sigma.pH(:,svi(p,i))=(sigma.pH(:,svi(p,i))*1e-9)./(log(10).*10.^(-outcalc(:,3))); % eH -> epH
        % dHin/dz -> dpHin/dz
        dcout(inpar(p,1),inpar(p,2),i,:)=(squeeze(dcout(inpar(p,1),inpar(p,2),i,:))*1e-9)./(-log(10).*10.^(-outcalc(:,3))); % first variable inpar(p,1), doutpar(p,i)/dinpar(p,1)
        dctit(inpar(p,1),inpar(p,2),i)= strrep(dctit(inpar(p,1),inpar(p,2),i),'dHin','dpHin');
        dcout(inpar(p,2),inpar(p,1),i,:)=(squeeze(dcout(inpar(p,2),inpar(p,1),i,:))*1e-9)./(-log(10).*10.^(-outcalc(:,3))); % second variable inpar(p,2), doutpar(p,i)/dinpar(p,2)
        dctit(inpar(p,2),inpar(p,1),i)= strrep(dctit(inpar(p,2),inpar(p,1),i),'dHin','dpHin');
    end
end
if inpar(p,1)==3, % check if first input variable is pH 
    % dz/dH -> dz/dpH
	for i=1:2 % first and second variable doutpar(p,:) / dinpar(p,1)
        dcout(inpar(p,1),inpar(p,2),i,:)=(squeeze(dcout(inpar(p,1),inpar(p,2),i,:))/1e-9).*(-log(10).*10.^(-outcalc(:,3))); 
        dctit(inpar(p,1),inpar(p,2),i)= strrep(dctit(inpar(p,1),inpar(p,2),i),'dH','dpH'); 
	end
elseif inpar(p,2)==3, % check if second input variable is pH 
    % dz/dH -> dz/dpH
    for i=1:2 % first and second variable doutpar(p,:) / dinpar(p,2)
        dcout(inpar(p,2),inpar(p,1),i,:)=(squeeze(dcout(inpar(p,2),inpar(p,1),i,:))/1e-9).*(-log(10).*10.^(-outcalc(:,3))); 
        dctit(inpar(p,2),inpar(p,1),i)= strrep(dctit(inpar(p,2),inpar(p,1),i),'dH','dpH'); 
    end
end
clear outcalc outerr headers
end

% build weighted mean covariance matrix for all parameters
for k=1:4
    cocov.(paramnames{k})=ones(4,4,nol)*NaN;
    for i=1:4, cocov.(paramnames{k})(i,i,:)=sigma.(paramnames{k})(:,i).^2; end
    for p=1:6
    if any(outpar(p,:)==k) % find calc (p) that calculated parameter k
    i=outpar(p,:)==k;
    % covariance from direct CANYON-B and calc (p): inpar(p,:)
    cocov.(paramnames{k})(1,svi(p,i),:)=...
        squeeze(1.*dcout(inpar(p,1),inpar(p,2),i,:)).*squeeze(cycov(k,inpar(p,1),:))+...
        squeeze(1.*dcout(inpar(p,2),inpar(p,1),i,:)).*squeeze(cycov(k,inpar(p,2),:));
    cocov.(paramnames{k})(svi(p,i),1,:)=cocov.(paramnames{k})(1,svi(p,outpar(p,:)==k),:); % and mirror
    % find second calc (o): inpar(o,:) for covariance term between calc (p) and calc (o)
        for o=p+1:6
        if any(outpar(o,:)==k) % find calc (o) that calculated parameter k
        j=outpar(o,:)==k;
        % covariance from calcs (1): 4,1 and (2): 3,1
        cocov.(paramnames{k})(svi(p,i),svi(o,j),:)=...
            squeeze(dcout(inpar(p,1),inpar(p,2),i,:).*dcout(inpar(o,1),inpar(o,2),j,:)).*squeeze(cycov(inpar(p,1),inpar(o,1),:))+...
            squeeze(dcout(inpar(p,1),inpar(p,2),i,:).*dcout(inpar(o,2),inpar(o,1),j,:)).*squeeze(cycov(inpar(p,1),inpar(o,2),:))+...
            squeeze(dcout(inpar(p,2),inpar(p,1),i,:).*dcout(inpar(o,1),inpar(o,2),j,:)).*squeeze(cycov(inpar(p,2),inpar(o,1),:))+...
            squeeze(dcout(inpar(p,2),inpar(p,1),i,:).*dcout(inpar(o,2),inpar(o,1),j,:)).*squeeze(cycov(inpar(p,2),inpar(o,2),:));
        cocov.(paramnames{k})(svi(o,j),svi(p,i),:)=cocov.(paramnames{k})(svi(p,i),svi(o,j),:); % and mirror
        end % outpar(o,:)
        end % for o
    end % outpar(p,:)
    end % for p
    % correlation matrix
    %cocov.(['r_' paramnames{k}])=cocov.(paramnames{k})*NaN; 
    %for i=1:nol
    %    cocov.(['r_' paramnames{k}])(:,:,i)=cocov.(paramnames{k})(:,:,i)./(ones(4,1)*sqrt(diag(cocov.(paramnames{k})(:,:,i)))')./(sqrt(diag(cocov.(paramnames{k})(:,:,i)))*ones(1,4)); 
    %end
end % all k params

% define weights
for i=1:4, 
    w.(paramnames{i})=1./sigma.(paramnames{i}).^2; % weights 
    w.([paramnames{i} 'sum'])=sum(w.(paramnames{i}),2); % and sum all to normalize weights to 1
end

% and output each variable
for i=1:4
% weighted mean
out.(paramnames{i})=reshape(sum(w.(paramnames{i}).*raw.(paramnames{i}),2)./w.([paramnames{i} 'sum']),size(pres)); % / [param unit]
% standard deviation about the mean (of weighted mean...; for 4 samples.)
sigma_delta=sqrt(sum(w.(paramnames{i}).*(raw.(paramnames{i})-(out.(paramnames{i})(:)*ones(1,4))).^2,2)./(w.([paramnames{i} 'sum']) - sum(w.(paramnames{i}).^2,2)./w.([paramnames{i} 'sum']))); % / [param unit]; is localized
% std propagation from correlated inputs for weighted mean
sigma_propagated=sqrt(squeeze(sum(sum(...
    repmat(permute((w.(paramnames{i})./repmat(w.([paramnames{i} 'sum']),1,4)) ,[3 2 1]),4,1) .*...
    repmat(permute((w.(paramnames{i})./repmat(w.([paramnames{i} 'sum']),1,4))',[1 3 2]),1,4) .*...
    cocov.(paramnames{i}),2),1))); % / [param unit]; is localized
out.([paramnames{i} '_sigma'])=reshape(sigma_delta+sigma_propagated,size(pres)); % / [param unit]
out.([paramnames{i} '_sigma_min'])=reshape(sigma_propagated,size(pres)); % % / [param unit]
clear sigma_delta sigma_propagated
end
% and raw calcs
for i=1:4, out.([paramnames{i} '_raw'])=raw.(paramnames{i}); end
out.sigma=sigma; % record of weight sigmas
out.cy=rmfield(cy,'inx'); % CANYON-B structure

%%{
% and add training density
if ~isempty(dir('GLODAPtraining_coverage_CONTENT_v4.mat'))
load('GLODAPtraining_coverage_CONTENT_v4.mat')
gvec=datevec(gtime(:));
inddoy=floor(datenum(gtime(:))-datenum(gvec(:,1),1,0)); % only full yearday used
inddoy(inddoy>365)=365;
lon(lon>180)=lon(lon>180)-360;
indlon=round(lon(:))+181;
indlat=round(lat(:))+91;
isubs=sub2ind([361,181,365],indlon,indlat,inddoy);
out.ncruise=ones(length(gtime(:)),6)*NaN;
out.nstation=ones(length(gtime(:)),6)*NaN;
inan=~isnan(isubs);
out.ncruise(inan,:) =[ ncruise_10_31d(isubs(inan))  ncruise_10_61d(isubs(inan))  ncruise_10_91d(isubs(inan))  ncruise_20_31d(isubs(inan))  ncruise_20_61d(isubs(inan))  ncruise_20_91d(isubs(inan))];
out.nstation(inan,:)=[nstation_10_31d(isubs(inan)) nstation_10_61d(isubs(inan)) nstation_10_91d(isubs(inan)) nstation_20_31d(isubs(inan)) nstation_20_61d(isubs(inan)) nstation_20_91d(isubs(inan))];
end
%%}
