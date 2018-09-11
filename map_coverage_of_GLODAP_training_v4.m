function map_coverage_of_GLODAP_training_v4
%function map_coverage_of_GLODAP_training
%
% load GLODAPv2 global data collection, re-do CANYON-B training QC selection,
% reduce data set to station level and then rearrange number of stations/
% locations/ dates/ number of cruises to a 1°x1°x365 days gridded look-up
% table to give the density of carbonate system training data as number of
% individual cruises and number of stations within
% - a box of +/-10° and +/-15 days (1 month window):  _10_31d
% - a box of +/-10° and +/-30 days (2 months window): _10_61d
% - a box of +/-10° and +/-45 days (3 months window): _10_91d
% - a box of +/-20° and +/-15 days (1 month window):  _20_31d
% - a box of +/-20° and +/-30 days (2 months window): _20_61d
% - a box of +/-20° and +/-45 days (3 months window): _20_91d
%
% Henry Bittig, LOV
% 20.09.2017

glodapdir='/Users/hbittig/Documents/Science/resources/Climatologies/GLODAPv2/';

% load month, location and CANYON inputs, as well as cruise refs
load([glodapdir 'GLODAPv2 Merged Master File.mat'],'G2day','G2month','G2year','G2oxygen','G2oxygenf','G2oxygenqc','G2latitude','G2longitude')
load([glodapdir 'GLODAPv2 Merged Master File.mat'],'G2pressure','G2temperature','G2salinity','G2salinityf','G2salinityqc','expocode','G2cruise','G2station')
% load carbon obs
load([glodapdir 'GLODAPv2 Merged Master File.mat'],'G2talk','G2talkf','G2talkqc','G2tco2','G2tco2f','G2tco2qc','G2phtsinsitutp','G2phtsinsitutpf','G2phtsqc')
% define Med Sea flag where 2nd QC criterion has been relaxed for training
flagmed=ismember(G2cruise,find(ismember(expocode,{'318M19771204';'06MT20011018';'06MT20110405'})));

% re-do quality flag selection of CANYON training
% inputs
flag=ismember(G2salinityf,[0 2]) & ismember(G2oxygenf,[0 2]) & ~isnan(G2salinity) & ~isnan(G2oxygen) & ~isnan(G2temperature) & ~isnan(G2pressure);
flag2ndqc=(G2salinityqc==1 & G2oxygenqc==1) | flagmed;
% at least 2 carbon parameters (for CONTENT)
flagcarb=(ismember(G2tco2f,[0 2]) & ~isnan(G2tco2)) + (ismember(G2talkf,[0 2]) & ~isnan(G2talk)) + (ismember(G2phtsinsitutpf,[0 2]) & ~isnan(G2phtsinsitutp))>=2;
flagcarb2ndqc=(G2tco2qc + G2talkqc + G2phtsqc)>=2 | flagmed;
allflag=flag & flag2ndqc & flagcarb & flagcarb2ndqc;
clear G2o* G2t* G2sa* G2p*
% apply flags
varnames=whos('G2*');varnames=cellstr(char(varnames.name));
for i=1:length(varnames),eval([varnames{i} '=' varnames{i} '(allflag);']), end
clear flag* allflag

% only stations with carbon data
[~,b]=unique(G2station+G2cruise/1e3); % b: index to unique stations

% now reduce to station level (no depth dimension any longer)
varnames=whos('G2*');varnames=cellstr(char(varnames.name));
for i=1:length(varnames),eval([varnames{i} '=' varnames{i} '(b);']), end
clear b varnames i

% create doy
doy=datetodoy(datenum([G2year G2month G2day])); % make explicit
doy(doy>365)=365; % and fudge day 366 to day 365 for 'nominal' year length
%{
% cycle each cruise to assign stations of each cruise to day; 
fdays=false(length(expocode),365); 
for i=1:length(expocode)
    ind=G2cruise==i; % ind: data of cruise i
    fdays(i,unique(doy(ind)))=1; % mark days on which there is a station
end
clear ind 
%}

% unique cruises
expocode=expocode(unique(G2cruise));
nocr=length(expocode);
[~,~,c]=unique(G2cruise); % c: reverse index to unique cruises
disp([num2str(length(G2cruise)) ' unique stations from ' num2str(nocr) ' carbon cruises'])

% now cycle 1°x1° grid to create +/-10° and +/-20° density
% number of cruises per day (1 cruise = 1)
% number of stations per day (1 station = 1)
[longrid,latgrid]=ndgrid(-180:180,-90:90); % duplicate at 180° longitude; ok for interpolation
ncruise_10=false(361,181,365,nocr);ncruise_20=false(361,181,365,nocr); % too large in memory
nstation_10=zeros(361,181,365);nstation_20=zeros(361,181,365);
% stupid for loops...
for j=1:size(longrid,2)
disp([num2str(j) '/181'])
    for i=1:size(longrid,1)
        % find stations within \pm10° / \pm20° box
        ind10=((G2longitude-longrid(i,j)<=10 & G2longitude-longrid(i,j)>=0) | mod(longrid(i,j)-G2longitude,360)<10) & (abs(G2latitude-latgrid(i,j))<10 | G2latitude-latgrid(i,j)==10);
        ind20=((G2longitude-longrid(i,j)<=20 & G2longitude-longrid(i,j)>=0) | mod(longrid(i,j)-G2longitude,360)<20) & (abs(G2latitude-latgrid(i,j))<20 | G2latitude-latgrid(i,j)==20);
        % only count those days with stations within the area
        nstation_10(i,j,:)=arrayfun(@(x)sum(doy(ind10)==x),1:365);
        nstation_20(i,j,:)=arrayfun(@(x)sum(doy(ind20)==x),1:365);
        % and only use those cruise days that have stations in the area;
        %%{
        ones10=ones(sum(ind10),1);ones20=ones(sum(ind20),1);
        % mark grid box and day of given cruise
        ind=sub2ind([361 181 365 nocr],i*ones10,j*ones10,doy(ind10),c(ind10));
        ncruise_10(ind)=1;
        ind=sub2ind([361 181 365 nocr],i*ones20,j*ones20,doy(ind20),c(ind20));
        ncruise_20(ind)=1;
        %%}        
        clear ind10 ind20 ind ones10 ones20
    end
end

doygrid=1:365;
% and sum #cruises within +/-15 days and +/-45 days
ncruise_10_31d=zeros(361,181,365);ncruise_20_31d=zeros(361,181,365);
ncruise_10_61d=zeros(361,181,365);ncruise_20_61d=zeros(361,181,365);
ncruise_10_91d=zeros(361,181,365);ncruise_20_91d=zeros(361,181,365);
for i=1:365
    if ~mod(i-1,30),disp([num2str(i) '/365']),end
    ind=abs(doygrid-i)<=15 | abs(doygrid-365-i)<=15;
    ncruise_10_31d(:,:,i)=sum(any(ncruise_10(:,:,ind,:),3),4);
    ncruise_20_31d(:,:,i)=sum(any(ncruise_20(:,:,ind,:),3),4);
end
for i=1:365
    if ~mod(i-1,30),disp([num2str(i) '/365']),end
    ind=abs(doygrid-i)<=30 | abs(doygrid-365-i)<=30;
    ncruise_10_61d(:,:,i)=sum(any(ncruise_10(:,:,ind,:),3),4);
    ncruise_20_61d(:,:,i)=sum(any(ncruise_20(:,:,ind,:),3),4);
end
for i=1:365
    if ~mod(i-1,30),disp([num2str(i) '/365']),end
    ind=abs(doygrid-i)<=45 | abs(doygrid-365-i)<=45;
    ncruise_10_91d(:,:,i)=sum(any(ncruise_10(:,:,ind,:),3),4);
    ncruise_20_91d(:,:,i)=sum(any(ncruise_20(:,:,ind,:),3),4);
end 

% and sum #stations within same ranges
nstation_10_31d=zeros(361,181,365);nstation_20_31d=zeros(361,181,365);
nstation_10_61d=zeros(361,181,365);nstation_20_61d=zeros(361,181,365);
nstation_10_91d=zeros(361,181,365);nstation_20_91d=zeros(361,181,365);
for i=1:365
    if ~mod(i-1,30),disp([num2str(i) '/365']),end
    ind=abs(doygrid-i)<=15 | abs(doygrid-365-i)<=15;
    nstation_10_31d(:,:,i)=sum(nstation_10(:,:,ind),3);
    nstation_20_31d(:,:,i)=sum(nstation_20(:,:,ind),3);
end
for i=1:365
    if ~mod(i-1,30),disp([num2str(i) '/365']),end
    ind=abs(doygrid-i)<=30 | abs(doygrid-365-i)<=30;
    nstation_10_61d(:,:,i)=sum(nstation_10(:,:,ind),3);
    nstation_20_61d(:,:,i)=sum(nstation_20(:,:,ind),3);
end
for i=1:365
    if ~mod(i-1,30),disp([num2str(i) '/365']),end
    ind=abs(doygrid-i)<=45 | abs(doygrid-365-i)<=45;
    nstation_10_91d(:,:,i)=sum(nstation_10(:,:,ind),3);
    nstation_20_91d(:,:,i)=sum(nstation_20(:,:,ind),3);
end 

%save('GLODAPtraining_coverage_CONTENT.mat','expocode','longrid','latgrid','ncruise_10_31d','ncruise_20_31d','ncruise_10_91d','ncruise_20_91d','nstation_10','nstation_20','-v7.3');
%save('GLODAPtraining_coverage_CONTENT.mat','expocode','longrid','latgrid','ncruise_10_31d','ncruise_20_31d','ncruise_10_91d','ncruise_20_91d','nstation_10','nstation_20');
save('GLODAPtraining_coverage_CONTENT_v4.mat','expocode','longrid','latgrid','ncruise_10_31d','ncruise_20_31d','ncruise_10_61d','ncruise_20_61d','ncruise_10_91d','ncruise_20_91d','nstation_10_31d','nstation_20_31d','nstation_10_61d','nstation_20_61d','nstation_10_91d','nstation_20_91d');
disp('done')