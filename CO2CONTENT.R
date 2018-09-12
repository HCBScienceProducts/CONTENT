CO2CONTENT <- function(date,lat,lon,pres,temp,psal,doxy,param,epres,etemp,epsal,edoxy){
# function out=CO2CONTENT(gtime,lat,lon,pres,temp,psal,doxy,epres,etemp,epsal,edoxy)
# 
# CONTENT calculation according to Bittig et al. (2018):
# Combination of CANYON-B Bayesian neural network mappings of AT, CT, pH,
# and pCO2 made consistent with carbonate chemistry constraints 
# for any set of water column P, T, S, O2, location data as an alternative
# to (spatial) climatological interpolation. Based on GLODAPv2/GO-SHIP
# bottle data.
#
# Please cite the paper when using CONTENT/CANYON-B.
#
# input:
# date  - date (UTC) as string ("yyyy-mm-dd HH:MM")
# lat   - latitude / °N  [-90 90]
# lon   - longitude / °E [-180 180] or [0 360]
# pres  - pressure / dbar
# temp  - in-situ temperature / °C
# psal  - salinity
# doxy  - dissolved oxygen / umol kg-1 (!)
# e(var)- input error for pres, temp, psal, doxy
#
# default values for input errors:
# epres:   0.5 dbar
# etemp:   0.005 ?C
# epsal:   0.005 psu
# edoxy:   1 % of doxy (<- This is a rather optimistic default, meant for
# GO-SHIP quality bottle data. A reasonable default sensor doxy error would
# be 3 % of doxy (or higher!).)
#
# output:
# out$(var)           - estimate of variable (same size as 'pres' input)
# out$(var)_sigma     - 1 sigma of variable (same size as 'pres' input)
# out$(var)_sigma_min - 1 sigma_min of variable (same size as 'pres' input)
# out$(var)_raw       - raw CANYON and 3 indirect estimates of variable (size numel(pres)x4)
# out$ncruise   [ONLY IN MATLAB VERSION] - number of cruises within...[see 1..6 below] (size numel(pres)x6)
# out$nstation  [ONLY IN MATLAB VERSION] - number of stations within...[see 1..6 below] (size numel(pres)x6)
#                       1: +/-10° and +/-15 d window (doy-based)
#                       2: +/-10° and +/-30 d window (doy-based)
#                       3: +/-10° and +/-45 d window (doy-based)
#                       4: +/-20° and +/-15 d window (doy-based)
#                       5: +/-20° and +/-30 d window (doy-based)
#                       6: +/-20° and +/-45 d window (doy-based)
#
# check values for 
# 09-Dec-2014 08:45, 17.6° N, -24.3° E, 180 dbar, 16 ?C, 36.1 psu, 104 umol O2 kg-1:
#                (var)    (var)_sigma  (var)_sigma_min      unit
#         AT:    2358.068     10.470        7.681      umol kg-1
#         CT:    2199.216     10.114        7.554      umol kg-1
#         pH:    7.869670      0.021982     0.017453   (insitu total scale, comparable to "calculated pH values" as described in Carter et al., 2018)
#         pCO2:  639.7380     36.1239      29.0812          uatm
#
# examples:
# out=CO2CONTENT(date="2014-12-09 08:45", lat=17.6, lon=-24.3, pres=180, temp=16, psal=36.1, doxy=104);  # all variables, default input errors
# out=CO2CONTENT(date="2014-12-09 08:45", lat=17.6, lon=-24.3, pres=180, temp=16, psal=36.1, doxy=104, edoxy=0.03*104);  # all variables, 3 % input error on doxy, otherwise default input errors
# out=CO2CONTENT(date=rep("2014-12-09 08:45",2), lat=rep(17.6,2), lon=rep(-24.3,2), pres=rep(180,2), temp=rep(16,2), psal=rep(36.1,2), doxy=rep(104,2));  # more than one single input; default input errors
#
# references:
#
# - CONTENT and CANYON-B method:
# Bittig et al. (2018). An alternative to static climatologies: Robust
# estimation of open ocean CO2 variables and nutrient concentrations from
# T, S and O2 data using Bayesian neural networks. Front. Mar. Sci. 5:328. doi:
# 10.3389/fmars.2018.00328 
#
# requires seacarb (available from CRAN):
# Gattuso et al. (2018). seacarb: Seawater Carbonate Chemistry. R package 
# version 3.2.8. http://CRAN.R-project.org/package=seacarb
#
# Henry Bittig, LOV
# v0.8, 02.03.2018, included correlations between carbon parameters;  uses CANYON-B
# v1.0, 11.09.2018, initial publication

#  The software is provided "as is", without warranty of any kind, 
#  express or implied, including but not limited to the warranties of
#  merchantability, fitness for a particular purpose and noninfringement.
#  In no event shall the authors or copyright holders be liable for any
#  claim, damages or other liability, whether in an action of contract,
#  tort or otherwise, arising from, out of or in connection with the
#  software or the use or other dealings in the software.      

# No input checks! Assumes informed use, e.g., same dimensions for all
# inputs, ...


## Nothing below here should need to be changed by the user ##

if (missing(date)) {date<-rep("2014-12-09 08:45",2); lat<-17.6*c(1,1); lon<--24.3*c(1,1); pres<-180*c(1,1); temp<-16*c(1,1); psal<-36.1*c(1,1);doxy<-104*c(1,1); } # test values
if (missing(epres)) epres=.5;
if (missing(etemp)) etemp=0.005;
if (missing(epsal)) epsal=0.005;
if (missing(edoxy)) edoxy=0.01*doxy;

# No input checks! Assumes informed use, e.g., same dimensions for all
# inputs, ...; only eX can either be scalar or same dimension as inputs

# get number of elements
nol=length(pres);
# and expand input errors if needed
if(length(epres)==1) epres=rep(epres,nol);
if(length(etemp)==1) etemp=rep(etemp,nol);
if(length(epsal)==1) epsal=rep(epsal,nol);
if(length(edoxy)==1) edoxy=rep(edoxy,nol);
# calculate (pure) carbonate system estimates (#1) and include estimates of
# nutrients from CANYON-B

# add CANYON-B carbonate estimates as well as phosphate and silicate
paramnames=c('AT','CT','pH','pCO2','SiOH4','PO4');
source('private/CANYONB_private.R')
cy=CANYONB_private(date=date, lat=lat, lon=lon, pres=pres, temp=temp, psal=psal, doxy=doxy, param=paramnames, epres=epres, etemp=etemp, epsal=epsal, edoxy=edoxy);

# copy to raw output and start mixing calculations
paramnames=c('AT','CT','pH','pCO2');
rawout=list();
for (i in (1:4)){ 
  rawout[[paramnames[i]]]=matrix(NaN,nol,4);
  rawout[[paramnames[i]]][,1]=cy[[paramnames[i]]];
} #end

# define arrangements of calculations
inpar=rbind(c(4,1),c(3,1),c(4,3),c(4,2),c(3,2),c(1,2)); # (1): 4,1; (2): 3,1; (3): 4,3; (4): 4,2; (5): 3,2; (6): 1,2
inflag=c(24, 8, 21, 25, 9, 15);
inscale=c(1e-6,1e-6,1,1);
outpar=matrix(NaN,6,2);
for (p in (1:6)) outpar[p,]=setdiff((1:4),inpar[p,]); # end
outparnames=c('ALK','DIC','pH','pCO2');
svi=rbind(c(2,4),c(3,3),c(4,4),c(2,3),c(3,4),c(2,2)); # save indices for outpar

# deal with weights / uncertainties
sigma=list();
for (i in (1:4)) sigma[[paramnames[i]]]=matrix(NaN,nol,4); # end
# nutrient uncertainty (var(inputs) + var(training data) + var(MLP))
cy$eSiOH4=cy$SiOH4_ci;cy$ePO4=cy$PO4_ci;
# covariance matrix of CANYON-B AT, CT, pH, pCO2 due to common inputs
cycov=array(NaN,c(4,4,nol));
for (i in (1:4)) {
  for (j in (i:4)) {
    cycov[i,j,]=rowSums(cy$inx[[paramnames[i]]]*cy$inx[[paramnames[j]]]*cbind(etemp,epsal,edoxy,epres)^2);cycov[j,i,]=cycov[i,j,];
  } # end
} # end

# correlation matrix of CANYON-B AT, CT, pH, pCO2 due to common inputs and
# CANYON-B estimation uncertainty
cyr=cycov*0;
# and full variance on diagonal: 
# var(inputs)[local] + var(training data)[global] + var(MLP)[local]
sigma$AT[,1]=cy$AT_ci; cycov[1,1,]=cy$AT_ci^2;
sigma$CT[,1]=cy$CT_ci; cycov[2,2,]=cy$CT_ci^2;
sigma$pH[,1]=cy$pH_ci; cycov[3,3,]=cy$pH_ci^2;
sigma$pCO2[,1]=cy$pCO2_ci; cycov[4,4,]=cy$pCO2_ci^2;
# correlation matrix with full variance 
for (i in (1:nol)) cyr[,,i]=cycov[,,i]/t(array(rep(sqrt(diag(cycov[,,i])),4),c(4,4)))/array(rep(sqrt(diag(cycov[,,i])),4),c(4,4)); # end

# calculate derivatives of carbonate system calculations (#2-#4) and store
# dx/dy's separately for each calculation (1)-(6)
dcout=array(NaN,c(4,4,2,nol));#dctit=cell(4,4,2);
for (p in (1:6)) {
  deriv=derivnum('var1',flag=inflag[p],var1=cy[[paramnames[inpar[p,1]]]]*inscale[inpar[p,1]], var2=cy[[paramnames[inpar[p,2]]]]*inscale[inpar[p,2]], S=psal, T=temp, P=pres/10, Patm=1, Pt=cy$PO4*1e-6, Sit=cy$SiOH4*1e-6, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="standard");
  dcout[inpar[p,1],inpar[p,2],1,]=deriv[[outparnames[outpar[p,1]]]]/inscale[outpar[p,1]]*inscale[inpar[p,1]]; 
  dcout[inpar[p,1],inpar[p,2],2,]=deriv[[outparnames[outpar[p,2]]]]/inscale[outpar[p,2]]*inscale[inpar[p,1]]; # dTCO2/dpCO2, dpHin/dpCO2
  deriv=derivnum('var2',flag=inflag[p],var1=cy[[paramnames[inpar[p,1]]]]*inscale[inpar[p,1]], var2=cy[[paramnames[inpar[p,2]]]]*inscale[inpar[p,2]], S=psal, T=temp, P=pres/10, Patm=1, Pt=cy$PO4*1e-6, Sit=cy$SiOH4*1e-6, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="standard");
  dcout[inpar[p,2],inpar[p,1],1,]=deriv[[outparnames[outpar[p,1]]]]/inscale[outpar[p,1]]*inscale[inpar[p,2]]; 
  dcout[inpar[p,2],inpar[p,1],2,]=deriv[[outparnames[outpar[p,2]]]]/inscale[outpar[p,2]]*inscale[inpar[p,2]];  # dTCO2/dALK, dpHin/dALK
  rm(deriv)
} # end


# do carbonate system calculations (#2-#4)
# using the K1K2 constants of Lueker et al, 2000, KSO4 of Dickson 1990 & TB of Uppstrom 1979
# add localized error calculations incl. parameter correlation due to common inputs
for (p in (1:6)) {
# outcalc=CO2SYS(cy.(paramnames{inpar(p,1)}),cy.(paramnames{inpar(p,2)}),inpar(p,1),inpar(p,2),psal(:),temp(:),NaN,pres(:),NaN,cy.SiOH4,cy.PO4,1,10,1); # (1): 4,1
  outcalc=carb(flag=inflag[p],var1=cy[[paramnames[inpar[p,1]]]]*inscale[inpar[p,1]], var2=cy[[paramnames[inpar[p,2]]]]*inscale[inpar[p,2]], S=psal, T=temp, P=pres/10, Patm=1, Pt=cy$PO4*1e-6, Sit=cy$SiOH4*1e-6, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="standard"); # (1): 4,1
#outerr= real( errors(cy.(paramnames{inpar(p,1)}),cy.(paramnames{inpar(p,2)}),inpar(p,1),inpar(p,2),psal(:),temp(:),NaN,pres(:),NaN,cy.SiOH4,cy.PO4,...
#  sqrt(cycov(inpar(p,1),inpar(p,1),:)),sqrt(cycov(inpar(p,2),inpar(p,2),:)),epsal,etemp,cy.eSiOH4,cy.ePO4,'',0,cyr(inpar(p,1),inpar(p,2),:),1,10,1) ) ; # (1): 4,1
  outerr=errors(flag=inflag[p],var1=cy[[paramnames[inpar[p,1]]]]*inscale[inpar[p,1]], var2=cy[[paramnames[inpar[p,2]]]]*inscale[inpar[p,2]], S=psal, T=temp, P=pres/10, Patm=1, Pt=cy$PO4*1e-6, Sit=cy$SiOH4*1e-6, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="standard",
    evar1=sqrt(cycov[inpar[p,1],inpar[p,1],])*inscale[inpar[p,1]],evar2=sqrt(cycov[inpar[p,2],inpar[p,2],])*inscale[inpar[p,2]],eS=epsal,eT=etemp,ePt=cy$ePO4*1e-6,eSit=cy$eSiOH4*1e-6,r=cyr[inpar[p,1],inpar[p,2],],epK=c(0.004,0.015,0.03,0.01,0.01,0.02,0.02),method="mo"); # (1): 4,1
  rawout[[paramnames[outpar[p,1]]]][,svi[p,1]]=outcalc[[outparnames[outpar[p,1]]]]/inscale[outpar[p,1]]; # outpar 1
  rawout[[paramnames[outpar[p,2]]]][,svi[p,2]]=outcalc[[outparnames[outpar[p,2]]]]/inscale[outpar[p,2]]; # outpar 2
  sigma[[paramnames[outpar[p,1]]]][,svi[p,1]]=outerr[[outparnames[outpar[p,1]]]]/inscale[outpar[p,1]]; # eoutpar 1
  sigma[[paramnames[outpar[p,2]]]][,svi[p,2]]=outerr[[outparnames[outpar[p,2]]]]/inscale[outpar[p,2]]; # eoutpar 2
  if (inpar[p,1]==3) { # check if first input variable is pH 
    # dz/dH -> dz/dpH
	  for (i in (1:2)) { # first and second variable doutpar(p,:) / dinpar(p,1)
      dcout[inpar[p,1],inpar[p,2],i,]=(dcout[inpar[p,1],inpar[p,2],i,])*(-log(10)*10^(-outcalc$pH)); 
    } # end for
  } # end if 
  if (inpar[p,2]==3) { # check if second input variable is pH 
    # dz/dH -> dz/dpH
    for (i in (1:2)) { # first and second variable doutpar(p,:) / dinpar(p,2)
      dcout[inpar[p,2],inpar[p,1],i,]=(dcout[inpar[p,2],inpar[p,1],i,])*(-log(10)*10^(-outcalc$pH)); 
    } # end for
  } # end if
  rm(outcalc,outerr)
} # end for p=1:6

# build weighted mean covariance matrix for all parameters
cocov=list();
for (k in (1:4)) {
  cocov[[paramnames[k]]]=array(NaN,c(4,4,nol));
  for (i in (1:4)) cocov[[paramnames[k]]][i,i,]=sigma[[paramnames[k]]][,i]^2; # end
  for (p in (1:6)) {
    if (any(outpar[p,]==k)) { # find calc (p) that calculated parameter k
      i=outpar[p,]==k;
      # covariance from direct CANYON-B and calc (p): inpar(p,:)
      cocov[[paramnames[k]]][1,svi[p,i],]=
        1*dcout[inpar[p,1],inpar[p,2],i,]*cycov[k,inpar[p,1],]+
        1*dcout[inpar[p,2],inpar[p,1],i,]*cycov[k,inpar[p,2],];
      cocov[[paramnames[k]]][svi[p,i],1,]=cocov[[paramnames[k]]][1,svi[p,outpar[p,]==k],]; # and mirror
      # find second calc (o): inpar(o,:) for covariance term between calc (p) and calc (o)
      if (p<6) { for (o in ((p+1):6)) {
        if (any(outpar[o,]==k)) { # find calc (o) that calculated parameter k
          j=outpar[o,]==k;
          # covariance from calcs (1): 4,1 and (2): 3,1
          cocov[[paramnames[k]]][svi[p,i],svi[o,j],]=
            dcout[inpar[p,1],inpar[p,2],i,]*dcout[inpar[o,1],inpar[o,2],j,]*cycov[inpar[p,1],inpar[o,1],]+
            dcout[inpar[p,1],inpar[p,2],i,]*dcout[inpar[o,2],inpar[o,1],j,]*cycov[inpar[p,1],inpar[o,2],]+
            dcout[inpar[p,2],inpar[p,1],i,]*dcout[inpar[o,1],inpar[o,2],j,]*cycov[inpar[p,2],inpar[o,1],]+
            dcout[inpar[p,2],inpar[p,1],i,]*dcout[inpar[o,2],inpar[o,1],j,]*cycov[inpar[p,2],inpar[o,2],];
          cocov[[paramnames[k]]][svi[o,j],svi[p,i],]=cocov[[paramnames[k]]][svi[p,i],svi[o,j],]; # and mirror
        } #end # outpar(o,:)
      } } # end # for o
    } # end # outpar(p,:)
  } # end # for p
} # end # all k params

# define weights
w=list();
for (i in (1:4)) {
  w[[paramnames[i]]]=1/sigma[[paramnames[i]]]^2; # weights 
  w[[paste0(paramnames[i],'sum')]]=rowSums(w[[paramnames[i]]]); # and sum all to normalize weights to 1
} # end

# and output each variable
out=list();
for (i in (1:4)) {
  # weighted mean
  out[[paramnames[i]]]=rowSums(w[[paramnames[i]]]*rawout[[paramnames[i]]])/w[[paste0(paramnames[i],'sum')]]; # / [param unit]
  # standard deviation about the mean (of weighted mean...; for 4 samples.)
  sigma_delta=sqrt(rowSums(w[[paramnames[i]]]*(rawout[[paramnames[i]]]-matrix(out[[paramnames[i]]],nol,4))^2)/(w[[paste0(paramnames[i],'sum')]] - rowSums(w[[paramnames[i]]]^2)/w[[paste0(paramnames[i],'sum')]])); # / [param unit]; is localized
  # std propagation from correlated inputs for weighted mean
  sigma_propagated=sqrt(apply(
    aperm(kronecker(array(1,c(1,1,4)),w[[paramnames[i]]]/matrix(w[[paste0(paramnames[i],'sum')]],nol,4)),c(3,2,1)) *
    aperm(kronecker(array(1,c(1,1,4)),w[[paramnames[i]]]/matrix(w[[paste0(paramnames[i],'sum')]],nol,4)),c(2,3,1)) *
    cocov[[paramnames[i]]],3,sum)); # / [param unit]; is localized
  out[[paste0(paramnames[i],'_sigma')]]=sigma_delta+sigma_propagated; # / [param unit]
  out[[paste0(paramnames[i],'_sigma_min')]]=sigma_propagated; # # / [param unit]
  rm(sigma_delta,sigma_propagated)
} # end
# and raw calcs
for (i in (1:4)) out[[paste0(paramnames[i],'_raw')]]=rawout[[paramnames[i]]]; # end
out$sigma=sigma; # record of weight sigmas
cy$inx=NULL;out$cy=cy; # CANYON-B structure


return(out)
} # end of function