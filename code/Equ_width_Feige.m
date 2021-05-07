% Input: 'obsspec.txt' the observed spectrum of Feige 110, 'FeigeA0_ALL' the synthetic spectrum of Feige 110 containing all IG elements, and 'Isolated_List1200.txt' a list of all isolated lines (with the corresponding transitions) with lambda > 1200A

% Output: 'ANS' a list of all isolated lines with their line properties taken from the obsspec and the synspec

% ---------------------- START -- CODE ----------------------

% reading the observed and synthetic spectrum, and the list of isolated lines

A = importdata('obsspec.txt');
obsspec = A.data;
A = importdata('FeigeA0_ALL');
synspec = A.data;
A = readcell('Isolated_List1200.txt');
linelist = cell2mat(A(:,1));
gflist = cell2mat(A(:,3));

% determining the continuum of synspec

[envelope, cont ] = FUNmakeCont(synspec,7.5,1467,1200,1230,2.5,1635,1645,0.5);

% Determination of equivalent widths in the synspec

% finding local maxima

[pks, locs] = findpeaks(synspec(:,2));
synspecval = synspec(:,2);

% finding corresponding isolated lines

l = length(linelist);

for k = 1 : l   
    
    clear linefit
    clear linelambda
    clear linevalue
    
dist = abs(linelist(k)-synspec(:,1));
minDist = min(dist);
idx = min(find(dist == minDist));

localcont = cont(synspec(idx,1));

% determining the next local maxima

dist = abs(linelist(k)-synspec(locs,1));
minDist = min(dist);
idloc1 = min(find(dist == minDist));

if linelist(k) > synspec(locs(idloc1))   
    idloc2 = idloc1+1;
else   
    idloc2 = idloc1;
    idloc1 = idloc2-1;
end

% determining the fitting interval

linelambda = (synspec(locs(idloc1) : locs(idloc2)));
linevalue  = (synspecval(locs(idloc1) : locs(idloc2)));

% local creation of the reduced spectrum

for n = 1 : length(linevalue)
    linevalue(n) = localcont - linevalue(n);
end

% fitting a Gaussian, determining line properties

linefit = fit(linelambda(:), linevalue(:), 'gauss1');

int = linefit.a1*sqrt(pi)*linefit.c1;
int = int*1000;
int = int/localcont;

Equwidth_syn(k) = int;
Sigma_syn(k) = linefit.c1;
Lambda_syn(k) = linefit.b1;

end

clear pks
clear locs

% Determination of equivalent widths in the obsspec

% finding local maxima

[pks, locs] = findpeaks(obsspec(:,2));
obsspecval = obsspec(:,2);

% finding corresponding isolated lines

for k = 1 : l    
    
    clear linefit
    clear linelambda
    clear linevalue
    
dist = abs(linelist(k)-obsspec(:,1));
minDist = min(dist);
idx = min(find(dist == minDist));

localcont = cont(obsspec(idx,1));

% determining the fitting interval: 95% of the integrated line flux (obtained from synspec) shall be inside, either an interval to the next local maxima is taken or at least pm 1.2 sigma (obtained from the synspec)

dist = abs(linelist(k)+1.2*Sigma_syn(k)-obsspec(:,1));
minDist = min(dist);
idobs2 = min(find(dist == minDist));

dist = abs(linelist(k)-1.2*Sigma_syn(k)-obsspec(:,1));
minDist = min(dist);
idobs1 = min(find(dist == minDist));

dist = abs(linelist(k)-obsspec(locs,1));
minDist = min(dist);
idloc1 = min(find(dist == minDist));

 if linelist(k) > obsspec(locs(idloc1))   
     idloc2 = idloc1+1;
 else   
     idloc2 = idloc1;
     idloc1 = idloc2-1;
 end

if locs(idloc1) < idobs1
    idobs1 = locs(idloc1);
end
if locs(idloc2) > idobs2
    idobs2 = locs(idloc2);
end

 linelambda = (obsspec(idobs1 : idobs2));
 linevalue  = (obsspecval(idobs1 : idobs2));

% local creation of the reduced flux

for n = 1 : length(linevalue)
    linevalue(n) = localcont - linevalue(n);
end

% fitting a Gaussian, determining line properties

options = fitoptions('gauss1', 'Lower', [-10 linelist(k)-1 0.0001], 'Upper', [100 linelist(k)+1 1]);
linefit = fit(linelambda(:), linevalue(:), 'gauss1', options);

int = linefit.a1*sqrt(pi)*linefit.c1;
int = int*1000;
int = int/localcont;

% determining W_obs,red

if int > 0
Sigma_obs(k) = linefit.c1;
Lambda_obs(k) = linefit.b1;
if Sigma_obs(k)>Sigma_syn(k) 
Equwidth_obsred(k) = int*Sigma_syn(k)/Sigma_obs(k)*(1-(abs(Lambda_syn(k)-Lambda_obs(k))/Sigma_obs(k))^2); 
else
Equwidth_obsred(k) = int*(1-(abs(Lambda_syn(k)-Lambda_obs(k))/Sigma_obs(k))^2);
end
Equwidth_obs(k) = int;
else   
Sigma_obs(k) = 0;
Lambda_obs(k) = Lambda_syn(k);
Equwidth_obs(k) = 0;
Equwidth_obsred(k) = 0;
end

if Equwidth_obsred(k)< 0
    Equwidth_obsred(k) = 0;
end

Sigma(k) = Sigma_obs(k)/Sigma_syn(k);
Del_Lambda(k) = abs(round(1000*(Lambda_syn(k)-Lambda_obs(k))));

end

ANS = round([ linelist, Equwidth_syn(:), Equwidth_obs(:),Equwidth_obsred(:), Sigma(:), Del_Lambda(:), gflist(:)],2);

