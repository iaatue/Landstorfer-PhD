% Input: 'Elem','Ion' element and ion, of which all lines are looked for, the synthetic spetrum of the corresponding ion 'FeigeA0_Elem_Ion', and 'idents.txt' a list containing the line centers of lighter-element lines to be excluded

% Output: 'linelist' containing the line centers of the specified ion's lines, 'equw' containing the associated equivalent width, and 'errorsum_final' containing the fit quality

% specifying the ion of which the lines are looked for 

Elem = 'Ni';
Ion = 'IV';

% ---------------------- START -- CODE ----------------------

% reading the input spectrum and the list of lighter-element lines

A = importdata(strcat(Elem,'/','FeigeA0_',Elem, Ion));
synspec = A.data;
idents = importdata('Idents.txt');

% determining the continuum of 'synspec'

 [envelope, cont ] = FUNmakeCont(synspec,7.5,1467,1200,1230,2.5,1635,1645,0.5);

% creating the reduced spectrum

for u = 1 : length(synspec(:,1))
    
    synspec(u,2) = cont(synspec(u,1))-synspec(u,2);
    
    if synspec(u,2) < cont(synspec(u,1))/100
        synspec(u,2) = 0;
    end
    
end

% determination of line properties

% local maxima of the reduced synspec = line centers

[pks, locs] = findpeaks(synspec(:,2));

% determining the fitting interval: either 1) pm 0.1A or 2) up to one third the distance to the next maximum
 
 for k = 1 : length(locs) 
      
dist = abs(synspec(locs(k),1)-0.1-synspec(:,1));
minDist = min(dist);
idloc1 = min(find(dist == minDist));
dist = abs(synspec(locs(k),1)+0.1-synspec(:,1));
minDist = min(dist);
idloc2 = min(find(dist == minDist));

if k > 1
testint1 = (synspec(locs(k),1) - synspec(locs(k-1),1));
else
testint1 = 0.2;
end

if testint1<0.2
dist = abs(synspec(locs(k),1)-testint1/3-synspec(:,1));
minDist = min(dist);
idloc1 = min(find(dist == minDist));
end

if k < length(locs)
testint2 = (synspec(locs(k+1),1) - synspec(locs(k),1));
else
testint2 = 0.2;
end

if testint2<0.2
dist = abs(synspec(locs(k),1)+testint2/3-synspec(:,1));
minDist = min(dist);
idloc2 = min(find(dist == minDist));
end

linelambda = synspec((idloc1 : idloc2),1);
linevalue  = synspec((idloc1 : idloc2),2);

if length(linelambda) > 5

% Gaussian fit, then determining the line properties

if k > 5 && k < length(locs)-5
options = fitoptions('gauss1', 'Lower', [-10 synspec(locs(k)-5,1) 0.0001], 'Upper', [100 synspec(locs(k)+5,1) 1]);
end
[linefit,gof] = fit(linelambda(:), linevalue(:), 'gauss1');

int = linefit.a1*sqrt(pi)*linefit.c1;
int = int*1000;
int = int/cont(synspec(locs(k),1));

 Equwidth_syn(k) = int;
 Sigma_syn(k) = linefit.c1;
 errorsum(k) = gof.sse;
 coefferr(k) = gof.rsquare;

else
 
 Equwidth_syn(k) = 0;
 Sigma_syn(k) = 0;
 errorsum(k) = 1;
 coefferr(k) = 1;
    
end  
 end

% excluding all lines which are closer than 0.03A to a lighter-element line center 

b = 1;
for m = 1:length(locs)
    
    a = 0;
 for n = 1:length(idents)   
 if  abs(synspec(locs(m),1) - idents(n)) < 0.03 
     a = a+1;
 
 end    
 end  
 
 if a == 0
  
    linelist(b) = synspec(locs(m));
    equw(b) = Equwidth_syn(m);
     
    errorsum_fin(b) = errorsum(m);
    coefferr_fin(b) = coefferr(m);
     
     b= b+1;
     
 end   
end

