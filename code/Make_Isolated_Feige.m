% Input: 'Isoprep_Lines.txt' a line list of all investigated ions, obtained from 'Prepare_Isolated_Feige' (containing the line properties), 'FeigeA0_RED_ALL.txt' the reduced synthetic spectrum of Feige 110 containing all investigated ions, and 'FeigeA0_RED_Elem,Ion.txt' reduced spectra for each of the investigated ions

% Output: 'ANS' a list of all isolated lines, containing the line center and the respective ion

% ---------------------- START -- CODE ----------------------

% reading the prepared line list and the reduced spectrum of all ions

A = importdata('Isoprep_Lines.txt');
spec = importdata('FeigeA0_RED_ALL.txt');

% somewhat complicated reading the reduced spectra of the different ions and saving in a cell array (LM denotes the lighter elements)

E(1) = "Ca";
E(2) = 'Sc';
E(3) = 'Ti';
E(4) = 'V_';
E(5) = 'Cr';
E(6) = 'Mn';
E(7) = 'Fe';
E(8) = 'Co';
E(9) = "Ni";
E(10) = 'LM';

I(1) = "III";
I(2) = 'IV';
I(3) = "V";
I(4) = "VI";
I(5) = 'III-VI';

z = 0;
                     
for x = 1:length(E)
    
    for y = 1:length(I)  
    
        if isfile(strcat(E(x),'/','FeigeA0_RED_',E(x),I(y),'.txt'))
           
            z = z+1;
            
    B = importdata(strcat(E(x),'/','FeigeA0_RED_',E(x),I(y),'.txt')) ; 
    
    B_cell{1,z} = B(:,1);
    B_cell{2,z} = B(:,2);
    B_cell{3,z} = strcat(E(x),I(y));

    clear B
    
        end
    end
end

% determining the isolated lines
 
 for k = 1 : length(A(:,1))

% determining the line centers and the interval of pm 1.82 sigma, so that about 99% of the integrated line flux is considered    
    
center = A(k,1);
sigma = A(k,3);
D(k,1) = center;
D(k,3) = sigma;

wlow = center-1.82*sigma;
wup  = center+1.82*sigma;

% determination of the integrated line flux

dist = abs(wlow-spec(:,1));
minDist = min(dist);
idloc1 = min(find(dist == minDist));
dist = abs(wup-spec(:,1));
minDist = min(dist);
idloc2 = min(find(dist == minDist));


linelambda = spec((idloc1 : idloc2),1);
linevalue  = spec((idloc1 : idloc2),2);

intflux_lineALL = trapz(linelambda,linevalue);

 clear linelambda
 clear linevalue

for l = 1:length(B_cell(1,:))
    
    C(:,1) = cell2mat(B_cell(1,l));
    C(:,2) = cell2mat(B_cell(2,l));
    
    dist = abs(wlow-C(:,1));
    minDist = min(dist);
    idloc1 = min(find(dist == minDist));
    dist = abs(wup-C(:,1));
    minDist = min(dist);
    idloc2 = min(find(dist == minDist));

    linelambda = C((idloc1 : idloc2),1);
    linevalue  = C((idloc1 : idloc2),2);

    intflux(l) = trapz(linelambda,linevalue);
   
    intflux_red(l) = intflux(l) / intflux_lineALL;
     
    clear linelambda
    clear linevalue
    
    clear C
end

    % test, if only one significant line is in the interval (>90% flux)

    if max(intflux_red) < 1.5 && max(intflux) > 0.90*sum(intflux)
        
        [maxval,idmax] = max(intflux_red);
        
        D(k,2) = idmax;
        
    else
        D(k,:) = 0;
    end
    
    % if a lighter-element line dominates, remove entry
    
    if D(k,2) == z
        
        D(k,:) = 0;
        
    end
end

D1 = D(any(D,2),:);

% if several transitions of the same ion are close, they were not removed yet. Here: for each found isolated line it is checked wether the integrated flux in the respective ion's reduced spectrum is matching the determined line flux

for n = 1:length(D1(:,1))
   
    center = D1(n,1);
    ion = D1(n,2);
    sigma = D1(n,3);
    
    % determining the fitting-interval
    
    wlow = center-2.5*sigma;
    wup  = center+2.5*sigma;
    
    F(:,1) = cell2mat(B_cell(1,ion));
    F(:,2) = cell2mat(B_cell(2,ion));
   
    dist = abs(wlow-F(:,1));
    minDist = min(dist);
    idloc1 = min(find(dist == minDist));
    dist = abs(wup-F(:,1));
    minDist = min(dist);
    idloc2 = min(find(dist == minDist));

    linelambda = F((idloc1 : idloc2),1);
    linevalue  = F((idloc1 : idloc2),2);

    % calculating the integrated flux
    
    int1 = trapz(linelambda,linevalue);
    
    options = fitoptions('gauss1', 'Lower', [-10 linelambda(1) 0.0001], 'Upper', [100 linelambda(length(linelambda)) 1]);
    [linefit,gof] = fit(linelambda(:), linevalue(:), 'gauss1', options);
    int2 = linefit.a1*sqrt(pi)*linefit.c1;
    
    % keeping the line if the deviation is < 0.05
    
    if 1.05 > int2/int1 > 0.95
    
    else
        
        D1(n,:) = 0;
        
    end
    
end

D2 = D1(any(D1,2),:);

% saving in a cell array

for m = 1:length(D2(:,1))
    
   ANS{m,1} = D2(m,1); 
   ANS{m,2} = B_cell{3,D2(m,2)};
   ANS{m,3} = D2(m,3);
    
end

