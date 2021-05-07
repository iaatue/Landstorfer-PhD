% Input: 'obsspec' and 'Feige_synspec', observed and synthetic spectra of Feige 110

% Output: 'corrfact' the correction factor for continuum calibration

% ---------------------- START -- CODE ----------------------

% reading the input spectra
           
    A = importdata('obsspec.txt') ; 
    obsspec = A.data;
   A = importdata('Feige_synspec');
    synspec = A.data;
    
% adding a low-pass filter (Savitzky-Golay filter)
    
    obsspec(:,2) = sgolayfilt(obsspec(:,2),3,11);
      
[envelope2,~] = FUNmakeCont(obsspec,2.5,0,1144,1180,0);

% S/N correction for the observed spectrum

for n = 1:length(envelope2)
    
    envelope2(n,2) = envelope2(n,2) - 1/80*envelope2(n,2);
    
end

% removing local maxima from 'envelope2'

[pks, locs] = findpeaks(envelope2(:,2));

for i = 1:length(locs)
   
    envelope2(locs(i),:) = 0;
    
end

envelope2 = envelope2(any(envelope2,2),:);

[pks, locs] = findpeaks(envelope2(:,2));

for i = 1:length(locs)
   
    envelope2(locs(i),:) = 0;
    
end

envelope2 = envelope2(any(envelope2,2),:);

% determining the continuum of 'synspec'

 [envelope, cont ] = FUNmakeCont(synspec,7.5,1467,1200,1230,2.5,1635,1645,0.5);

% finding 'minfact' to minimize Chi^2 when minfact*syncont = obscont applies

for i = 1:length(envelope2)
    
    envelope2(i,3) = cont(envelope2(i,1));
     
end

fact = (0.6:0.001:1.5);

for i = 1:length(fact)
    
    for j = 1:length(envelope2)
        
        difenvchi(j) = (abs(fact(i)*envelope2(j,2)-envelope2(j,3)))^2/envelope2(j,3);
               
    end
    
    difenv(i) = sqrt(sum(difenvchi));
    
end

mindifenv = min(difenv);

minfact = fact(min(find(difenv == mindifenv)));
corrfact = 1/minfact;

