% Input: 'spectrum' of which the local continuum is to be determined, 'maxint' the proper interval to find local maxima in, 'outlier1' can be used to exclude a maximum from the fit, 'wvllow', 'wvlup' and 'int' denote lower and upper wavelengths and the interval to find maxima for segments which need to be treated seperately

% Output: 'envelope' including the found maxima, and 'cont' which includes the spline-fit of the estimated continuum

function [envelope, cont] = FUNmakeCont(spectrum,maxint,outlier1,wvllow1,wvlup1,int1,wvllow2,wvlup2,int2,wvllow3,wvlup3,int3,wvllow4,wvlup4,int4)

% set 'int = 0' to exclude a segment from the fit

% ---------------------- START -- CODE ----------------------

% general continuum with 'maxint'

minspecwvl = spectrum(1,1);
maxspecwvl = spectrum(length(spectrum),1);

a = (minspecwvl:maxint:maxspecwvl+maxint/2);

for i = 1:length(a)
    
    dist1 = abs(a(i)-maxint/1.9-spectrum(:,1));
    minDist = min(dist1);
    idx1 = min(find(dist1 == minDist));
    dist2 = abs(a(i)+maxint/1.9-spectrum(:,1));
    minDist = min(dist2);
    idx2 = max(find(dist2 == minDist));
      
    for j = 1 : (idx2-idx1+1)        
    b(j) = spectrum(j+idx1-1,2);           
    end
    
    maxVal = max(b);
    idMax = max(find(b == maxVal));
    
    envelope(i,1) = spectrum(idMax + idx1-1,1);
    envelope(i,2) = spectrum(idMax + idx1-1,2);
    
    clear b
       
end

% segments treated seperately with smaller 'int'

if exist('wvllow1','var')
    
    for i = 1:length(envelope)
    if envelope(i,1) > wvllow1 && envelope(i,1) < wvlup1
        
        envelope(i,1) = 0;
    end
    end
       
a1 = (wvllow1:int1:wvlup1);

for i = 1 : length(a1)
    
    dist1 = abs(a1(i)-int1/2-spectrum(:,1));
    minDist = min(dist1);
    idx1 = min(find(dist1 == minDist));
    dist2 = abs(a1(i)+int1/2-spectrum(:,1));
    minDist = min(dist2);
    idx2 = max(find(dist2 == minDist));
      
    for j = 1 : (idx2-idx1+1)        
    b(j) = spectrum(j+idx1-1,2);           
    end
    
    maxVal = max(b);
    idMax = max(find(b == maxVal));
    
    envelope(i+length(a),1) = spectrum(idMax + idx1-1,1);
    envelope(i+length(a),2) = spectrum(idMax + idx1-1,2);
    
    clear b
           
end   
end

if exist('wvllow2','var')
    
    for i = 1:length(envelope)
    if envelope(i,1) > wvllow2 && envelope(i,1) < wvlup2
        
        envelope(i,1) = 0;
    end
    end
       
a2 = (wvllow2:int2:wvlup2);

for i = 1 : length(a2)
    
    dist1 = abs(a2(i)-int2/2-spectrum(:,1));
    minDist = min(dist1);
    idx1 = min(find(dist1 == minDist));
    dist2 = abs(a2(i)+int2/2-spectrum(:,1));
    minDist = min(dist2);
    idx2 = max(find(dist2 == minDist));
   
    for j = 1 : (idx2-idx1+1)        
    b(j) = spectrum(j+idx1-1,2);           
    end
    
    maxVal = max(b);
    idMax = max(find(b == maxVal));
    
    envelope(i+length(a)+length(a1),1) = spectrum(idMax + idx1-1,1);
    envelope(i+length(a)+length(a1),2) = spectrum(idMax + idx1-1,2);
    
    clear b
           
end   
end

if exist('wvllow3','var')
    
    for i = 1:length(envelope)
    if envelope(i,1) > wvllow3 && envelope(i,1) < wvlup3
        
        envelope(i,1) = 0;
    end
    end
        
a3 = (wvllow3:int3:wvlup3);

for i = 1 : length(a3)
    
    dist1 = abs(a3(i)-int3/2-spectrum(:,1));
    minDist = min(dist1);
    idx1 = min(find(dist1 == minDist));
    dist2 = abs(a3(i)+int3/2-spectrum(:,1));
    minDist = min(dist2);
    idx2 = max(find(dist2 == minDist));
   
    for j = 1 : (idx2-idx1+1)        
    b(j) = spectrum(j+idx1-1,2);           
    end
    
    maxVal = max(b);
    idMax = max(find(b == maxVal));
    
    envelope(i+length(a)+length(a1)+length(a2),1) = spectrum(idMax + idx1-1,1);
    envelope(i+length(a)+length(a1)+length(a2),2) = spectrum(idMax + idx1-1,2);
    
    clear b
          
end   
end

if exist('wvllow4','var')
    
    for i = 1:length(envelope)
    if envelope(i,1) > wvllow4 && envelope(i,1) < wvlup4
        
        envelope(i,1) = 0;
    end
    end
    
a4 = (wvllow4:int4:wvlup4);

for i = 1 : length(a4)
    
    dist1 = abs(a4(i)-int4/2-spectrum(:,1));
    minDist = min(dist1);
    idx1 = min(find(dist1 == minDist));
    dist2 = abs(a4(i)+int4/2-spectrum(:,1));
    minDist = min(dist2);
    idx2 = max(find(dist2 == minDist));
    
    for j = 1 : (idx2-idx1+1)        
    b(j) = spectrum(j+idx1-1,2);           
    end
    
    maxVal = max(b);
    idMax = max(find(b == maxVal));
    
    envelope(i+length(a)+length(a1)+length(a2)+length(a3),1) = spectrum(idMax + idx1-1,1);
    envelope(i+length(a)+length(a1)+length(a2)+length(a3),2) = spectrum(idMax + idx1-1,2);
    
    clear b
            
end   
end

% treating outliers 

if exist('outlier1','var')
    if outlier1 > 0
dist3 = abs(outlier1 - envelope(:,1));
minDist = min(dist3);
idx3 = find(dist3 == minDist);

for i = idx3 : length(envelope)-1   
    envelope(i,1) = envelope(i+1,1);
    envelope(i,2) = envelope(i+1,2);    
end
    end
end

% creating Output

envelope = sortrows(envelope);
envelope = envelope(any(envelope(:,1),2),:);

cont = fit(envelope(:,1), envelope(:,2), 'smoothingspline');
end

