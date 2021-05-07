% Input: 'Elem,Ion_Lines.txt' the prepared line list of a chosen ion (containing the line properties), 'Elem,Ion.POS' the Kurucz POS list of the ion, and 'Elem,Ion_Errors_Final.txt' a list containing lines with a goodness < 0.96

% Output: 'G' a matrix containing all identified transitions of the specified ion, with the applied reduction criteria

% specifying the ion of which the lines are looked for 

Elem = 'Ni';
Ion = 'IV';

% ---------------------- START -- CODE ----------------------

% reading the line list, the Kurucz POS list and the error list

A = importdata(strcat(Elem,'/',Elem,Ion,'_Lines.txt'));
B = importdata(strcat(Elem,'/',Elem,Ion,'.POS'));

if isfile(strcat(Elem,'/',Elem,Ion,'_Errors_final.txt'))      
C = importdata(strcat(Elem,'/',Elem,Ion,'_Errors_final.txt'));
end

% reducing the POS list, removing transitions with log(gf)< -2.5, and just keeping the interval 1150A < lamda < 1710A

q=1;

for i = 1:length(B(:,1))
    
    if q < length(B(:,1))+1
    
    if B(i,1) > 1150 && B(i,1) < 1710
    if B(i,2) > -2.5
        
    linelist(q) = round(B(i,1),2);
    gflist(q) = B(i,2);
    E_low(q) = B(i,4);
    E_up(q) = B(i,5);
    
    q = q+1;
    
    end
    end
    end  
end

% removing lines with equw < 1.5mA 

q=1;

for i = 1:length(A(:,1))
    
    if q < length(A(:,1))+1
    
    if A(i,1) > 1150 && A(i,1) < 1710
    if A(i,2) > 1.5
        
    linessyn(q) = A(i,1);
    equw(q) = A(i,2);
    
    q = q+1;
    
    end
    end
    end    
end

% A transition in the Kurucz POS list is taken when a corresponding line exists in the prepared list (line center is closer than 0.05A)

a = 1;
for m = 1:length(linelist)
    
 for n = 1:length(linessyn)   
 if  abs(linelist(m) - linessyn(n)) < 0.05 
     F(a,1) = linelist(m);
     F(a,2) = gflist(m);
     F(a,3) = E_low(m);
     F(a,4) = E_up(m);
     
     a = a+1;
 
 end    
 end  
end

% for all line centers in the error list, transitions from the POS list will be taken in a certain interval (0.05A to 0.15A) depending on the goodness-of-fit

 if isfile(strcat(Elem,'/',Elem,Ion,'_Errors_final.txt'))

m = length(F(:,1))+1;

for j = 1:length(C(:,1))
    
    errcorr = 1/8 * sqrt(C(j,2)) + 1/40;
    
for k = 1:length(linelist)
   
    if abs(linelist(k)-C(j,1))< errcorr
        
    F(m,1) = linelist(k);
    F(m,2) = gflist(k);
    F(m,3) = E_low(k);
    F(m,4) = E_up(k);  
    
    m = m+1;
        
    end       
end    
end
 end

% Arranging F to wavelength, deleting double entries

G = sortrows(F);
clear F

for a = 1:length(G(:,1))-1
    
    if G(a,1) == G(a+1,1) && G(a,2) == G(a+1,2) && G(a,3) == G(a+1,3) && G(a,4) == G(a+1,4)
    
    G(a,:) = 0; 
        
    end    
end    

F = sortrows(G);
clear G

G = F(any(F,2),:);

% reduction: for lines with log gf > -1, all nearby (< 0.03A) weaker transitions (delta log gf > 1.5) are removed

for i = 1:length(G(:,1))
    
  if G(i,2) > -1.0  && G(i,1) > 0
      
     for j = 1:length(G(:,1))
         
     if abs(G(i,1)-G(j,1))<0.03 && (G(i,2)-G(j,2))>1.5
         
         G(j,:)=0;
     end    
         
     end    
      
  end    
    
end    

F = sortrows(G);
clear G

G = F(any(F,2),:);

