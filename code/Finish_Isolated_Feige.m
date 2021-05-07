% Input: 'Isolated_Lines1.txt' the list of isolated lines, and 'Elem,Ion_FERTIG.txt' for all IG elements, the line lists created with 'Make_Linelists'

% Output: 'ANS' corresponding transitions for all isolated lines

% ---------------------- START -- CODE ----------------------

% reading the isolated list 

A = readcell('Isolated_Lines1.txt');
B(:,1) = cell2mat(A(:,1));
B(:,2) = cell2mat(A(:,3));

% for each isolated line, reading the corresponding list of transitions

z = 1;

for i = 1:length(B)
    
    center = B(i,1);
    sigma = B(i,2);
    ion = char(A(i,2));
    elem = ion(1:2);
   
    C = importdata(strcat(elem,'/',ion,'_FERTIG.txt'));
    
    wl = C(:,1);
    
    % a transition is considered isolated if only one occurs in pm 1 sigma around the line center
    
    a = 0;
   
    for j = 1:length(wl)
            
     if abs(center - wl(j)) < sigma 
         
         a = a+1;
         k = j;
         
     end     
    end
    
    if a == 1
          
        ANS{z,1} = C(k,1);
        ANS{z,2} = ion;
        ANS{z,3} = C(k,2);
        ANS{z,4} = C(k,3);
        ANS{z,5} = C(k,4);
        
        z = z+1;
               
    end
    clear wl   
end

