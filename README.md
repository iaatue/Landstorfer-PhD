Doktorarbeit Alexander Landstorfer, Tübingen, 2021


Inhalt



code/

The code/ files are written as Matlab Code (tested on R2020b), comments are indicated 
by a '%' at the beginning of a line. The basic build is the following:

% code description, input and output description  
additional input arguments  
% -start-code-  
input data  
code  
output data  

The file names correspond to the files mentioned in the thesis's appendix.


Equ_width_Feige.m,
Feige_make_Contfact.m,
Finish_Isolated_Feige.m,
FUNmakeCont.m,
Make_Isolated_Feige.m,
Make_Linelists_Feige.m,
Prepare_Linelists_Feige.m



tables/

In the tables/ files, comments are indicated by a '#'. Table columns are separated 
by ','. The column headers are given before the table entries start. The basic build 
is the following:

#Table description  
# 
#column headers  
# 
table  

Tables 25-30 correspond to the thesis's tables. Feige_autom_ALL and PG_autom_ALL 
comprise all found transitions in the respective stellar spectrum.


Feige_autom_ALL.txt,
PG_autom_ALL.txt,
Table_25_rel_lines.txt,
Table_26_Fei_iso_autom.txt,
Table_27_PG_iso_autom.txt,
Table_28_Fei_iso_manual.txt,
Table_29_Cu.txt,
Table_30_Zn.txt

