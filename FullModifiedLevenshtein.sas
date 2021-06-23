libname ocr "C:\Users\jamyle\Documents\OCR";

/*Set path to folder containing .txt files as a macro variable*/
%let dirname = C:\Users\jamyle\Documents\OCR\txt\all_docs;

/*In Windows OS use DOS command "dir" with the /B option to list only file names*/
filename all_list pipe "dir /B &dirname\*.txt";


/*create list of files*/
data all_list;
  length file_name $256;
  infile all_list length=reclen;
  input file_name $varying256. reclen;
run;

data scan_list gt_list;
	set all_list;
	if scan(file_name,2,"_")="gt.txt" then output gt_list;
	else output scan_list;
run;

/*import scanned text files*/
data temp_scan(keep=file_name text encoding=utf8);
    length text line $32767;
	retain text;
	set scan_list;
	by file_name;
	if first.file_name then do;
		line='';
		text='';
	end;
  	filepath = "&dirname\"||file_name;
  	infile dummy filevar = filepath
                  end=done
                  LRECL=32767
                  ENCODING='utf-8' 
                  TRUNCOVER
                  ;
	do while (not done); 
    	input line $32767.;
		text=catx(' ',text,line); 
	end;
	if done and last.file_name then output;
run; 

/*remove unprintables*/
data temp_scan(drop=test);
  length quality $50;
  set temp_scan;

  /* Verify each byte in string until no more non-printables are found */
  do until(test=0);

    /* For SAS 9.0 and above, NOTPRINT */
    test=notprint(text);

    /* If a non-printable is found...replace it with a space */
    if test>0 then do;
      substr(text,test,1)=' ';
      end;
    end;

	/*create id variables*/
	short_file_name=scan(file_name,1,'.');
	doc=scan(short_file_name,1,'_');
	font=scan(short_file_name,2,'_');
	dpi=scan(short_file_name,3,'_');
	if dpi="300" then quality="high";
	else if dpi="150" then quality="low";
run;

/*import gt text file*/
data temp_gt;
  length  text $32767;
  set gt_list;
  filepath = "&dirname\"||file_name;
  infile dummy filevar = filepath 
  		end=done
   		LRECL=32767
        ENCODING='wlatin1' 
		TRUNCOVER
		; 
  do while(not done);
    input text $32767.; 
  end;

  	/*create id var*/
	short_file_name=scan(file_name,1,'.');
	doc=scan(short_file_name,1,'_');
run;

/*combine gt & scanned*/
data temp_all;
	length category $50.;
	set temp_scan temp_gt;
	length=length(text);
	if missing(font)then category="gt" ;
	else category=catx("_",font,quality);
run;

proc sort data=temp_all;
	by doc;
run;

/*transpose to wide for edit distances*/
proc transpose data=temp_all out=all_transpose(drop=_name_) prefix=text_;
	id category;
	by doc;
	var text;
run;



/* Imports All_transpose dataset into IML */

/*parse each word in gt and predicted into vectors*/
proc iml worksize=4000000;

use Work.All_transpose;

read all var {doc text_arial_low text_arial_high text_tnr_low text_tnr_high text_gt};
/* Concats the vectors into a matrix that can be used in IML */
docs = text_arial_low || text_arial_high || text_tnr_low || text_tnr_high || text_gt;


/* Just using random GT and Predicted cells to test accuracy on entire document
gt = docs[1,5];

pd = docs[1,4]; */


gt = Hello World!;

pd = Helo Worlo!;


/* Makes a module called chop. Module is just a function*/
start chop(s);
  return ( substr(s, 1:nleng(s), 1) );   /* SAS/IML 12.1 */
finish;

/* Creates vectors with one char in each element */
pd_chopped = chop(pd);
gt_chopped = chop(gt);

/* Extracts number of characters from each document */
n1 = ncol(pd_chopped);
n2 = ncol(gt_chopped);

/* Creates matrix with n1+1 rows, n2+1 columns, with all zeros */
d = j(n1 + 1, n2 + 1, 0);

/* Fills the matrix, d, with the positions */
do i = 1 to n1;
	d[i, 1] = i;
	end;
do j = 1 to n2;
	d[1, j] = j;
	end;
/* Calculates the cost of each edit and stores them */
do i = 1 to n1;
	do j = 1 to n2;
		if pd_chopped[i] = gt_chopped[j] then cost = 0;
		else cost = 2;
		d[i + 1, j + 1] = min(d[i,j + 1] + 1, /*insert*/
		d[i + 1,j] + 1, /*delete*/
		d[i,j] + cost); /*replace*/
	end;
end;

/* */
shape = dimension(d);

/* i = rows and j = cols */
i = shape[1];
j = shape[2];
list_edits = {};
string_1_loc = {};
string_2_loc = {};
iterations = 1;

/* Translates cost (numeric) into an action*/
do until(i = 1 | j = 1 | iterations > 500000);
iterations = iterations + 1;
index_tar = d[i-1,j-1] || d[i,j-1] || d[i-1,j];
index = loc(index_tar = min(index_tar))[1];
	if index = 1 then do;
		if d[i,j] > d[i-1,j-1] then do;
			list_edits =  'replace' || list_edits ;
			string_1_loc = i - 1 || string_1_loc;
			string_2_loc = j - 1 || string_2_loc;
		end;
		i = i - 1;
		j = j - 1;
		end;
	else if index = 2 then do; 
			list_edits =  'insert' || list_edits ;
			string_1_loc = i - 1 || string_1_loc;
			string_2_loc = j - 1 || string_2_loc;
			j = j - 1;
			end;
	else if index = 3 then do;
			list_edits =  'delete' || list_edits ;
			string_1_loc = i - 1 || string_1_loc;
			string_2_loc = j - 1 || string_2_loc;
			i = i - 1;
			end;
end;

print i j iterations index;
print list_edits string_1_loc string_2_loc;
print pd_chopped, gt_chopped;
pd_shift = 0;
gt_shift = 0;

/* Executes edits */
do i = 1 to ncol(list_edits);
	edit = list_edits[i];
	pd_loc = string_1_loc[i];
	gt_loc = string_2_loc[i];
	if edit = "delete" then do;
		gt_chopped = gt_chopped[1:(gt_loc + gt_shift)]` || '¥' || gt_chopped[(gt_loc+gt_shift+1):ncol(gt_chopped)]`;
		gt_shift = gt_shift + 1;
	end;
	if edit = "insert" then do; 
		pd_chopped = pd_chopped[1:(pd_loc + pd_shift)]` || '¥' || pd_chopped[(pd_loc+pd_shift+1):ncol(pd_chopped)]`;
		pd_shift = pd_shift + 1;
	end;
end;

print pd_chopped, gt_chopped;

/* Transpose them to add to dataset */
pd_chopped = pd_chopped`;
gt_chopped = gt_chopped`;

/* Creates dataset in temp that  */
create temp var {pd_chopped gt_chopped};    /* a literal matrix of names      */
append;
close temp;

proc freq data=work.temp;
tables gt_chopped*pd_chopped / norow nocol nopct;
run;

quit;


/*  This works, but I just wanted to comment it out for now.

proc sql ;
	CREATE TABLE should_be_counts AS
	SELECT gt_chopped, count(*) AS GT_CT FROM work.Temp 
	GROUP BY gt_chopped;
quit;

proc sql ;
	CREATE TABLE actually_are_counts AS
	SELECT  gt_chopped, count(*)  FROM temp WHERE gt_chopped = pd_chopped 
	GROUP BY gt_chopped;
	quit;

*/
