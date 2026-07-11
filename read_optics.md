# refactor of read_optics.py

I want to refactor this script. 

I need to: 

1. Read the data-abs.txt output from molcas from the current directory. This gives an  list of intial and final states and oscillator strengths (MD, ED, Total). 

2. Provide a pair of initial and a pair of final states for the printing/plotting range. 
Also a name for the output .txt and .pdf files. 

3. Note that I would like to be able to plot emission, where the final states are lower than initial. 
In that case, use the absolute value of the energy. 

Use /home/users/mfr24/calculations/molcas/y2o3/er-site-y203-C2-Z as the test directory, with file data-abs.txt

Small changes: 
1. make output file .out not .txt

2. Start indexing at 1, as this is what we use for all output. 

3. Run tests using 
1,2 to 17,39 
17,18 to 1,16. 

4. On the stick plots the ED and MD are obscured by the total. Can we put at dot at the top of each stick? This dot would also be useful when the oscillator strength is zero, as we would still see it. 


