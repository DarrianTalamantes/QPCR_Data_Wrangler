The Idea of this project is to get data straight from the QPCR abs quant analysis and turn it into something usefull to me. 

Step 0: Run the program called multiply_List.py found in the List_Multiplyer directory to turn your ordered list of samples into triplicates that can easily be pasted on to label all triplicate samples.
Step 1: Run program QPCR_Data_Wrangler1.sh It will  ask for your data file pathway. Input it. Hit enter.
Step 2: Program will ask you to input qualatative data into the new four columns that are blank. The concentration column only needs to be filled out for the standards. Standards are labeled "stdX" X being a unique identifier for each standard. Recomend to use excel. Save this in a tab or comma delimiated file. File must have same name as it originally had.
Step 3: Run program QPCR_Data_Wrangler2.sh which will create graphs and place them in the output folder
Step 4: Rejoice your high quality R graphs are made!
