The Idea of this project is to get data straight from the QPCR abs quant analysis and turn it into something usefull to me. 

Step 0: Run the program called multiply_List.py found in the List_Multiplyer directory to turn your ordered list of samples into triplicates that can easily be pasted on to label all triplicate samples.
Step 1: Run program QPCR_Data_Wrangler1.sh It will  ask for your data file pathway. Input it. Hit enter.
Step 2: Program will ask you to input qualatative data into the new four columns that are blank. The concentration column only needs to be filled out for the standards. Standards are labeled "stdX" X being a unique identifier for each standard. Recomend to use excel. Save this in a tab or comma delimiated file. File must have same name as it originally had. Water sample must be labled "water" if it has no cp value input 0 (this will be removed automatically in the next step."
Step 3: Run program QPCR_Data_Wrangler2.sh which will create graphs and place them in the output folder
Step 4: Rejoice your high quality R graphs are made!
Step 5: Rename the graphs and data that you want to save and then run steps 1-4 again with the endophyte of the plant data. The data is called Sample_Data.csv
Step 6: Run the program called QPCR_Data_Wrangler_Relative_Biomass.sh and input the data for the endophyte and the plant into it.
Step 7: Rename the output (Data_for_project/Biomass_Data/Biomass_Data.csv) to something else
Step 8: Rejoice once more!
