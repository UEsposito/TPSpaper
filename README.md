# Working with the Temporal Population Structure (TPS) tool

TPS is a method to date ancient genomes solely from their genomic sequence.
In the following we provide the code to reproduce the results reported by Esposito et al.'s study

There are 3 folders:
  - data/ contains an  excel file with the input data for TPS.
  - src/ contains the Matlab scripts necessary to run TPS and to reproduce the results of the study.
  - results/ contains the output of TPS as published in the study. The folder incloudes an excel table, png figures, and a mat file.
  
  To run TPS do the following:
    1. Copy the folder structure to your computer.
    2. Execute the function Paper_results_main.m. This is the main function which loads the excel file in data/ and outputs the results to results/.
    3. The running time is 2-3 hours.
  
  Wednesday 30th October 2019
  Umberto Esposito
  umberto.espos@gmail.com
