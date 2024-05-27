# cell-survival-project
Code for Zhenlun's summer 2023 & summer 2024 cell survival rate project.

playground.py is the work done in summer 2023, and updated in summer 2024. 
newattempt.py is a new, recoded version of playground.py, done in summer 2024.
pideplot.py is a simple data visualizer, working with data in Particle Irradiation Data Ensemble Version 3.2.
The rest are Medras MC files found here: https://github.com/sjmcmahon/Medras-MC

files in SDDeditor folder are files used to process and edit the SDD data used. This includes:
  sddeditor.py: header (delete wording such as 'variable' and 'non specified') and first data entry (first data start with 2 instead of 1) editing.
  sddeditor2.py: changing SDD data field 4(chromosome position) from basepair to fractional position along the chromosome.
  sddeditor3.py: header editing (Dose or fluence).
  sddeditor4.py: header editing (Data entries).
  empty remover.py: removes files with nothing in it (no header, no data, no nothing).
  _nonsdd remover.py: removes non SDD files.
