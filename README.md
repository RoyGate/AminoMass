# AminoMass
Define which amino acids might compose a lipopeptide, given a single, simple fatty acid. Although somewhat universally built, it was initially designed to evaluate surfactins.

Written and tested with Python 3.9, although might run in 2.X as well.

In a computer without a Python interpreter, it is recommended to convert/compile the .py file in .exe (e.g. PyInstaller, py2exe).

Another quick-fix is to copy-paste the source code into a online compiler (e.g. Programiz, Online-Python).

IMPORTANT NOTE: 

The comments below are also in the .py file, but as they might get hidden after compilation, they are also added here.

The software was developed to analyze lipopeptides, mostly surfactins. Thus, its default settings revolve around them, i.e.:
  Total amino acids (a.a.): 7
  Fixed a.a.: 4 (glutamic acid, 2 x leucine, aspartic acid) --> only 03 a.a. are left to deduce.
  Search tolerance error: 0.5 Da (arbitrary)
  
Also, some presumptions were made:
  - There is ONLY ONE fatty acid moeity in the structure;
  - No derivatization (e.g. esterification) exists in the a.a.
