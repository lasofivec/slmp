Requirements:
============

  * python (< 3.)
  * numpy
  * scipy
  * prettyplotlib
  * matplotlib
  * igakit (See https://bitbucket.org/dalcinl/igakit)
  * caid (See https://github.com/ratnania/caid)
  * You will also need to execute the following command:
    	$ cp thesisstyle.mplstyle ~/.matplotlib/stylelib/
    (you might need to create this path)

To execute main program:
====
```
  cd main_program
  python main.py <$INPUT_FILE>
```

How to create a new Geometry and test
----------------------------------------

1. Create Geometry using CAID;
2. Export the created geometry:
   a. Click the geometry;
   b. Click "Export";
   c. Choose path to the domains directory;
   d. Save with xml extension;
4. Edit the input_options in main_programs/input_files:
   a. Go to the domain section;
   b. Choose a code name for the geometry;
   c. Associate it to a unique key
5. Edit the get_geometry function in reading_geometry;
6. Create your input file
