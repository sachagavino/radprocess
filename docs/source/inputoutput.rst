.. _chap-input-files:

The input and output files of radprocess
***************************************


Good practices
=======================
The RAMSES output folder should always be output_XXXXX and the file names should be consistent with it.
since it is installed with pip install -e ., a simple git pull will update the package. 

we need to commit things with pymses, because the octree.py needed some change to make it work with mf:

line 421 has to be if amr_file_type in ("grav", "mf"):

Some in any way we should keep this version somewhere:
cd /data/pebbles/software/pymses/pymses3
git add pymses/sources/ramses/octree.py
git commit -m "Fix mf file reading: treat mf files like gravity (is_hydro=0)"