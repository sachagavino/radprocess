#!/usr/bin/env bash
#   This file is part of PyMSES.
#
#   PyMSES is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PyMSES is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.

# cd pymses base directory
# update contact informations in ./doc/intro.rst (emails) and external web links
# update all doc text !
# update ./Release_notes.txt
# update screenshots in ./doc/_images/?
make all_tests #(nosetest needed)
# change version number in ./setup.py
# update pymses/doc/downloads.rst
python update_version_number.py
# commit and push release change to bitbucket !
make tarball
NEW_TAR_FILES=$(find . -maxdepth 1 -iname "pymses_*")
# move the new file
for NEW_TAR_FILE in $NEW_TAR_FILES; do
	echo "Moving $NEW_TAR_FILE --> ./doc/$NEW_TAR_FILE"
	mv "$NEW_TAR_FILE" "./doc/$NEW_TAR_FILE" || (echo "Error found while moving tar files ! "; exit 1)
done
cd doc
make html
python convert_png_to_pdf.py
make html
make latex (a fine latex installation is needed)
cd ..
make doc_pdf
# cp ./doc/_build/latex/PyMSES.pdf ./PyMSES_v4.1.0.pdf
# cp ./doc/PyMSES_v4.1.0.pdf ./doc/_build/html/_downloads/PyMSES_v4.1.0.pdf
# cp ./doc/pymses_4.1.0.tar.gz ./doc/_build/html/_downloads/pymses_4.1.0.tar.gz
# Open the pymses/doc/_build/html/intro.html file
# check that the website is working neatly : images, links, downloads, scripts, pdf from matplotlib examples !
# check http://irfu.cea.fr/Projets/PYMSES/ug_visu_raytrace.html#density-map webpage
# check directories names and scripts!
