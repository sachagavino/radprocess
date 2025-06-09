#!/bin/bash
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

# Retrieve the mercurial tag and hex id
HG_TAG=$( hg id -t )
HG_ID=$( hg id --id | sed -e 's/[^0-9a-f]//' )
HG_N=$( hg id -n )


HG_LABEL=$HG_N
echo "Using revision : $HG_LABEL"


# Create a temp dir
SCRATCH=$(mktemp -d)
echo "Using scratch dir: $SCRATCH"
PY_VER=$( grep "__version__" pymses/__init__.py | cut -d\' -f 2 )
TAR_NAME="pymses_"$PY_VER
TAR_DIR="$SCRATCH/$TAR_NAME"
mkdir "$TAR_DIR"

# Create an unversioned copy in the tarball dir.
# Make sure we exclude maketarball.sh,
#					   TODO,
#					   .hgignore,
#					   doc/*, assets/, .hg_archival.txt
echo "Exporting unversioned copy"
hg archive "$TAR_DIR" \
	-X "maketarball.sh" \
	-X "TODO" \
	-X ".hgignore" \
	-X ".hgtags"\
	-X ".pylintrc"\
	-X "how_to_do_a_new_PyMSES_release.sh"
rm -rf $TAR_DIR/doc/*
rm -rf $TAR_DIR/.hg_archival.txt
head -32 Makefile > $TAR_DIR/Makefile

# update version number:
CURRENT_DIR=$(pwd)
echo "cd $TAR_DIR"
cd $TAR_DIR
mkdir .hg
if [ ! -e $CURRENT_DIR/.hg/tags.cache ]; then
	cp $CURRENT_DIR/.hg/cache/tags .hg/tags.cache
else
	cp $CURRENT_DIR/.hg/tags.cache .hg/tags.cache
fi
python update_version_number.py
rm -f update_version_number.py
rm -rf .hg

echo "cd $CURRENT_DIR"
cd $CURRENT_DIR

# Handling Cython source files
#echo "Building Cython "
#(find . -name "*.pyx" | xargs --verbose cython) > /dev/null || (echo "make cython failed"; exit 1)

# List all the C filenames produced from Cython
C_FILES=$(find . -iname "*.pyx" | sed -e 's/\.pyx$/\.c/')
echo -e "Will pack Cython output:\n$C_FILES"

# Copy the Cython C files in there
for C_FILE in $C_FILES; do
	echo "Copying $C_FILE --> $TAR_DIR/$C_FILE"
	cp "$C_FILE" "$TAR_DIR/$C_FILE" || (echo "Error found while copying C source files ! "; exit 1)
done


## Generate the docs (PDF)
#echo "Building documentation..."
#make doc_pdf > /dev/null || (echo "make doc_pdf failed"; exit 1)
#
## Copy the built PDF doc
#echo "Copying PDF manual in '$TAR_DIR/doc' directory..."
#cp doc/_build/latex/PyMSES.pdf "$TAR_DIR/doc"
rmdir $TAR_DIR/doc

# Make the tarball and purge
TAR_FILE="$TAR_NAME.tar.gz"
tar  -C "$SCRATCH" -czf "$TAR_FILE" "$TAR_NAME"
rm -rf "$SCRATCH"

echo
echo "Produced $TAR_FILE"
