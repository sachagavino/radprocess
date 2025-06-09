PYSETUP_FLAGS=
PYMSES_BASE=pymses
UNIT_TEST_DIR=unit_tests
LONG_TEST_DIR=long_tests

# This compiles all modules for an inplace use
inplace_build:
	python setup.py build_ext --inplace $(PYSETUP_FLAGS)

# This generates the C sources from Cython sources
cython:
	find . -name "*.pyx" | xargs cython

clean:
	find . -type f -name "*.pyc" -exec rm {} \;
	find . -type f -name "*.so" -exec rm {} \;
	rm -rf build

# This runs long integration tests
long_tests:
	cd $(LONG_TEST_DIR) && nosetests -v -w . && cd -


# This runs fast unit tests
unit_tests:
	cd $(PYMSES_BASE)/$(UNIT_TEST_DIR) && nosetests -v -w . && cd -

# This runs all tests
all_tests: unit_tests long_tests


################################################
#           Developer targets                  #
################################################
# This builds a tarball suitable for distribution
tarball: 
	./maketarball.sh

# This builds the documentation (html and latex)
doc_pdf:
	cd doc && make latexpdf

doc_html:
	cd doc && make html

docs: doc_pdf doc_html

# This runs pylint
pylint:
	pylint --rcfile=./.pylintrc $(PYMSES_BASE)

