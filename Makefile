.PHONY: test
test: testrun clean

.PHONY: testrun
testrun: cdtw.so
	python -Bc "import evaluate; evaluate.test_cdtw()"

cdtw.so: setup.py cdtwmodule.c
	# Just throw the cdtw.so file in the current directory
	python setup.py install --install-lib="."

.PHONY: clean
clean:
	rm -rf build cDTW*egg-info
