.PHONY: all check clean

objects = $(wildcard *.in)
outputs := $(objects:.in=.txt)

all: $(outputs)

%.txt: %.in
	python check-python-version.py
	pip-compile -v --output-file $@ $<

test.txt: base.txt

dev.txt: test.txt

check:
	@which pip-compile > /dev/null

clean:
	- rm -f *.txt
