##
# ##
# stereo
#
# @file
# @version 0.0.1

# find python3
PYTHON=`/usr/bin/which python3`

# our testing targets
.PHONY: tests flake black isort all

all: mypy isort black flake tests

tests:
	${PYTHON} -m pytest --cov=stereo tests

flake:
	${PYTHON} -m flake8 stereo

black:
	${PYTHON} -m black -t py37 stereo tests

isort:
	${PYTHON} -m isort --atomic stereo tests

mypy:
	${PYTHON} -m mypy stereo

# end
