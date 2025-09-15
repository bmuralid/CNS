#!/usr/bin/bash
#
SHELL :=/bin/bash

all: 
	$(error please pick a target)

env:
	# source the local .env file
	# @if [ -f .env ]; then \
	# 	export $(cat .env | xargs) ; \
	# fi
	source .env

clean:
	rm -rf build
