#!/bin/bash
set -ex

# Run this script to open the project in codepod
# To install codepod: pip install --upgrade codepod
# You must also have docker installed
# Once in codepod, you can, for exampe, open vscode: code .

OPTS=""

eval "codepod -g -w $PWD $OPTS $@"
