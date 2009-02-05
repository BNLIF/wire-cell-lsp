#!/bin/sh

svn stat -v $1 | awk '{print $3}'
