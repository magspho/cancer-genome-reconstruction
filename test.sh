#!/bin/bash

greeting () {
  echo "Hello $1"
  echo $#
  echo "$*"
  echo $@
}

greeting John Joe Doe