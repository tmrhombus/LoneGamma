#!/bin/bash

grep -v -F -x -f ../data.h ../mc.h > xx_mc.txt
grep -F -x -f ../data.h ../mc.h > xx_common.txt
grep -v -F -x -f ../mc.h ../data.h > xx_data.txt

