#!/bin/sh

CC=gcc

$CC hog.c main.c -o hogc -lpthread -lm 

