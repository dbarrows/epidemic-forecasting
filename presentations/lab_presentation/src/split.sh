#!/bin/bash

echo Splitting $1
pdftk $1 cat odd output slides.pdf
pdftk $1 cat even output notes.pdf