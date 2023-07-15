#!/bin/zsh

# Run Ising.c 
make ising
./ising

# Run Python to show plot 
python ising_show.py
