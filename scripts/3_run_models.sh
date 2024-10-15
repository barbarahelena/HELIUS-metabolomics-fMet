#!/bin/bash
./XGBeast.py -path SBP/Fem -name fem_sbp -x reg -n 200 -param metabolite_grid.json
./XGBeast.py -path SBP/Male -name male_sbp -x reg -n 200 -param metabolite_grid.json
./XGBeast.py -path SBP/Fem -name fem_sbp -x reg -n 200 -param metabolite_grid.json -permute
./XGBeast.py -path SBP/Male -name male_sbp -x reg -n 200 -param metabolite_grid.json -permute

./XGBeast.py -path DBP/Fem -name fem_dbp -x reg -n 200 -param metabolite_grid.json
./XGBeast.py -path DBP/Male -name male_dbp -x reg -n 200 -param metabolite_grid.json
./XGBeast.py -path DBP/Fem -name fem_dbp -x reg -n 200 -param metabolite_grid.json -permute
./XGBeast.py -path DBP/Male -name male_dbp -x reg -n 200 -param metabolite_grid.json -permute
