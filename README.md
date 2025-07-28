# ChargedSelfAssembly
## 1. Raw Data Required For Calculations 

Several types of raw data are required to perform the calculations. The methods for obtaining these raw data files are described in the ChiralSelfAssembly repository. You will need: 

- `BulkResData_Type*_|A*|A*|A*_Mammmals.mat`, specified in the variable `FPathRes`.
- `EffDVecs_Mammals_Type*|A*|A*|A*.mat`, specified in the variable `FPathS`.

## 2. Raw Data Produced By Calculations 

The pairwise molecular interactions produced by these scripts are different from those in the ChiralSelfAssembly repository. New interaction data is produced by running the script `getIntCurves_pH_pKa.m`, remembering to specify the following variables: 

- `ColType`, collagen type input as a roman numeral.
- `TrimerType`, input in the format `"AX|AY|AZ"`, where X, Y and Z are positive integers.
- `Species`, the species that you want to run the calculation for. For example, for bovine collagen, set `Species = "BosTaurus"`.
- `FPathS`, `FPathRes`, discussed above.
- Specify the path to the directory `.../scripts/methods/...`  

Following the instructuons above, you will get two new outputs

- `IntCurves_Mammals_pH_Type*|A*|A*|A*.mat`, which contain the interaction energies $\textrm{E}_{i\textrm{-}j}(\Delta z)$ for different pH and pKa values.
- `SVecs_Mammals_pH_Type*|A*|A*|A*.mat`, which contain the discretised stagger vector, corresponding to discretised axial stagger $\Delta z$.
