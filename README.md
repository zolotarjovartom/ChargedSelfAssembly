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

Following the instructuons above, you will get the output

- `Species|Type*|A*|A*|A*|pKaGrid.m`, where the "Species" part is replaced in accordance with how you specified the variable `Species` above.

The output above contains both the interaction energies $\textrm{E}_{i\textrm{-}j}(\Delta z)$ as well as the discretisation of $\Delta z$. To access these data at a given pKa value, load the file `Species|Type*|A*|A*|A*|pKaGrid.m` to a variable, say, `out`. Enter `out.S.muMeasure_*.Energy` or `out.S.muMeasure_*.Energy.Staggers` to access the pairwise interaction energy or the discretised axial stagger respectively. The star in `out.S.muMeasure_*.Energy.Staggers` is replaced by the appropriate value of the overlap parameter $\mu(W_{\textrm{R}_{(1)}},W_{\textrm{R}_{(2)}})$. The convention for using the numerical values of the overlap parameter in construction of structure field names is as follows $-0.5\rightarrow \textrm{minus}0\textunderscore5$ for negative values and $0.2\rightarrow 0\textunderscore2$ for non-negative values. The specified value of the overlap parameter applies for both acidic and basic residues. 

