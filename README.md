# ChargedSelfAssembly
## 1. Raw Data Required For Calculations 

Several types of raw data are required to perform the calculations. The methods for obtaining these raw data files are described in the ChiralSelfAssembly repository. You will need: 

- `BulkResData_Type*_|A*|A*|A*_Mammmals.mat`, specified in the variable `FPathRes`.
- `EffDVecs_Mammals_Type*|A*|A*|A*.mat`, specified in the variable `FPathS`.

## 2. Calculation of Pairwise Interaction Energies (pKa/pH dependent)

The pairwise molecular interactions produced by these scripts are different from those in the ChiralSelfAssembly repository. New interaction data is produced by running the script `getIntCurves_pH_pKa.m`, remembering to specify the following variables: 

- `ColType`, collagen type input as a roman numeral.
- `TrimerType`, input in the format `"AX|AY|AZ"`, where X, Y and Z are positive integers.
- `Species`, the species that you want to run the calculation for. For example, for bovine collagen, set `Species = "BosTaurus"`.
- `FPathS`, `FPathRes`, discussed above.
- The path to the directory `.../scripts/methods/...`.
- `FPathpKa`, which corresponds to the directory where you want to store the output file (see below).

You may want to run the calculation above in parallel. To do this, make appropriate amends in L58-61 in the file `pKaGrid.m`, which is stored in  `.../scripts/methods/...`.

Following the instructuons above, you will get the output

- `Species|Type*|A*|A*|A*|pKaGrid.mat`, where the "Species" part is replaced in accordance with how you specified the variable `Species` above.

The output above contains both the interaction energies $`\textrm{U}_{i\textrm{-}j}(\Delta z)`$ as well as the discretisation of $`\Delta z`$. To access these data at a given set of pKa values, load the file `Species|Type*|A*|A*|A*|pKaGrid.mat` to a variable, say, `out`. Enter `out.S.muMeasure_*.Energy` or `out.S.muMeasure_*.Energy.Staggers` to access the pairwise interaction energy or the discretised axial stagger respectively. The star in `out.S.muMeasure_*.Energy.Staggers` is replaced by the appropriate value of the overlap parameter $`\mu(W_{\textrm{R}_{(1)}},W_{\textrm{R}_{(2)}})`$. The convention for using the numerical values of the overlap parameter in the structure field name is as follows:

- $-0.5\rightarrow \textrm{minus}0\textunderscore5$ for negative values.
- $0.2\rightarrow 0\textunderscore2$ for non-negative values.

The specified value of the overlap parameter applies to both acidic and basic residues. This means that `out.S.muMeasure_minus0_5.Energy` will output a data structure containing $`\textrm{U}_{i\textrm{-}j}(\Delta z)`$ for $`\mu(W_{\textrm{Glu}},W_{\textrm{Asp}}) = -0.5`$ at acidic pH values and for $`\mu(W_{\textrm{Arg}},W_{\textrm{Lys}}) = -0.5`$ at basic pH values. 

## 3. Calculating Equilibrium Statistics 

To calculate the equilibrium statistics, which are required for reproducing the main Figures in the paper, you will need to run the script `runProbsOnpKaGrid.m`, which is located in the directory `.../scripts/eqbmStatistics/...`. You will need to specify the following:

- Variables `ColType`, `TrimerType` and `Species`.
- `FPathpKa` and the path to the directory `.../scripts/methods/...`.
- Variable `FPathEqbm`, which specifies where you will store the output file (see below).

If you want to run this calculation in parallel (which I would highly recommend), the make appropriate changes on L35-39

This script produces a new type of output file: 

- `EqbmProbspKaGrid__Species|Type*|A*|A*|A*.mat`

The output contains the equilibrium probability $`\mathcal{P}_{\textrm{period}}(g)`$ and microfibril energy $\textrm{E}_{\textrm{M}}(g)$.


## 4. Generating Main Figures 

All of the relevant scripts are found in the directory `.../scripts/mainTextFigures/...`. The scripts are all in the format `genMainFig*.m`. You will need to specify the following:

- File paths in variables `FPathEqbm`, `FPathRes` and the file path to the directory `.../scripts/methods/...`.
- The variable `speciesSplit`, which is the same as `Species` variable with the following change: $`\underbrace{\textrm{BosTaurus}}_{\textrm{Species}}\rightarrow \underbrace{\textrm{Bos taurus}}_{\textrm{speciesSplit}}`$.

## 5. Generating Appendix Figures 

## 6. Computing Environment 

- **OS:** macOS Monterey 12.7.6
- **MATLAB Release:** R2024b (**note, for `genMainFig5.m` I strongly advise using an earlier release of MATLAB, such as R2022a**)  
- **Toolbox dependencies:** Bioinformatics Toolbox, Image Processing Toolbox, Optimization Toolbox, Parallel Computing Toolbox (optional, but recommended), Signal Processing Toolbox, Statistics and Machine Learning Toolbox, Symbolic Math Toolbox.

## 7. References

[1] `distinguishable_colors.m` script obtained from https://github.com/cortex-lab/MATLAB-tools/blob/master/distinguishable_colors.m.

[2] `fillout.m` script obtained from https://uk.mathworks.com/matlabcentral/fileexchange/10515-fillout.

[3] `freezecolors.m` script obtained from https://uk.mathworks.com/matlabcentral/fileexchange/7943-freezecolors-unfreezecolors

[4] S Kawashima, M Kanehisa, Aaindex: amino acid index database. Nucleic Acids Res. 28, 374–374 (2000).

[5] All collagen residue sequences are obtained from the NCBI RefSeq database. 

[6] K Rainey, MC Goh, A statistically derived parameterization for the collagen triple-helix. Protein Sci. 11, 2748–2754 (2002).
