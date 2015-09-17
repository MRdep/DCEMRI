# Joint-fitting
IDL code for joint fitting to variable flip angle and dynamic-contrast enhanced magnetic resonance imaging (DCE-MRI) data
# What the software does
The folder 'Joint Fitting' contains source files and synthetic data to replicate experiments described in 'Improved accuracy and precision of tracer kinetic parameters by joint fitting to variable flip angle and dynamic contrast enchanced data".

#How to run the software
1. clone the repository to your local machine, or simply download. Save to a suitable directory.
2. Download IDL virtual machine from http://www.exelisvis.com/Support/HelpArticlesDetail/TabId/219/ArtMID/900/ArticleID/12395/The-IDL-Virtual-Machine.aspx (not needed if you have an IDL development licence)
2. Open IDL and set current working directory to the location of the files.
3. Run batch_run.pro from the IDL console


#Understanding what happened
The software replicates a Monte Carlo experiment presented in the above paper. In short, it fits signal models to synthetic variable flip angle and DCE-MRI images under 3 experimental conditions (accurate arterial input function (AIF), underestimated AIF and overestimated AIF).It does this using standard sequential fitting (variable flip angle and dynamic signal models are fit one after the other) and jointly (models are fit simultaneously). The variable flip angle and dynamic signal models are based on the spoiled gradient recalled echo equations, and the dynamic signal model describes contrast kinetics using the two-compartment exchange model. The software saves the estimated parameter maps in the current working directory.   



