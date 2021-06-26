# Benthic-species-model
Brückner MZM, Schwarz C, Coco G, Baar A, Boechat Albernaz M, Kleinhans MG. 
Benthic species as mud patrol - modelled effects of bioturbators and biofilms
on large-scale estuarine mud and morphology. 
Earth Surf. Process. Landforms. 2021;1–17. 
https://doi.org/10.1002/esp.5080

This respository contains the model code used in the above publication. The Delft3D files of the domain can be downloaded here: https://surfdrive.surf.nl/files/index.php/s/BjFEN5S7za57bkL. The model consists of several modules that describe the interactions between MATLAB and Delft3D, the post-processing of the output data from the Delft3D model and the species computations. More detailed information can be found in the attached pdf. Please reference the paper when using (parts of) the model.
The model is tested in MATLAB 2017a and Delft3D FLOW 4.03.01. 

Please send any questions or remarks to m.bruckner@exeter.ac.uk.

## What has changed since publication

The updated version removes bugs in the competition and grazing computations. The macrobenthos fractions were updated for the grazing computations and the grazing calculations were corrected. For the grazing one can now choose from two different: 1) a direct average critical bed shear stress between the values induced by MPB and CV or 2) a linear reduction in MPB fraction depending on the fraction of CV and a subsequently computed mean value of crit. bed shear stress. To select the preferred method, a flag was added to Start_model.mat.
