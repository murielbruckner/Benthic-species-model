# Benthic-species-model
Brückner MZM, Schwarz C, Coco G, Baar A, Boechat Albernaz M, Kleinhans MG. 
Benthic species as mud patrol - modelled effects of bioturbators and biofilms
on large-scale estuarine mud and morphology. 
Earth Surf. Process. Landforms. 2021;1–17. 
https://doi.org/10.1002/esp.5080

This respository contains the model code used in the above publication with corrected competition. The competition was computed the wrong way around. The Delft3D input and setting files can be downloaded here: https://public.yoda.uu.nl/geo/UU01/GWLKT8.html. The model consists of several modules that describe the interactions between MATLAB and Delft3D, the post-processing of the output data from the Delft3D model and the species computations. More detailed information can be found in the attached pdf. Please reference the paper when using (parts of) the model.
The model is tested in MATLAB 2017a and Delft3D FLOW 4.03.01. 

Please send any questions or remarks to m.bruckner@exeter.ac.uk.

## What has changed since publication

The updated version removes minor bugs in the competition and grazing computations. The macrobenthos fractions were updated for the grazing computations and the grazing calculations were corrected. For the grazing one can now choose from two different options: 1) a direct average critical bed shear stress between the values induced by microphytobenthos and corophium volutator or 2) a linear reduction in microphytobenthos fraction depending on the fraction of corophium volutator and a subsequently computed mean value of crit. bed shear stress. To select the preferred method, a flag was added to Start_model.mat.
