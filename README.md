## Juno-M-Shell
Juno M-Shell calculated with Python 3.9.12 using JRM33+CON2020 between PJ1 (2016-08-27) and PJ56 (2023-11-22) with a 1s-resolution, limited to M < 30. Data are provided by Jonas RABIA (IRAP-CNRS, contact : jonas.rabia@irap.omp.eu) and are also available at https://amda.cdpp.eu. 

### Models and calculation parameters
- Jupiter's internal magnetic field is calculated using JRM33 with the 13th-order coefficients. The model is described in Connerney et al.,(2022).
- The magnetic field induced by the magnetodisc is calculated using the analytical version of the CON2020 model. The model is described in Connerney et al.,(2020).
- JRM33 and CON2020 models have been implemented in IDL, Matlab and Python as part of Magnetospheres of the Outer Planets Group Community Code and PSH: Planetary Spherical Harmonics community code project. They are publicly available at https://github.com/rjwilson-LASP/PSH (JRM33), https://github.com/mattkjames7/JupiterMag (JRM33) and https://github.com/gabbyprovan/con2020 (CON2020). See Wilson et al.,(2023) for further description of the numerical models.
- Magnetic field lines are traced from Juno towards the magnetic equator with a constant step size of 1/250 R<sub>J</sub>.
- Outer boundary of the calculation is set at R=30 R<sub>J</sub> because of model uncertainties beyond this distance.
  
### Data description 

Data are organized according to perijoves. Information on the dates corresponding to each perijove and on Juno's trajectory can be found at https://lasp.colorado.edu/mop/missions/juno/trajectory-information/. The data files contain 6 columns. A description of the contents of these columns is provided in the following table. 
| Column name  | Unit | Description |
| ------------ | ----------| ------------- |
| TimeUT       | / | Date, format is YYYY-MM-DDThh:mm:ss
| Orb          | / | Perijove number
| M            | / | M-Shell, radial distance of the intersection between the magnetic field line and the magnetic equator, normalized to 1 R<sub>J</sub> (1 R<sub>J</sub> = 71,492 km) 
| t            |deg| Colatitude (S3RH) of the intersection between the magnetic field line and the magnetic equator
| p            |deg| Longitude (S3RH) of the intersection between the magnetic field line and the magnetic equator
| S            |R<sub>J</sub>| Curvilinear distance along the magnetic field line between Juno and the magnetic equator


### Plot description 
Representation of the results can be found for each perijove. M-Shell values are represented along Juno's trajectory, plotted in cylindrical magnetic coordinates. CON2020 coordinate system is used, which has a z-axis tilted by 9.3° toward longitude λ<sub>S3RH</sub> = 155.8°. Large and small dots highlight Juno's position every day and hour, respectively. The corresponding day of year (DOY) is indicated for each new day. 
Magnetic field lines connected to the orbit of Io (R = 5.89 R<sub>J</sub>), Europa (R = 9.38 R<sub>J</sub>), Ganymede (R = 14.97 R<sub>J</sub>), Callisto (R = 26.3 R<sub>J</sub>) and the outer limit of the calculation (R = 30 R<sub>J</sub>) are displayed as gray lines. 

### References 
- Connerney, J. E. P., Timmins, S., Oliversen, R. J., Espley, J. R., Joergensen, J. L., Kotsiaros, S., et al. (2022). A new model of Jupiter's magnetic field at the completion of Juno's Prime Mission. Journal of Geophysical Research: Planets, 127, e2021JE007055. https://doi.org/10.1029/2021JE007055
- Connerney, J. E. P., Timmins, S., Herceg, M., & Joergensen, J. L. (2020). A Jovian magnetodisc model for the Juno era. Journal of Geophysical Research: Space Physics, 125, e2020JA028138. https://doi.org/10.1029/2020JA028138
- Wilson, R.J., Vogt, M.F., Provan, G. et al. Internal and External Jovian Magnetic Fields: Community Code to Serve the Magnetospheres of the Outer Planets Community. Space Sci Rev 219, 15 (2023). https://doi.org/10.1007/s11214-023-00961-3

 
