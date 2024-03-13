## Juno-M-Shell
Juno M-Shell calculated with Python 3.9.12 using JRM33+CON2020 between PJ1 (2016-08-27) and PJ56 (2023-11-22) with a 1s-resolution. Data, plots and numerical models used for these calculations are provided by Jonas RABIA (IRAP-CNRS, contact : jonas.rabia@irap.omp.eu). Data are also available at https://amda.cdpp.eu. 

**Dataset DOI:** [10.5281/zenodo.10813572 ](https://zenodo.org/doi/10.5281/zenodo.10813572)

### Models and calculation parameters
- Jupiter's internal magnetic field is calculated using JRM33 with the 13th-order coefficients. The model is described in Connerney et al.,(2022). Use of a 18th-order model may provide a better match with Juno observations, as reported by Moirano et al., (2024, doi: 10.1029/2023JE008130) but would increase computation time. 
- The magnetic field induced by the current sheet is calculated using an analytical version of the CON2020 model. The model is described in Connerney et al.,(2020). B<sub>ρ</sub> and B<sub>z</sub> components are calculated using the analytical equations provided in the Appendix of Connerney et al.,(1981). B<sub>φ</sub> calculation is done using a routine from the CON2020 Community code (https://github.com/gabbyprovan/con2020). The magnetodisc parameters that best fit the measurements obtained during the 24 first Juno orbits are used, i.e., R<sub>0</sub>, R<sub>1</sub>, D, μ<sub>0</sub>I/2, μ<sub>0</sub>I<sub>R</sub>/2π, θ<sub>d</sub>, φ<sub>d<sub>S3RH</sub></sub> = 7.8 R<sub>J</sub>, 51.4 R<sub>J</sub>, 3.6 R<sub>J</sub>, 139.6 nT, 16.7 MA, 9.3°, 155.8°.  
- The JRM33 model used for these calculations is adapted from a publicly available model that can be found at https://github.com/rjwilson-LASP/PSH (Wilson et al., 2023). The JRM33-order 13 model, the analytical CON2020 model and the routine used for the M-Shell calculations are provided in the Codes directory.
- M-Shell is found by iteratively tracing the magnetic field lines towards the magnetic equator,  with a constant step size of 1/250 R<sub>J</sub>, until reaching the minimum magnetic field strength.
- Outer boundary of the calculation is set at R=30 R<sub>J</sub> because of model accuracy beyond this distance (see Connerney et al., 2020 and Rabia et al., 2024 for further details).
- Due to the limitation in the spatial resolution of the computation performed (constant size steps of 1/250 R<sub>J</sub>), small stalls in colatitude, longitude, and curvilinear distance may exceptionally occur close to the magnetic equator when displaying the data. 


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
- Wiggle plot : M-Shell values are represented along Juno's trajectory, plotted in cylindrical magnetic coordinates, for each perijove. CON2020 coordinate system is used, which has a z-axis tilted by 9.3° toward longitude λ<sub>S3RH</sub> = 155.8°. Large and small dots highlight Juno's position every day and hour, respectively. The corresponding day of year (DOY) is indicated for each new day. 
Magnetic field lines connected to the orbit of Io (R = 5.89 R<sub>J</sub>), Europa (R = 9.38 R<sub>J</sub>), Ganymede (R = 14.97 R<sub>J</sub>), Callisto (R = 26.3 R<sub>J</sub>) and the outer limit of the calculation (R = 30 R<sub>J</sub>) are displayed as gray lines.

- Time series : M-Shell values are represented as a function of time, for each perijove. 

### References 
- Connerney, J. E. P., M. H. Acuña, and N. F. Ness (1981), Modeling the Jovian current sheet and inner magnetosphere, J. Geophys. Res., 86(A10), 8370–8384, https://doi.org/10.1029/JA086iA10p08370
- Connerney, J. E. P., Timmins, S., Herceg, M., & Joergensen, J. L. (2020). A Jovian magnetodisc model for the Juno era. Journal of Geophysical Research: Space Physics, 125, e2020JA028138. https://doi.org/10.1029/2020JA028138
- Connerney, J. E. P., Timmins, S., Oliversen, R. J., Espley, J. R., Joergensen, J. L., Kotsiaros, S., et al. (2022). A new model of Jupiter's magnetic field at the completion of Juno's Prime Mission. Journal of Geophysical Research: Planets, 127, e2021JE007055. https://doi.org/10.1029/2021JE007055
- Wilson, R.J., Vogt, M.F., Provan, G. et al. Internal and External Jovian Magnetic Fields: Community Code to Serve the Magnetospheres of the Outer Planets Community. Space Sci Rev 219, 15 (2023). https://doi.org/10.1007/s11214-023-00961-3
- Rabia, J., Nénon, Q., André, N., Hue, V., Santos-Costa, D., Kamran, A., & Blanc, M. (2024). Influence of the Jovian current sheet models on the mapping of the UV auroral footprints of Io, Europa, and Ganymede. Journal of Geophysical Research: Space Physics, 129, e2023JA032041. https://doi.org/10.1029/2023JA032041 
 
