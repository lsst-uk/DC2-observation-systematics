# DC2-observation-systematics

This repository contains codes and notebooks to look at the correlation between observational systematics and galaxy photometry and mean redshifts.

A more detailed documentation as well as various results can be found in this Overleaf note: https://www.overleaf.com/read/vjtwzsrmshzm

---

The [codes](https://github.com/lsst-uk/DC2-observation-systematics/tree/main/codes) folder contains some functions used for the analysis. 
- [measure_properties_with_systematics.py](https://github.com/lsst-uk/DC2-observation-systematics/blob/main/codes/measure_properties_with_systematics.py) contains some customised functions to loop over the tracts and selet data depending on the selected pixels given the systematics map. There is also a function to de-redden the galaxy magnitudes.
- [run_prop.py](https://github.com/lsst-uk/DC2-observation-systematics/blob/main/codes/run_prop.py) can be run to generate Fig.6 - 12 in the documentation.
- Dependences: Systematic maps: You can find the MAF maps on NERSC: `/global/cscratch1/sd/qhang/minion_1016/MAF-[1/5/10]year/`. Alternatively, these maps can be generated following the notebook: [DESC_DC2_minion_1016.ipynb](https://github.com/lsst-uk/DC2-observation-systematics/blob/main/notebooks/DESC_DC2_minion_1016.ipynb). Notice that to run this the `rubin` package needs to be installed: [rubin_sim](https://github.com/lsst/rubin_sim).

- Dependences: DC2 objects: For the DC2 catalogues with some basic cuts and `mag_i_cModel<25.3`, the files can be found here: `/global/cscratch1/sd/qhang/DESC_DC2_obs-dr[2/6]/`. Alternatively you can simply access through `GCRCatalogs` (this is included in the `DESC-python` environment), following e.g. [DC2_download_obj_with_pz-dr6.ipynb](https://github.com/lsst-uk/DC2-observation-systematics/blob/main/notebooks/DC2_download_obj_with_pz-dr6.ipynb).

---

The [notebooks](https://github.com/lsst-uk/DC2-observation-systematics/tree/main/notebooks) folder contains the Jupyter notebooks used to generate some of the result plots.
- [DESC_DC2_minion_1016.ipynb](https://github.com/lsst-uk/DC2-observation-systematics/blob/main/notebooks/DESC_DC2_minion_1016.ipynb) generates the MAF maps.
- [DC2_download_obj_with_pz-dr6.ipynb](https://github.com/lsst-uk/DC2-observation-systematics/blob/main/notebooks/DC2_download_obj_with_pz-dr6.ipynb) downloads catalogues of DC2 objects with simple selection and `mag_i_cModel<25.3`.

- [DC2_match_OpSim-dr6-multibands.ipynb](https://github.com/lsst-uk/DC2-observation-systematics/blob/main/notebooks/DC2_match_OpSim-dr6-multibands.ipynb) visualise the systematic maps and shows example corelations with magnitude error and mean photo-z.
- [DC2_galnum_density-dr6-multibands.ipynb](https://github.com/lsst-uk/DC2-observation-systematics/blob/main/notebooks/DC2_galnum_density-dr6-multibands.ipynb) looks at the correlation between systematics and galaxy number density near the $i$ band detection limit.

- [Compare_MAF_supreme.ipynb](https://github.com/lsst-uk/DC2-observation-systematics/blob/main/notebooks/Compare_MAF_supreme.ipynb) compares the MAF and [supreme maps](https://confluence.slac.stanford.edu/display/LSSTDESC/DC2+Run2.2i+DR6+Survey+Property+Maps).

- [test_BPZ_run22i.ipynb](https://github.com/lsst-uk/DC2-observation-systematics/blob/main/notebooks/test_BPZ_run22i.ipynb) is an example notebook to estimate photo-z on a subsample of DC2 objects using `BPZ_lite`. This is used for the dr2 data.

For more questions contact e.hang@ucl.ac.uk.

