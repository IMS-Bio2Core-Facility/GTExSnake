# Changelog

<!--next-version-placeholder-->

## v1.2.0 (2021-10-27)
### Feature
* **workflow:** Support snakemake report ([`2cb199f`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/2cb199fcb827737b745a5b4f71ccdb5c9bc5119f))
* **repo:** Comply to standardised usage ([`a860777`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/a86077797d23c61abc086f4c764b6f48cc5cacec))

### Fix
* **scripts:** Correct import ([`170d0a3`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/170d0a3a0b4b0af8ec0e4b28a6735ef63482d7c0))

## v1.1.2 (2021-07-21)
### Fix
* **deps:** Migrate to gtexquery v2.0.0 ([`d5ff215`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/d5ff21532aa01cf26c4141ed0f5fab974fa8d82c))

## v1.1.1 (2021-07-20)


## v1.1.0 (2021-07-15)
### Feature
* **workflow:** Add schema validation ([`0b35f95`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/0b35f95e8155fc7cea6b95223f66ec3ca262ace6))

## v1.0.0 (2021-07-14)
### Feature
* **workflow:** Provides resources in repository ([`d5ff4f2`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/d5ff4f2dcdc4c22ec1dbebf702000692f23c6c94))

### Fix
* **workflow:** Migrate to gtexquery v1.0.0 ([`82d5324`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/82d5324fc77ce198cd2ee375da6e662ce1388e51))

### Breaking
* While the output of the pipeline remains unchanged, users who do not already have the resources locally must either generate them or pull them from the repository. As this changes the input requirements of the pipeline, I have decided to class this as Breaking.  ([`d5ff4f2`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/d5ff4f2dcdc4c22ec1dbebf702000692f23c6c94))
* The output is now a collection of csv files, one per gene. As mentioned above, this eliminates the need to try to test excel file. It also greatly simplifies multithreading. Closes #1.  ([`82d5324`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/82d5324fc77ce198cd2ee375da6e662ce1388e51))

### Documentation
* **CONTRIBUTING:** Document snakemake tests ([`4fd39e7`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/4fd39e7313494b1ff88c974d5bb9eb5099a73f4f))
* **docs:** Update link targets ([`c467b3f`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/c467b3fb3038eeb0f956dc773918bba66cca4244))

## v0.1.0 (2021-07-13)
### Feature
* **scripts:** Add analysis scripts to workflow ([`ead7055`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/ead7055b846373feab8b3bd39443f11eef459773))
* **repo:** Initialise repository ([`a020dda`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/a020ddabd214b4739a9a1f50803f25512fccbfdb))

### Documentation
* **images:** Add rulegraph for pipeline ([`5938b20`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/5938b20aa258429972739756cab4d6b58434671f))
* **scripts:** Add documentation for scripts ([`8a2815f`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/8a2815fab3ef6463629caa5524bf1171b4528409))
* **repo:** Conform to commonmark standards ([`40425e8`](https://github.com/IMS-Bio2Core-Facility/GTExSnake/commit/40425e8216c1b00d017e05e77c9a2f9a7c626952))
