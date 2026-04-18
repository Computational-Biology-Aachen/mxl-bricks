<p align="center">
    <img src="https://raw.githubusercontent.com/Computational-Biology-Aachen/mxl-bricks/refs/heads/main/docs/assets/logo.png" width="400px" alt='mxlbricks-logo'>
</p>

# MxlBricks

[![pypi](https://img.shields.io/pypi/v/mxlbricks.svg)](https://pypi.python.org/pypi/mxlbricks)
[![docs][docs-badge]][docs]
![License](https://img.shields.io/badge/license-MIT-blue?style=flat-square)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![security: bandit](https://img.shields.io/badge/security-bandit-yellow.svg)](https://github.com/PyCQA/bandit)
[![PyPI Downloads](https://static.pepy.tech/badge/mxlbricks)](https://pepy.tech/projects/mxlbricks)

[docs-badge]: https://img.shields.io/badge/docs-main-green.svg?style=flat-square
[docs]: https://computational-biology-aachen.github.io/mxl-bricks/

MxlBricks is a Python package to build mechanistic models composed of pre-defined reactions (bricks). This facilitates re-use and interoperability between different models by sharing common parts.

## Installation


You can install mxlpy using pip: `pip install mxlbricks`.


If you want access to the sundials solver suite via the [assimulo](https://jmodelica.org/assimulo/) package, we recommend setting up a virtual environment via [pixi](https://pixi.sh/) or [mamba / conda](https://mamba.readthedocs.io/en/latest/) using the [conda-forge](https://conda-forge.org/) channel.

```bash
pixi init
pixi add python assimulo
pixi add --pypi mxlbricks
```


## Development setup

Install pixi [as described in the docs](https://pixi.sh/latest/#installation).

Run

```bash
pixi install
```


## Models

| Name                 | Description                                                                 |
| -------------------- | --------------------------------------------------------------------------- |
| Ebenhöh 2011         | PSII & two-state quencher & ATP synthase                                    |
| Ebenhöh 2014         | PETC & state transitions & ATP synthase from Ebenhoeh 2011                  |
| Matuszyńska 2016 NPQ | 2011 + PSII & four-state quencher                                           |
| Matuszyńska 2016 PhD | ?                                                                           |
| Matuszyńska 2019     | Merges PETC (Ebenhöh 2014), NPQ (Matuszynska 2016) and CBB (Poolman 2000)   |
| Saadat 2021          | 2019 + Mehler (Valero ?) & Thioredoxin & extendend PSI states & consumption |
| van Aalst 2023       | Saadat 2021 & Yokota 1985 & Witzel 2010                                     |


## References

| Name         | Description                                           |
| ------------ | ----------------------------------------------------- |
| Poolman 2000 | CBB cycle, based on Pettersson & Ryde-Pettersson 1988 |
| Yokota 1985  | Photorespiration                                      |
| Valero ?     |                                                       |

## Tool family 🏠

`MxlBricks` is part of a larger family of tools that are designed with a similar set of abstractions. Check them out!

- [MxlPy](https://github.com/Computational-Biology-Aachen/MxlPy) is a Python package for mechanistic learning (Mxl)
- [MxlWeb](https://github.com/Computational-Biology-Aachen/mxl-web) brings simulation of mechanistic models to the browser!
- [MxlModels](https://github.com/Computational-Biology-Aachen/mxl-models) supplies you with flat, easy-to-use and inspect models.
