# A Simple Pharmacophore-Toolkit

![Static Badge](https://img.shields.io/badge/Pharmacophore--Toolset-v0.0.3-blue)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/py50?style=flat&logo=python&logoColor=white)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

The Pharmacophore-Toolkit is built on RDKit and allows for building simple pharmacophore models. The
Pharmacophore-Toolset can generate models from crystal structures, docking poses, or SMILES string. To generate a 3D
model, a .pml file will be generated. This files contains scripts to generate spheres with color and XYZ coordinates 
defined. The final 3D image can be rendered in PyMOL. Additional information can be found under the
[Tutorials](tutorials/) section. The Pharmacophore-Toolkit can result in two types:

### 3D Model

<figure>
    <img src="img/3d_example.png" width="400">
    <figcaption>3D conformation of molecules and alignment was performed using RDKit. Spheres were generated based on 
                Serotonin. Blue spheres represent hydrogen donors and gold spheres represent aromatic rings. The 3D 
                images were generated in PyMOL.
    </figcaption>
</figure>

### 2D Model
<figure>
    <img src="img/similarity_exmaple.png" width="400">
    <figcaption>Similarity map of molecules. All molecules were compared to Serotonin. More information can be seen 
                at RDKit documentation. 
    </figcaption>
</figure>

## Tutorials

Tutorials are written as JupyterNotebooks and can be found [here](tutorials/).

## Install
You can install the Pharmacophore-Toolkit using pip:
```
pip install pharmacophore-toolkit
```

Alternatively, the environemnt can be created by cloning the repository and then running the following conda script:
```
conda env create -f environment.yaml
```

> [!NOTE]
> **Note:** The Pharmacophore-Toolkit relies on cairosvg to create images before being converted into .png format. Depending on your 
workstation/machine, the [CairoSVG](https://github.com/Kozea/CairoSVG) package will need to be installed manually. 
Installation instructions can be found [here](https://cairosvg.org). If it is already installed on your machine globally, 
on my machine, I used conda to install cairosvg:
```
conda install conda-forge::cairosvg
```

This solved the problem on my machine, but may differ from yours.
