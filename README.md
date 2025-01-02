# A Simple Pharmacophore-Toolkit

![Static Badge](https://img.shields.io/badge/Pharmacophore--Toolset-v0.1.0-blue)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

The Pharmacophore-Toolkit is built on RDKit and allows for building simple pharmacophore models. The
Pharmacophore-Toolset can generate models from crystal structures, docking poses, or SMILES string. To generate a 3D
model, a .pml file will be generated. This files contains scripts to generate spheres with color and XYZ coordinates 
defined. The final 3D image can be rendered in PyMOL. Additional information can be found under the
[Tutorials](tutorial/) section. The Pharmacophore-Toolkit can result in two types:

### 3D Model

<figure>
    <img src="img/3d_example.png" width="400">
    <figcaption>3D conformation of molecules and alignment was performed using RDKit. Spheres were generated based on 
                Serotonin. Blue spheres represent hydrogen donors and gold spheres represent aromatic rings. 
    </figcaption>
</figure>

The 3D images were generated in PyMOL. 

### 2D Model
<figure>
    <img src="img/similarity_exmaple.png" width="400">
    <figcaption>Similarity map of molecules. All molecules were compared to Serotonin. </figcaption>
</figure>

## Tutorials

Tutorials are written as JupyterNotebooks and can be found [here](tutorial/).

## Install
You can install the Pharmacophore-Toolkit using pip:
```python
pip install pharmacophore-toolkit
```

or directly from the github repository as follows:
```python
pip install git+https://github.com/tlint101/pharmacophore-toolkit.git
```

