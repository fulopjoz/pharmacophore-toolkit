"""
Script to draw pharmacophore in 2D figure
"""

import os
import matplotlib.image as img
import matplotlib.pyplot as plt
from collections import defaultdict
from cairosvg import svg2png
from IPython.display import SVG
from pharmacophore.constants import FEATURES
from pharmacophore import Pharmacophore
from rdkit.Chem import rdDepictor, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw.MolDrawing import DrawingOptions

class Draw:
    def __init__(self,):
        pass

    def pharm(self):
        Pharmacophore.calc_pharm()

    def draw_pharm(rdkit_mol, features, savepath="pharm.jpg"):
        atom_highlights = defaultdict(list)
        highlight_rads = {}
        for feature in features:
            if feature[0] in FEATURES:
                color = FEATURES[feature[0]][2]
                for atom_id in feature[1]:
                    atom_highlights[atom_id].append(color)
                    highlight_rads[atom_id] = 0.5

        rdDepictor.Compute2DCoords(rdkit_mol)
        rdDepictor.SetPreferCoordGen(True)
        drawer = rdMolDraw2D.MolDraw2DSVG(800, 800)
        # Use black for all elements
        drawer.drawOptions().updateAtomPalette(
            {k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()}
        )
        drawer.SetLineWidth(2)
        drawer.SetFontSize(6.0)
        drawer.drawOptions().continuousHighlight = False
        drawer.drawOptions().splitBonds = False
        drawer.drawOptions().fillHighlights = True

        for atom in rdkit_mol.GetAtoms():
            atom.SetProp("atomLabel", atom.GetSymbol())
        drawer.DrawMoleculeWithHighlights(
            rdkit_mol, "", dict(atom_highlights), {}, highlight_rads, {}
        )
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")
        SVG(svg)
        with open(f"pharm.svg", "w") as f:
            f.write(svg)

        svg2png(bytestring=svg, write_to=f"image.png")

        fig, (ax, picture) = plt.subplots(
            nrows=2,
            figsize=(4, 4),
            gridspec_kw={"height_ratios": [1, 5]},
        )

        mol_image = img.imread(f"image.png")
        picture.imshow(mol_image)
        picture.axis("off")
        os.remove(f"image.png")
        os.remove(f"pharm.svg")

        # Data for the circles
        circle_radii = [0, 50, 100, 150, 200, 250]
        feature_values = list(FEATURES.values())
        circle_colors = [i[2] for i in feature_values]
        circle_annotations = [
            "Donor",
            "Acceptor",
            "Cation",
            "Anion",
            "Ring",
            "Hydrophobe",
        ]
        # Draw the circles and annotations
        for radius, color, annotation in zip(
            circle_radii, circle_colors, circle_annotations
        ):
            x = radius
            circle = plt.Circle((x, -5), 5, color=color)  # , alpha=0.5)
            ax.add_patch(circle)
            ax.annotate(
                annotation,
                (x, 10),
                va="center",
                ha="center",
                fontsize=6,
                fontweight="bold",
            )

        # Set axis limits
        ax.set_xlim(-10, 270)
        ax.set_ylim(-20, 20)
        ax.axis("off")

        # Set aspect ratio to equal
        ax.set_aspect("equal", adjustable="box")
        plt.savefig(f"{savepath}", dpi=300)

    # support function to draw molecule with atom index
    def atom_number(mol, label):
        # Check if 2D coordinates are missing, and compute them if necessary
        if not mol.GetNumConformers() or mol.GetConformer().Is3D():
            AllChem.Compute2DCoords(mol)

        for atom in mol.GetAtoms():
            if label == 'atomNote' or label == 'atomLabel' or label == 'molAtomMapNumber':
                atom.SetProp(label, str(atom.GetIdx() + 1))
            else:
                raise ValueError("Only 'atomNote', 'atomLabel', and 'molAtomMapNumber' accepted!")
        return mol
