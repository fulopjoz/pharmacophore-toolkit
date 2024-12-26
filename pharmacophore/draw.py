"""
Script to draw pharmacophore in 2D figure
"""

import os
import matplotlib.image as img
import matplotlib.pyplot as plt
from typing import Optional
from collections import defaultdict
from cairosvg import svg2png
from IPython.display import SVG
from pharmacophore.constants import FEATURES, FEATURE_COLORS
from pharmacophore import Pharmacophore
from rdkit import Chem
from rdkit.Chem import rdDepictor, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw.MolDrawing import DrawingOptions


class Draw:
    def __init__(self, mol: Optional[Chem.Mol] = None):
        self.mol = mol

    def draw_pharm(self, mol: Chem.Mol, features: str = 'default', savepath: str = None):
        """
        Draw a 2D pharmacophore image
        :param mol: Chem.Mol
            A molecule in ROMol format.
        :param features: str
            Determines which feature algorithm to use. Can use 'default' or 'rdkit'.
        :param savepath: str
            Set the savepath to save the iamge.
        :return:
        """
        # set instance variables
        if mol is None:
            mol = self.mol

        # calculate pharmacophore features
        pharm = Pharmacophore()
        if features == 'default' or features == 'rdkit':
            calc_features = pharm.calc_pharm(mol=mol, features=features)
        else:
            raise ValueError('Unknown features! Only "default" or "rdkit" are supported!')

        # dictionaries to highlight atoms and highlight radius
        atom_highlights = defaultdict(list)
        highlight_rads = {}

        # calc features and append results to atom_highlights and highlight_rads
        for feature in calc_features:
            if feature[0] in FEATURE_COLORS:
                # extract colors from constants
                color = FEATURE_COLORS[feature[0]]
                # # for troubleshooting
                # print(color)
                for atom_id in feature[1]:
                    atom_highlights[atom_id].append(color)
                    highlight_rads[atom_id] = 0.5

        # flatten molecule into 2D
        rdDepictor.Compute2DCoords(mol)
        rdDepictor.SetPreferCoordGen(True)
        drawer = rdMolDraw2D.MolDraw2DSVG(800, 800)

        # set drawing optinos
        # Use black for all elements
        drawer.drawOptions().updateAtomPalette(
            {k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()}
        )
        drawer.SetLineWidth(2)
        drawer.SetFontSize(6.0)
        drawer.drawOptions().continuousHighlight = False
        drawer.drawOptions().splitBonds = False
        drawer.drawOptions().fillHighlights = True

        # get atom label
        for atom in mol.GetAtoms():
            atom.SetProp("atomLabel", atom.GetSymbol())
        drawer.DrawMoleculeWithHighlights(
            mol, "", dict(atom_highlights), {}, highlight_rads, {}
        )
        drawer.FinishDrawing()

        # draw molecule and save a temporary file
        svg = drawer.GetDrawingText().replace("svg:", "")
        SVG(svg)
        with open(f"pharm.svg", "w") as f:
            f.write(svg)

        # convert svg into png
        svg2png(bytestring=svg, write_to=f"image.png")

        # Set figure legend for feature and color type
        fig, (ax, picture) = plt.subplots(
            nrows=2,
            figsize=(4, 4),
            gridspec_kw={"height_ratios": [1, 5]},
        )

        # draw image and remove temporary file
        mol_image = img.imread(f"image.png")
        picture.imshow(mol_image)
        picture.axis("off")
        os.remove(f"image.png")
        os.remove(f"pharm.svg")

        # Data for the circles
        circle_radii = [0, 50, 100, 150, 200, 250]
        feature_values = list(FEATURE_COLORS.values()) # extract color values
        circle_colors = [i for i in feature_values] # match circle to color values
        circle_annotations = [
            "Donor",
            "Acceptor",
            "Aromatic",
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
        if savepath:
            plt.savefig(f"{savepath}", dpi=300)

    # support function to draw molecule with atom index
    def atom_number(self, mol: Chem.Mol = None, label: str = "atomNote", size: tuple = (300, 300)):
        """
        Draw query molecule with labeled RDKit atom indices.
        :param mol: Chem.Mol
            A molecule in ROMol format.
        :param label: str
            Determines which style to label the molecule. Defaults to atomNote. In total, can use 'atomNote',
            'atomLabel', and 'molAtomMapNumber'.
        :param size: tuple
            Determine the size of the molecule to draw.
        :return:
        """
        # Check if 2D coordinates are missing, and compute them if necessary. This will also strip away 3D coordinates.
        if not mol.GetNumConformers() or mol.GetConformer().Is3D():
            AllChem.Compute2DCoords(mol)

        # get atom indices
        for atom in mol.GetAtoms():
            if label == 'atomNote' or label == 'atomLabel' or label == 'molAtomMapNumber':
                atom.SetProp(label, str(atom.GetIdx()))
            else:
                raise ValueError("Only 'atomNote', 'atomLabel', and 'molAtomMapNumber' accepted!")

        # draw molecule
        img = Chem.Draw.MolToImage(mol, size=size)

        return img

    def similarity_maps(self, mol: Chem.Mol = None):
        pass
