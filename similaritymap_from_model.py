"""
Get similarity map from ML model
There is an issue with the code. Copied from here:
https://github.com/rdkit/rdkit/issues/5087
"""

from rdkit import Chem
from rdkit.Chem import Draw
import numpy
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem.Draw import rdMolDraw2D
import xgboost as xgb
import io
from PIL import Image


def show_png(data):
    bio = io.BytesIO(data)
    img = Image.open(bio)
    return img

def get_pred(fp, pred_function):
    fp = numpy.array([list(fp)])
    return pred_function(fp)[0]


def plot_similarity_map(mol, model):
    d = Draw.MolDraw2DCairo(900, 900)
    fig , maxweight = SimilarityMaps.GetSimilarityMapForModel(mol, SimilarityMaps.GetMorganFingerprint, lambda x : get_pred(x, model.predict), draw2d=d)
    d.FinishDrawing()
    return d

model = xgb.XGBClassifier()
model.load_model('cyp1a2.h5')

test_smiles = 'Nc5ncc(c2nc(N1CCOCC1)nc3c2CCC3c4ccncc4)cn5'
mol = Chem.MolFromSmiles(test_smiles)


test_fps = numpy.array(AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=2048))
test_fps = test_fps.reshape(-1, 2048)

activity_prediction = model.predict(test_fps)[0]
print(f'The activity class is : {activity_prediction}')

res = plot_similarity_map(mol, model)
activity_sim_map = show_png(res.GetDrawingText())
activity_sim_map.save('SimMap.png')

