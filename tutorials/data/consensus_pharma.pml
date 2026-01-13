set_color Donor_color, (0.2549019607843137, 0.4117647058823529, 0.8823529411764706)
set_color Acceptor_color, (1.0, 0.27058823529411763, 0.0)
set_color Aromatic_color, (0.8549019607843137, 0.6470588235294118, 0.12549019607843137)
set_color Hydrophobe_color, (0.1803921568627451, 0.5450980392156862, 0.3411764705882353)
set_color LumpedHydrophobe_color, (0.1803921568627451, 0.5450980392156862, 0.3411764705882353)
set_color PosIonizable_color, (0.0, 0.7490196078431373, 1.0)
pseudoatom Donor_1, pos=[-4.615622269567313, 0.1451455783503386, 0.5863201578506856]
pseudoatom Donor_2, pos=[-2.9688883127180237, 2.1191857624459898, 1.4756023911236902]
pseudoatom Donor_3, pos=[2.939045487574166, -0.815396252511792, 0.7394092666299179]
pseudoatom Aromatic_1, pos=[-1.88430386157689, -0.014552394821863308, 0.06132409850125146]

show spheres, Acceptor_*
color acceptor_color, Acceptor_*

show spheres, Donor_*
color donor_color, Donor_*

show spheres, Hydrophobe_*
color hydrophobe_color, Hydrophobe_*

show spheres, Aromatic_*
color aromatic_color, Aromatic_*

show spheres, LumpedHydrophobe_*
color lumpedhydrophobe, LumpedHydrophobe_*

show spheres, PosIonizable*
color posionizable, PosIonizable*

set sphere_scale, 0.8
