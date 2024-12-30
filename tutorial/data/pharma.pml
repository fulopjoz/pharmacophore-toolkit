set_color Donor_color, (0.2549019607843137, 0.4117647058823529, 0.8823529411764706)
set_color Acceptor_color, (1.0, 0.27058823529411763, 0.0)
set_color Aromatic_color, (0.8549019607843137, 0.6470588235294118, 0.12549019607843137)
set_color Hydrophobe_color, (0.1803921568627451, 0.5450980392156862, 0.3411764705882353)
set_color LumpedHydrophobe_color, (0.1803921568627451, 0.5450980392156862, 0.3411764705882353)
set_color PosIonizable_color, (0.0, 0.7490196078431373, 1.0)
pseudoatom Donor_1, pos=[-0.9046174520168838, -2.235699615827842, 0.13645934200480056]
pseudoatom Donor_2, pos=[2.949925200279259, -0.15648413497184238, 1.2254205073233995]
pseudoatom Donor_3, pos=[-2.5995394084260526, 2.931207478597561, -0.5905278269049035]
pseudoatom Aromatic_1, pos=[-1.8583650464904398, 0.2926900071440297, -0.1871054571985934]
pseudoatom Aromatic_2, pos=[-0.384216675975982, -1.2654866328983938, -0.21240652662909393]
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
set sphere_scale, 0.5
