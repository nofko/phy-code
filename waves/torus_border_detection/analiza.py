

fname = "/data/torus_mi/6_10_2022/sweep_67hz_A_500mv.mp4"


sin = torus(fname,1,500,800)

sin.openVideo()

sin.seuil = 20
sin.showVideo = 1

sin.runVideo()

sin.saveData()



