import glob


for i in range(15,17):

    files = ["/data/torus_solitons/solitons_run_6_10/"+str(i)+".mp4"]
    
    files = files+sorted(glob.glob("/data/torus_solitons/solitons_run_6_10/"+str(i)+"_*.mp4"))
    
    for f in files:

        name = f.split("/")[-1]

        fnames = getFilenames("solitons_run_6_10","torus_solitons",name,"1,jpg","1.jpg")
        
        sin = torus(fnames,1,100,700)

        sin.seuil = 25
        
        sin.showVideo=0
        
        sin.openVideo()
        
        sin.runVideo()

        sin.saveData()
