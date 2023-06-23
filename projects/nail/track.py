import matplotlib.pyplot as plt
import numpy as np
import cv2
import sys

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                            TRACKING OF THE NAIL                           ##
##                                                                           ##
## ============================================================================

"""
Detection of a white point on dark background and tracking of its 
position in time. The input is a mp4 video, and the result is an 
array of the form [[x0,y0],...]
"""



############################      DEFINITIONS      ############################


showVideo = input("Show video? (answer 0 or 1) \n")
showVideo = int(showVideo)

threshold = 90           # Threshold for grayscale
rc = 550                 # Crop radius

############################        LOAD FILE      ############################


path = "/data/nail/12_06/2_Hz_A6Vpp.mp4"


video = cv2.VideoCapture(path,0)

numFrames = int(video.get(cv2.CAP_PROP_FRAME_COUNT))
fps = int(video.get(cv2.CAP_PROP_FPS))


############################        TRACKING      ############################


x = np.zeros(numFrames)
y = np.zeros(numFrames)

N = numFrames

for i in range(N):

    ret, frame = video.read()

    height, width = frame.shape[:2]

    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)                # Convert to grayscale

    mask = np.zeros_like(gray)

    cv2.circle(mask,                                              # Circular crop with radius rc
               (int(width/2),
                int(height/2)),
               rc,1,-1,8,0)

    gray = np.where(gray>threshold,255,0)
    gray = gray.astype(np.uint8)
    gray *= mask

    cnts, hierarchy = cv2.findContours(gray,                      # Contour detection
                                       cv2.RETR_TREE,
                                       cv2.CHAIN_APPROX_NONE)

    cMax = max(cnts,key=cv2.contourArea)                          # Get the contour with maximal area

    shape = np.vstack(cMax).squeeze()
    
    x[i] = np.mean(shape[:,0])
    y[i] = np.mean(shape[:,1])


    if showVideo == 1:
    
        blank = np.ones((height,width,3),np.uint8)
        blank = cv2.cvtColor(gray,cv2.COLOR_GRAY2RGB)

        cv2.drawContours(blank, [cMax], 0, (255,0,0),5)

        cv2.namedWindow('shown',cv2.WINDOW_NORMAL)
        cv2.imshow('shown',blank)
        cv2.resizeWindow('shown',800,600)
        cv2.waitKey(10)

    sys.stdout.flush()
    sys.stdout.write("Progress: "+str(int(100*i/numFrames))+" %\r")

video.release()


############################        SAVE FILE      ############################



coord = np.array(list(zip(x,y)))

name = path.split(".mp4")
name = name[0]

np.save(name,coord)


############################       PLOT COORDS     ############################



# t = np.arange(0,N)/fps

# plt.plot(t,y,".")
# plt.plot(t,x,".")

# plt.show()
